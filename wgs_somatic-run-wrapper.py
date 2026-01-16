"""
Wrapper to be used by cron for automatic start of wgs_somatic
"""

import sys
import argparse
import os
import re
import glob
from datetime import datetime
import json
import threading

from definitions import WRAPPER_CONFIG_PATH, ROOT_DIR, LAUNCHER_CONFIG_PATH #, INSILICO_CONFIG, INSILICO_PANELS_ROOT
from tools.context import RunContext, SampleContext
from tools.helpers import setup_logger, read_config
from tools.slims import get_sample_slims_info, find_or_download_fastqs, get_pair_dict, link_fastqs_to_outputdir
from tools.custom_email import start_email, end_email, error_email, error_admin_qc_email, error_setup_email
from launch_snakemake import analysis_main, yearly_stats, copy_results, get_timestamp
from tools.wgs_admin_summary.combine_wgsadmin_qc_summary import combine_qc_stats


def look_for_runs(config, instrument):
    """Look for runs in demultiplexdir"""
    instrument_root_path = config[instrument]["demultiplex_path"]
    found_paths = glob.glob(os.path.join(instrument_root_path, "*"))
    regex = config[instrument]["seq_name_regex"]
    return [path for path in found_paths if re.search(regex, os.path.basename(path))]


def generate_context_objects(Rctx, logger):
    """Create Rctx and Sctx for a demultiplexed run"""

    # Read demultiplex stats file for sample names, fastq paths, and nr reads
    with open(Rctx.demultiplex_summary_path, "r") as inp:
        demuxer_info = json.load(inp)

    for sample_id, sample_info in demuxer_info["samples"].items():
        logger.info(f"Setting up context for {sample_id}.")

        # Setup Sample context class and add listed fastq paths
        Sctx = SampleContext(sample_id)
        Sctx.add_fastq(sample_info["fastq_paths"])

        # Query Slims for clinical information and add to sample context
        logger.info("Fetching SLIMS info.")
        Sctx.slims_info = get_sample_slims_info(Sctx, run_tag=Rctx.run_tag)

        if not Sctx.slims_info:
            # If no slims info is found/available we need to be notified
            # The sample has been added to the runlist without running potential wgs-somatic samples
            logger.error(f"No SLIMS info available for {sample_id}!")
            raise ValueError(f"No SLIMS info available for {sample_id}!")

        # NOTE: 54 is Slims internal primary key for wgs_somatic
        if 54 not in Sctx.slims_info["secondary_analysis"]:
            logger.info("Sample not set for wgs_somatic.")
            continue

        # Add sample context to run context list
        logger.info("Sample set for wgs_somatic, adding to RunContext.")
        Rctx.add_sample_context(Sctx)

    return Rctx


def return_first_new_run(config, instrument, logger):
    """Return first new run for given instrument, which has files for wgs_somatic analysis."""
    local_run_paths = look_for_runs(config, instrument)

    previous_runs_file = config["previous_runs_file_path"]
    previous_runs_file_path = os.path.join(ROOT_DIR, previous_runs_file)
    if not os.path.exists(previous_runs_file_path):
        raise FileNotFoundError(f"Runlist file not found: {previous_runs_file_path}")

    with open(previous_runs_file_path, "r") as prev:
        previous_runs = [line.rstrip() for line in prev]

    for run_path in local_run_paths:
        Rctx = RunContext(run_path)
        # Check if demultiplexing is completed
        if not Rctx.demultiplex_complete:
            continue

        # Check if run has been previously analysed
        if Rctx.run_name in previous_runs:
            continue

        # Write run name to previously analysed list to ensure no double-running
        with open(previous_runs_file_path, "a") as prev:
            logger.info(f"Writing {Rctx.run_name} to previous runs list.")
            print(Rctx.run_name, file=prev)

        Rctx_run = generate_context_objects(Rctx, logger)

        if not Rctx_run.sample_contexts:
            logger.info(f"No samples for wgs_somatic found in run {Rctx.run_name}.")
            continue

        return Rctx_run


def call_script(**kwargs):
    '''Function to call main function from launch_snakemake.py'''
    args = argparse.Namespace(**kwargs)
    analysis_main(args, **kwargs)


def check_ok(outputdir):
    '''Function to check if analysis has finished correctly'''

    if os.path.isfile(f"{outputdir}/workflow_finished.txt"):
        return True
    else:
        return False


def analysis_end(outputdir, tumorsample=None, normalsample=None):
    '''Function to check if analysis has finished correctly and add to yearly stats and copy results'''

    if check_ok(outputdir):
        if tumorsample:
            if normalsample:
                # these functions are only executed if snakemake workflow has finished successfully
                yearly_stats(tumorsample, normalsample)
            else:
                yearly_stats(tumorsample, 'None')
        else:
            yearly_stats('None', normalsample)
        copy_results(outputdir)
    else:
        pass


def submit_pipeline(tumorsample, normalsample, outpath, config, logger, threads):
    timestamp = get_timestamp()
    if tumorsample and normalsample:
        logger.info(f'Preparing run: Tumor {tumorsample} and Normal {normalsample}')
        fastq_dict_tumor = find_or_download_fastqs(tumorsample, logger)
        fastq_dict_normal = find_or_download_fastqs(normalsample, logger)
        tumorid = list(fastq_dict_tumor.keys())[0]  # E.g. DNA123456_250101_AHJLJHBGXF
        outputdir = os.path.join(outpath, f"{tumorid}_{timestamp}")
        os.makedirs(outputdir, exist_ok=False)  # Make sure a new outputdir is created, not overwriting old results
        tumor_fastq_dir = link_fastqs_to_outputdir(fastq_dict_tumor, outputdir, logger)
        normal_fastq_dir = link_fastqs_to_outputdir(fastq_dict_normal, outputdir, logger)
        pipeline_args = {'outputdir': f'{outputdir}',
                         'normalname': f'{normalsample}',
                         'normalfastqs': f'{normal_fastq_dir}',
                         'tumorname': f'{tumorsample}',
                         'tumorfastqs': f'{tumor_fastq_dir}'}

    elif tumorsample:
        logger.info(f'Preparing run: Tumor-only {tumorsample}')
        fastq_dict_tumor = find_or_download_fastqs(tumorsample, logger)
        outpath = os.path.join(outpath, "tumor_only")
        os.makedirs(outpath, exist_ok=True)
        tumorid = list(fastq_dict_tumor.keys())[0]
        outputdir = os.path.join(outpath, f"{tumorid}_{timestamp}")
        os.makedirs(outputdir, exist_ok=False)
        tumor_fastq_dir = link_fastqs_to_outputdir(fastq_dict_tumor, outputdir, logger)
        pipeline_args = {'outputdir': f'{outputdir}',
                         'tumorname': f'{tumorsample}',
                         'tumorfastqs': f'{tumor_fastq_dir}'}

    elif normalsample:
        logger.info(f'Preparing run: Normal-only {normalsample}')
        fastq_dict_normal = find_or_download_fastqs(normalsample, logger)
        outpath = os.path.join(outpath, "normal_only")
        os.makedirs(outpath, exist_ok=True)
        normalid = list(fastq_dict_normal.keys())[0]
        outputdir = os.path.join(outpath, f"{normalid}_{timestamp}")
        os.makedirs(outputdir, exist_ok=False)
        normal_fastq_dir = link_fastqs_to_outputdir(fastq_dict_normal, outputdir, logger)
        pipeline_args = {'outputdir': f'{outputdir}',
                         'normalname': f'{normalsample}',
                         'normalfastqs': f'{normal_fastq_dir}'}

    threads.append(threading.Thread(target=call_script, kwargs=pipeline_args))
    logger.info(f'Starting wgs_somatic with arguments {pipeline_args}')
    return outputdir


def wrapper(instrument=None, outpath=None):
    '''Automatic wrapper function'''
    
    # === Setup run ===
    try:
        config = read_config(WRAPPER_CONFIG_PATH)
        wrapper_log_path = config["wrapper_log_path"]
        logger = setup_logger('wrapper', os.path.join(wrapper_log_path, f'{instrument}_WS_wrapper.log'))

        # Empty dict, will update later with T/N pair info
        pair_dict_all_pairs = {}

        # prepare hcp download directory
        hcptmp = config["hcp_download_dir"]
        if not os.path.isdir(hcptmp):
            try:
                os.makedirs(hcptmp)
            except Exception as e:
                logger.error(f"outputdirectory: {hcptmp} does not exist and could not be created")
                raise e

        # If outputpath is not specified, get from config
        if not outpath:
            try:
                outpath = config['cron_outpath']
            except KeyError:
                logger.error('Output path for cron job not specified in the configuration.')
                raise ValueError('Output path for cron job not specified in the configuration.')

        # NOTE: we only process one run at a time. The next run will start upon the next cron execution
        Rctx = return_first_new_run(config, instrument, logger)

        if Rctx is None:
            sys.exit(0)

        logger.info(f'Found {len(Rctx.sample_contexts)} samples for wgs_somatic in run {Rctx.run_name}.')

        # Register start time
        start_time = datetime.now()
        logger.info(f'Started {Rctx.run_name} at {start_time}')

        # Get T/N pair info in a dict for samples and link additional fastqs from other runs
        for sctx in Rctx.sample_contexts:
            pair_dict = get_pair_dict(sctx, Rctx, logger)
            pair_dict_all_pairs.update(pair_dict)
            logger.info(f'Sample {sctx.sample_name} info: {pair_dict[sctx.sample_name]}')

        # Uses the dictionary of T/N samples to put the correct pairs together and finds the correct input arguments to the pipeline
        threads = []
        outputdirs = []
        end_threads = []
        final_pairs = []
        paired_samples = []

        tumor_samples = {}
        normal_samples = {}

        # Separate tumor and normal samples
        for key, value in pair_dict_all_pairs.items():
            if value[0] == 'tumor':
                tumor_samples[key] = value
            elif value[0] == 'normal':
                normal_samples[key] = value

        logger.info(f"tumor_samples: {tumor_samples}")
        logger.info(f"normal_samples: {normal_samples}\n")
        # Pair samples based on tumorNormalID
        for t_key, t_value in tumor_samples.items():
            t_ID = t_value[1]  # tumorNormalID
            paired = False
            for n_key, n_value in normal_samples.items():
                n_ID = n_value[1]  # tumorNormalID
                if t_ID == n_ID or t_ID == n_key.split("DNA")[1] or n_ID == t_key.split("DNA")[1]:
                    paired_samples.append((t_key, n_key))
                    paired = True
                    outputdir = submit_pipeline(t_key, n_key, outpath, config, logger, threads)
                    outputdirs.append(outputdir)
                    end_threads.append(threading.Thread(target=analysis_end, args=(outputdir, t_key, n_key)))
                    final_pairs.append(f'{t_key} (T) {n_key} (N), {n_value[2]} {["prio" if (n_value[3] or t_value[3]) else ""][0]}')
                    break

            if not paired:
                outputdir = submit_pipeline(t_key, None, outpath, config, logger, threads)
                outputdirs.append(outputdir)
                end_threads.append(threading.Thread(target=analysis_end, args=(outputdir, t_key, None)))
                final_pairs.append(f'{t_key} (T), {t_value[2]} {["prio" if t_value[3] else ""][0]}')

        # We are currently not running the normal-only samples
        for n_key, n_value in normal_samples.items():
            if not any(n_key == pair[1] for pair in paired_samples):
                logger.info(f'Skipping normal-only sample {n_key} as it is not paired with a tumor sample.')

        # If there are no samples to process, skip the rest of the code
        if not threads:
            logger.info("No samples to process for this run. Skipping emails and further processing.")
            sys.exit(0)
    except Exception as e:
        print(f"Error during setup: {e}")
        error_setup_email(instrument)
        raise e

    # === Start analysis ===
    start_email(Rctx.run_name, final_pairs)

    # Start several samples at the same time
    for t in threads:
        t.start()
        logger.info(f'Thread {t} has started')

    for u in threads:
        u.join()
        logger.info(f'Thread {u} is finished')

    ok_samples = []
    bad_samples = []
    # Check if all samples in run have finished successfully. If not, exit script and send error email.
    for outputdir, sample_info in zip(outputdirs, final_pairs):
        if check_ok(outputdir):
            ok_samples.append(sample_info)
            logger.info(f'Finished correctly: {sample_info}')
        else:
            logger.info(f'Not finished correctly: {sample_info}')
            bad_samples.append(sample_info)

    if bad_samples:
        # send emails about which samples ok and which not ok
        error_email(Rctx.run_name, ok_samples, bad_samples)
    else:
        logger.info('All jobs have finished successfully')
        end_email(Rctx.run_name, final_pairs)

    # Run analysis_end for all samples in run, will check again which (if any) are ok
    # Will add ok samples to yearly stats and copy results
    for t in end_threads:
        t.start()
    for u in end_threads:
        u.join()
        
    # Combine all qc stats for the samples in the same run
    # the defaults for base_directory and output_directory are defined in the launcher config file
    # and don't need to be added as arguments
    try:
        logger.info(f'Combining qc stats for run {Rctx.run_name}')
        combine_qc_stats(launcher_config = LAUNCHER_CONFIG_PATH, outputdirs=outputdirs, runname=Rctx.run_name, logger=logger)
        logger.info(f'Done with combining qc stats for run {Rctx.run_name}')
    except Exception as e:
        logger.error(f"Error combining qc stats: {e}")
        error_admin_qc_email(Rctx.run_name)
 

def manual(tumorsample=None, normalsample=None, outpath=None, copyresults=False, qcsummary=False):
    '''Manual pipeline submission'''
    config = read_config(WRAPPER_CONFIG_PATH)
    wrapper_log_path = config["wrapper_log_path"]
    logger = setup_logger('wrapper', os.path.join(wrapper_log_path, 'Manual_WS_wrapper.log'))

    # If outputpath is not specified, get from config
    if not outpath:
        try:
            outpath = config['manual_outpath']
        except KeyError:
            logger.error('Output path for manual submission not specified in the configuration.')
            raise ValueError('Output path for manual submission not specified in the configuration.')

    threads = []
    outputdir = submit_pipeline(tumorsample, normalsample, outpath, config, logger, threads)
    threads[0].start()  # For manual runs we only have one thread

    threads[0].join()  # Wait for the thread to finish

    if copyresults and check_ok(outputdir):
        copy_results(outputdir)

    if qcsummary and check_ok(outputdir):
        try:
            logger.info(f'Combining qc stats for manual run with outputdir {outputdir}')
            combine_qc_stats(launcher_config = LAUNCHER_CONFIG_PATH, outputdirs=[outputdir], runname='manual_run', logger=logger)
            logger.info(f'Done with combining qc stats for manual run with outputdir {outputdir}')
        except Exception as e:
            logger.error(f"Error combining qc stats: {e}")
    return


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--instrument', help='For example novaseq_687_gc or novaseq_A01736', required=False)
    parser.add_argument('-t','--tumorsample', help='Specify the name of the tumor sample (e.g. DNA123456)', required=False)
    parser.add_argument('-n','--normalsample', help='Specify the name of the normal sample (e.g. DNA123456)', required=False)
    parser.add_argument('-o', '--outpath', help='Manually specify the path where the outputdir will go', required=False)
    parser.add_argument('-c', '--copyresults', help='Copy the results from a manual run to webstore', required=False, action='store_true', default=False)
    parser.add_argument('-q', '--qcsummary', help='Create combined qc summary for the run', required=False, action='store_true', default=False)
    args = parser.parse_args()

    if args.instrument:
        if args.tumorsample or args.normalsample or args.copyresults:
            parser.warning("When specifying --instrument, --tumorsample, --normalsample and --copyresults are ignored.")
        wrapper(args.instrument, args.outpath)
    elif args.tumorsample or args.normalsample:
        manual(args.tumorsample, args.normalsample, args.outpath, args.copyresults, args.qcsummary)
    else:
        parser.error("You must specify either --instrument or --tumorsample/--normalsample.")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
