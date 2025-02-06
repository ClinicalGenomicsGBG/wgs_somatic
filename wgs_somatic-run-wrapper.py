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
from itertools import chain
import traceback
import subprocess
import threading

from definitions import WRAPPER_CONFIG_PATH, ROOT_DIR #, INSILICO_CONFIG, INSILICO_PANELS_ROOT
from tools.context import RunContext, SampleContext
from tools.helpers import setup_logger, read_config
from tools.slims import get_sample_slims_info, SlimsSample, find_or_download_fastqs, get_pair_dict, link_fastqs_to_workingdir
from tools.email import start_email, end_email, error_email
from launch_snakemake import analysis_main, yearly_stats, copy_results, get_timestamp


# Store info about samples to use for sending report emails
sample_status = {'missing_slims': [],
                 'unset_WS': [],
                 'approved': []}


def look_for_runs(config, instrument):
    '''Look for runs in demultiplexdir'''
    instrument_root_path = config[instrument]['demultiplex_path']
    found_paths = glob.glob(os.path.join(instrument_root_path, '*'))
    regex = config[instrument]['seq_name_regex']
    return [path for path in found_paths if re.search(regex, os.path.basename(path))]


def generate_context_objects(Rctx, logger):
    '''Create Rctx and Sctx for a demultiplexed run'''

    # Read demultiplex stats file for sample names, fastq paths, and nr reads
    with open(Rctx.demultiplex_summary_path, 'r') as inp:
        demuxer_info = json.load(inp)


    for sample_id, sample_info in demuxer_info['samples'].items():
        logger.info(f'Setting up context for {sample_id}.')

        # Setup Sample context class and add listed fastq paths
        Sctx = SampleContext(sample_id)
        Sctx.add_fastq(sample_info['fastq_paths'])

        # Query Slims for clinical information and add to sample context
        logger.info(f'Fetching SLIMS info.')
        Sctx.slims_info = get_sample_slims_info(Sctx, run_tag = Rctx.run_tag)

        if not Sctx.slims_info:
            logger.warning(f'No SLIMS info available!')
            logger.warning(f'Sample will not be analysed.')
            sample_status['missing_slims'].append(Sctx)
            continue

        # NOTE: 54 is Slims internal primary key for wgs_somatic
        if 54 not in Sctx.slims_info['secondary_analysis']:
            sample_status['unset_WS'].append(Sctx)
            continue

        # Add sample context to run context list
        logger.info(f'Adding sample context to list of contexts.')
        Rctx.add_sample_context(Sctx)

    if not Rctx.sample_contexts:
        # If no samples set for wgs_somatic
        # doesn't skip continuing with the rest of the code - need to fix this
        logger.info('No samples set for wgs_somatic. Skipping run.')
        # exit script if no samples set for wgs_somatic in run
        sys.exit()

    for Sctx in Rctx.sample_contexts:
        sample_status['approved'].append(Sctx)

    return Rctx


def call_script(**kwargs):
    '''Function to call main function from launch_snakemake.py'''
    args = argparse.Namespace(**kwargs)
    analysis_main(args, **kwargs)


def check_ok(workingdir):
    '''Function to check if analysis has finished correctly'''

    if os.path.isfile(f"{workingdir}/reporting/workflow_finished.txt"):
        return True
    else:
        return False


def analysis_end(workingdir, tumorsample=None, normalsample=None):
    '''Function to check if analysis has finished correctly and add to yearly stats and copy results'''

    if os.path.isfile(f"{workingdir}/reporting/workflow_finished.txt"):
        if tumorsample:
            if normalsample:
                # these functions are only executed if snakemake workflow has finished successfully
                yearly_stats(tumorsample, normalsample)
                copy_results(workingdir)
            else:
                yearly_stats(tumorsample, 'None')
                copy_results(workingdir)
        else:
            yearly_stats('None', normalsample)
            copy_results(workingdir)
    else:
        pass


def wrapper(instrument):
    '''Wrapper function'''
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
            error_list.append(f"workingdirectory: {hcptmp} does not exist and could not be created")


    # Grab all available local run paths
    local_run_paths = look_for_runs(config, instrument)
    # Read all previously analysed runs
    previous_runs_file = config['previous_runs_file_path']
    previous_runs_file_path = os.path.join(ROOT_DIR, previous_runs_file)

    if not os.path.exists(previous_runs_file_path):
        raise FileNotFoundError(f"Runlist file not found: {previous_runs_file_path}")

    with open(previous_runs_file_path, 'r') as prev:
        previous_runs = [line.rstrip() for line in prev]

    # Loop through each run path and setup Run context class
    for run_path in local_run_paths:
        Rctx = RunContext(run_path)
        # Check if demultiplexing is completed
        if not Rctx.demultiplex_complete:
            continue

        # Check if run has been previously analysed
        if Rctx.run_name in previous_runs:
            continue

        # Register start time
        start_time = datetime.now()
        logger.info(f'Started {Rctx.run_name} at {start_time}')

        # Write run name to previously analysed list to ensure no double-running
        with open(previous_runs_file_path, 'a') as prev:
            logger.info(f'Writing {Rctx.run_name} to previous runs list.')
            print(Rctx.run_name, file=prev)

        # get Rctx and Sctx for current run
        Rctx_run = generate_context_objects(Rctx, logger)

        # Get T/N pair info in a dict for samples and link additional fastqs from other runs
        for sctx in Rctx_run.sample_contexts:
            pair_dict = get_pair_dict(sctx, Rctx, logger)
            pair_dict_all_pairs.update(pair_dict)


        # Uses the dictionary of T/N samples to put the correct pairs together and finds the correct input arguments to the pipeline
        threads = []
        check_ok_outdirs = []
        end_threads = []
        final_pairs = []
        paired_samples = []

        def submit_pipeline(tumorsample, normalsample):
            hg38ref = config['hg38ref']['GMS-BT']
            timestamp = get_timestamp()
            workingdir = None
            if tumorsample and normalsample:
                fastq_dict_tumor = find_or_download_fastqs(tumorsample, logger)
                fastq_dict_normal = find_or_download_fastqs(normalsample, logger)
                tumorid = list(fastq_dict_tumor.keys())[0]  # E.g. DNA123456_250101_AHJLJHBGXF
                workingdir = os.path.join(config['workingdir'], f"{tumorid}_{timestamp}")
                os.makedirs(workingdir, exist_ok=False)  # Make sure a new workingdir is created, not overwriting old results
                tumor_fastq_dir = link_fastqs_to_workingdir(fastq_dict_tumor, workingdir, logger)
                normal_fastq_dir = link_fastqs_to_workingdir(fastq_dict_normal, workingdir, logger)
                pipeline_args = {'workingdir': f'{workingdir}', 
                                 'normalname': f'{normalsample}', 
                                 'normalfastqs': f'{normal_fastq_dir}', 
                                 'tumorname': f'{tumorsample}', 
                                 'tumorfastqs': f'{tumor_fastq_dir}', 
                                 'hg38ref': f'{hg38ref}'}

            elif tumorsample:
                fastq_dict_tumor = find_or_download_fastqs(tumorsample, logger)
                workingdir = os.path.join(config['workingdir'], "tumor_only")
                os.makedirs(workingdir, exist_ok=True)
                tumorid = list(fastq_dict_tumor.keys())[0]
                workingdir = os.path.join(workingdir, f"{tumorid}_{timestamp}")
                os.makedirs(workingdir, exist_ok=False)
                tumor_fastq_dir = link_fastqs_to_workingdir(fastq_dict_tumor, workingdir, logger)
                pipeline_args = {'workingdir': f'{workingdir}', 
                                 'tumorname': f'{tumorsample}', 
                                 'tumorfastqs': f'{tumor_fastq_dir}', 
                                 'hg38ref': f'{hg38ref}'}
                
            elif normalsample:
                fastq_dict_normal = find_or_download_fastqs(normalsample, logger)
                workingdir = os.path.join(config['workingdir'], "normal_only")
                os.makedirs(workingdir, exist_ok=True)
                normalid = list(fastq_dict_normal.keys())[0]
                workingdir = os.path.join(workingdir, f"{normalid}_{timestamp}")
                os.makedirs(workingdir, exist_ok=False)
                normal_fastq_dir = link_fastqs_to_workingdir(fastq_dict_normal, workingdir, logger)
                pipeline_args = {'workingdir': f'{workingdir}', 
                                 'normalname': f'{normalsample}', 
                                 'normalfastqs': f'{normal_fastq_dir}', 
                                 'hg38ref': f'{hg38ref}'}
            
            threads.append(threading.Thread(target=call_script, kwargs=pipeline_args))
            logger.info(f'Starting wgs_somatic with arguments {pipeline_args}')
            check_ok_outdirs.append(workingdir)
            end_threads.append(threading.Thread(target=analysis_end, args=(workingdir, tumorsample, normalsample)))

        tumor_samples = {}
        normal_samples = {}

        # Separate tumor and normal samples
        for key, value in pair_dict_all_pairs.items():
            if value[0] == 'tumor':
                tumor_samples[key] = value
            elif value[0] == 'normal':
                normal_samples[key] = value

        # Pair samples based on tumorNormalID
        for t_key, t_value in tumor_samples.items():
            t_ID = t_value[1]  # tumorNormalID
            paired = False
            for n_key, n_value in normal_samples.items():
                n_ID = n_value[1]  # tumorNormalID
                if t_ID == n_ID or t_ID == n_key.split("DNA")[1] or n_ID == t_key.split("DNA")[1]:
                    paired_samples.append((t_key, n_key))
                    paired = True
                    submit_pipeline(t_key, n_key)
                    final_pairs.append(f'{t_key} (T) {n_key} (N), {n_value[2]} {["prio" if (n_value[3] or t_value[3]) else ""][0]}')
                    break
            if not paired:
                submit_pipeline(t_key, None)
                final_pairs.append(f'{t_key} (T), {t_value[2]} {["prio" if t_value[3] else ""][0]}')

        for n_key in normal_samples:
            if not any(n_key == pair[1] for pair in paired_samples):
                submit_pipeline(None, n_key)
                final_pairs.append(f'{n_key} (N), {n_value[2]} {["prio" if n_value[3] else ""][0]}')

        
        # Start several samples at the same time
        for t in threads:
            t.start()
            logger.info(f'Thread {t} has started')

        # Send start email
        start_email(Rctx_run.run_name, final_pairs)

        for u in threads:
            u.join()
            logger.info(f'Thread {u} is finished')

        ok_samples = []
        bad_samples = []
        # Check if all samples in run have finished successfully. If not, exit script and send error email.
        for outdir, sample_info in zip(check_ok_outdirs, final_pairs):
            if check_ok(outdir) == True:
                ok_samples.append(sample_info)
                logger.info(f'Finished correctly: {sample_info}')
            else:
                logger.info(f'Not finished correctly: {sample_info}')
                bad_samples.append(sample_info)
        if bad_samples:
            # send emails about which samples ok and which not ok
            error_email(Rctx_run.run_name, ok_samples, bad_samples)
            if ok_samples:
                # yearly stats ok samples
                # even though thread starts for all samples, function checks if sample ok
                # so it will only do yearly stats for ok samples
                for t in end_threads:
                    t.start()
                for u in end_threads:
                    u.join()
            sys.exit()

        logger.info('All jobs have finished successfully')
        end_email(Rctx_run.run_name, final_pairs)

        # If all jobs have finished successfully - add to yearly stats
        for t in end_threads:
            t.start()
        for u in end_threads:
            u.join()

        # break out of for loop to avoid starting pipeline for a possible other run that was done sequenced at the same time.
        # if cron runs every 30 mins it will find other runs at the next cron instance and run from there instead (and add to novaseq_runlist)
        break

    # DON'T FORGET TO UNCOMMENT PETAGENE COMPRESSION AND CHANGE BACK TO CORRECT workingdir AND IGVUSER!!!!!!

    # some arguments are hardcoded right now, need to fix this. 
    # only considers barncancer hg38 (GMS-AL + GMS-BT samples) right now. 

    # the arguments runtumor and runnormal could be "wrong" by doing it like this 
    # since they use run name of current run 
    # but if for example normal is in older run it has the wrong argument for "runnormal". 
    # this argument is not that important, it is only used to create a unique sample name. 
    # maybe it would be better discard/modify this argument than spending time on getting the correct value of it for samples from different runs. 
    # also, if fastqs come from more than one run, what will the value of this argument be then to be "correct"?...


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--instrument', help='For example novaseq_687_gc or novaseq_A01736', required=True)
    args = parser.parse_args()

    wrapper(args.instrument)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
    #except Exception:
        #format_exc = traceback.format_exc()
        #logger.error(format_exc)
        # TODO: add send email about error here
