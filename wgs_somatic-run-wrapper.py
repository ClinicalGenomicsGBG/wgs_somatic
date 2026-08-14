"""
Wrapper to be used by cron for automatic start of wgs_somatic
"""

import argparse
import os
import sys
import threading
import time
from pathlib import Path

from tools.slims import return_pending_patients, add_merge_samples, add_matched_samples
from tools.pre_pipeline import pre_pipeline, Run, Batch
from definitions import WRAPPER_CONFIG_PATH, ROOT_DIR, LAUNCHER_CONFIG_PATH 
from tools.custom_email import start_email, end_email, error_email, error_admin_qc_email, error_setup_email, manual_start_email, manual_end_email
from tools.helpers import setup_logger, read_config
from launch_snakemake import analysis_main, yearly_stats, copy_results
from tools.wgs_admin_summary.combine_wgsadmin_qc_summary import combine_qc_stats


def get_timestamp():
    return time.strftime("%y%m%d-%H%M%S")


def call_script(**kwargs):
    """Function to call main function from launch_snakemake.py."""
    args = argparse.Namespace(**kwargs)
    analysis_main(args, **kwargs)


def analysis_end(run):
    """Check whether analysis finished successfully and process the results."""
    outputdir = run.run_work_dir
    tumorsample = run.tumor_sample
    normalsample = run.normal_sample

    if run.run_ok:
        if tumorsample:
            if normalsample:
                # these functions are only executed if snakemake workflow has finished successfully
                yearly_stats(tumorsample, normalsample)
            else:
                yearly_stats(tumorsample, 'None')
        else:
            yearly_stats('None', normalsample)
        copy_results(outputdir)


def submit_pipeline(batch, run, config, logger):
    timestamp = get_timestamp()

    tumorsample = run.tumor_sample.id
    normalsample = run.normal_sample.id

    if tumorsample and normalsample:
        logger.info(f'Preparing run: Tumor {tumorsample} and Normal {normalsample}')
        #fastq_dict_tumor = find_or_download_fastqs(tumorsample, logger)
        #fastq_dict_normal = find_or_download_fastqs(normalsample, logger)
        #tumorid = list(fastq_dict_tumor.keys())[0]  # E.g. DNA123456_250101_AHJLJHBGXF
        #outputdir = os.path.join(run.run_work_dir, f"{tumorsample}_{timestamp}")
        #os.makedirs(outputdir, exist_ok=False)  # Make sure a new outputdir is created, not overwriting old results
        tumor_fastq_dir = run.prepared_fastq_dir
        normal_fastq_dir = run.prepared_fastq_dir
        pipeline_args = {'outputdir': f'{run.run_work_dir}',
                         'normalname': f'{normalsample}',
                         'normalfastqs': f'{normal_fastq_dir}',
                         'tumorname': f'{tumorsample}',
                         'tumorfastqs': f'{tumor_fastq_dir}',
                         'wrapper_logger': logger,}

    elif tumorsample:
        logger.info(f'Preparing run: Tumor-only {tumorsample}')
        #fastq_dict_tumor = find_or_download_fastqs(tumorsample, logger)
        #outpath = os.path.join(outpath, "tumor_only")
        #os.makedirs(outpath, exist_ok=True)
        #tumorid = list(fastq_dict_tumor.keys())[0]
        outputdir = os.path.join(run.run_work_dir, f"{run.tumor_sample}_{timestamp}")
        os.makedirs(outputdir, exist_ok=False)
        tumor_fastq_dir = run.prepared_fastq_dir
        pipeline_args = {'outputdir': f'{run.run_work_dir}',
                            'tumorname': f'{tumorsample}',
                            'tumorfastqs': f'{tumor_fastq_dir}',
                            'wrapper_logger': logger,}

    elif normalsample:
        logger.info(f'Preparing run: Normal-only {normalsample}')
        #fastq_dict_normal = find_or_download_fastqs(normalsample, logger)
        outpath = os.path.join(outpath, "normal_only")
        os.makedirs(outpath, exist_ok=True)
        #normalid = list(fastq_dict_normal.keys())[0]
        #outputdir = os.path.join(outpath, f"{run.normal_sample}_{timestamp}")
        #os.makedirs(outputdir, exist_ok=False)
        normal_fastq_dir = run.prepared_fastq_dir
        pipeline_args = {'outputdir': f'{run.run_work_dir}',
                         'normalname': f'{normalsample}',
                         'normalfastqs': f'{normal_fastq_dir}',
                         'wrapper_logger': logger,}


    batch.threads.append(threading.Thread(target=call_script, kwargs=pipeline_args))
    logger.info(f'Starting wgs_somatic with arguments {pipeline_args}')
    batch.outputdirs.append(run.run_work_dir)


def wrapper(outpath=None):
    '''Automatic wrapper function'''
    # === Setup run ===
    try:
        devmode = os.environ.get("DEVMODE") == "true"
        config = read_config(WRAPPER_CONFIG_PATH)

        if devmode:
            wrapper_log_path = config["develop_mode"]["log_path"]
            outpath = config["develop_mode"]["outpath"]
        else:
            wrapper_log_path = config["wrapper_log_path"]
        logger = setup_logger('Wrapper', os.path.join(wrapper_log_path, f'wgs-somatic-run-wrapper.log'))
        
        # Find pending samples
        if not devmode:
            patients = return_pending_patients(logger)
        else:
            patients = return_pending_patients(logger) # TODO Add test samples

        if not patients:
            logger.info('No pending samples found for wgs_somatic. Exiting wrapper.')
            sys.exit(0)

        # prepare hcp download directory
        hcptmp = config["hcp_download_dir"]
        if not os.path.isdir(hcptmp):
            try:
                os.makedirs(hcptmp)
            except Exception as e:
                logger.error(f'outputdirectory: {hcptmp} does not exist and could not be created')
                raise e

        # If outputpath is not specified, get from config
        if not outpath:
            try:
                outpath = config['cron_outpath']
            except KeyError:
                logger.error('Output path for cron job not specified in the configuration.')
                raise ValueError('Output path for cron job not specified in the configuration.')

        add_matched_samples(patients, logger)
        add_merge_samples(patients, logger)

        # Make a batch with all runs
        batch = Batch()

        # Add all patients to this batch of runs
        for patient in patients:
            batch.runs.append(Run(logger=logger, patient=patient, run_root_dir=Path(outpath)))

        # Prepare all samples for analysis
        pre_pipeline(batch, config, logger)

    except Exception as e:
        error_setup_email(f'Error during setup: {e}')
        raise

    for run in batch.runs:
        if not run.ready_for_pipeline:
            logger.info(f'Run {run.run_work_dir} is not ready for pipeline, skipping.')
            continue

        run.pipeline_args["wrapper_logger"] = logger
        
        submit_pipeline(batch, run, config, logger)
        #logger.info(f"Starting wgs_somatic with arguments {run.pipeline_args}")
        #batch.threads.append(threading.Thread(target=call_script,kwargs=run.pipeline_args,))

        batch.end_threads.append(threading.Thread(target=analysis_end, args=(run,)))

        if run.patient.has_normal:
            batch.run_summary.append(f"{run.tumor_sample.id} (T) {run.normal_sample.id} (N)")
        else:
            batch.run_summary.append(f"{run.tumor_sample.id} (T)")

    if not batch.threads:
        error_message = f"Failed to start analysis for all samples in {batch.batch_name}."
        error_setup_email(error_message)
        sys.exit(0)

    for t in batch.threads:
        t.start()
        logger.info(f"Thread {t} has started")

    for u in batch.threads:
        u.join()
        logger.info(f"Thread {u} has finished")

    # === Start analysis ===
    start_email(batch.batch_name, batch.run_summary)

    # Check if all samples in run have finished successfully. If not, exit script and send error email.
    if batch.runs_ok:
        logger.info('All jobs have finished successfully')
        end_email(batch.batch_name, batch.run_summary) 
    else:
        error_email(batch.batch_name, batch.ok_samples, batch.bad_samples)

    # Run analysis_end for all samples in run, will check again which (if any) are ok
    # Will add ok samples to yearly stats and copy results
    if not devmode:
        for t in batch.end_threads:
            t.start()
        for u in batch.end_threads:
            u.join()
        
    # Combine all qc stats for the samples in the same run
    # the defaults for base_directory and output_directory are defined in the launcher config file
    # and don't need to be added as arguments
    try:
        logger.info(f'Combining qc stats for run {batch.batch_name}')
        combine_qc_stats(launcher_config=LAUNCHER_CONFIG_PATH, outputdirs=batch.outputdirs, runname=batch.batch_name, logger=logger)
        logger.info(f'Done with combining qc stats for run {batch.batch_name}')
    except Exception as e:
        logger.error(f"Error combining qc stats: {e}")
        error_admin_qc_email(batch.batch_name)
 

def manual(tumorsample=None, normalsample=None, outpath=None, copyresults=False, qcsummary=False):
    '''Manual pipeline submission'''
    devmode = os.environ.get("DEVMODE")
    config = read_config(WRAPPER_CONFIG_PATH)

    if devmode:
        wrapper_log_path = config["develop_mode"]["log_path"]
        outpath = config["develop_mode"]["outpath"]
    else:
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
    
    if email:
        manual_start_email(tumorsample, normalsample)

    outputdir = submit_pipeline(tumor_name, normal_name, outpath, config, logger, threads)

    threads[0].start()  # For manual runs we only have one thread
    threads[0].join()  # Wait for the thread to finish

    if check_ok(outputdir):
        sucess_run = True

        if copyresults:
            copy_results(outputdir)

        if qcsummary:
            try:
                logger.info(f'Combining qc stats for manual run with outputdir {outputdir}')
                combine_qc_stats(launcher_config = LAUNCHER_CONFIG_PATH, outputdirs=[outputdir], runname='manual_run', logger=logger)
                logger.info(f'Done with combining qc stats for manual run with outputdir {outputdir}')
            except Exception as e:
                logger.error(f"Error combining qc stats: {e}")

    if email:
        manual_end_email(sucess_run, tumorsample, normalsample)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--integrated', help='Integrated run querying for novel samples in slims', required=False, action='store_true', default=False)
    parser.add_argument('-t', '--tumorsample', help='Specify the name of the tumor sample (e.g. DNA123456)', required=False)
    parser.add_argument('-n', '--normalsample', help='Specify the name of the normal sample (e.g. DNA123456)', required=False)
    parser.add_argument('-o', '--outpath', help='Manually specify the path where the outputdir will go', required=False)
    parser.add_argument('-c', '--copyresults', help='Copy the results from a manual run to webstore', required=False, action='store_true', default=False)
    parser.add_argument('-q', '--qcsummary', help='Create combined qc summary for the run', required=False, action='store_true', default=False)
    parser.add_argument('-d', '--develop-mode', help='Run wrapper in "develop mode"', required=False, action='store_true', default=False)
    parser.add_argument('-e', '--email', help='Send start, end and crash emails from a manual run and in development mode', required=False, action='store_true', default=False)

    args = parser.parse_args()

    os.environ["DEVMODE"] = str(args.develop_mode).lower()
    os.environ["EMAIL"] = str(args.email).lower()
    
    if args.integrated:
        if args.tumorsample or args.normalsample or args.copyresults:
            parser.warning("When specifying --integrated, --tumorsample, --normalsample and --copyresults are ignored.")
        wrapper(args.outpath)
    elif args.tumorsample or args.normalsample:
        manual(args.tumorsample, args.normalsample, args.outpath, args.copyresults, args.qcsummary, args.email)
    else:
        parser.error("You must specify either --integrated or --tumorsample/--normalsample.")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
