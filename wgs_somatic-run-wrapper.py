"""
Wrapper to be used by cron for automatic start of wgs_somatic
"""

import sys
import argparse
import os
import re
import glob
import json
import threading
from collections import defaultdict
from pathlib import Path

from definitions import WRAPPER_CONFIG_PATH, ROOT_DIR, LAUNCHER_CONFIG_PATH #, INSILICO_CONFIG, INSILICO_PANELS_ROOT
from tools.custom_email import( start_email, end_email, error_email, error_admin_qc_email, error_setup_email, manual_start_email, manual_end_email
from tools.helpers import setup_logger, read_config, make_long_term_storage_info
from tools.slims import return_pending_patients, add_merge_samples, add_matched_samples, Sample
from tools.pre_pipeline import pre_pipeline, Batch, Run
from launch_snakemake import analysis_main, yearly_stats, copy_results, get_timestamp
from tools.wgs_admin_summary.combine_wgsadmin_qc_summary import combine_qc_stats


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


def start_analysis(batch, run, logger):
    '''Function to start the pipeline'''
    run.pipeline_args["wrapper_logger"] = logger # Pass the wrapper logger to analysis_main(), so the process can be tracked in the termeinal, mainly during development

    if run.patient.has_normal:
        run_summary.append(f"{run.tumor_sample.id} (T) {run.normal_sample.id} (N)")
    else:
        run_summary.append(f"{run.tumor_sample.id} (T)")

    logger.info(f"Preparing wgs_somatic with arguments {run.pipeline_args}")
    batch.threads.append(threading.Thread(target=call_script,kwargs=run.pipeline_args,))
    batch.end_threads.append(threading.Thread(target=analysis_end, args=(batch, run, logger)))


def analysis_end(batch, run, logger):
    '''Function to check if analysis has finished correctly and add to yearly stats and copy results'''

    if run.check_ok:
        if run.has_tumor:
            if run.has_normal:
                # these functions are only executed if snakemake workflow has finished successfully
                yearly_stats(run.tumor_sample, normal_sample)
            else:
                yearly_stats(tumor_sample, 'None')
        elif run.has_normal:
            yearly_stats('None', normal_sample)
        else:
            logger.warning("No samples in run {run.main_id}")
        logger.info(f"Copying results from {run.run_work_dir}")
#TODO CONTINUE HERE
        copy_results(run_work_dir)
    else:
        pass


def wrapper(outpath=None):
    '''Automatic wrapper function'''
    # Fetch the enivoromental variables.
    devmode = os.environ.get("DEVMODE")
    silence = os.environ.get("SILENCE")

    # === Setup runs ===
    try:
        config = read_config(WRAPPER_CONFIG_PATH)

        if devmode:
            wrapper_log_path = config["develop_mode"]["log_path"]
            outpath = config["develop_mode"]["outpath"]
        else:
            wrapper_log_path = config["wrapper_log_path"]
        logger = setup_logger('wrapper', os.path.join(wrapper_log_path, f'wgs-somatic-run-wrapper.log'))

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

        # Find pending samples
        if not devmode:
            patients = return_pending_patients(config, logger)
        else:
            run_type = 'success' #For future testing of predefined run types, e.g., successfull, failed due to, quality, slims error, et.c.,
            patients = return_pending_patients(config, logger) #TODO Add test samples
            #patients = generate_develop_patients(config, run_type, logger)
        
        if not patients:
            logger.info("No pending samples found for wgs_somatic. Exiting wrapper.")
            sys.exit(0)
        
        add_matched_samples(patients, logger)        
        add_merge_samples(patients, logger)

        # Make a batch with all runs
        batch = Batch()

        # Iterate over all patients found with samples to be analysed and add them to this batch of runs
        for patient in patients:
            batch.runs.append(Run(logger=logger, patient=patient, run_root_dir=Path(outpath)))

            for sample in patient.samples:
                if not sample.long_term_storage_info:
                    logger.debug("This is a temporary solution due to missing slims info")
                    bucket = config["hcp"]["download_locations"]["sg1_illumina"]["bucket"]
                    endpoint = "https://sg1.vgregion.se"
                    credentials = config["hcp"]["download_locations"]["sg1_illumina"]["credentials_file"]
                    sample.long_term_storage_info = make_long_term_storage_info(sample.id, bucket, endpoint, credentials)
                    sample.remote_keys, sample.bucket = sample._parse_long_term_storage_info()

        # Prepare all the samples for analysis
        pre_pipeline(runs, config, logger)

    batch.determine_batch_id

    except Exception as e:
        error_message = f"Error during setup: {e}"
        logger.warning(error_message)
        error_setup_email(error_message)
        raise e
 

    # === Start runs ===
    for run in batch.runs:
        if not run.ready_for_pipeline:
            logger.warning(f"Run {run.patient} is not ready for pipeline, skipping.")
            continue
        
        start_analysis(batch, run, logger)

    if not threads:
        error_message = f"Failed to start analysis for all samples in {batch.batch_name}."
        logger.warning(error_message)
        error_setup_email(error_message)
        sys.exit(0)

    if not silence: start_email(batch.batch_name, run_summary)
    #
    # === Start analysis ===
    start_email(Rctx.run_name, final_pairs)
    
    # Start several samples at the same time
    for thread in threads:
        thread.start()
        logger.info(f"Thread {thread} has started")

    for thread in threads:
        thread.join()
        logger.info(f"Thread {thread} has finished")

    # Check if all samples in run have finished successfully. If not, exit script and send error email.
    for outputdir, sample_info in zip(outputdirs, final_pairs):
        if check_ok(outputdir):
            batch.ok_samples.append(sample_info)
            logger.info(f'Finished correctly: {sample_info}')
        else:
            logger.info(f'Not finished correctly: {sample_info}')
            batch.bad_samples.append(sample_info)

    if bad_samples:
        # send emails about which samples ok and which not ok
        error_email(batch.batch_name, ok_samples, bad_samples)
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
 

def manual(tumorsample=None, normalsample=None, outpath=None, copyresults=False, qcsummary=False, email=False):
    '''Manual pipeline submission'''
    devmode = os.environ.get("DEVMODE") == "true"
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

    outputdir = submit_pipeline(tumorsample, normalsample, outpath, config, logger, threads)

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
    parser.add_argument('-e', '--email', help='Send start, end and crash emails from a manual run', required=False, action='store_true', default=False)
    parser.add_argument('-s', '--silence-email', help='Do not send any emails', required=False, action='store_true', default=False)
    parser.add_argument(-'r', '--snakemake-launcher', help='Restart snakemake in specified directory', required=False)

    args = parser.parse_args()

    os.environ["DEVMODE"] = str(args.develop_mode).lower()
    os.environ["SILENCE"] = str(args.silence_email).lower()

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
