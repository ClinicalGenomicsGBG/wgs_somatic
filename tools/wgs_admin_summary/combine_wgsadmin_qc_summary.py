import os
import pandas as pd
import click
import json
import logging
from collections import defaultdict
import glob
from launch_snakemake import get_timestamp


def setup_logging(log_dir):
    os.makedirs(log_dir, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(os.path.join(log_dir, "combine_wgsadmin_qc_summary.log")),
            logging.StreamHandler()
        ]
    )


def collect_wgsadmin_qc(outputdirs, logger):
    """
    For each outputdir, find the first *_qc_stats_wgsadmin.xlsx in qc_report or subdirs,
    read into a DataFrame, and return a dict: {outputdir: df or None}
    """
    qc_data = {}
    for outputdir in outputdirs:
        if not os.path.exists(outputdir):
            logger.warning(f"Output directory does not exist: {outputdir}. Skipping...")
            continue
        found_files = glob.glob(os.path.join(outputdir, "qc_report", "*_qc_stats_wgsadmin.xlsx"))
        if not found_files:
            logger.info(f"No qc_report directory or no matching files in {outputdir}. Searching all subdirectories...")
            found_files = glob.glob(os.path.join(outputdir, "*", "*_qc_stats_wgsadmin.xlsx"))
            if not found_files:
                logger.warning(f"No matching files found in any subdirectory of {outputdir}. Searching directly in {outputdir}...")
                found_files = glob.glob(os.path.join(outputdir, "*_qc_stats_wgsadmin.xlsx"))
                if not found_files:
                    logger.error(f"No qc_admin files found in {outputdir}. Skipping...")
                    qc_data[outputdir] = None
                    continue
        logger.info(f"Found qc_admin file(s) in {outputdir}: {found_files}")
        if len(found_files) > 1:
            logger.warning(f"Multiple qc_admin files found in {outputdir}. Only {found_files[0]} will be processed.")
        try:
            df = pd.read_excel(found_files[0])
            qc_data[outputdir] = df
        except Exception as e:
            logger.error(f"Error reading {found_files[0]}: {e}")
            qc_data[outputdir] = None
    return qc_data


def combine_qc_stats(launcher_config, outputdirs, runname=None, qc_summary_directory=None, logger=None):
    """
    Combine the first *_qc_stats_wgsadmin.xlsx file found in each output directory into a single summary Excel and TSV file.
    """
    config_data = {}
    if launcher_config:
        try:
            with open(launcher_config, 'r') as launcher_config_file:
                config_data = json.load(launcher_config_file)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            raise FileNotFoundError(f"Launcher config not found or invalid: {e}")

    if not logger:
        logdir = config_data.get('logdir', os.getcwd())
        setup_logging(logdir)
        logger = logging.getLogger(__name__)

    logger.info("Starting WGS admin QC summary per run")

    # Prepare output directories and files
    if qc_summary_directory is None:
        qc_summary_directory = config_data.get('wgsadmin_dir', os.getcwd())
    if runname:
        outputfile_prefix = os.path.join(qc_summary_directory, f"{runname}_wgs_somatic_qc_summary_{get_timestamp()}")
    else:
        outputfile_prefix = os.path.join(qc_summary_directory, f"wgs_somatic_qc_summary_{get_timestamp()}")
    output_file_xlsx = f"{outputfile_prefix}.xlsx"
    output_file_tsv = f"{outputfile_prefix}.tsv"
    if os.path.exists(output_file_xlsx) or os.path.exists(output_file_tsv):
        raise FileExistsError(f"Output file(s) already exist: {output_file_xlsx} or {output_file_tsv}")
    logger.info(f"QC summary directory: {qc_summary_directory}")
    os.makedirs(qc_summary_directory, exist_ok=True)
    logger.info(f"Output files will be: {output_file_xlsx} and {output_file_tsv}")

    # Collect QC data from each outputdir
    qc_data = collect_wgsadmin_qc(outputdirs, logger)

    # Combine all DataFrames into one and save
    if not qc_data:
        logger.warning("No QC data collected from any output directories.")
        return
    combined_df = pd.concat(qc_data.values(), ignore_index=True)
    combined_df.to_excel(output_file_xlsx, index=False)
    combined_df.to_csv(output_file_tsv, sep='\t', index=False)
    logger.info(f"Combined QC summary saved to {output_file_xlsx} and {output_file_tsv}")

@click.command()
@click.option('--launcher-config', required=False, type=click.Path(exists=True, dir_okay=False), help='Path to launcher config JSON file (optional).')
@click.option('-o', '--outputdirs', required=True, multiple=True, type=click.Path(exists=True, file_okay=False), help='List of output directories to search for QC files. Can be specified multiple times.')
@click.option('--runname', required=False, help='Optional run name for output file naming.')
@click.option('--qc-summary-directory', required=False, type=click.Path(file_okay=False), help='Directory to save combined summary files. If not set, uses wgsadmin_dir from launcher config or current directory.')
def cli(launcher_config, outputdirs, runname, qc_summary_directory):
    """
    Combine the first *_qc_stats_wgsadmin.xlsx file found in each output directory into a single summary Excel and TSV file.

    For each output directory, searches for qc_report/*_qc_stats_wgsadmin.xlsx, or elsewhere if not found.
    Reads the first matching file into a DataFrame and combines all found DataFrames into a single summary file.
    Output files are written to the specified QC summary directory, the launcher config's wgsadmin_dir, or the current directory.
    """
    combine_qc_stats(
        launcher_config=launcher_config,
        outputdirs=outputdirs,
        runname=runname,
        qc_summary_directory=qc_summary_directory
    )

if __name__ == "__main__":
    cli()