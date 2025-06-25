import os
import pandas as pd
import click
import json
import logging
from collections import defaultdict
import re

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

def combine_qc_stats(launcher_config, runtag_results, base_directory=None, output_directory=None, logger=None,regex=None):
    """
    Command-line tool to combine '_qc_stats_wgsadmin.xlsx' files for each runtag in the given BASE_DIRECTORY.
    """
    # Load launcher_config once and parse it
    try:
        with open(launcher_config, 'r') as launcher_config_file:
            config_data = json.load(launcher_config_file)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        raise FileNotFoundError(f"Launcher config not found or invalid: {e}")
    
    if not logger:
        setup_logging(config_data.get('logdir'))
        logger = logging.getLogger(__name__)
    
    logger.info("Starting WGS admin QC summary per run")
    
    # If output_directory is not provided, read it from the launcher_config
    # If base_directory is not provided, read it from the launcher_config
    try:
        if output_directory is None:
            output_directory = config_data.get('wgsadmin_dir')
            if output_directory is None:
                raise KeyError("Key 'wgsadmin_dir' not found in launcher_config.")
        if base_directory is None:
            base_directory = config_data.get('resultdir_hg38')
            if base_directory is None:
                raise KeyError("Key 'resultdir_hg38' not found in launcher_config.")
    except KeyError as e:
        raise FileNotFoundError(f"Launcher config missing required key: {e}")
        
    logger.info(f"Base directory: {base_directory}")
    logger.info(f"Output directory: {output_directory}")
    
    # Dictionary to store dataframes for each runtag
    runtag_dataframes = defaultdict(list)

    if runtag_results:
        runtag_parts = runtag_results.split("+")[0].split("_")
        input_runtag = f"{runtag_parts[0]}_{runtag_parts[3]}"
    else:
        input_runtag = None

    # Walk through all subdirectories in the base directory
    for root, dirs, files in os.walk(base_directory):
        # Extract the runtag (second and third parts of the current directory name)
        parent_dir = os.path.basename(root)
        if regex and re.match(regex, parent_dir):
            parts = parent_dir.split('_')
            current_runtag = f"{parts[1]}_{parts[2]}"
            # If a specific runtag is provided, skip directories that don't match exactly
            if input_runtag and current_runtag != input_runtag:
                # logger.info(f"Skipping directory: {root} as it does not match the specified runtag: {input_runtag}")
                continue

            logger.info(f"Processing directory: {root} for runtag: {current_runtag}")

            # Find all matching files in the current directory
            matching_files = [file for file in files if file.endswith("_qc_stats_wgsadmin.xlsx")]

            # Skip the directory if no matching files are found
            if not matching_files:
                logger.warning(f"No '_qc_stats_wgsadmin.xlsx' files found in {root}. Skipping...")
                continue

            # Process each matching file
            for file in matching_files:
                file_path = os.path.join(root, file)
                try:
                    df = pd.read_excel(file_path)
                    runtag_dataframes[current_runtag].append((file_path, df))
                    logger.info(f"Successfully read file: {file_path}")
                except Exception as e:
                    logger.error(f"Error reading {file_path}: {e}")

    if not runtag_dataframes:
        raise RuntimeError(f"No QC summary files found for runtag {input_runtag}")
    
    os.makedirs(output_directory, exist_ok=True)

    # Combine and save the dataframes for each runtag
    for runtag, file_dataframes in runtag_dataframes.items():
        logger.info(f"Processing runtag: {runtag} with {len(file_dataframes)} files")
        
        output_file_xlsx = os.path.join(output_directory, f"{runtag}_wgs_somatic_qc_summary.xlsx")
        output_file_tsv = os.path.join(output_directory, f"{runtag}_wgs_somatic_qc_summary.tsv")
        
        if os.path.exists(output_file_xlsx):
            logger.info(f"Existing combined file found for runtag {runtag}: {output_file_xlsx}")
            # Read the existing combined file
            existing_df = pd.read_excel(output_file_xlsx)
            combined_df = existing_df
            # newer_file_detected = False
            new_content_detected = False
            for file_path, new_df in file_dataframes:
                # check if the content of the files in the directory is already in the combined file
                if not new_df.isin(existing_df.to_dict(orient='list')).all().all():
                    logger.info(f"New content detected in file: {file_path}. Combining with existing summary.")
                    combined_df = pd.concat([combined_df, new_df], ignore_index=True).drop_duplicates()
                    new_content_detected = True
                else:
                    logger.info(f"Skipping {file_path} as its content is already present in the existing combined file.")
            
            if not new_content_detected:
                logger.info(f"No new content detected for runtag {runtag}. Skipping...")
                continue
        else:
            logger.info(f"No existing combined file found for runtag {runtag}. Creating a new one.")
            # Combine all new dataframes if no existing file
            combined_df = pd.concat([df for _, df in file_dataframes], ignore_index=True).drop_duplicates()
        
        # Save the combined file
        combined_df.to_excel(output_file_xlsx, index=False)
        combined_df.to_csv(output_file_tsv, sep='\t', index=False)
        logger.info(f"Combined files saved for runtag {runtag}: {output_file_xlsx} and {output_file_tsv}")

    logger.info("Finished WGS admin QC summary summary per run.")

@click.command()
@click.option('--launcher_config', required=True, help='Path to launcher config JSON file.')
@click.option('--runtag_results', required=False, help='Runtag of the directories to process.')
@click.option('--base_directory', required=False, help='Base directory to search for QC files.')
@click.option('--output_directory', required=False, help='Directory to save combined summary files.')
@click.option('--regex', required=False, help='Regex to match directory names.')
def cli(launcher_config, runtag_results, base_directory, output_directory, regex):
    combine_qc_stats(
        launcher_config=launcher_config,
        runtag_results=runtag_results,
        base_directory=base_directory,
        output_directory=output_directory,
        regex=regex
    )

if __name__ == "__main__":
    cli()