import os
import pandas as pd
import click
import json
import logging

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

@click.command()
@click.option('--base_directory', 
              default="/clinical/data/wgs_somatic/test_output/test_qc_admin/test_combined", 
              type=click.Path(file_okay=False, dir_okay=True, exists=True),
              help="Base directory containing subdirectories for each sample/sample T-N pairs. Each subdirectory should contain '_qc_stats_wgsadmin.xlsx' files to be combined.", 
              required=True, show_default=True)
@click.option('--output_directory', 
              default=None, 
              type=click.Path(file_okay=False, dir_okay=True),
              help="Directory to save the combined output files. If not provided, the directory is determined from the launcher config file.",
              show_default=True)
@click.option('--launcher_config', 
              type=click.Path(exists=True), 
              required=True, 
              default=os.path.join(os.path.dirname(__file__), "../../configs/launcher_config.json"), 
              help="Path to the launcher config file. This file should contain the 'wgsadmin_dir' key specifying the default output directory.", 
              show_default=True)
def combine_qc_stats(base_directory, output_directory, launcher_config):
    """
    Command-line tool to combine '_qc_stats_wgsadmin.xlsx' files for each runtag in the given BASE_DIRECTORY.
    """
    # set up logger. log dir is taken from launcher_config
    with open(launcher_config, 'r') as launcher_config_file:
        setup_logging(json.load(launcher_config_file).get('logdir'))
    
    logger = logging.getLogger(__name__)
    logger.info("Starting WGS admin QC summary summary per run")
    
    if output_directory is None:
        try:
            with open(launcher_config, 'r') as launcher_config_file:
                output_directory = json.load(launcher_config_file).get('wgsadmin_dir')
                if output_directory is None:
                    raise KeyError("Key 'wgsadmin_dir' not found in launcher_config.")
        except (FileNotFoundError, json.JSONDecodeError, KeyError) as e:
            raise FileNotFoundError(f"Launcher config not found: {e}")
        
    logger.info(f"Base directory: {base_directory}")
    logger.info(f"Output directory: {output_directory}")
    
    # Dictionary to store dataframes for each runtag
    runtag_dataframes = {}

    # Iterate through directories in the base directory
    for directory in os.listdir(base_directory):
        dir_path = os.path.join(base_directory, directory)
        
        # Check if it's a directory and matches the naming pattern
        if os.path.isdir(dir_path) and len(directory.split('_')) >= 4:
            # Extract the runtag (second and third parts of the directory name)
            parts = directory.split('_')
            runtag = f"{parts[1]}_{parts[2]}"
            
            logger.info(f"Processing directory: {dir_path} for runtag: {runtag}")

            # Find all matching files in the directory
            matching_files = [file for file in os.listdir(dir_path) if file.endswith("_qc_stats_wgsadmin.xlsx")]
            
            # Skip the directory if no matching files are found
            if not matching_files:
                logger.warning(f"No '_qc_stats_wgsadmin.xlsx' files found in {dir_path}. Skipping...")
                continue
            
            # Initialize a list to store dataframes for this runtag
            if runtag not in runtag_dataframes:
                runtag_dataframes[runtag] = []
            
            # Process each matching file
            for file in matching_files:
                file_path = os.path.join(dir_path, file)
                try:
                    df = pd.read_excel(file_path)
                    runtag_dataframes[runtag].append((file_path, df))
                    logger.info(f"Successfully read file: {file_path}")
                except Exception as e:
                    logger.error(f"Error reading {file_path}: {e}")

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
            newer_file_detected = False
            for file_path, new_df in file_dataframes:
                # Compare modification times
                if os.path.getmtime(file_path) > os.path.getmtime(output_file_xlsx):
                    logger.info(f"Newer file detected: {file_path}. Combining with existing report.")
                    combined_df = pd.concat([combined_df, new_df], ignore_index=True).drop_duplicates()
                    newer_file_detected = True
                else:
                    logger.info(f"Skipping {file_path} as it is older than the existing combined file.")
            
            # If no newer file was detected, skip saving
            if not newer_file_detected:
                logger.info(f"No newer files detected for runtag {runtag}. Skipping...")
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

if __name__ == "__main__":
    combine_qc_stats()