import os
import pandas as pd
# import click
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
    
# @click.command()
# @click.option('--base_directory', 
#               default=None,
#               type=click.Path(file_okay=False, dir_okay=True, exists=True),
#               help="Base directory containing subdirectories for each sample/sample T-N pairs. Each subdirectory should contain '_qc_stats_wgsadmin.xlsx' files to be combined. If not provided, the directory is determined from the launcher config file.", 
#               show_default=True)
# @click.option('--output_directory', 
#               default=None,
#               type=click.Path(file_okay=False, dir_okay=True),
#               help="Directory to save the combined output files. If not provided, the directory is determined from the launcher config file.",
#               show_default=True)
# @click.option('--launcher_config', 
#               type=click.Path(exists=True), 
#               required=True, 
#               default=os.path.join(os.path.dirname(__file__), "../../configs/launcher_config.json"), 
#               help="Path to the launcher config file. This file should contain the 'wgsadmin_dir' key specifying the default output directory.", 
#               show_default=True)
# @click.option('--runtag_results',
#               type=str,
#               default=None,
#               help="Specific runtag to process. If provided, only this runtag will be processed. Otherwise, all runtags in the base directory will be processed.",
#               show_default=True)
def combine_qc_stats(base_directory, output_directory, launcher_config, runtag_results):
    """
    Command-line tool to combine '_qc_stats_wgsadmin.xlsx' files for each runtag in the given BASE_DIRECTORY.
    """
    # set up logger. log dir is taken from launcher_config
    with open(launcher_config, 'r') as launcher_config_file:
        setup_logging(json.load(launcher_config_file).get('logdir'))
    
    logger = logging.getLogger(__name__)
    logger.info("Starting WGS admin QC summary per run")
    
    # if output_directory is not provided, read it from the launcher_config
    # if base_directory is not provided, read it from the launcher_config
    try:
        with open(launcher_config, 'r') as launcher_config_file:
            config_data = json.load(launcher_config_file)
            if output_directory is None:
                output_directory = config_data.get('wgsadmin_dir')
                if output_directory is None:
                    raise KeyError("Key 'wgsadmin_dir' not found in launcher_config.")
            if base_directory is None:
                base_directory = config_data.get('resultdir_hg38')
                if base_directory is None:
                    raise KeyError("Key 'resultdir_hg38' not found in launcher_config.")
    except (FileNotFoundError, json.JSONDecodeError, KeyError) as e:
        raise FileNotFoundError(f"Launcher config not found or invalid: {e}")
        
    logger.info(f"Base directory: {base_directory}")
    logger.info(f"Output directory: {output_directory}")
    
    # Dictionary to store dataframes for each runtag
    runtag_dataframes = {}
    if runtag_results:
        runtag_parts = runtag_results.split("+")[0].split("_")
        input_runtag = f"{runtag_parts[0]}_{runtag_parts[3]}"
    else:
        input_runtag = None

    # Iterate through directories in the base directory. Also checks subdirectories such as tumor_only and normal_only
    for root, dirs, files in os.walk(base_directory):
        for directory in dirs:
            dir_path = os.path.join(root, directory)  # Full path of the directory
            
            # Check if it's a directory and matches the naming pattern
            if len(directory.split('_')) >= 4:
                # Extract the runtag (second and third parts of the directory name)
                parts = directory.split('_')
                current_runtag = f"{parts[1]}_{parts[2]}"
                
                # If a specific runtag is provided, skip directories that don't match
                if input_runtag and current_runtag != input_runtag:
                    logger.info(f"Skipping directory: {dir_path} as it does not match the specified runtag: {input_runtag}")
                    continue
                
                logger.info(f"Processing directory: {dir_path} for runtag: {current_runtag}")

                # Find all matching files in the directory
                matching_files = [file for file in os.listdir(dir_path) if file.endswith("_qc_stats_wgsadmin.xlsx")]
                
                # Skip the directory if no matching files are found
                if not matching_files:
                    logger.warning(f"No '_qc_stats_wgsadmin.xlsx' files found in {dir_path}. Skipping...")
                    continue
                
                # Initialize a list to store dataframes for this runtag
                if current_runtag not in runtag_dataframes:
                    runtag_dataframes[current_runtag] = []
                
                # Process each matching file
                for file in matching_files:
                    file_path = os.path.join(dir_path, file)
                    try:
                        df = pd.read_excel(file_path)
                        runtag_dataframes[current_runtag].append((file_path, df))
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

if __name__ == "__main__":
    combine_qc_stats()