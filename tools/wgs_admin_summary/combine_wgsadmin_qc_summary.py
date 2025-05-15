import os
import pandas as pd
import click
import json

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
    
    if output_directory is None:
        try:
            with open(launcher_config, 'r') as launcher_config_file:
                output_directory = json.load(launcher_config_file).get('wgsadmin_dir')
                if output_directory is None:
                    raise KeyError("Key 'wgsadmin_dir' not found in launcher_config.")
        except (FileNotFoundError, json.JSONDecodeError, KeyError) as e:
            raise FileNotFoundError(f"Launcher config not found: {e}")
    
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
            
            # Find all matching files in the directory
            matching_files = [file for file in os.listdir(dir_path) if file.endswith("_qc_stats_wgsadmin.xlsx")]
            
            # Skip the directory if no matching files are found
            if not matching_files:
                print(f"No '_qc_stats_wgsadmin.xlsx' files found in {dir_path}. Skipping...")
                continue
            
            # Initialize a list to store dataframes for this runtag
            if runtag not in runtag_dataframes:
                runtag_dataframes[runtag] = []
            
            # Process each matching file
            for file in matching_files:
                file_path = os.path.join(dir_path, file)
                
                # Read the Excel file into a dataframe
                try:
                    df = pd.read_excel(file_path)
                    runtag_dataframes[runtag].append(df)
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")

    os.makedirs(output_directory, exist_ok=True)

    # Combine and save the dataframes for each runtag
    for runtag, dataframes in runtag_dataframes.items():
        print(f"Processing runtag: {runtag} with {len(dataframes)} files")
        if dataframes:
            combined_df = pd.concat(dataframes, ignore_index=True)
            # Remove identical rows
            combined_df = combined_df.drop_duplicates()
            
            output_file_xlsx = os.path.join(output_directory, f"{runtag}_combined.xlsx")
            output_file_tsv = os.path.join(output_directory, f"{runtag}_combined.tsv")
            
            combined_df.to_excel(output_file_xlsx, index=False)
            combined_df.to_csv(output_file_tsv, sep='\t', index=False)
            print(f"Combined files saved for runtag {runtag}: {output_file_xlsx} and {output_file_tsv}")

if __name__ == "__main__":
    combine_qc_stats()