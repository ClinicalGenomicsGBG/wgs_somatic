import os
import sys
import time
import shutil
import logging
import subprocess
import click
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from helpers import read_config
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configure logging to write to a file
log_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "logs")
os.makedirs(log_dir, exist_ok=True)  # Ensure the logs directory exists

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(os.path.join(log_dir, "replace_bam_w_cram.log")),  # Log file name
        logging.StreamHandler()  # Also log to console
    ]
)

@click.command(help="Replace BAM files with CRAM files in webstore")
@click.option('--webstore_dir', type=click.Path(exists=True), required=True, default="/webstore/clinical/routine/wgs_somatic/current", help='Base webstore directory', show_default=True)
@click.option('--workdir', type=click.Path(), required=True, default="/clinical/data/wgs_somatic/bam_cleanup_workdir", help='Base execution directory where snakemake will create cram and crai', show_default=True)
@click.option('--age_threshold', type=int, default=50, help='Age threshold in days', show_default=True)
@click.option('--dry_run', is_flag=True, help='Perform a dry run without making any changes', show_default=True)
@click.option('--extra_snakemake_args', type=str, default="", help='Extra arguments for snakemake command', show_default=True)
@click.option('--launcher_config', type=click.Path(exists=True), required=True, default=os.path.join(os.path.dirname(__file__), "../../configs/launcher_config.json"), help='Path to the launcher config file', show_default=True)
@click.option('--snakemake_config', type=click.Path(exists=True), required=True, default=os.path.join(os.path.dirname(__file__), "../../configs/cluster.yaml"), help='Path to the snakemake config file', show_default=True)
@click.option('--keep_bam', is_flag=True, help='Do not delete BAM files after conversion', show_default=True)
def cli(webstore_dir, workdir, age_threshold, dry_run, extra_snakemake_args, launcher_config, snakemake_config, keep_bam):
    main(webstore_dir, workdir, age_threshold, dry_run, extra_snakemake_args, launcher_config, snakemake_config, keep_bam)

def is_older_than(file_path, age_threshold):
    """Check if a file is older than the specified number of days."""
    file_age = (time.time() - os.path.getmtime(file_path)) / (60 * 60 * 24)  # file_age is in days
    return file_age > age_threshold

def construct_singularity_args(binddirs, additional_binds):
    """Construct Singularity arguments for binding directories."""
    singularity_binds = ",".join(
        f"{bind['source']}:{bind['destination']}"
        for bind in binddirs.values()
    )
    additional_binds_str = ",".join(additional_binds)
    return f"-e --bind {singularity_binds},{additional_binds_str}"

def setup_snakemake(dir_to_process, workdir, launcher_config, snakemake_config, extra_snakemake_args, dry_run):
    """
    Set up and execute the Snakemake command for BAM to CRAM conversion.
    """
    # Load configuration
    config = read_config(launcher_config)  # Update with the correct path
    
    # Prepare bind directories for Singularity
    binddirs = config["singularitybinddirs"]
    additional_binds = [
        "/medstore",
        "/seqstore",
        "/apps",
        "/clinical",
        "/oldseqstore",
        "/webstore"
    ]
    singularity_args = construct_singularity_args(binddirs, additional_binds)
  
    cluster_args = [
        "qsub",
        "-S", "/bin/bash",
        "-pe", "mpi", "{cluster.threads}",
        "-q", "{cluster.queue}",
        "-N", "{cluster.name}",
        "-o", f"{workdir}/logs/{{cluster.output}}",
        "-e", f"{workdir}/logs/{{cluster.error}}",
        "-l", "{cluster.excl}"
    ]
    cluster_args_str = " ".join(cluster_args)

    snakemake_command = (
        "snakemake -s tools/bam_cleanup/convert_bam_to_cram.smk"
        f" --configfile {snakemake_config}"
        f" --use-singularity --singularity-args \"{singularity_args}\" "
        " --cluster-config configs/cluster.yaml"
        f" --cluster \"{cluster_args_str}\""
        " --jobs 999"
        " --latency-wait 60"
        f" --config dir_to_process={dir_to_process} launcher_config_path={launcher_config}"
        f" --directory {workdir} "
        f"{extra_snakemake_args if extra_snakemake_args else ''}"
    )

    if dry_run:
        snakemake_command += " --dry-run -p"
    
    return(snakemake_command)

def postprocess_directory(processed_directory, processed_complete_workdir, process, lock_file, keep_bam):
    """
    Perform postprocessing for a single directory after Snakemake execution.
    """
    process.wait()
    if process.returncode == 0:
        logging.info(f"Pipeline for directory {processed_directory} completed successfully.")
        # Transfer created files from complete_workdir to webstore_dir in parallel
        with ThreadPoolExecutor() as executor:
            futures = [
            executor.submit(
                shutil.copy,
                os.path.join(processed_complete_workdir, file_name),
                processed_directory
            )
            for file_name in os.listdir(processed_complete_workdir)
            if ".cram" in file_name
            ]

            for future in as_completed(futures):
                try:
                    source_file = future.result()  # shutil.copy returns the destination path
                    logging.info(f"Moved {source_file} to {processed_directory}.")
                except Exception as e:
                    logging.error(f"Error moving file: {e}")
        # Delete BAM files from the webstore processed_directory
        if not keep_bam:
            for file_name in os.listdir(processed_directory): #get bam files in the webstore directory
                if file_name.endswith(".bam"):
                    bam_file = os.path.join(processed_directory, file_name)
                    corresponding_cram = bam_file.replace(".bam", ".cram")
                    if os.path.exists(corresponding_cram): # if the corresponding cram file exists, delete the bam file
                        os.remove(bam_file)
                        logging.info(f"Deleted BAM file: {bam_file}.")
                    else:
                        logging.warning(f"CRAM file not found for BAM file: {bam_file}. Skipping deletion.")

                if file_name.endswith(".bai"):
                    bai_file = os.path.join(processed_directory, file_name)
                    corresponding_crai = bai_file.replace(".bam.bai", ".cram.crai")
                    if os.path.exists(corresponding_crai): # if the corresponding crai file exists, delete the bai file
                        os.remove(bai_file)
                        logging.info(f"Deleted BAI file: {bai_file}.")
                    else:
                        logging.warning(f"CRAM file not found for BAI file: {bai_file}. Skipping deletion.")
        else:
            logging.info(f"Keeping BAM files in directory: {processed_directory}.")
        
        # Delete the complete_workdir if the transfer and deletion are successful
        try:
            shutil.rmtree(processed_complete_workdir)
            logging.info(f"Deleted working directory: {processed_complete_workdir}.")
        except Exception as e:
            logging.error(f"Failed to delete working directory {processed_complete_workdir}: {e}")
    else:
        logging.error(f"Pipeline for directory {processed_directory} failed. Skipping file transfer and BAM deletion.")

    # Remove the lock file from the webstore directory
    if os.path.exists(lock_file):
        os.remove(lock_file)
        logging.info(f"Removed lock file for directory {processed_directory}.")

def main(webstore_dir, workdir, age_threshold, dry_run, extra_snakemake_args, launcher_config, snakemake_config, keep_bam):
    """Scan directories and process specific subdirectories."""

    if dry_run:
        logging.info("\n##### Dry run mode: No changes will be made. #####")
    else:
        logging.info("\n##### Processing directories. #####")

    directories_to_process = []

    for root_webstore, _, _ in os.walk(webstore_dir):
        logging.info(f"Scanning directory: {root_webstore}")
        for subdir, _, files in os.walk(root_webstore):
            bam_files = [os.path.join(subdir, file) for file in files if file.endswith(".bam")]
            if bam_files:  # If BAM files are found in the directory
                lock_file = os.path.join(subdir, ".bam2cram.processing.lock")
                
                # Check if the lock file exists and is not stale
                if os.path.exists(lock_file):
                    lock_age = time.time() - os.path.getmtime(lock_file)
                    if lock_age > 24 * 60 * 60:  # 24 hours
                        logging.warning(f"Stale lock file found for {subdir}. Removing it.")
                        os.remove(lock_file)
                    else:
                        logging.info(f"Skipping directory {subdir} because it is already being processed.")
                        continue
                
                # Check if CRAM files already exist for all BAM files in the directory
                crams_exist = any(
                    os.path.exists(os.path.join(subdir, bam_file.replace(".bam", ".cram")))
                    for bam_file in bam_files
                )
                if crams_exist:
                    logging.info(f"CRAM files already exist for at least one of the BAM files in {subdir}. Skipping processing. The existing BAMs will not be deleted.")
                    continue  # Skip processing for this directory
                
                # Check if all BAM files are older than the age threshold
                if all(is_older_than(bam_file, age_threshold) for bam_file in bam_files):
                    directories_to_process.append(subdir)

    directories_to_process = list(set(directories_to_process))  # Remove duplicates

    if not directories_to_process:
        logging.info("No new directories to process.")
        return
    logging.info("Directories to process:")
    for directory in directories_to_process:
        logging.info(f"- {directory}")

    processes = [] # List to track processes
    dry_run_directories = []  # Track directories created during dry run

    # submit snakemake jobs for each directory
    for directory in directories_to_process:
        # Skip lock file creation during dry run
        if not dry_run:
            # Create a hidden lock file in the directory to avoid that two runs of the script process the same directory at the same time
            lock_file = os.path.join(directory, ".bam2cram.processing.lock")

            # Create the lock file
            with open(lock_file, "w") as f:
                f.write(f"Processing started at {time.strftime('%Y-%m-%d %H:%M:%S')}")
            os.chmod(lock_file, 0o666)  # Set permissions so that everyone can read and remove the file

        random_id = uuid.uuid4().hex[:8]
        complete_workdir = os.path.join(workdir, random_id) # Create a unique workdir for each snakemake run/webstore directory to be processed

        # Track directories created during dry run
        if dry_run:
            dry_run_directories.append(complete_workdir)

        snakemake_command = setup_snakemake(
            dir_to_process=directory, # webstore directory to process. The bam files are in this directory and the cram/crai will be copied to here
            workdir=complete_workdir, # directory where snakemake will create cram and crai
            launcher_config=launcher_config,
            snakemake_config=snakemake_config,
            extra_snakemake_args=extra_snakemake_args,
            dry_run=dry_run
        )
        if dry_run:
            logging.info(f"Dry run: {snakemake_command}")
            try:
                result = subprocess.run(snakemake_command, shell=True, check=True, capture_output=True, text=True)
                logging.info(f"Command output: {result.stdout}")
                if result.stderr:
                    logging.error(f"Command error output: {result.stderr}")
            except subprocess.CalledProcessError as e:
                logging.error(f"Dry run command failed with exit code {e.returncode}: {e.stderr}")
        else:
            logging.info(f"Executing Snakemake command: {snakemake_command}")
            process = subprocess.Popen(snakemake_command, shell=True) # the snakemake pipeline for each directory is run in parallel and submitted as a job
            processes.append((directory, complete_workdir, process, lock_file))
            time.sleep(1)

    # Clean up directories created during dry run (the dry run is only used to print out which directories will be processed)
    # and to check if the snakemake command is set up correctly
    if dry_run:
        for dry_run_dir in dry_run_directories:
            if os.path.exists(dry_run_dir):
                try:
                    shutil.rmtree(dry_run_dir)
                    logging.info(f"Deleted dry run directory: {dry_run_dir}")
                except Exception as e:
                    logging.error(f"Failed to delete dry run directory {dry_run_dir}: {e}")
        logging.info("##### Dry run completed. #####")
        return

    # Process results after snakemake execution in parallel
    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(
                postprocess_directory,
                processed_directory,
                processed_complete_workdir,
                process,
                lock_file,
                keep_bam
            )
            for processed_directory, processed_complete_workdir, process, lock_file in processes
        ]

        for future in as_completed(futures):
            try:
                future.result()  # This will raise any exception that occurred during execution
            except Exception as e:
                logging.error(f"Error during postprocessing: {e}")

    logging.info("##### Processing completed. #####")

if __name__ == "__main__":
    cli()
