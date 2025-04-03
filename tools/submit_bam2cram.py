import os
import time
import shutil
import logging
import subprocess
import click
from helpers import read_config
import uuid

# Configure logging to write to a file
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(os.path.join(os.path.dirname(os.path.abspath(__file__)), "replace_bam_w_cram.log")),  # Log file name
        logging.StreamHandler()  # Also log to console
    ]
)

@click.command(help="Replace BAM files with CRAM files in webstore")
@click.option('--webstore_dir', type=click.Path(exists=True), required=True, default="/webstore/clinical/routine/wgs_somatic/current", help='Base webstore directory', show_default=True)
@click.option('--workdir', type=click.Path(), required=True, default="/clinical/data/wgs_somatic/bam_cleanup_workdir", help='Base execution directory where snakemake will create cram and crai', show_default=True)
@click.option('--age_threshold', type=int, default=50, help='Age threshold in days', show_default=True)
@click.option('--dry_run', is_flag=True, help='Perform a dry run without making any changes', show_default=True)
@click.option('--extra_snakemake_args', type=str, default="", help='Extra arguments for snakemake command', show_default=True)
@click.option('--launcher_config', type=click.Path(exists=True), required=True, default=os.path.join(os.path.dirname(__file__), "../configs/launcher_config.json"), help='Path to the launcher config file', show_default=True)
@click.option('--snakemake_config', type=click.Path(exists=True), required=True, default=os.path.join(os.path.dirname(__file__), "../configs/cluster.yaml"), help='Path to the snakemake config file', show_default=True)
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
        "/oldseqstore"
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
        "snakemake -s tools/convert_bam_to_cram.smk"
        f" --configfile {snakemake_config}"
        f" --use-singularity --singularity-args \"{singularity_args}\" "
        " --cluster-config configs/cluster.yaml"
        f" --cluster \"{cluster_args_str}\""
        " --jobs 999"
        " --latency-wait 60"
        f" --config dir_to_process={dir_to_process}"
        f" --directory {workdir} "
        f"{extra_snakemake_args if extra_snakemake_args else ''}"
    )

    if dry_run:
        snakemake_command += " --dry-run -p"
    
    return(snakemake_command)


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
            if any(file.endswith(".bam") for file in files):
                logging.info(f"Found BAM files in directory: {subdir}")
                directories_to_process.append(subdir)
    directories_to_process = list(set(directories_to_process))  # Remove duplicates
    
    # Filter directories based on age threshold
    directories_to_process = [
        d for d in directories_to_process if is_older_than(d, age_threshold)
    ]
    if not directories_to_process:
        logging.info("No new directories to process.")
        return
    logging.info("Directories to process:")
    for directory in directories_to_process:
        logging.info(f"- {directory}")

    processes = []
    dry_run_directories = []  # Track directories created during dry run

    for directory in directories_to_process:
        # Skip lock file creation during dry run
        if not dry_run:
            # Create a hidden lock file in the directory
            lock_file = os.path.join(directory, ".bam2cram.processing.lock")
            if os.path.exists(lock_file):
                lock_age = time.time() - os.path.getmtime(lock_file)
                if lock_age > 24 * 60 * 60:  # 24 hours
                    logging.warning(f"Stale lock file found for {directory}. Removing it.")
                    os.remove(lock_file)
                else:
                    logging.info(f"Skipping directory {directory} because it is already being processed.")
                    continue

            # Create the lock file
            with open(lock_file, "w") as f:
                f.write(f"Processing started at {time.strftime('%Y-%m-%d %H:%M:%S')}")

        random_id = uuid.uuid4().hex[:8]
        complete_workdir = os.path.join(workdir, random_id)

        # Track directories created during dry run
        if dry_run:
            dry_run_directories.append(complete_workdir)

        snakemake_command = setup_snakemake(
            dir_to_process=directory,
            workdir=complete_workdir,
            launcher_config=launcher_config,
            snakemake_config=snakemake_config,
            extra_snakemake_args=extra_snakemake_args,
            dry_run=dry_run
        )
        if dry_run:
            logging.info(f"Dry run: {snakemake_command}")
            os.system(snakemake_command)
        else:
            logging.info(f"Executing Snakemake command: {snakemake_command}")
            process = subprocess.Popen(snakemake_command, shell=True)
            processes.append((directory, complete_workdir, process, lock_file))
            time.sleep(1)

    # Clean up directories created during dry run
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

    # Process results after actual execution
    for processed_directory, processed_complete_workdir, process, lock_file in processes:
        process.wait()
        if process.returncode == 0:
            logging.info(f"Pipeline for directory {processed_directory} completed successfully.")
            # Transfer created files from complete_workdir to webstore_dir
            for file_name in os.listdir(processed_complete_workdir):
                if ".cram" in file_name:
                    source_file = os.path.join(processed_complete_workdir, file_name)
                    shutil.copy(source_file, processed_directory)
                    logging.info(f"Moved {source_file} to {processed_directory}.")
            
            # Delete BAM files from the webstore processed_directory
            if not keep_bam:
                for file_name in os.listdir(processed_directory):
                    if file_name.endswith(".bam"):
                        bam_file = os.path.join(processed_directory, file_name)
                        os.remove(bam_file)
                        logging.info(f"Deleted BAM file: {bam_file}.")
                        
                    if file_name.endswith(".bai"):
                        bai_file = os.path.join(processed_directory, file_name)
                        os.remove(bai_file)
                        logging.info(f"Deleted BAI file: {bai_file}.")
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

        # Remove the lock file
        if not dry_run and os.path.exists(lock_file):
            os.remove(lock_file)
            logging.info(f"Removed lock file for directory {processed_directory}.")

    logging.info("##### Processing completed. #####")

if __name__ == "__main__":
    cli()
