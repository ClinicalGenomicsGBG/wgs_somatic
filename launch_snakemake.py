#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python
import json
import argparse
import os
import re
import glob
from tools.helpers import read_config
import sys
import time
import traceback
from shutil import copyfile, copy
import subprocess
import stat
import yaml
import random
import string
from definitions import LAUNCHER_CONFIG_PATH


def get_time():
    nowtime = time.strftime("%Y-%m-%d-%H-%M-%S")
    return nowtime


def get_timestamp():
    return time.strftime("%y%m%d-%H%M%S")


def logger(message, logfile=False):
    config = read_config(LAUNCHER_CONFIG_PATH)
    logdir = config["logdir"]
    current_date = time.strftime("%Y-%m-%d")
    if not logfile:
        logname = f"{logdir}/{current_date}.log"
        logfile = open(logname, "a+")
        logfile.write(f"{get_time()}: {message}" + "\n")
    else:
        logfile = open(logfile, "a+")
        logfile.write(f"{get_time()}: {message}" + "\n")
    print(message)


def get_normalid_tumorid(normalfastqs=None, normalname=None, tumorfastqs=None, tumorname=None):
    '''Get tumorid and normalid based on normalfastq+normalname and tumorfastq+tumorname'''
    def get_id_from_fastq(fastq_dir, name):
        for filename in os.listdir(fastq_dir):
            if filename.startswith(name):
                parts = filename.split("_")
                if len(parts) >= 3:
                    return "_".join(parts[:3])  # Returns an id based on the first hit
                else:
                    # Add filler if the filename cannot be split into at least 4 parts
                    logger(f"Filename {filename} could not be split into at least 4 parts")
                    logger(f"Returning filler id for {name}")
                    return "_".join(parts + ["filler"] * (3 - len(parts)))
        raise ValueError(f"No fastq found for {name} in {fastq_dir}")
    # End of get_id_from_fastq

    normalid = None
    tumorid = None

    if normalfastqs and normalname:
        normalid = get_id_from_fastq(normalfastqs, normalname)

    if tumorfastqs and tumorname:
        tumorid = get_id_from_fastq(tumorfastqs, tumorname)

    return normalid, tumorid


def yearly_stats(tumorname, normalname):
    '''Update yearly stats file with sample that has finished running in pipeline correctly'''
    yearly_stats = "yearly_stats.txt"
    try:
        if not os.path.exists(yearly_stats):
            os.mknod(yearly_stats)
            os.chmod(yearly_stats, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
    except FileExistsError:
        # Ignore error when multiple threads try to create the file at the same time
        pass

    with open(yearly_stats, "a") as yearly_stats_file:
        date_time = time.strftime("%Y-%m-%d-%H-%M-%S")
        yearly_stats_file.write(f"Tumor ID: {tumorname} Normal ID: {normalname} Date and Time: {date_time}\n")


def copy_results(outputdir, runconfig=None):
    '''Rsync result files from outputdir to resultdir'''
    try:
        if not runconfig:
            # Automatic detection of runconfig
            config_dir = os.path.join(outputdir, 'configs')
            config_pattern = re.compile(r'DNA[\dA-Za-z]+_.+_.+_config\.json')

            if os.path.isdir(config_dir):
                for f in os.listdir(config_dir):
                    if config_pattern.match(f):
                        runconfig = os.path.join(config_dir, f)
                        logger(f"Found runconfig: {runconfig}")
                        break
                else:
                    logger(f"Automatic detection of config file failed. No matching configuration file found in {config_dir}")
                    logger(f"Pattern used to match the config file: {config_pattern.pattern}")
                    raise ValueError("Automatic detection of config file failed. No matching configuration file found.")
        else:
            if not os.path.isfile(runconfig):
                logger(f"Provided runconfig file does not exist: {runconfig}")
                raise FileNotFoundError(f"Runconfig file not found: {runconfig}")

        try:
            # Read the runconfig file to get resultdir and resultsconf
            with open(runconfig, 'r') as cf:
                config_data = json.load(cf)
                try:
                    resultdir = config_data.get('resultdir')
                    logger(f"Resultdir found in config file: {resultdir}")
                except KeyError:
                    logger(f"Key 'resultdir' not found in config file {runconfig}")
                    raise KeyError("Key 'resultdir' not found in config file")
                try:
                    resultsconf = config_data.get('resultfilesconf')
                    logger(f"Results configuration file found: {resultsconf}")
                except KeyError:
                    logger(f"Key 'resultfilesconf' not found in config file {runconfig}")
                    raise KeyError("Key 'resultfilesconf' not found in config file")
        except Exception as e:
            logger(f"Error reading config file {runconfig}: {e}")
            raise

        try:
            results = read_config(resultsconf)
        except Exception as e:
            logger(f"Failed to read results config from {resultsconf}: {e}")
            raise

        # Read the results into the subdirectories unless they are marked toplevel
        for category, relpaths in results.items():
            dest_dir = resultdir if category == "toplevel" else os.path.join(resultdir, category)

            try:
                os.makedirs(dest_dir, exist_ok=True)
            except Exception as e:
                logger(f"Error creating result directory {dest_dir}: {e}")
                raise

            for relpath in relpaths:
                src_path = os.path.join(outputdir, relpath)
                dest_path = os.path.join(dest_dir, os.path.basename(relpath))

                if os.path.exists(src_path):
                    try:
                        copyfile(src_path, dest_path)
                        logger(f"Copied {src_path} to {dest_path}")

                        # Remove .bam and .bai files from outputdir once copied
                        if src_path.endswith('.bam') or src_path.endswith('.bai'):
                            try:
                                logger(f"Removing {src_path} from outputdir.")
                                os.remove(src_path)
                            except Exception as e:
                                logger(f"Error occurred while removing {src_path}: {e}")
                    except Exception as e:
                        logger(f"Error copying {src_path} to {dest_path}: {e}")
                else:
                    logger(f"Source file {src_path} does not exist, skipping copy.")

    except Exception as e:
        logger(f"Unhandled error in copy_results: {e}")
        raise


def analysis_main(args, outputdir, normalname=False, normalfastqs=False, tumorname=False, tumorfastqs=False, hg38ref=False, starttype=False, notemp=False):
    try:
        ################################################################
        # Write InputArgs to logfile
        ################################################################
        config = read_config(LAUNCHER_CONFIG_PATH)
        commandlogs = config["commandlogs"]
        os.makedirs(commandlogs, exist_ok=True)
        command = f"{sys.argv[0]}"
        current_date = time.strftime("%Y-%m-%d")
        commandlog = f"{commandlogs}/commands_{current_date}.log"
        for arg in vars(args):
            command = f"{command} --{arg} {getattr(args, arg)}"
        commandlogfile = open(commandlog, "a+")
        commandlogfile.write(f"{get_time()}" + "\n")
        commandlogfile.write(command + "\n")

        #################################################################
        # Validate Inputs
        ################################################################
        if outputdir.endswith("/"):
            outputdir = outputdir[:-1]
        if normalfastqs:
            if normalfastqs.endswith("/"):
                normalfastqs = normalfastqs[:-1]
        if tumorfastqs:
            if tumorfastqs.endswith("/"):
                tumorfastqs = tumorfastqs[:-1]

        error_list = []

        if hg38ref:
            logger(f"hg38 argument given with value: {hg38ref}")
            if hg38ref != "yes":
                logger("argument is not yes, if you want hg19 simply dont provide hg38 argument, exiting")
                error_list.append(f"invalid hg38 argument value: {hg38ref}")

        if hg38ref == "yes":
            mainconf = "hg38conf"
        else:
            mainconf = "hg19conf"
        configdir = config["configdir"]
        mainconf_name = config[mainconf]

        # validate fastqdirs
        if starttype == "force":
            f_tumorfastqs = ""
            f_normalfastqs = ""

        else:
            if normalfastqs:
                if not os.path.isdir(normalfastqs):
                    error_list.append(f"{normalfastqs} does not appear to be a directory")
                else:
                    f_normalfastqs = glob.glob(f"{normalfastqs}/*{normalname}*fastq.gz")
                    if not f_normalfastqs:
                        logger("Warning: No fastqs found in normaldir")
                        f_normalfastqs = glob.glob(f"{normalfastqs}/*{normalname}*fasterq")
                        if not f_normalfastqs:
                            error_list.append("No fastqs or fasterqs found in normaldir")

            if tumorfastqs:
                if not os.path.isdir(tumorfastqs):
                    error_list.append(f"{tumorfastqs} does not appear to be a directory")
                else:
                    f_tumorfastqs = glob.glob(f"{tumorfastqs}/*{tumorname}*fastq.gz")
                    if not f_tumorfastqs:
                        logger("Warning: No fastqs found in tumordir")
                        f_tumorfastqs = glob.glob(f"{tumorfastqs}/*{tumorname}*fasterq")
                        if not f_tumorfastqs:
                            error_list.append("No fastqs or fasterqs found in tumordir")

        # prepare outputdirectory
        if not os.path.isdir(outputdir):
            try:
                os.mkdir(outputdir)
            except Exception as e:
                error_list.append(f"outputdirectory: {outputdir} does not exist and could not be created: {e}")

        if error_list:
            logger("Errors found in arguments to script:")
            for arg in vars(args):
                logger(f"{arg} = {getattr(args, arg)}")
            for error in error_list:
                logger(error)
            logger("Exiting")
            sys.exit()

        #################################################################
        # Prepare AnalysisFolder
        #################################################################
        normalid, tumorid = get_normalid_tumorid(normalfastqs, normalname, tumorfastqs, tumorname)
        samplelogs = f"{outputdir}/logs"
        if not os.path.isdir(samplelogs):
            os.mkdir(samplelogs)
        runconfigs = os.path.join(outputdir, configdir)
        if not os.path.isdir(runconfigs):
            os.mkdir(runconfigs)

        # copying configfiles to analysisdir
        clusterconf = config["clusterconf"]
        filterconf = config["filterconf"]
        copyfile(os.path.join(configdir, clusterconf), os.path.join(runconfigs, clusterconf))
        copyfile(os.path.join(configdir, filterconf), os.path.join(runconfigs, filterconf))
        copyfile(os.path.join(configdir, mainconf_name), os.path.join(runconfigs, mainconf_name))

        # Use wildcards to fill in the resultfiles templates
        wildcards = {
                "tumor": "tumor",
                "normal": "normal",
                "DNAtumor": tumorid,
                "DNAnormal": normalid
            }

        # Read result_files templates from config
        resultsconf = config["resultfilesconf"]
        resultsconf_path = os.path.join(configdir, resultsconf)
        if tumorname and normalname:
            results_headers = read_config(resultsconf_path)["tumor-normal"]
        elif tumorname:
            results_headers = read_config(resultsconf_path)["tumor-only"]
        elif normalname:
            results_headers = read_config(resultsconf_path)["normal-only"]

        results = {
            cat: [ pattern.format(**wildcards) for pattern in patterns ]
            for cat, patterns in results_headers.items()
        }

        # Write the result_files configuration to the working directory
        with open(os.path.join(runconfigs, resultsconf), 'w') as results_file:
            yaml.dump(results, results_file)

        # Prepare log
        if tumorname:
            samplelog = f"{samplelogs}/{tumorid}.log"
        else:
            samplelog = f"{samplelogs}/{normalid}.log"
        logger("Input validated:", samplelog)
        logger(f"{command}", samplelog)
        if normalname:
            logger("Fastqs found for normal:", samplelog)
            logger(f"{f_normalfastqs}", samplelog)
        if tumorname:
            logger("Fastqs found for tumor:", samplelog)
            logger(f"{f_tumorfastqs}", samplelog)

        ##################################################################
        # Create AnalysisConfigfile
        ##################################################################
        analysisdict = {}
        analysisdict["normalname"] = normalname
        analysisdict["normalid"] = normalid
        analysisdict["normalfastqs"] = [normalfastqs]
        analysisdict["tumorname"] = tumorname
        analysisdict["tumorid"] = tumorid
        analysisdict["tumorfastqs"] = [tumorfastqs]

        # configs
        analysisdict["filterconfig"] = os.path.join(configdir, filterconf)
        analysisdict["clusterconfig"] = os.path.join(configdir, clusterconf)
        analysisdict["pipeconfig"] = os.path.join(configdir, mainconf_name)
        analysisdict["resultfilesconf"] = os.path.join(runconfigs, resultsconf)

        # insilico
        analysisdict["insilico"] = config["insilicopanels"]

        basename_outputdir = os.path.basename(outputdir)

        if hg38ref == "yes":
            analysisdict["reference"] = "hg38"
            if tumorname:
                if normalname:
                    analysisdict["resultdir"] = f'{config["resultdir_hg38"]}/{basename_outputdir}'  # Use f'{config["testresultdir"]}/{tumorname}'for testing
                else:
                    analysisdict["resultdir"] = f'{config["resultdir_hg38"]}/tumor_only/{basename_outputdir}'
            else:
                analysisdict["resultdir"] = f'{config["resultdir_hg38"]}/normal_only/{basename_outputdir}'
        else:
            analysisdict["reference"] = "hg19"
            if tumorname:
                if normalname:
                    analysisdict["resultdir"] = f'{config["resultdir_hg19"]}/{basename_outputdir}'
                else:
                    analysisdict["resultdir"] = f'{config["resultdir_hg19"]}/tumor_only/{basename_outputdir}'
            else:
                analysisdict["resultdir"] = f'{config["resultdir_hg19"]}/normal_only/{basename_outputdir}'

        if tumorname:
            snakemake_config = f"{runconfigs}/{tumorid}_config.json"
        else:
            snakemake_config = f"{runconfigs}/{normalid}_config.json"

        with open(snakemake_config, 'w') as analysisconf:
            json.dump(analysisdict, analysisconf, ensure_ascii=False, indent=4)

        ###################################################################
        # Prepare Singularity Binddirs
        ###################################################################
        binddirs = config["singularitybinddirs"]
        binddir_string = ""
        for binddir in binddirs:
            source = binddirs[binddir]["source"]
            if not analysisdict["reference"] in source:
                if "sentieon" not in source:
                    continue
            destination = binddirs[binddir]["destination"]
            logger(f"preparing binddir variable {binddir} source: {source} destination: {destination}")
            binddir_string = f"{binddir_string}{source}:{destination},"
            if normalname:
                for normalfastqdir in analysisdict["normalfastqs"]:
                    binddir_string = f"{binddir_string}{normalfastqdir},"
            if tumorname:
                for tumorfastqdir in analysisdict["tumorfastqs"]:
                    binddir_string = f"{binddir_string}{tumorfastqdir},"
        binddir_string = f"{binddir_string}{outputdir}"
        print(binddir_string)

    except Exception as e:
        tb = traceback.format_exc()
        logger("Error in setting up the snakemake run:")
        logger(f"{e} Traceback: {tb}")
        sys.exit(1)

    try:
        ###################################################################
        # Start SnakeMake pipeline
        ###################################################################
        # Generate random hash for shadow directory
        letters = string.ascii_lowercase
        letters = ''.join(random.choice(letters) for i in range(10))
        shadow_dir = os.path.join('/tmp', letters)

        snakemake_path = config["snakemake_env"]
        os.environ["PATH"] += os.pathsep + snakemake_path
        my_env = os.environ.copy()

        # Set development arguments
        notemp_arg = ["--notemp"] if notemp else []

        # Construct Snakemake command

        singularity_args = [
            "-e",
            "--bind", binddir_string,
            "--bind", "/medstore",
            "--bind", "/seqstore",
            "--bind", "/apps",
            "--bind", "/clinical",
            "--bind", "/webstore",
        ]

        cluster_args = [
            "qsub",
            "-S", "/bin/bash",
            "-pe", "mpi", "{cluster.threads}",
            "-q", "{cluster.queue}",
            "-N", "{cluster.name}",
            "-o", f"{samplelogs}/{{cluster.output}}",
            "-e", f"{samplelogs}/{{cluster.error}}",
            "-l", "{cluster.excl}"
        ]

        # Create DAG of pipeline
        snakemake_args_DAG = [
            "snakemake", "-s", "Snakefile",
            "--configfile", snakemake_config,
            "--directory", outputdir,
            "--dag", "|", "dot", "-Tsvg", ">", f"{samplelogs}/dag_{current_date}.svg"
        ]
        snakemake_args_DAG_str = " ".join(snakemake_args_DAG)
        subprocess.run(snakemake_args_DAG_str, shell=True, env=my_env)

        snakemake_args = [
            "snakemake", "-s", "Snakefile",
            "--configfile", snakemake_config,
            "--use-singularity", "--singularity-args", " ".join(singularity_args),
            "--cluster-config", "configs/cluster.yaml",
            "--cluster", " ".join(cluster_args),
            "--jobs", "999",
            "--latency-wait", "60",
            "--directory", outputdir,
            "--shadow-prefix", shadow_dir,
            "--rerun-incomplete",
            "--stats", f"{samplelogs}/stats_{current_date}.json"
        ] + notemp_arg

        # Execute Snakemake command with outputdir redirection
        with open(samplelog, "a") as log_file:
            subprocess.run(snakemake_args, env=my_env, check=True, stdout=log_file, stderr=log_file)

    except subprocess.CalledProcessError as e:
        tb = traceback.format_exc()
        logger(f"Error running Snakemake: {e}")
        logger(f"Traceback: {tb}")
        sys.exit(1)
    except Exception as e:
        tb = traceback.format_exc()
        logger(f"An error occurred: {e}\n{tb}")
        sys.exit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outputdir', nargs='?', help='output directory, where to put results', required=True)
    parser.add_argument('-ns', '--normalsample', nargs='?', help='normal samplename', required=False)
    parser.add_argument('-nf', '--normalfastqs', nargs='?', help='path to directory containing normal fastqs', required=False)
    parser.add_argument('-ts', '--tumorsample', nargs='?', help='tumor samplename', required=False)
    parser.add_argument('-tf', '--tumorfastqs', nargs='?', help='path to directory containing tumor fastqs', required=False)
    parser.add_argument('-hg38', '--hg38ref', nargs='?', help='run analysis on hg38 reference (write yes if you want this option)', required=False)
    parser.add_argument('-stype', '--starttype', nargs='?', help='write forcestart if you want to ignore fastqs', required=False)
    parser.add_argument('-cr', '--copyresults', action="store_true", help='Copy results to resultdir on webstore', required=False)
    parser.add_argument('--notemp', action="store_true", help='Run the pipeline in notemp mode (all intermediate files kept)', required=False)
    parser.add_argument('-onlycopy', '--onlycopyresults', action="store_true", help='Only run the copy_results function', required=False)
    args = parser.parse_args()

    if not args.outputdir.startswith("/"):
        args.outputdir = os.path.abspath(args.outputdir)
        logger(f"Adjusted outputdir to {args.outputdir}")
    if args.onlycopyresults:
        copy_results(args.outputdir)
    else:
        if args.tumorfastqs:
            if not args.tumorfastqs.startswith("/"):
                args.tumorfastqs = os.path.abspath(args.tumorfastqs)
                logger(f"Adjusted tumorfastqs to {args.tumorfastqs}")
        if args.normalfastqs:
            if not args.normalfastqs.startswith("/"):
                args.normalfastqs = os.path.abspath(args.normalfastqs)
                logger(f"Adjusted normalfastqs to {args.normalfastqs}")
        analysis_main(args, args.outputdir, args.normalsample, args.normalfastqs, args.tumorsample, args.tumorfastqs, args.hg38ref, args.starttype, args.notemp)

        if os.path.isfile(f"{args.outputdir}/reporting/workflow_finished.txt"):
            if args.tumorsample:
                if args.normalsample:
                    # these functions are only executed if snakemake workflow has finished successfully
                    yearly_stats(args.tumorsample, args.normalsample)
                    if args.copyresults:
                        copy_results(args.outputdir)
                else:
                    yearly_stats(args.tumorsample, 'None')
                    if args.copyresults:
                        copy_results(args.outputdir)
            else:
                yearly_stats('None', args.normalsample)
                if args.copyresults:
                    copy_results(args.outputdir)
