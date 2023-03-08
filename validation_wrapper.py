#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python
import json
import argparse
import os
import glob
import helpers
import sys
import time
import traceback
from shutil import copyfile
import subprocess
from definitions import ROOT_DIR

def read_wrapperconf():
    with open(f"{ROOT_DIR}/configs/wrapper_conf.json", 'r') as configfile:
        config_data = json.load(configfile)
        return config_data

def get_time():
    nowtime = time.strftime("%Y-%m-%d-%H-%M-%S")
    return nowtime

def logger(message, logfile=False):
    config = read_wrapperconf()
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


config = read_wrapperconf()
valconf = f"{ROOT_DIR}/validation_scripts/validation_config.json"
valconf = helpers.read_config(valconf)
current_date = time.strftime("%Y-%m-%d")

for trange in valconf["fractions"]["hg38"]:
    
    ##################################################################
    # Create AnalysisConfigfile
    ##################################################################
    analysisdict = {}
    analysisdict["data"] = valconf["fractions"]["hg38"][trange]
    analysisdict["workingdir"] = valconf["outputdir"] 
    analysisdict["rtg"] = {}
    analysisdict["rtg"]["tools"] = valconf["rtgtools"]["rtg"] 
    analysisdict["rtg"]["sdf"]  = valconf["rtgtools"]["hg38"]["sdf"]


    output = valconf["outputdir"]
    if not os.path.isdir(f"{output}/{trange}"):
        os.mkdir(f"{output}/{trange}")
    
#    for tnsetting in valconf["tnscope_settings"]:
#        if not os.path.isdir(f"{output}/{trange}/{tnsetting}"):
#            os.mkdir(f"{output}/{trange}/{tnsetting}")
        
    runnormal = valconf["fractions"]["hg38"][trange]["runname"]
    runtumor = valconf["fractions"]["hg38"][trange]["runname"]
    normalname = valconf["fractions"]["hg38"][trange]["normalname"]
    tumorname = valconf["fractions"]["hg38"][trange]["tumorname"]

    date, _, _, chip, *_ = runnormal.split('_')
    normalid= '_'.join([normalname, date, chip])
    date, _, _, chip, *_ = runtumor.split('_')
    tumorid = '_'.join([tumorname, date, chip])       

    # AnalysisConfig Continued...
    output = f"{output}/{trange}"
    analysisdict["tnscope"] = valconf["tnscope_settings"]
    analysisdict["sample"] = {}
    analysisdict["sample"]["normalname"] = normalname
    analysisdict["sample"]["tumorname"] = tumorname
    analysisdict["sample"]["tumorid"] = tumorid
    analysisdict["sample"]["normalid"] = normalid
    analysisdict["tumorfastqs"] = [valconf["fractions"]["hg38"][trange]["tumor"]]
    analysisdict["normalfastqs"] = [valconf["fractions"]["hg38"][trange]["normal"] ]
    analysisdict["workingdir"] = output

    #################################################################
    # Prepare AnalysisFolder
    #################################################################
    date, _, _, chip, *_ = runnormal.split('_')
    normalid= '_'.join([normalname, date, chip])
    date, _, _, chip, *_ = runtumor.split('_')
    tumorid = '_'.join([tumorname, date, chip])

    
    mainconf = "hg38conf"
    configdir = config["configdir"]
    mainconf_name = config[mainconf]
    mainconf_path = f"{configdir}/{mainconf_name}"

    samplelogs = f"{output}/logs"
    if not os.path.isdir(samplelogs):
        os.mkdir(samplelogs)
    runconfigs = f"{output}/configs"
    if not os.path.isdir(runconfigs):
        os.mkdir(runconfigs)

    # copying configfiles to analysisdir
    clusterconf = config["clusterconf"]
    copyfile(f"{configdir}/{clusterconf}", f"{runconfigs}/{clusterconf}")
    copyfile(f"{configdir}/{mainconf_name}", f"{runconfigs}/{mainconf_name}")

    samplelog = f"{samplelogs}/{tumorid}.log"
    logger("Starting validation analysis with information:", samplelog)
    logger(f"{analysisdict}", samplelog)
    
    with open(f"{runconfigs}/{tumorid}_config.json", 'w') as analysisconf:
        json.dump(analysisdict, analysisconf, ensure_ascii=False, indent=4)

    ###################################################################
    # Prepare Singularity Binddirs
    binddirs = config["singularitybinddirs"]
    binddir_string = ""
    for binddir in binddirs:
        source = binddirs[binddir]["source"]
        destination = binddirs[binddir]["destination"]
        logger(f"preparing binddir variable {binddir} source: {source} destination: {destination}")
        binddir_string = f"{binddir_string}{source}:{destination},"
        for normalfastqdir in analysisdict["normalfastqs"]:
             binddir_string = f"{binddir_string}{normalfastqdir},"
        for tumorfastqdir in analysisdict["tumorfastqs"]:
            binddir_string = f"{binddir_string}{tumorfastqdir},"
    binddir_string = f"{binddir_string}{output}"


    ###################################################################
    # Start SnakeMake pipeline
    ###################################################################
    scriptdir = os.path.dirname(os.path.realpath(__file__)) # find current dir

    snakemake_path = config["snakemake_env"]
    os.environ["PATH"] += os.pathsep + snakemake_path
    my_env = os.environ.copy()
    snakemake_args = f"snakemake -s tnscope_eval.snakefile --configfile {runconfigs}/{tumorid}_config.json --dag | dot -Tsvg > {samplelogs}/dag_{current_date}.svg"
    # >>>>>>>>>>>> Create Dag of pipeline
    subprocess.run(snakemake_args, shell=True, env=my_env) # CREATE DAG
    snakemake_args = f"snakemake -s tnscope_eval.snakefile --configfile {runconfigs}/{tumorid}_config.json --use-singularity --singularity-args '--bind {binddir_string}' --cluster-config {ROOT_DIR}/configs/cluster.yaml --cluster \"qsub -S /bin/bash -pe mpi {{cluster.threads}} -q {{cluster.queue}} -N {{cluster.name}} -o {samplelogs}/{{cluster.output}} -e {samplelogs}/{{cluster.error}} -l {{cluster.excl}}\" --jobs 999 --latency-wait 60 --directory {output} &>> {samplelog}"
    # >>>>>>>>>>>> Start pipeline
    subprocess.run(snakemake_args, shell=True, env=my_env) # Shellscript pipeline
