#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python
import json
import argparse
import os
import glob 
import helpers
import sys
import time
import traceback
from shutil import copyfile, copy
import subprocess
import stat
import requests
import random
import string
from definitions import ROOT_DIR

def read_wrapperconf():
    with open(f"{ROOT_DIR}/configs/wrapper_conf.json", 'r') as configfile:
        config_data = json.load(configfile)
        return config_data


def get_time():
    nowtime = time.strftime("%Y-%m-%d-%H-%M-%S")
    return nowtime

def get_timestamp():
    return time.strftime("%y%m%d-%H%M%S")

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


def get_normalid_tumorid(runnormal=None, normalname=None, runtumor=None, tumorname=None):
    '''Get tumorid and normalid based on runnormal/tumor and normal/tumorname '''
    if runnormal:
        date, _, _, chip, *_ = runnormal.split('+')[0].split('_')
    if normalname:
        normalid= '_'.join([normalname, date, chip])
    else:
        normalid = None
    if runtumor:
        date, _, _, chip, *_ = runtumor.split('+')[0].split('_')
    if tumorname:
        tumorid = '_'.join([tumorname, date, chip])
    else:
        tumorid = None
    return normalid, tumorid

def yearly_stats(tumorname, normalname):
    '''Update yearly stats file with sample that has finished running in pipeline correctly'''
    #config_data = read_wrapperconf()
    #yearly_stats = open(config_data["yearly_stats"], "a")
    yearly_stats = "yearly_stats.txt"
    if not os.path.exists(yearly_stats):
        os.mknod(yearly_stats)
        os.chmod(yearly_stats, stat.S_IRWXU|stat.S_IRWXG|stat.S_IRWXO)
    yearly_stats = open(yearly_stats, "a")
    date_time = time.strftime("%Y-%m-%d-%H-%M-%S")
    yearly_stats.write("Tumor ID: " + tumorname + " Normal ID: " + normalname + " Date and Time: " + date_time + "\n")
    yearly_stats.close()


def alissa_upload(outputdir, normalname, runnormal, ref=False):
    '''Upload germline SNV_CNV vcf to Alissa'''
    ref = 'hg38' # reference genome, change so it can also be hg19. but probably shouldn't upload vcf for hg19 samples
    size = '240_000_000' # needed as argument for alissa upload. if vcf is larger than this, it is split
    date, _, _, chip, *_ = runnormal.split('_')
    normalid = f"{normalname}_{date}_{chip}"
    vcfpath = f"{outputdir}/{normalid}_{ref}_SNV_CNV_germline.vcf.gz"
    config = read_wrapperconf()
    logger(f"Uploading vcf to Alissa for {normalname}")
    queue = config["alissa"]["queue"]
    threads = config["alissa"]["threads"]
    qsub_script = config["alissa"]["qsub_script"]
    standardout = f"{outputdir}/logs/{normalid}_alissa_upload_standardout.txt"
    standarderr = f"{outputdir}/logs/{normalid}_alissa_upload_standarderr.txt"
    qsub_args = ["qsub", "-N", f"WGSSomatic-{normalname}_alissa_upload", "-q", queue, "-o", standardout, "-e", standarderr, qsub_script, normalname, vcfpath, size, ref]
    subprocess.call(qsub_args, shell=False)

def copy_results(outputdir, runnormal=None, normalname=None, runtumor=None, tumorname=None):
    '''Rsync result files from workingdir to resultdir'''

    config = read_wrapperconf()

    # Find correct resultdir on webstore from sample config in workingdir
    normalid, tumorid = get_normalid_tumorid(runnormal, normalname, runtumor, tumorname)
    if tumorid:
        with open(f"{outputdir}/configs/{tumorid}_config.json", "r") as analysisdict:
            analysisdict = json.load(analysisdict)
    else:
        with open(f"{outputdir}/configs/{normalid}_config.json", "r") as analysisdict:
            analysisdict = json.load(analysisdict)
    resultdir = analysisdict["resultdir"]
    workdir = analysisdict["workingdir"]
    os.makedirs(resultdir, exist_ok=True)
    igv_dir = os.path.join(resultdir, 'igv_files')
    os.makedirs(igv_dir, exist_ok=True)
    

    # Find resultfiles to copy to resultdir on webstore
    copy_files = []
    files_match = ['.xlsx', 'CNV_SNV_germline.vcf.gz', 'somatic.vcf.gz', 'refseq3kfilt.vcf.gz']
    for files in files_match:
        copy_files = copy_files + glob.glob(os.path.join(workdir, f'*{files}*'))
    copy_files = set(copy_files)
    for f in os.listdir(workdir):
        f = os.path.join(workdir, f)
        if os.path.isdir(f):
            if 'configs' in f: 
                copy(os.path.join(workdir, 'configs', 'config_hg38.json'), resultdir)
                # copy config file to resultdir
                logger(f"Run configuration file copied successfully")
        if os.path.isfile(f):
            if f in copy_files:
                try:
                    copy(f, resultdir)
                    logger(f"{f} copied successfully")
                except:
                    logger(f"Error occurred while copying {f}")
            else:
                # Copy all files that are not cram and that do not match the files_match pattern to webstore igv_dir
                if not f.endswith('.cram') and not f.endswith('.crai'):
                    try:
                        copy(f, igv_dir)
                        logger(f"{f} copied successfully")
                        if f.endswith('.bam') or f.endswith('.bai'):
                            try:
                                # Remove the bam files from the workingdir
                                logger(f"Removing {f} from workdir.")
                                os.remove(f)
                            except:
                                logger(f"Error occurred while removing {f}")
                    except:
                        logger(f"Error occurred while copying {f}")

    # Make webstore portal API call to make path searchable
    webstore_api_url = config["webstore_api_url"]
    json_payload = {'path': resultdir}
    response = requests.post(webstore_api_url, json=json_payload)
    if response.status_code != 200:
        raise WebstoreError(f'Webstore api call returned a non 200 return: {response.text}')


def analysis_main(args, output, runnormal=False, normalname=False, normalfastqs=False, runtumor=False, tumorname=False, tumorfastqs=False, hg38ref=False, starttype=False):
    try:
        ################################################################
        # Write InputArgs to logfile
        config = read_wrapperconf()
        commandlogs = config["commandlogs"]
        #if not os.path.exists(commandlogs):
        #    os.makedirs(commandlogs)
        command = f"{sys.argv[0]}"
        current_date = time.strftime("%Y-%m-%d")
        commandlog = f"{commandlogs}/commands_{current_date}.log"
        for arg in vars(args):
            command = f"{command} --{arg} {getattr(args, arg)}"
        commandlogfile = open(commandlog, "a+")
        commandlogfile.write(f"{get_time()}" + "\n")
        commandlogfile.write(command + "\n")
        ################################################################

        if output.endswith("/"):
            output = output[:-1]
        if normalfastqs:
            if normalfastqs.endswith("/"):
                normalfastqs = normalfastqs[:-1]
        if tumorfastqs:
            if tumorfastqs.endswith("/"):
                tumorfastqs = tumorfastqs[:-1]


        #################################################################
        # Validate Inputs
        ################################################################
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
        mainconf_path = f"{configdir}/{mainconf_name}"

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
                        logger(f"Warning: No fastqs found in normaldir")
                        f_normalfastqs = glob.glob(f"{normalfastqs}/*{normalname}*fasterq")
                        if not f_normalfastqs:
                            error_list.append(f"No fastqs or fasterqs found in normaldir")

            if tumorfastqs:
                if not os.path.isdir(tumorfastqs):
                    error_list.append(f"{tumorfastqs} does not appear to be a directory")
                else:
                    f_tumorfastqs = glob.glob(f"{tumorfastqs}/*{tumorname}*fastq.gz")
                    if not f_tumorfastqs:
                        logger(f"Warning: No fastqs found in tumordir")
                        f_tumorfastqs = glob.glob(f"{tumorfastqs}/*{tumorname}*fasterq")
                        if not f_tumorfastqs:
                            error_list.append(f"No fastqs or fasterqs found in tumordir")
       
        # prepare outputdirectory
        if not os.path.isdir(output):
            try:
                os.mkdir(output)
            except Exception as e:
                error_list.append(f"outputdirectory: {output} does not exist and could not be created")

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
        normalid, tumorid = get_normalid_tumorid(runnormal, normalname, runtumor, tumorname)
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
        analysisdict["workingdir"] = output
        #insilico
        analysisdict["insilico"] = config["insilicopanels"]

        basename_output = os.path.basename(output)

        if hg38ref == "yes":
            analysisdict["reference"] = "hg38"
            if tumorname:
                if normalname:
                    analysisdict["resultdir"] =  f'{config["resultdir_hg38"]}/{basename_output}' #Use f'{config["testresultdir"]}/{tumorname}'for testing
                else:
                    analysisdict["resultdir"] = f'{config["resultdir_hg38"]}/tumor_only/{basename_output}'
            else:
                analysisdict["resultdir"] = f'{config["resultdir_hg38"]}/normal_only/{basename_output}'
        else:
            analysisdict["reference"] = "hg19"
            if tumorname:
                if normalname:
                    analysisdict["resultdir"] = f'{config["resultdir_hg19"]}/{basename_output}'
                else:
                    analysisdict["resultdir"] = f'{config["resultdir_hg19"]}/tumor_only/{basename_output}'
            else:
                analysisdict["resultdir"] = f'{config["resultdir_hg19"]}/normal_only/{basename_output}'
        if tumorname:
            with open(f"{runconfigs}/{tumorid}_config.json", 'w') as analysisconf:
                json.dump(analysisdict, analysisconf, ensure_ascii=False, indent=4)
        else:
            with open(f"{runconfigs}/{normalid}_config.json", 'w') as analysisconf:
                json.dump(analysisdict, analysisconf, ensure_ascii=False, indent=4)
        ###################################################################
        # Prepare Singularity Binddirs
        binddirs = config["singularitybinddirs"]
        binddir_string = ""
        for binddir in binddirs:
            source = binddirs[binddir]["source"]
            if not analysisdict["reference"] in source:
                if not "petagene" in source:
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
        binddir_string = f"{binddir_string}{output}"
        print(binddir_string)


        ###################################################################
        # Start SnakeMake pipeline
        ###################################################################
        scriptdir = os.path.dirname(os.path.realpath(__file__)) # find current dir

        # Generate random hash for shadow directory
        letters = string.ascii_lowercase
        letters = ''.join(random.choice(letters) for i in range(10))
        shadow_dir = os.path.join('/tmp', letters)

        snakemake_path = config["snakemake_env"]
        os.environ["PATH"] += os.pathsep + snakemake_path
        my_env = os.environ.copy() 
        if tumorname:
            snakemake_args = f"snakemake -s pipeline.snakefile --configfile {runconfigs}/{tumorid}_config.json --dag | dot -Tsvg > {samplelogs}/dag_{current_date}.svg"
        else:
            snakemake_args = f"snakemake -s pipeline.snakefile --configfile {runconfigs}/{normalid}_config.json --dag | dot -Tsvg > {samplelogs}/dag_{current_date}.svg"
        # >>>>>>>>>>>> Create Dag of pipeline
        subprocess.run(snakemake_args, shell=True, env=my_env) # CREATE DAG
        if tumorname:
            snakemake_args = f"snakemake -s pipeline.snakefile --configfile {runconfigs}/{tumorid}_config.json --use-singularity --singularity-args '-e --bind {binddir_string} --bind /medstore --bind /seqstore --bind /apps' --cluster-config configs/cluster.yaml --cluster \"qsub -S /bin/bash -pe mpi {{cluster.threads}} -q {{cluster.queue}} -N {{cluster.name}} -o {samplelogs}/{{cluster.output}} -e {samplelogs}/{{cluster.error}} -l {{cluster.excl}}\" --jobs 999 --latency-wait 60 --directory {output} --shadow-prefix {shadow_dir} &>> {samplelog}"
        else:
            snakemake_args = f"snakemake -s pipeline.snakefile --configfile {runconfigs}/{normalid}_config.json --use-singularity --singularity-args '-e --bind {binddir_string} --bind /medstore --bind /seqstore --bind /apps' --cluster-config configs/cluster.yaml --cluster \"qsub -S /bin/bash -pe mpi {{cluster.threads}} -q {{cluster.queue}} -N {{cluster.name}} -o {samplelogs}/{{cluster.output}} -e {samplelogs}/{{cluster.error}} -l {{cluster.excl}}\" --jobs 999 --latency-wait 60 --directory {output} --shadow-prefix {shadow_dir} &>> {samplelog}"
        # >>>>>>>>>>>> Start pipeline
        subprocess.run(snakemake_args, shell=True, env=my_env) # Shellscript pipeline

    except Exception as e:
        tb = traceback.format_exc()
        logger(f"Error in script:")
        logger(f"{e} Traceback: {tb}")
        sys.exit(1)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-rn', '--runnormal', nargs='?', help='the sequencing run the normalsample was sequenced in', required=False)
    parser.add_argument('-o', '--outputdir', nargs='?', help='output directory, where to put results', required=True)
    parser.add_argument('-ns', '--normalsample', nargs='?', help='normal samplename', required=False)
    parser.add_argument('-nf', '--normalfastqs', nargs='?', help='path to directory containing normal fastqs', required=False)
    parser.add_argument('-rt', '--runtumor', nargs='?', help='the sequencing run the tumorsample was sequenced in',     required=False)
    parser.add_argument('-tn', '--tumorsample', nargs='?', help='tumor samplename', required=False)
    parser.add_argument('-tf', '--tumorfastqs', nargs='?', help='path to directory containing tumor fastqs', required=False)
    parser.add_argument('-hg38', '--hg38ref', nargs='?', help='run analysis on hg38 reference (write yes if you want this option)', required=False)
    parser.add_argument('-stype', '--starttype', nargs='?', help='write forcestart if you want to ignore fastqs', required=False)
    parser.add_argument('-na', '--noalissa', action="store_true", help='Disables Alissa upload', required=False)
    parser.add_argument('-cr', '--copyresults', action="store_true", help='Copy results to resultdir on seqstore', required=False)
    args = parser.parse_args()
    timestamp = get_timestamp()
    args.outputdir = f'{args.outputdir}_{timestamp}'
    analysis_main(args, args.outputdir, args.runnormal, args.normalsample, args.normalfastqs, args.runtumor, args.tumorsample, args.tumorfastqs, args.hg38ref, args.starttype)

    if os.path.isfile(f"{args.outputdir}/reporting/workflow_finished.txt"):
        if args.tumorsample:
            if args.normalsample:
            # these functions are only executed if snakemake workflow has finished successfully
                yearly_stats(args.tumorsample, args.normalsample)
                if args.copyresults:
                    copy_results(args.outputdir, args.runnormal, args.normalsample, args.runtumor, args.tumorsample)
                if args.hg38ref and not args.noalissa:
                    alissa_upload(args.outputdir, args.normalsample, args.runnormal, args.hg38ref)
            else:
                yearly_stats(args.tumorsample, 'None')
                if args.copyresults:
                    copy_results(args.outputdir, runtumor=args.runtumor, tumorname=args.tumorsample)
        else:
            yearly_stats('None', args.normalsample)
            if args.copyresults:
                copy_results(args.outputdir, runnormal=args.runnormal, normalname=args.normalsample)
            if args.hg38ref and not args.noalissa:
                alissa_upload(args.outputdir, args.normalsample, args.runnormal, args.hg38ref)
