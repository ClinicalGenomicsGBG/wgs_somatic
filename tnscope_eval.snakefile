# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
from os.path import join
import glob
import time
from pathlib import Path
import yaml
import helpers
import os
import json
from definitions import ROOT_DIR

__author__ = "Rickard 'Ricksy' Rickardsson"

normaldata = config["normalfastqs"]
normalfastqdirs = normaldata
normalname = config["sample"]["normalname"]
normalid = config["sample"]["normalid"]

tumordata = config["tumorfastqs"]
tumorfastqdirs = tumordata
tumorname = config["sample"]["tumorname"]
tumorid = config["sample"]["tumorid"]

rtg = config["rtg"]["tools"]
rtgsdf = config["rtg"]["sdf"]
bedfile = config["data"]["bed"]
truthset = config["data"]["tset"]
tnscopesetting = config["tnscope"]

tnscopesetting_list = []
for setting in tnscopesetting:
    tnscopesetting_list.append(setting)

##################################################
# Chose Config based on Reference
# ---------------------------------------------
reference = "hg38"
configfilepath = f"{ROOT_DIR}/configs/config_hg38.json"

#----------------------------------------------


pipeconfig = helpers.read_config(configfilepath)
clusterconf = helpers.read_config(f"{ROOT_DIR}/configs/cluster.yaml")

shell.executable("/bin/bash")

analysistime = time.strftime("%Y-%m-%d-%H-%M-%S")

sampleconfig = {}
sampleconfig[normalname] = {}
sampleconfig[normalname]["stype"] = "normal"
sampleconfig[tumorname] = {}
sampleconfig[tumorname]["stype"] = "tumor"
sampleconfig["normal"] = normalid
sampleconfig["tumor"] = tumorid
sampleconfig["normalname"] = normalname
sampleconfig["tumorname"] = tumorname

####################################################
# Prepare Fastq Variables 
# -------------------------------------------------

fwdpatterns = ["_1.fastq.gz", "_R1_001.fastq.gz", "_1.fasterq", "_R1_001.fasterq"] 
revpatterns = ["_2.fastq.gz", "_R2_001.fastq.gz", "_2.fasterq", "_R2_001.fasterq"]

fastq_dict = {}
fastq_dict["normal"] = {}
fastq_dict["normal"]["fastqpair_patterns"] = {}

fastq_dict["tumor"] = {}
fastq_dict["tumor"]["fastqpair_patterns"] = {}

# Prepare Normal Fastq Variables
for normalfastqdir in normalfastqdirs:
    for fwdpattern in fwdpatterns:
        normal_fwd_fastqs = glob.glob(f"{normalfastqdir}/*{fwdpattern}")
        if normal_fwd_fastqs:
            for normal_fwd_fastq in normal_fwd_fastqs:
                fastqpair_pattern = os.path.basename(normal_fwd_fastq).replace(fwdpattern, "")
                fastq_dict["normal"]["fastqpair_patterns"][fastqpair_pattern] = {}
                fastq_dict["normal"]["fastqpair_patterns"][fastqpair_pattern]["fwd"] = normal_fwd_fastq
    for revpattern in revpatterns:
        normal_rev_fastqs = glob.glob(f"{normalfastqdir}/*{revpattern}")
        if normal_rev_fastqs:
            for normal_rev_fastq in normal_rev_fastqs:
                fastqpair_pattern = os.path.basename(normal_rev_fastq).replace(revpattern, "")
                fastq_dict["normal"]["fastqpair_patterns"][fastqpair_pattern]["rev"] = normal_rev_fastq

# Prepare Tumor Fastq Variables
for tumorfastqdir in tumorfastqdirs:
    for fwdpattern in fwdpatterns:
        tumor_fwd_fastqs = glob.glob(f"{tumorfastqdir}/*{fwdpattern}")
        if tumor_fwd_fastqs:
            for tumor_fwd_fastq in tumor_fwd_fastqs:
                fastqpair_pattern = os.path.basename(tumor_fwd_fastq).replace(fwdpattern, "")
                fastq_dict["tumor"]["fastqpair_patterns"][fastqpair_pattern] = {}
                fastq_dict["tumor"]["fastqpair_patterns"][fastqpair_pattern]["fwd"] = tumor_fwd_fastq
    for revpattern in revpatterns:
        tumor_rev_fastqs = glob.glob(f"{tumorfastqdir}/*{revpattern}")
        if tumor_rev_fastqs:
            for tumor_rev_fastq in tumor_rev_fastqs:
                fastqpair_pattern = os.path.basename(tumor_rev_fastq).replace(revpattern, "")
                fastq_dict["tumor"]["fastqpair_patterns"][fastqpair_pattern]["rev"] = tumor_rev_fastq 

wildcard_constraints:
    sname="[^_]*_[^_]*_[^_]*"

###########################################################
# Defining Non Cluster Rules
localrules: all, validation_wf, filter_truthset, share_eval_to_igv
###########################################################

########################################
# Workflows
include:        "validation_scripts/validation_workflow.smk"

########################################
# Mapping

#########################################
# VariantCalling
include:        "validation_scripts/tnscope_eval.smk"

#########################################
# ResultSharing:


if reference == "hg38":
    ###########################################################
    # HG38 rules
    ###########################################################
    # Mapping
    include:        "workflows/rules/mapping/mapping_hg38.smk"
else:
    ###########################################################
    # HG19 rules
    ###########################################################
    # Mapping
    include:        "workflows/rules/mapping/mapping.smk"

rule all:
    input:
        "reporting/workflow_finished.txt"
