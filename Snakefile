# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
from os.path import join
import glob
import time
from pathlib import Path
import yaml
from tools.helpers import read_config
import os
from definitions import ROOT_DIR

normalfastqdirs = config["normalfastqs"]
normalname = config["normalname"]
normalid = config["normalid"]

tumorfastqdirs = config["tumorfastqs"]
tumorname = config["tumorname"]
tumorid = config["tumorid"]

reference = config["reference"]

# It uses the following configs from the working directory
pipeconfig = read_config(config["pipeconfig"])  # In launch_snakemake.py the pipeconfig is adjusted to the genome (hg19/hg38)
clusterconf = read_config(config["clusterconfig"])
filterconfig = read_config(config["filterconfig"])
resultsconf = read_config(config["resultfilesconf"])

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
if normalfastqdirs:
    for normalfastqdir in normalfastqdirs:
        for fwdpattern in fwdpatterns:
            normal_fwd_fastqs = glob.glob(f"{normalfastqdir}/{normalname}*{fwdpattern}")
            if normal_fwd_fastqs:
                for normal_fwd_fastq in normal_fwd_fastqs:
                    fastqpair_pattern = os.path.basename(normal_fwd_fastq).replace(fwdpattern, "")
                    fastq_dict["normal"]["fastqpair_patterns"][fastqpair_pattern] = {}
                    fastq_dict["normal"]["fastqpair_patterns"][fastqpair_pattern]["fwd"] = normal_fwd_fastq
        for revpattern in revpatterns:
            normal_rev_fastqs = glob.glob(f"{normalfastqdir}/{normalname}*{revpattern}")
            if normal_rev_fastqs:
                for normal_rev_fastq in normal_rev_fastqs:
                    fastqpair_pattern = os.path.basename(normal_rev_fastq).replace(revpattern, "")
                    fastq_dict["normal"]["fastqpair_patterns"][fastqpair_pattern]["rev"] = normal_rev_fastq

# Prepare Tumor Fastq Variables
if tumorfastqdirs:
    for tumorfastqdir in tumorfastqdirs:
        for fwdpattern in fwdpatterns:
            tumor_fwd_fastqs = glob.glob(f"{tumorfastqdir}/{tumorname}*{fwdpattern}")
            if tumor_fwd_fastqs:
                for tumor_fwd_fastq in tumor_fwd_fastqs:
                    fastqpair_pattern = os.path.basename(tumor_fwd_fastq).replace(fwdpattern, "")
                    fastq_dict["tumor"]["fastqpair_patterns"][fastqpair_pattern] = {}
                    fastq_dict["tumor"]["fastqpair_patterns"][fastqpair_pattern]["fwd"] = tumor_fwd_fastq
        for revpattern in revpatterns:
            tumor_rev_fastqs = glob.glob(f"{tumorfastqdir}/{tumorname}*{revpattern}")
            if tumor_rev_fastqs:
                for tumor_rev_fastq in tumor_rev_fastqs:
                    fastqpair_pattern = os.path.basename(tumor_rev_fastq).replace(revpattern, "")
                    fastq_dict["tumor"]["fastqpair_patterns"][fastqpair_pattern]["rev"] = tumor_rev_fastq 

#wildcard_constraints:
#    sname="[^_]*_[^_]*_[^_]*"

###########################################################
# Defining Non Cluster Rules
if tumorid:
    if normalid:
        # Runs tn_workflow / paired if tumorid and normalid
        localrules: all, tn_workflow, share_to_resultdir, excel_qc, tmb_calculation
    else:
        # Runs tumoronly_workflow if tumorid but not normalid
        localrules: all, tumoronly_workflow, share_to_resultdir, excel_qc, tmb_calculation
else: 
    # Runs normalonly_workflow if normalid but not tumorid
    localrules: all, normalonly_workflow, share_to_resultdir, excel_qc
###########################################################

#########################################
# Mapping
include:    "workflows/rules/mapping/mapping.smk"
include:    "workflows/rules/mapping/cram.smk"
include:    "workflows/rules/mapping/generate_tdf.smk"

#########################################
# VariantCalling
if tumorid:
    include:        "workflows/rules/variantcalling/tnscope.smk"
    include:        "workflows/rules/variantcalling/pindel.smk"
    include:        "workflows/rules/variantcalling/control-freec.smk"
    include:        "workflows/rules/variantcalling/ascat.smk"
include:        "workflows/rules/variantcalling/dnascope.smk"
include:        "workflows/rules/variantcalling/canvas.smk"
include:        "workflows/rules/variantcalling/manta.smk"

#########################################
# Small Tools
if tumorid:
    include:        "workflows/rules/small_tools/tmb_calculation.smk"
    include:        "workflows/rules/small_tools/msi.smk"
include:        "workflows/rules/small_tools/ballele.smk"
include:        "workflows/rules/small_tools/bgzip.smk"

#########################################
# QC
include:       "workflows/rules/qc/aggregate_qc.smk"
include:       "workflows/rules/qc/coverage.smk"


ruleorder: merge_snvs_cnvs > dnascope_vcffilter
ruleorder: canvas_germline > bgzip_vcf

if tumorid and normalid:
    ruleorder: canvas_somatic > bgzip_vcf


all_result_files = []
for result in resultsconf.values():
    all_result_files.extend(result)

rule all:
    input: 
        all_result_files
