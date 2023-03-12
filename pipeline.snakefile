# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
from os.path import join
import glob
import time
from pathlib import Path
import yaml
import helpers
import os
from definitions import ROOT_DIR

__author__ = "Rickard 'Ricksy' Rickardsson"

normalfastqdirs = config["normalfastqs"]
normalname = config["normalname"]
normalid = config["normalid"]

tumorfastqdirs = config["tumorfastqs"]
tumorname = config["tumorname"]
tumorid = config["tumorid"]

reference = config["reference"]
workingdir = config["workingdir"]
insilico_panels = config["insilico"]

##################################################
# Chose Config based on Reference
# ---------------------------------------------

if reference == "hg38":
    configfilepath = f"{ROOT_DIR}/configs/config_hg38.json"
else:
    configfilepath = f"{ROOT_DIR}/configs/config_hg19.json"



pipeconfig = helpers.read_config(configfilepath)
clusterconf = helpers.read_clusterconf()

shell.executable("/bin/bash")

####################################################
# Prepare Fastq Variables 
# -------------------------------------------------

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
                    fastq_dict["normal"]["fastqpair_patterns"][fastqpair_pattern] = {"fwd": normal_fwd_fastq}
        for revpattern in revpatterns:
            normal_rev_fastqs = glob.glob(f"{normalfastqdir}/{normalname}*{revpattern}")
            if normal_rev_fastqs:
                for normal_rev_fastq in normal_rev_fastqs:
                    fastqpair_pattern = os.path.basename(normal_rev_fastq).replace(revpattern, "")
                    fastq_dict["normal"]["fastqpair_patterns"][fastqpair_pattern] |= {"rev": normal_rev_fastq}

# Prepare Tumor Fastq Variables
if tumorfastqdirs:
    for tumorfastqdir in tumorfastqdirs:
        for fwdpattern in fwdpatterns:
            tumor_fwd_fastqs = glob.glob(f"{tumorfastqdir}/{tumorname}*{fwdpattern}")
            if tumor_fwd_fastqs:
                for tumor_fwd_fastq in tumor_fwd_fastqs:
                    fastqpair_pattern = os.path.basename(tumor_fwd_fastq).replace(fwdpattern, "")
                    fastq_dict["tumor"]["fastqpair_patterns"][fastqpair_pattern] = {"fwd": tumor_fwd_fastq}
        for revpattern in revpatterns:
            tumor_rev_fastqs = glob.glob(f"{tumorfastqdir}/{tumorname}*{revpattern}")
            if tumor_rev_fastqs:
                for tumor_rev_fastq in tumor_rev_fastqs:
                    fastqpair_pattern = os.path.basename(tumor_rev_fastq).replace(revpattern, "")
                    fastq_dict["tumor"]["fastqpair_patterns"][fastqpair_pattern] |= {"rev": tumor_rev_fastq} 
# -------------------------------------------------

####################################################
# Include Rules 
# -------------------------------------------------

# Mapping
include: "workflows/rules/mapping/mapping.smk"
include: "workflows/rules/mapping/generate_tdf.smk"
# VariantCalling
include: "workflows/rules/variantcalling/manta.smk"
include: "workflows/rules/variantcalling/tnscope.smk"
include: "workflows/rules/variantcalling/pindel.smk"
include: "workflows/rules/variantcalling/dnascope.smk"
include: "workflows/rules/small_tools/ballele.smk"
include: "workflows/rules/variantcalling/canvas.smk"
include: "workflows/rules/small_tools/bgzip.smk"
include: "workflows/rules/small_tools/filter_bed.smk"
# QC
include: "workflows/rules/qc/aggregate_qc.smk"
include: "workflows/rules/qc/insilico_coverage.smk"
include: "workflows/rules/qc/coverage.smk"
# ResultSharing:
include: "workflows/rules/results_sharing/share_to_resultdir.smk"

# Non-cluster rules
localrules: all, share_to_resultdir, excel_qc

# Priority for ambiguous rules
ruleorder: merge_snvs_cnvs > dnascope_vcffilter
ruleorder: canvas > bgzip_vcf

wildcard_constraints:
    vartype=r"(germline)|(somatic)"


####################################################
# Result definitions
# -------------------------------------------------
def get_result_input(wildcards):
    inputs = [
        # Mapping (Tumor > Normal)
        f"{'tumor' if tumorid else 'normal'}/reports/{tumorid or normalid}_baf.igv",
        # Canvas Germline (Normal > Tumor)
        f"{'normal' if normalid else 'tumor'}/canvas/{normalid or tumorid}_CNV_germline.vcf.xlsx",
        *multiext(f"{'normal' if normalid else 'tumor'}/canvas/{normalid or tumorid}_CNV_germline", ".vcf.gz", ".vcf.gz.csi"),
        f"{'normal' if normalid else 'tumor'}/canvas/{normalid or tumorid}_germline_CNV_observed.seg",
        f"{'normal' if normalid else 'tumor'}/canvas/{normalid or tumorid}_germline_CNV_called.seg",
        # Insilico (Tumor > Normal)
        f"qc_report/{tumorname or normalname}_qc_stats.xlsx",
    ]

    if normalid:
        inputs += [
            # Mapping (Normal)
            *multiext(f"normal/realign/{normalid}_REALIGNED", ".bam", ".bam.bai"),
            f"normal/reports/{normalid}_REALIGNED.bam.tdf",
            # DNAscope Germline (Normal)
            *multiext(f"normal/dnascope/{normalid}_germline_refseq3kfilt", ".vcf.gz", ".vcf.gz.csi"),
            *multiext(f"normal/dnascope/{normalid}_{reference}_SNV_CNV_germline", ".vcf.gz", ".vcf.gz.csi"),
            # Manta Germline (Normal)
            *multiext(f"normal/manta/{normalid}_germline_MantaBNDs", ".vcf.gz", ".vcf.gz.csi"),
            *multiext(f"normal/manta/{normalid}_germline_MantaNOBNDs", ".vcf.gz", ".vcf.gz.csi"),
        ]

    if tumorid:
        inputs += [
            # Mapping (Tumor)
            *multiext(f"tumor/realign/{tumorid}_REALIGNED", ".bam", ".bam.bai"),
            f"tumor/reports/{tumorid}_REALIGNED.bam.tdf",
            # TNscope Somatic (Tumor)
            *multiext(f"tumor/tnscope/{tumorid}_somatic_refseq3kfilt", ".vcf.gz", ".vcf.gz.csi"),
            # Manta Somatic (Tumor)
            f"tumor/manta/{tumorid}_somatic_mantaSV.vcf.xlsx",
            f"tumor/manta/{tumorid}_somatic_mantaSV_Summary.xlsx",
            *multiext(f"tumor/manta/{tumorid}_somatic_MantaBNDs", ".vcf.gz", ".vcf.gz.csi"),
            *multiext(f"tumor/manta/{tumorid}_somatic_MantaNOBNDs", ".vcf.gz", ".vcf.gz.csi"),
        ]

    if tumorid and reference == "hg38":
        # Pindel Somatic (Tumor, hg38 only)
        inputs += [f"tumor/pindel/{tumorid}_pindel.xlsx"]

    if tumorid and normalid:
        inputs += [
            # Canvas (Tumor + Normal)
            *multiext(f"tumor/canvas/{tumorid}_CNV_somatic", ".vcf.gz", ".vcf.gz.csi"),
            f"tumor/canvas/{tumorid}_somatic_CNV_observed.seg",
            f"tumor/canvas/{tumorid}_somatic_CNV_called.seg",
        ]

    return inputs

rule all:
    input:
        get_result_input,
        "reporting/workflow_finished.txt"
