import sys
sys.path.append("..") 
from os.path import join
import glob
import time
from pathlib import Path
import yaml
from helpers import read_config
import os
# sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# from definitions import ROOT_DIR
# ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT_DIR = os.path.join(workflow.basedir, "..")
config_path = os.path.join(ROOT_DIR, 'configs', 'launcher_config.json')
launcher_config = read_config(config_path)
# print(launcher_config)
pipeconfig = read_config(os.path.join(ROOT_DIR, "configs",launcher_config["hg38conf"]))  # In launch_snakemake.py the pipeconfig is adjusted to the genome (hg19/hg38)
clusterconf = read_config(os.path.join(ROOT_DIR, "configs",launcher_config["clusterconf"]))
filterconfig = read_config(os.path.join(ROOT_DIR, "configs",launcher_config["filterconf"]))

dir_with_bams = config.get("dir_to_process")

sname, = glob_wildcards(f"{dir_with_bams}/{{sname}}.bam") # variable is called sname and not something more descriptive so that the default settings in cluster.yaml apply
rule all:
    input:
        expand("{sname}.cram", sname=sname),
        expand("{sname}.cram.crai", sname=sname),

rule cram_crai:
    input:
        bam = f"{dir_with_bams}/{{sname}}.bam",
        bai = f"{dir_with_bams}/{{sname}}.bam.bai"
    singularity:
        pipeconfig["singularities"]["samtools"]["sing"]
    params:
        threads = clusterconf["cram"]["threads"],
        referencegenome = pipeconfig["referencegenome"]
    output:
        cram = "{sname}.cram",
        crai = "{sname}.cram.crai"
    shell:
        "samtools view -C --threads {params.threads} -T {params.referencegenome} -o {output.cram} {input.bam} ; "
        "samtools index {output.cram}"
