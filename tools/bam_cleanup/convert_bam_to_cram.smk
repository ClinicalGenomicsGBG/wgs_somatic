import sys
# print(sys.path)
# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from os.path import join
import glob
import time
from pathlib import Path
import yaml
ROOT_DIR = os.path.join(workflow.basedir, "..")
sys.path.append(ROOT_DIR)
from helpers import read_config
import os

# config_path = os.path.join(ROOT_DIR, 'configs', 'launcher_config.json')
launcher_config = read_config(config.get("launcher_config_path"))
launcher_config_parentdir = os.path.dirname(config.get("launcher_config_path"))

pipeconfig = read_config(os.path.join(launcher_config_parentdir, launcher_config["hg38conf"]))  # In launch_snakemake.py the pipeconfig is adjusted to the genome (hg19/hg38)
# clusterconf = read_config(os.path.join(launcher_config_parentdir, launcher_config["clusterconf"]))

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
        threads = config["cram"]["threads"], # config is passed as snakemake argument
        referencegenome = pipeconfig["referencegenome"]
    output:
        cram = "{sname}.cram",
        crai = "{sname}.cram.crai"
    shell:
        "samtools view -C --threads {params.threads} -T {params.referencegenome} -o {output.cram} {input.bam} ; "
        "samtools index {output.cram}"
