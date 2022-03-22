# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os

def get_bedfile_path(wcs):
    return config["insilico"][f"{wcs.insiliconame}"]['bedfile']

def get_bedfile_version(wcs):
    return config["insilico"][f"{wcs.insiliconame}"]['version']

rule insilico_coverage:
    input:
        bamfile = "{workingdir}/dedup/{sampleid}_DEDUP.bam"
        #bedinfo = get_insilico_info
    params:
        insilico_wrapper = pipeconfig["rules"]["insilico"]["insilico_main"],
        samtools = pipeconfig["rules"]["insilico"]["samtools"],
        pythonversion = pipeconfig["rules"]["insilico"]["python"],
        bedfile_path = get_bedfile_path,
        bedfile_version = get_bedfile_version
    output:
        "{workingdir}/insilico/{insiliconame}/{sampleid}_{insiliconame}_10x.xlsx",
        "{workingdir}/insilico/{insiliconame}/{sampleid}_{insiliconame}_20x.xlsx",
        "{workingdir}/insilico/{insiliconame}/{sampleid}_{insiliconame}_genes_below10x.xlsx",
        "{workingdir}/insilico/{insiliconame}/{sampleid}_{insiliconame}.csv",
        "{workingdir}/insilico/{insiliconame}/{sampleid}_{insiliconame}_cov.tsv"
    run:
        os.makedirs(f'{wildcards.workingdir}/insilico', exist_ok=True)
        os.makedirs(f'{wildcards.workingdir}/insilico/{wildcards.insiliconame}', exist_ok=True)

        insilico_level = config["insilico"][f"{wildcards.insiliconame}"]["levels"]
        shell(f"{params.pythonversion} {params.insilico_wrapper} --bamfile {input.bamfile} --bedfile {params.bedfile_path} --version {params.bedfile_version} --outputdir {wildcards.workingdir}/insilico/{wildcards.insiliconame} --annotationlevel {insilico_level} --samplename {wildcards.sampleid} --samtools {params.samtools}")
