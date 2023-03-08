# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os

def get_bedfile_path(wcs):
    return config["insilico"][f"{wcs.insiliconame}"]['bed']

def get_bedfile_version(wcs):
    return config["insilico"][f"{wcs.insiliconame}"]['version']

def get_insiliconame(wcs):
    return config["insilico"]


# If normal sample is available, we will get the insilico coverage values from there
if normalid:
    rule insilico_coverage:
        input:
            bamfile = "{normalid}/dedup/{sname}_DEDUP.bam"
            #bedinfo = get_insilico_info
        params:
            insilico_wrapper = pipeconfig["rules"]["insilico"].get("insilico_main", f"{ROOT_DIR}/workflows/scripts/insilico_coverage/insilico_main.py"),
            samtools = pipeconfig["rules"]["insilico"]["samtools"],
            pythonversion = pipeconfig["rules"]["insilico"]["python"],
            bedfile_path = get_bedfile_path,
            bedfile_version = get_bedfile_version
        output:
            "{normalid}/insilico/{insiliconame}/{sname}_{insiliconame}_10x.xlsx",
            "{normalid}/insilico/{insiliconame}/{sname}_{insiliconame}_20x.xlsx",
            "{normalid}/insilico/{insiliconame}/{sname}_{insiliconame}_genes_below10x.xlsx",
            "{normalid}/insilico/{insiliconame}/{sname}_{insiliconame}.csv",
            "{normalid}/insilico/{insiliconame}/{sname}_{insiliconame}_cov.tsv",
        run:
            os.makedirs(f'{wildcards.normalid}/insilico', exist_ok=True)
            os.makedirs(f'{wildcards.normalid}/insilico/{wildcards.insiliconame}', exist_ok=True)

            insilico_level = config["insilico"][f"{wildcards.insiliconame}"]["levels"]
            shell(f"{params.pythonversion} {params.insilico_wrapper} --bamfile {input.bamfile} --bedfile {params.bedfile_path} --version {params.bedfile_version} --outputdir {wildcards.normalid}/insilico/{wildcards.insiliconame} --annotationlevel {insilico_level} --samplename {wildcards.sname} --samtools {params.samtools}")

# Otherwise (only when running tumour only) we will try go get coverage values from the tumour sample
else:
     rule insilico_coverage:
        input:
            bamfile = "{tumorid}/dedup/{sname}_DEDUP.bam"
            #bedinfo = get_insilico_info
        params:
            insilico_wrapper = pipeconfig["rules"]["insilico"].get("insilico_main", f"{ROOT_DIR}/workflows/scripts/insilico_coverage/insilico_main.py"),
            samtools = pipeconfig["rules"]["insilico"]["samtools"],
            pythonversion = pipeconfig["rules"]["insilico"]["python"],
            bedfile_path = get_bedfile_path,
            bedfile_version = get_bedfile_version
        output:
            "{tumorid}/insilico/{insiliconame}/{sname}_{insiliconame}_10x.xlsx",
            "{tumorid}/insilico/{insiliconame}/{sname}_{insiliconame}_20x.xlsx",
            "{tumorid}/insilico/{insiliconame}/{sname}_{insiliconame}_genes_below10x.xlsx",
            "{tumorid}/insilico/{insiliconame}/{sname}_{insiliconame}.csv",
            "{tumorid}/insilico/{insiliconame}/{sname}_{insiliconame}_cov.tsv",
        run:
            os.makedirs(f'{wildcards.tumorid}/insilico', exist_ok=True)
            os.makedirs(f'{wildcards.tumorid}/insilico/{wildcards.insiliconame}', exist_ok=True)

            insilico_level = config["insilico"][f"{wildcards.insiliconame}"]["levels"]
            shell(f"{params.pythonversion} {params.insilico_wrapper} --bamfile {input.bamfile} --bedfile {params.bedfile_path} --version {params.bedfile_version} --outputdir {wildcards.tumorid}/insilico/{wildcards.insiliconame} --annotationlevel {insilico_level} --samplename {wildcards.sname} --samtools {params.samtools}")


