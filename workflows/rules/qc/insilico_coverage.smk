# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

def get_bedfile_path(wcs):
    return config["insilico"][f"{wcs.insiliconame}"]["bed"]


def get_bedfile_version(wcs):
    return config["insilico"][f"{wcs.insiliconame}"]["version"]


def get_insiliconame(wcs):
    return config["insilico"]


rule insilico_coverage:
    input:
        bam="{stype}/dedup/{sname}_DEDUP.bam",
        bai="{stype}/dedup/{sname}_DEDUP.bam.bai",
    params:
        insilico_level = lambda wildcards: config["insilico"][wildcards.insiliconame]["levels"],
        insilico_wrapper = pipeconfig["rules"]["insilico"].get("insilico_main", f"{ROOT_DIR}/workflows/scripts/insilico_coverage/insilico_main.py",), 
        pythonversion = pipeconfig["rules"]["insilico"]["python"],    
        samtools = pipeconfig["rules"]["insilico"]["samtools"],
        bedfile_path = get_bedfile_path,
        bedfile_version = get_bedfile_version,
    output:
        temp("{stype}/insilico/{insiliconame}/{sname}_{insiliconame}_10x.xlsx"),
        temp("{stype}/insilico/{insiliconame}/{sname}_{insiliconame}_20x.xlsx"),
        temp("{stype}/insilico/{insiliconame}/{sname}_{insiliconame}_genes_below10x.xlsx"),
        temp("{stype}/insilico/{insiliconame}/{sname}_{insiliconame}.csv"),
        temp("{stype}/insilico/{insiliconame}/{sname}_{insiliconame}_cov.tsv"),
        touch("{stype}/insilico/{insiliconame}/{sname}.complete"),
    shadow:
        pipeconfig["rules"].get("insilico", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "{params.pythonversion} {params.insilico_wrapper} "
            "--bamfile {input.bam} "
            "--bedfile {params.bedfile_path} "
            "--version {params.bedfile_version} "
            "--outputdir {wildcards.stype}/insilico/{wildcards.insiliconame} "
            "--annotationlevel {params.insilico_level} "
            "--samplename {wildcards.sname} "
            "--samtools {params.samtools}"
