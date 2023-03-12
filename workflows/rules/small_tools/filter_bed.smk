# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os

rule filter_variants_in_bed:
    input:
        "{stype}/{caller}/{sname}_{vcftype}.vcf"    
    params:
        rtg = pipeconfig["rules"]["filter_variants_in_bed"]["rtg"],
        bedfile = pipeconfig["rules"]["filter_variants_in_bed"]["bedfile"],
    output:
        "{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.vcf"
    shell:
        "{params.rtg} vcffilter "
            "--include-bed={params.bedfile} "
            "--output={output} "
            "--input={input}; "
        "gunzip {output}.gz"