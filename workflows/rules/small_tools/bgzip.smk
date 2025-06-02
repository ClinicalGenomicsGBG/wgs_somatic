# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
from tools.helpers import conditional_temp

rule bgzip_vcf:
    input:
        "{path}/{file}.vcf"
    params:
        bgzip = pipeconfig["rules"]["bgzip"]["bgzip"],
    output:
        conditional_temp("{path}/{file}.vcf.gz", keepfiles)
    wildcard_constraints:
        file = r"[^/]+"
    shadow:
        pipeconfig["rules"].get("bgzip_vcf", {}).get("shadow", pipeconfig.get("shadow", False))    
    shell:
        "{params.bgzip} -c {input} > {output}"

rule index_vcf:
    input:
        "{path}/{file}.vcf.gz"
    params:
        bcftools = pipeconfig["rules"]["bgzip"]["bcftools"],
    output:
        conditional_temp("{path}/{file}.vcf.gz.csi", keepfiles)
    wildcard_constraints:
        file = r"[^/]+"
    shadow:
        pipeconfig["rules"].get("bgzip", {}).get("shadow", pipeconfig.get("shadow", False))    
    shell:
        "{params.bcftools} index {input}"