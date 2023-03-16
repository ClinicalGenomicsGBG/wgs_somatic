# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule bgzip_vcf:
    input:
        "{path}/{file}.vcf"
    params:
        bgzip = pipeconfig["rules"]["bgzip"]["bgzip"],
    output:
        temp("{path}/{file}.vcf.gz")
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
        temp("{path}/{file}.vcf.gz.csi")
    wildcard_constraints:
        file = r"[^/]+"
    shadow:
        pipeconfig["rules"].get("bgzip", {}).get("shadow", pipeconfig.get("shadow", False))    
    shell:
        "{params.bcftools} index {input}"