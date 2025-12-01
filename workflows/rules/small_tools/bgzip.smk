# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule bgzip_vcf:
    input:
        "{path}/{file}.vcf"
    params:
        bgzip = pipeconfig["rules"]["bgzip"]["bgzip"],
        vstamp = f"{VDIR}/bgzip_vcf.txt"
    output:
        "{path}/{file}.vcf.gz",
    wildcard_constraints:
        file = r"[^/]+"
    shadow:
        pipeconfig["rules"].get("bgzip_vcf", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        """
        {params.bgzip} --version | head -n 1 > {params.vstamp}
        {params.bgzip} -c {input} > {output}
        """

rule index_vcf:
    input:
        "{path}/{file}.vcf.gz"
    params:
        bcftools = pipeconfig["rules"]["bgzip"]["bcftools"],
        vstamp = f"{VDIR}/index_vcf.txt"
    output:
        "{path}/{file}.vcf.gz.csi",
    wildcard_constraints:
        file = r"[^/]+"
    shadow:
        pipeconfig["rules"].get("bgzip", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        """
        {params.bcftools} --version | head -n 2 > {params.vstamp}
        {params.bcftools} index {input}
        """
