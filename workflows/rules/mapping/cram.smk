# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
from tools.helpers import conditional_temp

rule cram:
    input:
        bam = "{stype}/realign/{sname}_REALIGNED.bam",
        bai = "{stype}/realign/{sname}_REALIGNED.bam.bai"
    singularity:
        pipeconfig["singularities"]["samtools"]["sing"]
    params:
        threads = clusterconf["cram"]["threads"],
        referencegenome = pipeconfig["referencegenome"]
    shadow:
        pipeconfig["rules"].get("cram", {}).get("shadow", pipeconfig.get("shadow", False))
    output:
        cram = conditional_temp("{stype}/realign/{sname}_REALIGNED.cram", keepfiles),
        crai = conditional_temp("{stype}/realign/{sname}_REALIGNED.cram.crai", keepfiles)
    shell:
        "samtools view -C --threads {params.threads} -T {params.referencegenome} -o {output.cram} {input.bam} ; "
        "samtools index {output.cram}"
