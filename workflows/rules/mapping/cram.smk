# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule cram:
    input:
        bam = "{stype}/realign/{sname}_REALIGNED.bam",
        bai = "{stype}/realign/{sname}_REALIGNED.bam.bai"
    singularity:
        pipeconfig["singularities"]["samtools"]["sing"]
    params:
        threads = clusterconf["mapping"]["threads"],
        referencegenome = pipeconfig["singularities"]["samtools"]["reference"]
    shadow:
        pipeconfig["rules"].get("cram", {}).get("shadow", pipeconfig.get("shadow", False))
    output:
        cram = temp("{stype}/realign/{sname}_REALIGNED.cram"),
        crai = temp("{stype}/realign/{sname}_REALIGNED.cram.crai")
    shell:
        "samtools view -C --threads 10 -T {params.referencegenome} -o {output.cram} {input.bam} ; "
        "samtools index {output.cram}"
