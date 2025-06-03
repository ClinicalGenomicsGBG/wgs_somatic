# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

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
        cram = "{stype}/realign/{sname}_REALIGNED.cram",
        crai = "{stype}/realign/{sname}_REALIGNED.cram.crai",
    shell:
        "samtools view -C --threads {params.threads} -T {params.referencegenome} -o {output.cram} {input.bam} ; "
        "samtools index {output.cram}"

# Convert cram to bam if the bam is removed or not present but necessary for downstream analysis.
rule cram_to_bam:
    input:
        cram = "{stype}/realign/{sname}_REALIGNED.cram",
        crai = "{stype}/realign/{sname}_REALIGNED.cram.crai"
    singularity:
        pipeconfig["singularities"]["samtools"]["sing"]
    params:
        threads = clusterconf["cram_to_bam"]["threads"],
        referencegenome = pipeconfig["referencegenome"]
    shadow:
        pipeconfig["rules"].get("cram_to_bam", {}).get("shadow", pipeconfig.get("shadow", False))
    output:
        bam = "{stype}/realign/{sname}_REALIGNED.bam",
        bai = "{stype}/realign/{sname}_REALIGNED.bam.bai",
    shell:
        "samtools view -b --threads {params.threads} -T {params.referencegenome} -o {output.bam} {input.cram} ; "
        "samtools index {output.bam}"