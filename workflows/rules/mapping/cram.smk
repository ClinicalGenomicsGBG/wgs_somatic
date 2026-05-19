# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule cram:
    input:
        bam = "{path}/{file}.bam",
        bai = "{path}/{file}.bam.bai"
    singularity:
        pipeconfig["singularities"]["samtools"]["sing"]
    params:
        referencegenome = pipeconfig["referencegenome"],
        vstamp = f"{VDIR}/cram.txt"
    threads:
        clusterconf["cram"]["threads"]
    shadow:
        pipeconfig["rules"].get("cram", {}).get("shadow", pipeconfig.get("shadow", False))
    output:
        cram = "{path}/{file}.cram",
        crai = "{path}/{file}.cram.crai",
    shell:
        """
        samtools --version | head -n 2 > {params.vstamp}
        samtools view -C -@ {threads} -T {params.referencegenome} -o {output.cram} {input.bam}
        samtools index -@ {threads} {output.cram}
        """
