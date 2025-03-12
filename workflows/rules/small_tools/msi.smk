# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule msi:
    input:
        tumor_bam = "tumor/realign/{sname}_REALIGNED.bam",
        tumor_bai = "tumor/realign/{sname}_REALIGNED.bam.bai",
        normal_bam = "normal/realign/{sname}_REALIGNED.bam",
        normal_bai = "normal/realign/{sname}_REALIGNED.bam.bai",
    params:
        msisensor = pipeconfig["rules"]["msi"]["msisensor-pro"],
        reference_list = pipeconfig["rules"]["msi"]["msi_list"],
        threads = clusterconf["msi"]["threads"],
    output:
        msi_out = temp("tumor/msi/{sname}_msi.txt"),
        msi_out_dis = temp("tumor/msi/{sname}_msi.txt_dis"),
        msi_out_germline = temp("tumor/msi/{sname}_msi.txt_germline"),
        msi_out_somatic = temp("tumor/msi/{sname}_msi.txt_somatic"),
    shell:
        """
        {params.msisensor} msi \
            -d {params.reference_list} \
            -n {input.normal_bam} \
            -t {input.tumor_bam} \
            -b {params.threads} \
            -z 1 \
            -o {output.msi_out}
        """

rule msi_filter_bam:
    input:
        bam = "{stype}/realign/{sname}_REALIGNED.bam",
    params:
        bed = pipeconfig["rules"]["msi"]["msi_bed"],
        bedtools = pipeconfig["rules"]["msi"]["bedtools"],
        samtools = pipeconfig["rules"]["msi"]["samtools"],
    output:
        filtered_bam = "{stype}/msi/{sname}_filtered.bam",
        filtered_bam_bai = "{stype}/msi/{sname}_filtered.bam.bai"
    shell:
        """
        {params.bedtools} intersect -abam {input.bam} -b {params.bed} > {output.filtered_bam}
        {params.samtools} index {output.filtered_bam}
        """

rule msi_reduced:
    input:
        tumor_bam = "tumor/msi/{sname}_filtered.bam",
        tumor_bai = "tumor/msi/{sname}_filtered.bam.bai",
        normal_bam = "normal/msi/{sname}_filtered.bam",
        normal_bai = "normal/msi/{sname}_filtered.bam.bai",
    params:
        msisensor = pipeconfig["rules"]["msi"]["msisensor-pro"],
        reference_list = pipeconfig["rules"]["msi"]["msi_list"],
        threads = clusterconf["msi"]["threads"],
    output:
        msi_out = temp("tumor/msi/{sname}_msi_reduced.txt"),
        msi_out_dis = temp("tumor/msi/{sname}_msi_reduced.txt_dis"),
        msi_out_germline = temp("tumor/msi/{sname}_msi_reduced.txt_germline"),
        msi_out_somatic = temp("tumor/msi/{sname}_msi_reduced.txt_somatic"),
    shell:
        """
        {params.msisensor} msi \
            -d {params.reference_list} \
            -n {input.normal_bam} \
            -t {input.tumor_bam} \
            -b {params.threads} \
            -z 1 \
            -o {output.msi_out}
        """
