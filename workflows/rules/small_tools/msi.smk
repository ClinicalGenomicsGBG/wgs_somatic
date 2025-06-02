# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
from tools.helpers import conditional_temp

rule msi:
    input:
        tumor_bam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        tumor_bai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        normal_bam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        normal_bai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=normalid, stype=sampleconfig[normalname]["stype"]),
    params:
        reference_list = pipeconfig["rules"]["msi"]["msi_list"],
    threads:
        clusterconf["msi"]["threads"]
    singularity:
        pipeconfig["singularities"]["msi"]["sing"]
    output:
        msi_out = conditional_temp("{stype}/msi/{sname}_msi.txt", keepfiles),
        msi_out_dis = conditional_temp("{stype}/msi/{sname}_msi.txt_dis", keepfiles),
        msi_out_all = conditional_temp("{stype}/msi/{sname}_msi.txt_all", keepfiles),
        msi_out_unstable = conditional_temp("{stype}/msi/{sname}_msi.txt_unstable", keepfiles),
    shell:
        """
        msisensor-pro msi \
            -d {params.reference_list} \
            -n {input.normal_bam} \
            -t {input.tumor_bam} \
            -b {threads} \
            -z 1 \
            -o {output.msi_out}
        """

rule msi_filter_bam:
    input:
        bam = "{stype}/realign/{sname}_REALIGNED.bam",
    params:
        bed = pipeconfig["rules"]["msi"]["msi_bed"],
    singularity:
        pipeconfig["singularities"]["msi"]["sing"]
    output:
        filtered_bam = conditional_temp("{stype}/msi/{sname}_filtered.bam", keepfiles),
        filtered_bam_bai = conditional_temp("{stype}/msi/{sname}_filtered.bam.bai", keepfiles)
    shell:
        """
        bedtools intersect -wa -abam {input.bam} -b {params.bed} > {output.filtered_bam}
        samtools index {output.filtered_bam}
        """

rule msi_filtered:
    input:
        tumor_bam = expand("{stype}/msi/{sname}_filtered.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        tumor_bai = expand("{stype}/msi/{sname}_filtered.bam.bai", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        normal_bam = expand("{stype}/msi/{sname}_filtered.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        normal_bai = expand("{stype}/msi/{sname}_filtered.bam.bai", sname=normalid, stype=sampleconfig[normalname]["stype"]),
    params:
        reference_list = pipeconfig["rules"]["msi"]["msi_list"],
    threads:
        clusterconf["msi"]["threads"]
    singularity:
        pipeconfig["singularities"]["msi"]["sing"]
    output:
        msi_out = conditional_temp("{stype}/msi/{sname}_msi_filtered.txt", keepfiles),
        msi_out_dis = conditional_temp("{stype}/msi/{sname}_msi_filtered.txt_dis", keepfiles),
        msi_out_all = conditional_temp("{stype}/msi/{sname}_msi_filtered.txt_all", keepfiles),
        msi_out_unstable = conditional_temp("{stype}/msi/{sname}_msi_filtered.txt_unstable", keepfiles),
    shell:
        """
        msisensor-pro msi \
            -d {params.reference_list} \
            -n {input.normal_bam} \
            -t {input.tumor_bam} \
            -b {threads} \
            -z 1 \
            -o {output.msi_out}
        """
