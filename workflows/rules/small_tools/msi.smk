# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule msi:
    input:
        tumor_bam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        tumor_bai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        normal_bam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        normal_bai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=normalid, stype=sampleconfig[normalname]["stype"]),
    params:
        reference_list = pipeconfig["rules"]["msi"]["msi_list"],
        threads = clusterconf["msi"]["threads"],
    singularity:
        pipeconfig["singularities"]["msi"]["sing"]
    output:
        msi_out = temp("{stype}/msi/{sname}_msi.txt"),
        msi_out_dis = temp("{stype}/msi/{sname}_msi.txt_dis"),
        msi_out_all = temp("{stype}/msi/{sname}_msi.txt_all"),
        msi_out_unstable = temp("{stype}/msi/{sname}_msi.txt_unstable"),
    shell:
        """
        msisensor-pro msi \
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
    singularity:
        pipeconfig["singularities"]["msi"]["sing"]
    output:
        filtered_bam = temp("{stype}/msi/{sname}_filtered.bam"),
        filtered_bam_bai = temp("{stype}/msi/{sname}_filtered.bam.bai")
    shell:
        """
        bedtools intersect -wa -abam {input.bam} -b {params.bed} > {output.filtered_bam}
        samtools index {output.filtered_bam}
        """

rule msi_reduced:
    input:
        tumor_bam = expand("{stype}/msi/{sname}_filtered.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        tumor_bai = expand("{stype}/msi/{sname}_filtered.bam.bai", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        normal_bam = expand("{stype}/msi/{sname}_filtered.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        normal_bai = expand("{stype}/msi/{sname}_filtered.bam.bai", sname=normalid, stype=sampleconfig[normalname]["stype"]),
    params:
        reference_list = pipeconfig["rules"]["msi"]["msi_list"],
        threads = clusterconf["msi"]["threads"],
    singularity:
        pipeconfig["singularities"]["msi"]["sing"]
    output:
        msi_out = temp("{stype}/msi/{sname}_msi_reduced.txt"),
        msi_out_dis = temp("{stype}/msi/{sname}_msi_reduced.txt_dis"),
        msi_out_all = temp("{stype}/msi/{sname}_msi_reduced.txt_all"),
        msi_out_unstable = temp("{stype}/msi/{sname}_msi_reduced.txt_unstable"),
    shell:
        """
        msisensor-pro msi \
            -d {params.reference_list} \
            -n {input.normal_bam} \
            -t {input.tumor_bam} \
            -b {params.threads} \
            -z 1 \
            -o {output.msi_out}
        """
