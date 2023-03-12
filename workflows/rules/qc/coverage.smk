# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


rule wgs_coverage:
    input:
        bam="{stype}/realign/{sname}_REALIGNED.bam",
        bai="{stype}/realign/{sname}_REALIGNED.bam.bai",
    params:
        threads=clusterconf["wgs_coverage"]["threads"],
        sentieon=pipeconfig["singularities"]["sentieon"]["tool_path"],
        reference=pipeconfig["singularities"]["sentieon"]["reference"],
        interval=(
            "1,2,3,4,5,"
            "6,7,8,9,10,"
            "11,12,13,14,15,"
            "16,17,18,19,20,"
            "21,22,X,Y,MT"
        )
        if reference == "hg19"
        else (
            "chr1,chr2,chr3,chr4,chr5,"
            "chr6,chr7,chr8,chr9,chr10,"
            "chr11,chr12,chr13,chr14,chr15,"
            "chr16,chr17,chr18,chr19,chr20,"
            "chr21,chr22,chrX,chrY,chrM"
        ),
    singularity:
        pipeconfig["singularities"]["sentieon"]["simg"]
    output:
        temp("{stype}/reports/{sname}_WGScov.tsv"),
    shadow:
        pipeconfig["rules"].get("wgs_coverage", {}).get(
            "shadow", pipeconfig.get("shadow", False)
        )
    shell:
        "{params.sentieon} driver "
            "-t {params.threads} "
            "-i {input.bam} "
            "-r {params.reference} "
            "--interval {params.interval} "
            "--algo WgsMetricsAlgo "
            "{output}"


rule y_coverage:
    input:
        bam="{stype}/realign/{sname}_REALIGNED.bam",
        bai="{stype}/realign/{sname}_REALIGNED.bam.bai",
    params:
        threads=clusterconf["y_coverage"]["threads"],
        sentieon=pipeconfig["singularities"]["sentieon"]["tool_path"],
        reference=pipeconfig["singularities"]["sentieon"]["reference"],
        interval="Y" if reference == "hg19" else "chrY",
    singularity:
        pipeconfig["singularities"]["sentieon"]["simg"]
    output:
        temp("{stype}/reports/{sname}_Ycov.tsv"),
    shadow:
        pipeconfig["rules"].get("y_coverage", {}).get(
            "shadow", pipeconfig.get("shadow", False)
        )
    shell:
        "{params.sentieon} driver "
            "-t {params.threads} "
            "-i {input.bam} "
            "-r {params.reference} "
            "--interval {params.interval} "
            "--algo WgsMetricsAlgo "
            "{output}"
