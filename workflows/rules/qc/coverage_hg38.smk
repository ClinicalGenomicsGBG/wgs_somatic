# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


rule wgs_coverage:
    input:
        "{stype}/realign/{sname}_REALIGNED.bam"
    params:
        threads = clusterconf["wgs_coverage"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        reference = pipeconfig["singularities"]["sentieon"]["reference"]
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        "{stype}/reports/{sname}_WGScov.tsv"
    shell:
        "{params.sentieon} driver -t {params.threads} -i {input} -r {params.reference} "
            "--interval chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM --algo WgsMetricsAlgo {output}"

rule y_coverage:
    input:
        "{stype}/realign/{sname}_REALIGNED.bam"
    params:
        threads = clusterconf["y_coverage"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        reference = pipeconfig["singularities"]["sentieon"]["reference"]
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        "{stype}/reports/{sname}_Ycov.tsv"
    shell:
        "{params.sentieon} driver -t {params.threads} -i {input} -r {params.reference} --interval chrY --algo WgsMetricsAlgo {output}"
