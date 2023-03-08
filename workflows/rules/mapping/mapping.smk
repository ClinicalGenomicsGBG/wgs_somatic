# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

def format_fwd(wcs):
    fastq = fastq_dict[f"{wcs.stype}"]["fastqpair_patterns"][f"{wcs.fastqpattern}"]["fwd"]
    if fastq.endswith(".fasterq"):
        fastq = fastq.replace(".fasterq", ".fastq.gz")
    return fastq

def format_rev(wcs):
    fastq = fastq_dict[f"{wcs.stype}"]["fastqpair_patterns"][f"{wcs.fastqpattern}"]["rev"]
    if fastq.endswith(".fasterq"):
        fastq = fastq.replace(".fasterq", ".fastq.gz")
    return fastq

def get_samplename(wcs):
    return sampleconfig[f"{wcs.stype}name"]

def get_mapping(wcs):
    fastqpatterns = []
    for fastqpattern in fastq_dict[f"{wcs.stype}"]["fastqpair_patterns"]:
        fastqpatterns.append(fastqpattern)
    return expand("{stype}/mapping/{fastqpattern}.bam", stype=f"{wcs.stype}", fastqpattern=fastqpatterns)

rule mapping:
    input:
        fwd_fmt = format_fwd,
        rev_fmt = format_rev
    params:
        threads = clusterconf["mapping"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        samplename = get_samplename,
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        temp("{stype}/mapping/{fastqpattern}.bam")
    shell:
        "echo $HOSTNAME;"
        "{params.sentieon} bwa mem "
            "-M -R '@RG\\tID:{wildcards.fastqpattern}\\tSM:{params.samplename}\\tPL:ILLUMINA' "
            "-t {params.threads} {params.referencegenome} {input.fwd_fmt} {input.rev_fmt} "
        "| {params.sentieon} util sort -o {output} -t {params.threads} --sam2bam -i -"

rule dedup:
    input:
        bamfiles = get_mapping
    output:
        temp("{stype}/dedup/{sname}_DEDUP.bam"),
        "{stype}/dedup/{sname}_DEDUP.txt"
    params:
        threads = clusterconf["dedup"]["threads"],
        samplename = get_samplename,
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    shell:
        "echo $HOSTNAME;"
        "shellbamfiles=$(echo {input.bamfiles} | sed 's/ / -i /g') ;"
        "{params.sentieon} driver -t {params.threads} "
            "-i $shellbamfiles "
            "--algo LocusCollector "
            "--fun score_info "
            "{wildcards.stype}/dedup/{wildcards.sname}_DEDUP_score.txt ;"
        "{params.sentieon} driver "
            "-t {params.threads} "
            "-i $shellbamfiles "
            "--algo Dedup "
            "--rmdup "
            "--score_info {wildcards.stype}/dedup/{wildcards.sname}_DEDUP_score.txt "
            "--metrics {wildcards.stype}/dedup/{wildcards.sname}_DEDUP.txt "
            "{wildcards.stype}/dedup/{wildcards.sname}_DEDUP.bam"

rule realign_mapping:
    input:
        "{stype}/dedup/{sname}_DEDUP.bam"
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        "{stype}/realign/{sname}_REALIGNED.bam"
    params:
        threads = clusterconf["realign_mapping"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
        mills = pipeconfig["singularities"]["sentieon"]["mills"],
        tgenomes = pipeconfig["singularities"]["sentieon"]["tgenomes"]
    shell:
        "echo $HOSTNAME;"
        "{params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo Realigner -k {params.mills} -k {params.tgenomes} {output}"

rule baserecal:
    input:
        "{stype}/realign/{sname}_REALIGNED.bam"
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    params:
        threads = clusterconf["baserecal"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
        dbsnp = pipeconfig["singularities"]["sentieon"]["dbsnp"],
        mills = pipeconfig["singularities"]["sentieon"]["mills"],
        tgenomes = pipeconfig["singularities"]["sentieon"]["tgenomes"]
    output:
        "{stype}/recal/{sname}_RECAL_DATA.TABLE"
    shell:
        "echo $HOSTNAME;"
        "{params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo QualCal -k {params.mills} -k {params.dbsnp} -k {params.tgenomes} {output}"
