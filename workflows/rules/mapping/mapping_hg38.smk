# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

def get_fwd_pattern(wcs):
    return fastq_dict[f"{wcs.stype}"]["fastqpair_patterns"][f"{wcs.fastqpattern}"]["fwd"] 

def get_rev_pattern(wcs):
    return fastq_dict[f"{wcs.stype}"]["fastqpair_patterns"][f"{wcs.fastqpattern}"]["rev"]

def format_fwd(wcs):
    fastq = fastq_dict[f"{wcs.stype}"]["fastqpair_patterns"][f"{wcs.fastqpattern}"]["rev"]
    if fastq.endswith(".fasterq"):
        fastq = fastq.replace(".fasterq", ".fastq.gz")
    return fastq

def format_rev(wcs):
    fastq = fastq_dict[f"{wcs.stype}"]["fastqpair_patterns"][f"{wcs.fastqpattern}"]["fwd"]
    if fastq.endswith(".fasterq"):
        fastq = fastq.replace(".fasterq", ".fastq.gz")
    return fastq

def get_samplename(wcs):
    return sampleconfig[f"{wcs.stype}name"]

def get_mapping(wcs):
    fastqpatterns = []
    for fastqpattern in fastq_dict[f"{wcs.stype}"]["fastqpair_patterns"]:
        fastqpatterns.append(fastqpattern)
    return expand("{workingdir}/{stype}/mapping/{fastqpattern}.bam", workingdir=f"{wcs.workingdir}", stype=f"{wcs.stype}", fastqpattern=fastqpatterns)

rule mapping:
    input:
        fwd = get_fwd_pattern,
        rev = get_rev_pattern
    params:
        threads = clusterconf["mapping"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        samplename = get_samplename,
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
        fwd_fmt = format_fwd,
        rev_fmt = format_rev
    shadow: "copy-minimal"
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        temp("{workingdir}/{stype}/mapping/{fastqpattern}.bam")
    shell:
        "echo $HOSTNAME;"
        "{params.sentieon} bwa mem "
            "-M -R '@RG\\tID:{wildcards.fastqpattern}\\tSM:{params.samplename}\\tPL:ILLUMINA' "
            "-t {params.threads} {params.referencegenome} {params.fwd_fmt} {params.rev_fmt} "
        "| {params.sentieon} util sort -o {output} -t {params.threads} --sam2bam -i -"

rule dedup:
    input:
        bamfiles = get_mapping
    output:
        temp("{workingdir}/{stype}/dedup/{sname}_DEDUP.bam"),
        "{workingdir}/{stype}/dedup/{sname}_DEDUP.txt"
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
            "{wildcards.workingdir}/{wildcards.stype}/dedup/{wildcards.sname}_DEDUP_score.txt ;"
        "{params.sentieon} driver "
            "-t {params.threads} "
            "-i $shellbamfiles "
            "--algo Dedup "
            "--rmdup "
            "--score_info {wildcards.workingdir}/{wildcards.stype}/dedup/{wildcards.sname}_DEDUP_score.txt "
            "--metrics {wildcards.workingdir}/{wildcards.stype}/dedup/{wildcards.sname}_DEDUP.txt "
            "{wildcards.workingdir}/{wildcards.stype}/dedup/{wildcards.sname}_DEDUP.bam"

rule realign_mapping:
    input:
        "{workingdir}/{stype}/dedup/{sname}_DEDUP.bam"
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        "{workingdir}/{stype}/realign/{sname}_REALIGNED.bam"
    params:
        threads = clusterconf["realign_mapping"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
        mills = pipeconfig["singularities"]["sentieon"]["mills"]
    shell:
        "echo $HOSTNAME;"
        "{params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo Realigner -k {params.mills} {output}"

rule baserecal:
    input:
        "{workingdir}/{stype}/realign/{sname}_REALIGNED.bam"
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    params:
        threads = clusterconf["baserecal"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
        dbsnp = pipeconfig["singularities"]["sentieon"]["dbsnp"],
        mills = pipeconfig["singularities"]["sentieon"]["mills"],
    output:
        "{workingdir}/{stype}/recal/{sname}_RECAL_DATA.TABLE"
    shell:
        "echo $HOSTNAME;"
        "{params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo QualCal -k {params.mills} -k {params.dbsnp} {output}"
