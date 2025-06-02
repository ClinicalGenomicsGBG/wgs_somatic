# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
from tools.helpers import conditional_temp

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
    return dict(
        bam = expand("{stype}/mapping/{fastqpattern}.bam", stype=f"{wcs.stype}", fastqpattern=fastqpatterns),
        bai = expand("{stype}/mapping/{fastqpattern}.bam.bai", stype=f"{wcs.stype}", fastqpattern=fastqpatterns)
    )

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
    shadow:
        pipeconfig["rules"].get("mapping", {}).get("shadow", pipeconfig.get("shadow", False))
    output:
        bam = conditional_temp("{stype}/mapping/{fastqpattern}.bam", keepfiles),
        bai = conditional_temp("{stype}/mapping/{fastqpattern}.bam.bai", keepfiles)
    shell:
        "echo $HOSTNAME;"
        "{params.sentieon} bwa mem "
            "-M -R '@RG\\tID:{wildcards.fastqpattern}\\tSM:{params.samplename}\\tPL:ILLUMINA' "
            "-t {params.threads} {params.referencegenome} {input.fwd_fmt} {input.rev_fmt} "
        "| {params.sentieon} util sort -o {output.bam} -t {params.threads} --sam2bam -i -"

rule dedup:
    input:
        unpack(get_mapping)
    output:
        bam = conditional_temp("{stype}/dedup/{sname}_DEDUP.bam", keepfiles),
        bai = conditional_temp("{stype}/dedup/{sname}_DEDUP.bam.bai", keepfiles),
        score = conditional_temp("{stype}/dedup/{sname}_DEDUP_score.txt", keepfiles),
        metrics = conditional_temp("{stype}/dedup/{sname}_DEDUP.txt", keepfiles)
    params:
        threads = clusterconf["dedup"]["threads"],
        samplename = get_samplename,
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        bamfiles = lambda wildcards, input: "-i " + " -i ".join(input.bam)
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    shadow:
        pipeconfig["rules"].get("dedup", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "echo $HOSTNAME;"
        "{params.sentieon} driver -t {params.threads} "
            "{params.bamfiles} "
            "--algo LocusCollector "
            "--fun score_info "
            "{output.score}; "
        "{params.sentieon} driver "
            "-t {params.threads} "
            "{params.bamfiles} "
            "--algo Dedup "
            "--rmdup "
            "--score_info {output.score} "
            "--metrics {output.metrics} "
            "{output.bam}"

rule realign_mapping:
    input:
        bam = "{stype}/dedup/{sname}_DEDUP.bam",
        bai = "{stype}/dedup/{sname}_DEDUP.bam.bai"
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        bam = conditional_temp("{stype}/realign/{sname}_REALIGNED.bam", keepfiles),
        bai = conditional_temp("{stype}/realign/{sname}_REALIGNED.bam.bai", keepfiles)
    params:
        threads = clusterconf["realign_mapping"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
        mills = pipeconfig["singularities"]["sentieon"]["mills"]
    shadow:
        pipeconfig["rules"].get("realign_mapping", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "echo $HOSTNAME;"
        "{params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input.bam} --algo Realigner -k {params.mills} {output.bam}"

rule baserecal:
    input:
        bam = "{stype}/realign/{sname}_REALIGNED.bam",
        bai = "{stype}/realign/{sname}_REALIGNED.bam.bai"
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    params:
        threads = clusterconf["baserecal"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
        dbsnp = pipeconfig["singularities"]["sentieon"]["dbsnp"],
        mills = pipeconfig["singularities"]["sentieon"]["mills"]
    output:
        conditional_temp("{stype}/recal/{sname}_RECAL_DATA.TABLE", keepfiles)
    shadow:
        pipeconfig["rules"].get("baserecal", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "echo $HOSTNAME;"
        "{params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input.bam} --algo QualCal -k {params.mills} -k {params.dbsnp} {output}"
