# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


def get_mapping_input(wcs):
    return {
        d: fastq_dict[wcs.stype]["fastqpair_patterns"][wcs.fastqpattern][d].replace(".fasterq", ".fastq.gz")
        for d in ["fwd", "rev"]
    }


def get_mapping(wcs):
    fastqpatterns = fastq_dict[wcs.stype]["fastqpair_patterns"]
    return dict(
        bam = expand(f"{wcs.stype}/mapping/{{fastqpattern}}.bam", fastqpattern=fastqpatterns),
        bai = expand(f"{stype}/mapping/{{fastqpattern}}.bam.bai", fastqpattern=fastqpatterns)
    )

rule mapping:
    input:
        unpack(get_mapping_input)
    params:
        threads = clusterconf["mapping"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        samplename = lambda wildcards: tumorname if wildcards.stype == "tumor" else normalname,
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
    singularity:
        pipeconfig["singularities"]["sentieon"]["simg"]
    shadow:
        pipeconfig["rules"].get("mapping", {}).get("shadow", pipeconfig.get("shadow", False))
    output:
        bam = temp("{stype}/mapping/{fastqpattern}.bam"),
        bai = temp("{stype}/mapping/{fastqpattern}.bam.bai")
    shell:
        "echo $HOSTNAME;"
        "{params.sentieon} bwa mem "
            "-M "
            "-R '@RG\\tID:{wildcards.fastqpattern}\\tSM:{params.samplename}\\tPL:ILLUMINA' "
            "-t {params.threads} "
            "{params.referencegenome} "
            "{input.fwd} "
            "{input.rev} | "
                "{params.sentieon} util sort "
                    "-o {output.bam} "
                    "-t {params.threads} "
                    "--sam2bam "
                    "-i -"

rule dedup:
    input:
        unpack(get_mapping)
    output:
        bam = temp("{stype}/dedup/{sname}_DEDUP.bam"),
        bai = temp("{stype}/dedup/{sname}_DEDUP.bam.bai"),
        score = temp("{stype}/dedup/{sname}_DEDUP_score.txt"),
        metrics = temp("{stype}/dedup/{sname}_DEDUP.txt")
    params:
        threads = clusterconf["dedup"]["threads"],
        samplename = get_samplename,
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
    singularity:
        pipeconfig["singularities"]["sentieon"]["simg"]
    shadow:
        pipeconfig["rules"].get("dedup", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "echo $HOSTNAME;"
        "{params.sentieon} driver -t {params.threads} "
            "-i {input.bam} "
            "--algo LocusCollector "
            "--fun score_info "
            "{output.score}; "
        "{params.sentieon} driver "
            "-t {params.threads} "
            "-i {input.bam} "
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
        pipeconfig["singularities"]["sentieon"]["simg"]
    output:
        bam = temp("{stype}/realign/{sname}_REALIGNED.bam"),
        bai = temp("{stype}/realign/{sname}_REALIGNED.bam.bai")
    params:
        threads = clusterconf["realign_mapping"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
        mills = pipeconfig["singularities"]["sentieon"]["mills"],
        hg19args = f"-k {pipeconfig['singularities']['sentieon']['tgenomes']}" if reference == "hg19" else ""
    shadow:
        pipeconfig["rules"].get("realign_mapping", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "echo $HOSTNAME;"
        "{params.sentieon} driver "
            "-t {params.threads} "
            "-r {params.referencegenome} "
            "-i {input.bam} "
            "--algo Realigner "
            "-k {params.mills} "
            "{params.hg19args} "
            "{output.bam}"

rule baserecal:
    input:
        bam = "{stype}/realign/{sname}_REALIGNED.bam",
        bai = "{stype}/realign/{sname}_REALIGNED.bam.bai"
    singularity:
        pipeconfig["singularities"]["sentieon"]["simg"]
    params:
        threads = clusterconf["baserecal"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
        dbsnp = pipeconfig["singularities"]["sentieon"]["dbsnp"],
        mills = pipeconfig["singularities"]["sentieon"]["mills"],
        hg19args = f"-k {pipeconfig['singularities']['sentieon']['tgenomes']}" if reference == "hg19" else ""
    output:
        temp("{stype}/recal/{sname}_RECAL_DATA.TABLE")
    shadow:
        pipeconfig["rules"].get("baserecal", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "echo $HOSTNAME;"
        "{params.sentieon} driver "
            "-t {params.threads} "
            "-r {params.referencegenome} "
            "-i {input.bam} "
            "--algo QualCal "
            "-k {params.mills} "
            "-k {params.dbsnp} "
            "{params.hg19args} "
            "{output}"
