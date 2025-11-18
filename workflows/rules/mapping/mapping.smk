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
        vstamp = f"{VDIR}/mapping.txt"
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    shadow:
        pipeconfig["rules"].get("mapping", {}).get("shadow", pipeconfig.get("shadow", False))
    output:
        bam = temp("{stype}/mapping/{fastqpattern}.bam"),
        bai = temp("{stype}/mapping/{fastqpattern}.bam.bai")
    shell:
        """
        # Version info
        echo bwa: $({params.sentieon} bwa 2>&1 | grep Version) > {params.vstamp}
        {params.sentieon} driver --version >> {params.vstamp}

        # Run bwa mem and sort
        echo $HOSTNAME
        {params.sentieon} bwa mem \
            -M -R '@RG\\tID:{wildcards.fastqpattern}\\tSM:{params.samplename}\\tPL:ILLUMINA' \
            -t {params.threads} {params.referencegenome} {input.fwd_fmt} {input.rev_fmt} \
            | {params.sentieon} util sort -o {output.bam} -t {params.threads} --sam2bam -i -
        """

rule dedup:
    input:
        unpack(get_mapping)
    params:
        threads = clusterconf["dedup"]["threads"],
        samplename = get_samplename,
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        bamfiles = lambda wildcards, input: "-i " + " -i ".join(input.bam),
        vstamp = f"{VDIR}/dedup.txt"
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    shadow:
        pipeconfig["rules"].get("dedup", {}).get("shadow", pipeconfig.get("shadow", False))
    output:
        bam = temp("{stype}/dedup/{sname}_DEDUP.bam"),
        bai = temp("{stype}/dedup/{sname}_DEDUP.bam.bai"),
        score = temp("{stype}/dedup/{sname}_DEDUP_score.txt"),
        metrics = temp("{stype}/dedup/{sname}_DEDUP.txt"),
    shell:
        """
        # Version info
        {params.sentieon} driver --version > {params.vstamp}

        # Run Dedup
        echo $HOSTNAME
        {params.sentieon} driver \
            -t {params.threads} \
            {params.bamfiles} \
            --algo LocusCollector \
            --fun score_info \
            {output.score}
        {params.sentieon} driver \
            -t {params.threads} \
            {params.bamfiles} \
            --algo Dedup \
            --rmdup \
            --score_info {output.score} \
            --metrics {output.metrics} \
            {output.bam}
        """

rule realign_mapping:
    input:
        bam = "{stype}/dedup/{sname}_DEDUP.bam",
        bai = "{stype}/dedup/{sname}_DEDUP.bam.bai"
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    params:
        threads = clusterconf["realign_mapping"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
        mills = pipeconfig["singularities"]["sentieon"]["mills"],
        vstamp = f"{VDIR}/realign_mapping.txt"
    shadow:
        pipeconfig["rules"].get("realign_mapping", {}).get("shadow", pipeconfig.get("shadow", False))
    output:
        bam = "{stype}/realign/{sname}_REALIGNED.bam",
        bai = "{stype}/realign/{sname}_REALIGNED.bam.bai"
    shell:
        """
        # Version info
        {params.sentieon} driver --version > {params.vstamp}

        # Run Realigner
        echo $HOSTNAME
        {params.sentieon} driver \
            -t {params.threads} \
            -r {params.referencegenome} \
            -i {input.bam} \
            --algo Realigner \
            -k {params.mills} \
            {output.bam}
        """

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
        mills = pipeconfig["singularities"]["sentieon"]["mills"],
        vstamp = f"{VDIR}/baserecal.txt"
    shadow:
        pipeconfig["rules"].get("baserecal", {}).get("shadow", pipeconfig.get("shadow", False))
    output:
        recal = temp("{stype}/recal/{sname}_RECAL_DATA.TABLE"),
    shell:
        """
        # Version info
        {params.sentieon} driver --version > {params.vstamp}
        basename {params.dbsnp} >> {params.vstamp}
        basename {params.mills} >> {params.vstamp}

        # Run QualCal
        echo $HOSTNAME
        {params.sentieon} driver \
            -t {params.threads} \
            -r {params.referencegenome} \
            -i {input.bam} \
            --algo QualCal \
            -k {params.mills} \
            -k {params.dbsnp} \
            {output}
        """
