# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule dnascope:
    input:
        bam = "{stype}/dedup/{sname}_REALIGNED.bam",
        bai = "{stype}/dedup/{sname}_REALIGNED.bam.bai"
    params:
        threads = clusterconf["dnascope"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        reference = pipeconfig["singularities"]["sentieon"]["reference"],
        dbsnp = pipeconfig["singularities"]["sentieon"]["dbsnp"],
        model = pipeconfig["singularities"]["sentieon"]["dnascope_m"],
        callsettings = pipeconfig["rules"]["dnascope"]["settings"]
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        vcf = temp("{stype}/dnascope/{sname}_DNAscope.vcf"), 
        idx = temp("{stype}/dnascope/{sname}_DNAscope.vcf.idx"),
    shadow:
        pipeconfig["rules"].get("dnascope", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "echo $HOSTNAME;"
        "{params.sentieon} driver -t {params.threads} -r {params.reference} "
            "-i {input.bam} --algo DNAscope -d {params.dbsnp} "
            "--var_type snp,indel --model {params.model} {params.callsettings} {output.vcf}"
        
rule dnascope_modelfilter:
    input:
        vcf = "{stype}/dnascope/{sname}_DNAscope.vcf",
        idx = "{stype}/dnascope/{sname}_DNAscope.vcf.idx"
    params:
        threads = clusterconf["dnascope_modelfilter"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        reference = pipeconfig["singularities"]["sentieon"]["reference"],
        model = pipeconfig["singularities"]["sentieon"]["dnascope_m"],
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        vcf = temp("{stype}/dnascope/{sname}_DNAscope_modelfiltered.vcf"),
        idx = temp("{stype}/dnascope/{sname}_DNAscope_modelfiltered.vcf.idx"),
    shadow:
        pipeconfig["rules"].get("dnascope_modelfilter", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "echo $HOSTNAME;"
        "{params.sentieon} driver -t {params.threads} -r {params.reference} --algo DNAModelApply --model {params.model} -v {input.vcf} {output.vcf}"

rule dnascope_vcffilter:
    input:
        vcf = "{stype}/dnascope/{sname}_DNAscope_modelfiltered.vcf",
        idx = "{stype}/dnascope/{sname}_DNAscope_modelfiltered.vcf.idx"
    params:
        threads = clusterconf["dnascope_vcffilter"]["threads"],
        bcftools = pipeconfig["rules"]["dnascope_vcffilter"]["bcftools"],
        vcftools = pipeconfig["rules"]["dnascope_vcffilter"]["vcftools"],
        passfilter = "'FILTER=\"PASS\"'"
    output:
        germline_vcf = "{stype}/dnascope/{sname}_germline.vcf", 
        germline_snv_vcf = temp("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf"),
    shadow:
        pipeconfig["rules"].get("dnascope_vcffilter", {}).get("shadow", pipeconfig.get("shadow", False))
    run:
        shell("{params.bcftools} filter -s 'ML_FAIL' -i 'INFO/ML_PROB <= 0.95' -m x {input.vcf} > {wildcards.stype}/dnascope/{wildcards.sname}_DNAscope_modelfiltered_0.95.vcf")
        shell("{params.bcftools} filter -i {params.passfilter} -m x {wildcards.stype}/dnascope/{wildcards.sname}_DNAscope_modelfiltered_0.95.vcf > {output.germline_vcf}")
        shell("{params.vcftools} --vcf {output.germline_vcf} --remove-indels --recode --stdout > {output.germline_snv_vcf}")
