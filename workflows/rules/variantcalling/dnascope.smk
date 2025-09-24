# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule dnascope:
    input:
        bam = "{stype}/dedup/{sname}_DEDUP.bam",
        bai = "{stype}/dedup/{sname}_DEDUP.bam.bai"
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
    output:
        germline_vcf = "{stype}/dnascope/{sname}_germline.vcf", 
        germline_snv_vcf = temp("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf"),
    #shadow:
        #pipeconfig["rules"].get("dnascope_vcffilter", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        """
        {params.bcftools} filter -s 'ML_FAIL' -e 'INFO/ML_PROB > 0.5' -m + {input.vcf} > {wildcards.stype}/dnascope/{wildcards.sname}_DNAscope_filtered1.vcf
        {params.bcftools} filter -s 'low_depth' -e 'FORMAT/DP < 10' -m + {wildcards.stype}/dnascope/{wildcards.sname}_DNAscope_filtered1.vcf> {wildcards.stype}/dnascope/{wildcards.sname}_DNAscope_filtered2.vcf
        {params.bcftools} filter -s 'low_genotype_quality' -e 'FORMAT/GQ < 20' -m + {wildcards.stype}/dnascope/{wildcards.sname}_DNAscope_filtered2.vcf > {wildcards.stype}/dnascope/{wildcards.sname}_DNAscope_filtered3.vcf
        {params.bcftools} filter -s 'low_qual_per_depth' -e 'INFO/QD < 4' -m + {wildcards.stype}/dnascope/{wildcards.sname}_DNAscope_filtered3.vcf > {wildcards.stype}/dnascope/{wildcards.sname}_DNAscope_filtered4.vcf
        {params.bcftools} view -f PASS {wildcards.stype}/dnascope/{wildcards.sname}_DNAscope_filtered4.vcf > {output.germline_vcf}
        {params.vcftools} --vcf {output.germline_vcf} --remove-indels --recode --stdout > {output.germline_snv_vcf}
        """