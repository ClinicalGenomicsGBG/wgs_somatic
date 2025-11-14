
rule tmb_calculation:
    input:
        somatic_vcf = "{stype}/tnscope/{sname}_TNscope_somatic_w_normal.vcf" if normalid else "{stype}/tnscope/{sname}_TNscope_somatic.vcf"
    params:
        calculate_tmb = os.path.join(workflow.basedir, "workflows/scripts/calculate_tmb.sh"),
        bcftools = pipeconfig["rules"]["tmb_calculation"]["bcftools"],
        effective_genome_size = pipeconfig["rules"]["tmb_calculation"]["effective_genome_size"],
        min_mutant_allele_fraction = filterconfig["tmb_filter"]["min_mutant_allele_fraction"],
        min_mutant_allele_reads = filterconfig["tmb_filter"]["min_mutant_allele_reads"],
        min_coverage = filterconfig["tmb_filter"]["min_coverage"],
        include_normal = ("True" if normalid else "False"),
        vstamp = f"{VDIR}/tmb_calculation.txt"
    output:
        tmb = temp("{stype}/reports/{sname}_tmb.txt"),
    shell:
        """
        bcftools --version | head -n 2 > {params.vstamp}
        bash {params.calculate_tmb} \
        {params.bcftools} \
        {input.somatic_vcf} \
        {params.effective_genome_size} \
        {params.min_mutant_allele_fraction} \
        {params.min_mutant_allele_reads} \
        {params.min_coverage} \
        {output.tmb} \
        {params.include_normal}
        """
