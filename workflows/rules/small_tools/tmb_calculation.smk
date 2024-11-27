if normalid:
    rule tmb_calculation:
        input:
            somatic_n = "{stype}/tnscope/{sname}_somatic_w_normal.vcf"
        params:
            bcftools = pipeconfig["rules"]["tmb_calculation"]["bcftools"],
            genome_size = pipeconfig["rules"]["tmb_calculation"]["genome_size"],
            min_mutant_allele_fraction = filterconfig["tmb_filter"]["min_mutant_allele_fraction"],
            min_mutant_allele_reads = filterconfig["tmb_filter"]["min_mutant_allele_reads"],
            min_coverage = filterconfig["tmb_filter"]["min_coverage"],
        output:
            tmb = temp("{stype}/reports/{sname}_tmb.txt"),
        run:
            shell_command = (
                f"{params.bcftools} filter -i 'FORMAT/AF[0]>{params.min_mutant_allele_fraction} "
                f"&& FORMAT/AD[0:1]>{params.min_mutant_allele_reads} "
                f"&& FORMAT/AFDP[0,1]>{params.min_coverage}' {input.somatic_n} | "
                f"{params.bcftools} stats | grep 'number of records:' | awk '{{{{print $NF}}}}' | "
                f"awk -v genome_size={params.genome_size} '{{{{ printf \"%.2f\", ($1 / genome_size * 1e6) }}}}' > {output.tmb}"
            )
            print(shell_command)
            shell(shell_command)

else:
    rule tmb_calculation:
        input:
            somatic = "{stype}/tnscope/{sname}_somatic.vcf"
        params:
            bcftools = pipeconfig["rules"]["tmb_calculation"]["bcftools"],
            genome_size = pipeconfig["rules"]["tmb_calculation"]["genome_size"],
            min_mutant_allele_fraction = filterconfig["tmb_filter"]["min_mutant_allele_fraction"],
            min_mutant_allele_reads = filterconfig["tmb_filter"]["min_mutant_allele_reads"],
            min_coverage = filterconfig["tmb_filter"]["min_coverage"],
        output:
            tmb = temp("{stype}/reports/{sname}_tmb.txt"),
        run:
            shell_command = (
                f"{params.bcftools} filter -i 'FORMAT/AF[0]>{params.min_mutant_allele_fraction} "
                f"&& FORMAT/AD[0:1]>{params.min_mutant_allele_reads} "
                f"&& FORMAT/AFDP[0]>{params.min_coverage}' {input.somatic} | "  # Note, here we only use the first sample
                f"{params.bcftools} stats | grep 'number of records:' | awk '{{{{print $NF}}}}' | "
                f"awk -v genome_size={params.genome_size} '{{{{ printf \"%.2f\", ($1 / genome_size * 1e6) }}}}' > {output.tmb}"
            )
            print(shell_command)
            shell(shell_command)