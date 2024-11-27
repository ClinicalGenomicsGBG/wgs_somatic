if normalid:
    rule tmb_calculation:
        input:
            somatic_n = "{stype}/tnscope/{sname}_somatic_w_normal.vcf"
        params:
            bcftools = pipeconfig["rules"]["tmb_calculation"]["bcftools"],
            effective_genome_size = pipeconfig["rules"]["tmb_calculation"]["effective_genome_size"],
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
                f"awk -v gsize={params.effective_genome_size} '{{{{ printf \"TMB\\t%.2f\\n\", ($1 / gsize * 1e6) }}}}' > {output.tmb}"
            )
            add_parameters = (
                f"echo -e 'min_mutant_allele_fraction\\t{params.min_mutant_allele_fraction}\\n"
                f"min_mutant_allele_reads\\t{params.min_mutant_allele_reads}\\n"
                f"min_coverage\\t{params.min_coverage}\\n"
                f"effective_genome_size\\t{params.effective_genome_size}' >> {output.tmb}"
            )
            print(shell_command)
            shell(shell_command)
            shell(add_parameters)

else:
    rule tmb_calculation:
        input:
            somatic = "{stype}/tnscope/{sname}_somatic.vcf"
        params:
            bcftools = pipeconfig["rules"]["tmb_calculation"]["bcftools"],
            effective_genome_size = pipeconfig["rules"]["tmb_calculation"]["effective_genome_size"],
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
                f"awk -v gsize={params.effective_genome_size} '{{{{ printf \"TMB\\t%.2f\\n\", ($1 / gsize * 1e6) }}}}' > {output.tmb}"
            )
            add_parameters = (
                f"echo -e 'min_mutant_allele_fraction\\t{params.min_mutant_allele_fraction}\\n"
                f"min_mutant_allele_reads\\t{params.min_mutant_allele_reads}\\n"
                f"min_coverage\\t{params.min_coverage}\\n"
                f"effective_genome_size\\t{params.effective_genome_size}' >> {output.tmb}"
            )
            print(shell_command)
            shell(shell_command)
            shell(add_parameters)