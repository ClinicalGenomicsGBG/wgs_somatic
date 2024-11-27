if normalid:
    rule tmb_calculation:
        input:
            somatic_n = "{stype}/tnscope/{sname}_somatic_w_normal.vcf"
        params:
            bcftools = pipeconfig["rules"]["tmb_calculation"]["bcftools"],
            genome_size = pipeconfig["rules"]["tmb_calculation"]["genome_size"]
        output:
            tmb = temp("{stype}/reports/{sname}_tmb.txt"),
        run:
            shell_command = (
                f"{params.bcftools} filter -i 'FORMAT/AF[0]>0.05 && FORMAT/AD[0:1]>3 && FORMAT/AFDP[0,1]>10' {input.somatic_n} | "
                f"{params.bcftools} query -f '%POS\\n' | wc -l | "
                f"awk -v genome_size={params.genome_size} '{{{{ printf \"%.2f\", ($0 / genome_size * 1e6) }}}}' > {output.tmb}"
            )
            print(shell_command)
            shell(shell_command)

else:
    rule tmb_calculation:
        input:
            somatic = "{stype}/tnscope/{sname}_somatic.vcf"
        params:
            bcftools = pipeconfig["rules"]["tmb_calculation"]["bcftools"],
            genome_size = pipeconfig["rules"]["tmb_calculation"]["genome_size"]
        output:
            tmb = temp("{stype}/reports/{sname}_tmb.txt"),
        run:
            shell_command = (
                f"{params.bcftools} filter -i 'FORMAT/AF[0]>0.05 && FORMAT/AD[0:1]>3 && FORMAT/AFDP[0]>10' {input.somatic} | "
                f"{params.bcftools} query -f '%POS\\n' | wc -l | "
                f"awk -v genome_size={params.genome_size} '{{{{ printf \"%.2f\", ($0 / genome_size * 1e6) }}}}' > {output.tmb}"
            )
            print(shell_command)
            shell(shell_command)
