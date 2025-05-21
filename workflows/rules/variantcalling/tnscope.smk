# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

if normalid:
    rule tnscope:
        input:
            tumorbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            tumorbai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            normalbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            normalbai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            tumortable = expand("{stype}/recal/{sname}_RECAL_DATA.TABLE", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            normaltable = expand("{stype}/recal/{sname}_RECAL_DATA.TABLE", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        params:
            threads = clusterconf["tnscope"]["threads"],
            tumorname = tumorname,
            normalname = normalname,
            sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
            reference = pipeconfig["singularities"]["sentieon"]["reference"],
            dbsnp = pipeconfig["singularities"]["sentieon"]["dbsnp"],
            modelpath = pipeconfig["singularities"]["sentieon"]["tnscope_m"],
        singularity:
            pipeconfig["singularities"]["sentieon"]["sing"]
        output:
            tnscope_vcf = temp("{stype}/tnscope/{sname}_TNscope_tn.vcf"),
            tnscope_idx = temp("{stype}/tnscope/{sname}_TNscope_tn.vcf.idx"),
        shadow:
            pipeconfig["rules"].get("tnscope", {}).get("shadow", pipeconfig.get("shadow", False))
        shell:
            """
            {params.sentieon} driver -r {params.reference} -t {params.threads} \
                -i {input.tumorbam} -q {input.tumortable} \
                -i {input.normalbam} -q {input.normaltable} \
                --algo TNscope --tumor_sample {params.tumorname} --normal_sample {params.normalname} \
                --model {params.modelpath} {output.tnscope_vcf} || \
                {{ echo 'TNscope failed'; exit 1; }}
            """
else:
    rule tnscope:
        input:
            tumorbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            tumorbai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            tumortable = expand("{stype}/recal/{sname}_RECAL_DATA.TABLE", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        params:
            threads = clusterconf["tnscope"]["threads"],
            tumorname = tumorname,
            sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
            reference = pipeconfig["singularities"]["sentieon"]["reference"],
            dbsnp = pipeconfig["singularities"]["sentieon"]["dbsnp"],
            modelpath = pipeconfig["singularities"]["sentieon"]["tnscope_m"],
            pon = pipeconfig["rules"]["tnscope"]["pon"],
        singularity:
            pipeconfig["singularities"]["sentieon"]["sing"]
        output:
            tnscope_vcf = temp("{stype}/tnscope/{sname}_TNscope_tn.vcf"),
            tnscope_idx = temp("{stype}/tnscope/{sname}_TNscope_tn.vcf.idx"),
        shadow:
            pipeconfig["rules"].get("tnscope", {}).get("shadow", pipeconfig.get("shadow", False))
        shell:
            """
            {params.sentieon} driver -r {params.reference} -t {params.threads} \
                -i {input.tumorbam} -q {input.tumortable} \
                --algo TNscope --tumor_sample {params.tumorname} --pon {params.pon} \
                --model {params.modelpath} {output.tnscope_vcf} || \
                {{ echo 'TNscope failed'; exit 1; }}
            """

if normalid:
    rule tnscope_modelfilter:
        input:
            tnscopevcf = "{stype}/tnscope/{sname}_TNscope_tn.vcf",
            tnscopeidx = "{stype}/tnscope/{sname}_TNscope_tn.vcf.idx"
        params:
            threads = clusterconf["tnscope_modelfilter"]["threads"],
            sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
            reference = pipeconfig["singularities"]["sentieon"]["reference"],
            modelpath = pipeconfig["singularities"]["sentieon"]["tnscope_m"],
        singularity:
            pipeconfig["singularities"]["sentieon"]["sing"]
        output:
            vcf = temp("{stype}/tnscope/{sname}_TNscope_tn_ML.vcf"),
            idx = temp("{stype}/tnscope/{sname}_TNscope_tn_ML.vcf.idx")
        shadow:
            pipeconfig["rules"].get("tnscope_modelfilter", {}).get("shadow", pipeconfig.get("shadow", False))
        shell:
            """
            {params.sentieon} driver -r {params.reference} -t {params.threads} --algo TNModelApply \
                --model {params.modelpath} -v {input.tnscopevcf} {output.vcf} || \
                {{ echo 'TNModelApply failed'; exit 1; }}
            """

# In the latest version the ML_PROB was inverted, with low scores meaning high confidence
if normalid:
    rule tnscope_vcffilter:
        input:
            tnscopevcf_ml = "{stype}/tnscope/{sname}_TNscope_tn_ML.vcf"
        params:
            threads = clusterconf["tnscope_vcffilter"]["threads"],
            outputdir = pipeconfig["rules"]["tnscope_vcffilter"]["outputdir"],
            bcftools = pipeconfig["rules"]["tnscope_vcffilter"]["bcftools"]
        output:
            somatic_n = temp("{stype}/tnscope/{sname}_somatic_w_normal.vcf"),
            somatic = "{stype}/tnscope/{sname}_somatic.vcf",
            tnscope_filterstats = "{stype}/tnscope/{sname}_tnscope_filterstats.txt",
        #shadow:
        #    pipeconfig["rules"].get("tnscope_vcffilter", {}).get("shadow", pipeconfig.get("shadow", False))
        run:
            vcfname = os.path.basename(f"{input.tnscopevcf_ml}")
            vcfname = vcfname.replace(".vcf", "")
            shell_command = [f"{params.bcftools} filter -s LowQual -e 'QUAL < 1' -m + {input.tnscopevcf_ml}", f"> {params.outputdir}/{vcfname}_lowqual1.vcf ;",
                            f"{params.bcftools} annotate -x FILTER/triallelic_site {params.outputdir}/{vcfname}_lowqual1.vcf", f"> {params.outputdir}/{vcfname}_triallelic2.vcf ;",
                            f"{params.bcftools} annotate -x FILTER/alt_allele_in_normal {params.outputdir}/{vcfname}_triallelic2.vcf", f"> {params.outputdir}/{vcfname}_altalleleinnormal4.vcf ;",
                            f"{params.bcftools} filter -s uncertainAF -e 'FORMAT/AF[0]<0.045 && FORMAT/AD[0:1]<4' -m + {params.outputdir}/{vcfname}_altalleleinnormal4.vcf", f"> {params.outputdir}/{vcfname}_uncertainaf6.vcf ;",
                            f"{params.bcftools} filter -s likely_artifact -e 'FORMAT/AF[0]<0.1 && FORMAT/AF[1]>0.06' -m + {params.outputdir}/{vcfname}_uncertainaf6.vcf", f"> {params.outputdir}/{vcfname}_likelyartifact7.vcf ;",
                            f"{params.bcftools} filter -s lowAD -e 'FORMAT/AD[0:1] < 3' {params.outputdir}/{vcfname}_likelyartifact7.vcf", f"> {params.outputdir}/{vcfname}_lowad8.vcf ;",
                            f"{params.bcftools} filter -s MLrejected -e 'INFO/ML_PROB>0.63' -m + {params.outputdir}/{vcfname}_lowad8.vcf", f"> {params.outputdir}/{vcfname}_mladjusted9.vcf ;",
                            f"{params.bcftools} filter -s 'orientation_bias' -e 'FMT/FOXOG[0] == 1' -m + {params.outputdir}/{vcfname}_mladjusted9.vcf", f"> {params.outputdir}/{vcfname}_foxogadj.vcf ;",
                            f"{params.bcftools} filter -s 'strand_bias' -e 'SOR > 3' -m + {params.outputdir}/{vcfname}_foxogadj.vcf", f"> {params.outputdir}/{vcfname}_strandbiasadj.vcf ;",
                            f"{params.bcftools} filter -i 'FILTER=\"PASS\"' {params.outputdir}/{vcfname}_strandbiasadj.vcf > {output.somatic_n} ;",
                            f"{params.bcftools} view -s {tumorname} {output.somatic_n} > {output.somatic}"]
            shell_command = " ".join(shell_command)
            print(shell_command)      
            shell(shell_command)
            filter_stats_command = ["awk -F'\t' '!/^#/ {{ split($7, filters, \";\");", 
                                    "for (i in filters) count[filters[i]]++; total++; }}",
                                    "END {{ for (f in count) {{ fraction = count[f] / total; printf \"%s\\t%d\\t%.4f\\n\", f, count[f], fraction; }} }}'",
                                    f"{params.outputdir}/{vcfname}_strandbiasadj.vcf | sort > {output.tnscope_filterstats}"]
            filter_stats_command = " ".join(filter_stats_command)
            print(filter_stats_command)
            shell(filter_stats_command)

else:
    rule tnscope_vcffilter:
        input:
            tnscopevcf_ml = "{stype}/tnscope/{sname}_TNscope_tn.vcf"
        params:
            threads = clusterconf["tnscope_vcffilter"]["threads"],
            outputdir = pipeconfig["rules"]["tnscope_vcffilter"]["outputdir"],
            bcftools = pipeconfig["rules"]["tnscope_vcffilter"]["bcftools"]
        output:
            somatic_n = temp("{stype}/tnscope/{sname}_somatic_w_normal.vcf"),
            somatic = "{stype}/tnscope/{sname}_somatic.vcf",
            tnscope_filterstats = "{stype}/tnscope/{sname}_tnscope_filterstats.txt",
        shadow:
            pipeconfig["rules"].get("tnscope_vcffilter", {}).get("shadow", pipeconfig.get("shadow", False))
        run:
            vcfname = os.path.basename(f"{input.tnscopevcf_ml}")
            vcfname = vcfname.replace(".vcf", "")
            shell_command = [f"{params.bcftools} filter -s LowQual -e 'QUAL < 1' -m + {input.tnscopevcf_ml}", f"> {params.outputdir}/{vcfname}_lowqual1.vcf ;",
                            f"{params.bcftools} annotate -x FILTER/triallelic_site {params.outputdir}/{vcfname}_lowqual1.vcf", f"> {params.outputdir}/{vcfname}_triallelic2.vcf ;",
                            f"{params.bcftools} annotate -x FILTER/alt_allele_in_normal {params.outputdir}/{vcfname}_triallelic2.vcf", f"> {params.outputdir}/{vcfname}_altalleleinnormal4.vcf ;",
                            f"{params.bcftools} filter -s uncertainAF -e 'FORMAT/AF[0]<0.045 && FORMAT/AD[0:1]<4' -m + {params.outputdir}/{vcfname}_altalleleinnormal4.vcf", f"> {params.outputdir}/{vcfname}_uncertainaf6.vcf ;",
                            #f"{params.bcftools} filter -s likely_artifact -e 'FORMAT/AF[0]<0.1 && FORMAT/AF[1]>0.06' -m + {params.outputdir}/{vcfname}_uncertainaf6.vcf", f"> {params.outputdir}/{vcfname}_likelyartifact7.vcf ;",
                            f"{params.bcftools} filter -s lowAD -e 'FORMAT/AD[0:1] < 3' {params.outputdir}/{vcfname}_uncertainaf6.vcf", f"> {params.outputdir}/{vcfname}_lowad8.vcf ;",
                            #f"{params.bcftools} filter -s MLrejected -e 'INFO/ML_PROB>0.63' -m + {params.outputdir}/{vcfname}_lowad8.vcf", f"> {params.outputdir}/{vcfname}_mladjusted9.vcf ;",
                            f"{params.bcftools} filter -s 'orientation_bias' -e 'FMT/FOXOG[0] == 1' -m + {params.outputdir}/{vcfname}_lowad8.vcf", f"> {params.outputdir}/{vcfname}_foxogadj.vcf ;",
                            f"{params.bcftools} filter -s 'strand_bias' -e 'SOR > 3' -m + {params.outputdir}/{vcfname}_foxogadj.vcf", f"> {params.outputdir}/{vcfname}_strandbiasadj.vcf ;",
                            f"{params.bcftools} filter -i 'FILTER=\"PASS\"' {params.outputdir}/{vcfname}_strandbiasadj.vcf > {output.somatic_n} ;",
                            f"{params.bcftools} view -s {tumorname} {output.somatic_n} > {output.somatic}"]
            shell_command = " ".join(shell_command)
            print(shell_command)
            shell(shell_command)
            filter_stats_command = ["awk -F'\t' '!/^#/ {{ split($7, filters, \";\");", 
                                    "for (i in filters) count[filters[i]]++; total++; }}",
                                    "END {{ for (f in count) {{ fraction = count[f] / total; printf \"%s\\t%d\\t%.4f\\n\", f, count[f], fraction; }} }}'",
                                    f"{params.outputdir}/{vcfname}_strandbiasadj.vcf | sort > {output.tnscope_filterstats}"]
            filter_stats_command = " ".join(filter_stats_command)
            print(filter_stats_command)
            shell(filter_stats_command)
