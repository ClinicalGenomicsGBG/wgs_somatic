# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from workflows.scripts.create_qc_excel import create_excel_main
import os

if tumorid:
    if normalid:
        # excel_qc rule for paired analysis
        rule excel_qc:
            input:
                tumorcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                ycov = expand("{stype}/reports/{sname}_Ycov.tsv", stype=sampleconfig[normalname]["stype"], sname=normalid),
                normalcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[normalname]["stype"], sname=normalid),
                tumordedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                normaldedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=sampleconfig[normalname]["stype"], sname=normalid),
                tumorvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                normalvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[normalname]["stype"], sname=normalid),
                canvasvcf = expand("{stype}/canvas/{sname}_CNV_somatic.vcf", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                insilicofile = expand("{stype}/insilico/{insiliconame}/{sname}_{insiliconame}_genes_below10x.xlsx", stype=sampleconfig[normalname]["stype"], insiliconame=sampleconfig["insilico"], sname=normalid),
                tmb = expand("{stype}/reports/{sname}_tmb.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                msi_filtered = expand("{stype}/msi/{sname}_msi_filtered.txt", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                msi = expand("{stype}/msi/{sname}_msi.txt", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                tumor_info_files = expand("{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_info.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid, ploidy=pipeconfig["rules"]["control-freec"]["ploidy"]),
            output:
                temp("qc_report/{tumorname}_qc_stats.xlsx")
            run:
                my_insilicofile = f"{input.insilicofile}".split(" ", 1)[0]
                insilicodir = os.path.dirname(f"{my_insilicofile}").rsplit("/", 1)[0] # Somehow this stuff works
                create_excel_main(tumorcov = f"{input.tumorcov}", ycov = f"{input.ycov}" , normalcov = f"{input.normalcov}", tumordedup = f"{input.tumordedup}", normaldedup = f"{input.normaldedup}", tumorvcf = f"{input.tumorvcf}", normalvcf = f"{input.normalvcf}", canvasvcf = f"{input.canvasvcf}", tmb = f"{input.tmb}", msi = f"{input.msi}", msi_red = f"{input.msi_filtered}", tumor_info_files=[f"{file}" for file in input.tumor_info_files], output = f"{output}", insilicodir = f"{insilicodir}") 
    else:
        # excel_qc rule for tumor only
        rule excel_qc:
            input:
                tumorcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                ycov = expand("{stype}/reports/{sname}_Ycov.tsv", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                tumordedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                tumorvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                insilicofile = expand("{stype}/insilico/{insiliconame}/{sname}_{insiliconame}_genes_below10x.xlsx", stype=sampleconfig[tumorname]["stype"], insiliconame=sampleconfig["insilico"], sname=tumorid),
                tmb = expand("{stype}/reports/{sname}_tmb.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                tumor_info_files = expand("{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_info.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid, ploidy=pipeconfig["rules"]["control-freec"]["ploidy"]),
            output:
                temp("qc_report/{tumorname}_qc_stats.xlsx")
            run:
                my_insilicofile = f"{input.insilicofile}".split(" ", 1)[0]
                insilicodir = os.path.dirname(f"{my_insilicofile}").rsplit("/", 1)[0] # Somehow this stuff works
                create_excel_main(tumorcov = f"{input.tumorcov}", ycov = f"{input.ycov}", tumordedup = f"{input.tumordedup}", tumorvcf = f"{input.tumorvcf}", tmb = f"{input.tmb}", tumor_info_files=[f"{file}" for file in input.tumor_info_files], output = f"{output}", insilicodir = f"{insilicodir}")

else:
    # excel_qc rule for normal only
    rule excel_qc:
        input:
            normalcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[normalname]["stype"], sname=normalid),
            ycov = expand("{stype}/reports/{sname}_Ycov.tsv", stype=sampleconfig[normalname]["stype"], sname=normalid),
            normaldedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=sampleconfig[normalname]["stype"], sname=normalid),
            normalvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[normalname]["stype"], sname=normalid),
            insilicofile = expand("{stype}/insilico/{insiliconame}/{sname}_{insiliconame}_genes_below10x.xlsx", stype=sampleconfig[normalname]["stype"], insiliconame=sampleconfig["insilico"], sname=normalid),
        output:
            temp("qc_report/{normalname}_qc_stats.xlsx")
        run:
            my_insilicofile = f"{input.insilicofile}".split(" ", 1)[0]
            insilicodir = os.path.dirname(f"{my_insilicofile}").rsplit("/", 1)[0] # Somehow this stuff works
            create_excel_main(normalcov = f"{input.normalcov}", ycov = f"{input.ycov}", normaldedup = f"{input.normaldedup}", normalvcf =f"{input.normalvcf}", output = f"{output}", insilicodir = f"{insilicodir}")
