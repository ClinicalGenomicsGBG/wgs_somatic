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
                tmb = expand("{stype}/reports/{sname}_tmb.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                msi_filtered = expand("{stype}/msi/{sname}_msi_filtered.txt", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                msi = expand("{stype}/msi/{sname}_msi.txt", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                tumor_info_files = expand("{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_info.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid, ploidy=pipeconfig["rules"]["control-freec"]["ploidy"]),
            output:
                "qc_report/{tumorname}_qc_stats.xlsx"
            run:
                create_excel_main(tumorcov = f"{input.tumorcov}", ycov = f"{input.ycov}" , normalcov = f"{input.normalcov}", tumordedup = f"{input.tumordedup}", normaldedup = f"{input.normaldedup}", tumorvcf = f"{input.tumorvcf}", normalvcf = f"{input.normalvcf}", canvasvcf = f"{input.canvasvcf}", tmb = f"{input.tmb}", msi = f"{input.msi}", msi_red = f"{input.msi_filtered}", tumor_info_files=[f"{file}" for file in input.tumor_info_files], output = f"{output}", insilicodir = f"{insilicodir}") 
    else:
        # excel_qc rule for tumor only
        rule excel_qc:
            input:
                tumorcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                ycov = expand("{stype}/reports/{sname}_Ycov.tsv", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                tumordedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                tumorvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                tmb = expand("{stype}/reports/{sname}_tmb.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                tumor_info_files = expand("{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_info.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid, ploidy=pipeconfig["rules"]["control-freec"]["ploidy"]),
            output:
                "qc_report/{tumorname}_qc_stats.xlsx"
            run:
                create_excel_main(tumorcov = f"{input.tumorcov}", ycov = f"{input.ycov}", tumordedup = f"{input.tumordedup}", tumorvcf = f"{input.tumorvcf}", tmb = f"{input.tmb}", tumor_info_files=[f"{file}" for file in input.tumor_info_files], output = f"{output}"")

else:
    # excel_qc rule for normal only
    rule excel_qc:
        input:
            normalcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[normalname]["stype"], sname=normalid),
            ycov = expand("{stype}/reports/{sname}_Ycov.tsv", stype=sampleconfig[normalname]["stype"], sname=normalid),
            normaldedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=sampleconfig[normalname]["stype"], sname=normalid),
            normalvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[normalname]["stype"], sname=normalid),
        output:
            "qc_report/{normalname}_qc_stats.xlsx"
        run:
            create_excel_main(normalcov = f"{input.normalcov}", ycov = f"{input.ycov}", normaldedup = f"{input.normaldedup}", normalvcf =f"{input.normalvcf}", output = f"{output}")
