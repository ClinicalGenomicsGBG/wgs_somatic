# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from workflows.scripts.create_qc_excel import create_excel_main, create_qc_toaggregate
import os

# set name of output qc files taking into account if tumorid or normalid exist
if tumorid:
    excel_qc_output = "qc_report/{tumorname}_qc_stats.xlsx"
    qcstats_wgs_admin_output = "qc_report/{tumorname}_qc_stats_wgsadmin.xlsx"
else:
    excel_qc_output = "qc_report/{normalname}_qc_stats.xlsx"
    qcstats_wgs_admin_output = "qc_report/{normalname}_qc_stats_wgsadmin.xlsx"


rule excel_qc:
    input: # the empty input is given as an empty list and not string because snakemake interprets "" as a missing file, instead of "no input"
        tumorcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid else [],
        ycov = expand("{stype}/reports/{sname}_Ycov.tsv", stype=sampleconfig[(normalname if normalid else tumorname)]["stype"], sname=(normalid if normalid else tumorid)),
        normalcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[normalname]["stype"], sname=normalid) if normalid else [],
        tumordedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid else [],
        normaldedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=sampleconfig[normalname]["stype"], sname=normalid) if normalid else [],
        tumorvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid else [],
        normalvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[normalname]["stype"], sname=normalid) if normalid else [],
        canvasvcf = expand("{stype}/canvas/{sname}_CNV_somatic.vcf", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid and normalid else [],
        tmb = expand("{stype}/reports/{sname}_tmb.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid else [],
        msi_filtered = expand("{stype}/msi/{sname}_msi_filtered.txt", sname=tumorid, stype=sampleconfig[tumorname]["stype"]) if tumorid and normalid else [],
        msi = expand("{stype}/msi/{sname}_msi.txt", sname=tumorid, stype=sampleconfig[tumorname]["stype"]) if tumorid and normalid else [],
        tumor_info_files = expand("{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_info.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid, ploidy=pipeconfig["rules"]["control-freec"]["ploidy"]),
    output:
        excel_qc_output
    run:
        # Call create_excel_main with all inputs
        create_excel_main(
            tumorcov=f"{input.tumorcov}",
            ycov=f"{input.ycov}",
            normalcov=f"{input.normalcov}",
            tumordedup=f"{input.tumordedup}",
            normaldedup=f"{input.normaldedup}",
            tumorvcf=f"{input.tumorvcf}",
            normalvcf=f"{input.normalvcf}",
            canvasvcf=f"{input.canvasvcf}",
            tmb=f"{input.tmb}",
            msi=f"{input.msi}",
            msi_red=f"{input.msi_filtered}",
            output=f"{output}",
            tumor_info_files=[file for file in input.tumor_info_files]
        )

rule qcstats_wgs_admin:
# Rule to create a QC report for WGS admin
# works for tumor+normal, tumor-only and normal-only
    input: 
        tumorcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid else [],
        ycov = expand("{stype}/reports/{sname}_Ycov.tsv", stype=sampleconfig[(normalname if normalid else tumorname)]["stype"], sname=(normalid if normalid else tumorid)),
        normalcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[normalname]["stype"], sname=normalid) if normalid else [],
        tumordedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid else [],
        normaldedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=sampleconfig[normalname]["stype"], sname=normalid) if normalid else [],
        tumorvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid else [],
        normalvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[normalname]["stype"], sname=normalid) if normalid else [],
        tmb = expand("{stype}/reports/{sname}_tmb.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid else [],
    output: 
        qcstats_wgs_admin_output
    params:
        tumorid = tumorid if tumorid else '',
        normalid = normalid if normalid else ''
    run: 
        create_qc_toaggregate(
        tumorcov=f"{input.tumorcov}",
        ycov=f"{input.ycov}",
        normalcov=f"{input.normalcov}" ,
        tumordedup=f"{input.tumordedup}" ,
        normaldedup=f"{input.normaldedup}" ,
        tumorvcf=f"{input.tumorvcf}" ,
        normalvcf=f"{input.normalvcf}" ,
        tmb=f"{input.tmb}" ,
        output=f"{output}",
        tumorid = f"{params.tumorid}",
        normalid = f"{params.normalid}"
        )
