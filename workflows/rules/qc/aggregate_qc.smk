# # vim: syntax=python tabstop=4 expandtab
# # coding: utf-8

from workflows.scripts.create_qc_excel import create_excel_main, create_qc_toaggregate
import os

# set name of output qc files taking into account if tumorid or normalid exist
if tumorid:
    excel_qc_output = temp("qc_report/{tumorname}_qc_stats.xlsx")
    qcstats_to_aggregate_output = temp("qc_report/{tumorname}_qc_stats_wgsadmin.xlsx")
else:
    excel_qc_output = temp("qc_report/{normalname}_qc_stats.xlsx")
    qcstats_to_aggregate_output = temp("qc_report/{normalname}_qc_stats_wgsadmin.xlsx")


rule excel_qc:
    input: # the empty input is given as an empty list and not string because snakemake interprets "" as a missing file, instead of "no input"
        tumorcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid else [],
        ycov = expand("{stype}/reports/{sname}_Ycov.tsv", stype=sampleconfig[normalname]["stype"], sname=(normalid if normalid else tumorid)),
        normalcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[normalname]["stype"], sname=normalid) if normalid else [],
        tumordedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid else [],
        normaldedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=sampleconfig[normalname]["stype"], sname=normalid) if normalid else [],
        tumorvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid else [],
        normalvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[normalname]["stype"], sname=normalid) if normalid else [],
        canvasvcf = expand("{stype}/canvas/{sname}_CNV_somatic.vcf", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid and normalid else [],
        insilicofile = expand("{stype}/insilico/{insiliconame}/{sname}_{insiliconame}_genes_below10x.xlsx", stype=sampleconfig[normalname]["stype"], insiliconame=sampleconfig["insilico"], sname=(normalid if normalid else tumorid)),
        tmb = expand("{stype}/reports/{sname}_tmb.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid else [],
        msi_filtered = expand("{stype}/msi/{sname}_msi_filtered.txt", sname=tumorid, stype=sampleconfig[tumorname]["stype"]) if tumorid and normalid else [],
        msi = expand("{stype}/msi/{sname}_msi.txt", sname=tumorid, stype=sampleconfig[tumorname]["stype"]) if tumorid and normalid else [],
        tumor_info_files = expand("{stype}/control-freec_{ploidy}/{sname}.pileup_info.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid, ploidy=pipeconfig["rules"]["control-freec"]["ploidy"]) if tumorid else [],
    output:
        excel_qc_output
    run:
        # Determine the insilico directory if applicable
        my_insilicofile = f"{input.insilicofile}".split(" ", 1)[0] if "insilicofile" in input else None
        insilicodir = os.path.dirname(f"{my_insilicofile}").rsplit("/", 1)[0] if my_insilicofile else None

        # Call create_excel_main with all inputs
        create_excel_main(
            tumorcov=f"{input.tumorcov}" if input.tumorcov else "",
            ycov=f"{input.ycov}" if input.ycov else "",
            normalcov=f"{input.normalcov}" if input.normalcov else "",
            tumordedup=f"{input.tumordedup}" if input.tumordedup else "",
            normaldedup=f"{input.normaldedup}" if input.normaldedup else "",
            tumorvcf=f"{input.tumorvcf}" if input.tumorvcf else "",
            normalvcf=f"{input.normalvcf}" if input.normalvcf else "",
            canvasvcf=f"{input.canvasvcf}" if input.canvasvcf else "",
            tmb=f"{input.tmb}" if input.tmb else "",
            msi=f"{input.msi}" if input.msi else "",
            msi_red=f"{input.msi_filtered}" if input.msi_filtered else "",
            output=f"{output}",
            insilicodir=f"{insilicodir}" if insilicodir else "",
            tumor_info_files=[f"{file}" for file in input.tumor_info_files]
        )

rule qcstats_to_aggregate:
    input: 
        tumorcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid else [],
        ycov = expand("{stype}/reports/{sname}_Ycov.tsv", stype=sampleconfig[normalname]["stype"], sname=(normalid if normalid else tumorid)),
        normalcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[normalname]["stype"], sname=normalid) if normalid else [],
        tumordedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid else [],
        normaldedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=sampleconfig[normalname]["stype"], sname=normalid) if normalid else [],
        tumorvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid else [],
        normalvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[normalname]["stype"], sname=normalid) if normalid else [],
        tmb = expand("{stype}/reports/{sname}_tmb.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid) if tumorid else [],
    output: 
        qcstats_to_aggregate_output
    run: 
        create_qc_toaggregate(
        tumorcov=f"{input.tumorcov}" if input.tumorcov else "",
        ycov=f"{input.ycov}" if input.ycov else "",
        normalcov=f"{input.normalcov}" if input.normalcov else "",
        tumordedup=f"{input.tumordedup}" if input.tumordedup else "",
        normaldedup=f"{input.normaldedup}" if input.normaldedup else "",
        tumorvcf=f"{input.tumorvcf}" if input.tumorvcf else "",
        normalvcf=f"{input.normalvcf}" if input.normalvcf else "",
        tmb=f"{input.tmb}" if input.tmb else "",
        output=f"{output}"
        )
