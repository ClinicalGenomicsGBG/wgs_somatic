# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from workflows.scripts.create_qc_excel import create_excel_main, create_qc_toaggregate
from workflows.scripts.parse_somalier import SomalierParser
import os

# set name of output qc files taking into account if tumorid or normalid exist
if tumorid:
    excel_qc_output = "qc_report/{tumorid}_qc_stats.xlsx"
    qcstats_wgs_admin_output = "qc_report/{tumorid}_qc_stats_wgsadmin.xlsx"
else:
    excel_qc_output = "qc_report/{normalid}_qc_stats.xlsx"
    qcstats_wgs_admin_output = "qc_report/{normalid}_qc_stats_wgsadmin.xlsx"


rule excel_qc:
    input: # the empty input is given as an empty list and not string because snakemake interprets "" as a missing file, instead of "no input"
        tumorcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=stype_tumor, sname=tumorid) if tumorid else [],
        normalcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=stype_normal, sname=normalid) if normalid else [],
        tumordedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=stype_tumor, sname=tumorid) if tumorid else [],
        normaldedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=stype_normal, sname=normalid) if normalid else [],
        somalier_pairs = expand("{stype}/somalier/somalier.pairs.tsv", stype=stype_normal if normalid else stype_tumor),
        somalier_samples = expand("{stype}/somalier/somalier.samples.tsv", stype=stype_normal if normalid else stype_tumor),
        canvasvcf = expand("{stype}/canvas/{sname}_canvas_somatic.vcf.gz", stype=stype_tumor, sname=tumorid) if tumorid and normalid else [],
        tmb = expand("{stype}/reports/{sname}_tmb.txt", stype=stype_tumor, sname=tumorid) if tumorid else [],
        msi_filtered = expand("{stype}/msi/{sname}_msi_filtered.txt", sname=tumorid, stype=stype_tumor) if tumorid and normalid else [],
        msi = expand("{stype}/msi/{sname}_msi.txt", sname=tumorid, stype=stype_tumor) if tumorid and normalid else [],
        ascatstats = expand("{stype}/ascat/{sname}_ascat_stats.tsv", sname=tumorid, stype=stype_tumor) if tumorid else [],
    output:
        excel_qc_output
    run:
        # Parse somalier data
        somalier = SomalierParser(
            pairs_file=f"{input.somalier_pairs}",
            samples_file=f"{input.somalier_samples}",
            tumorstring=stype_tumor,
            normalstring=stype_normal,
            match_cutoff=filterconfig["somalier_filter"]["min_relatedness"]
        )
        # Call create_excel_main with all inputs
        create_excel_main(
            tumorcov=f"{input.tumorcov}",
            normalcov=f"{input.normalcov}",
            tumordedup=f"{input.tumordedup}",
            normaldedup=f"{input.normaldedup}",
            somalier_obj=somalier,
            canvasvcf=f"{input.canvasvcf}",
            tmb=f"{input.tmb}",
            msi=f"{input.msi}",
            msi_red=f"{input.msi_filtered}",
            output=f"{output}",
            ascatstats=f"{input.ascatstats}"
        )

rule qcstats_wgs_admin:
# Rule to create a QC report for WGS admin
# works for tumor+normal, tumor-only and normal-only
    input: 
        tumorcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=stype_tumor, sname=tumorid) if tumorid else [],
        normalcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=stype_normal, sname=normalid) if normalid else [],
        tumordedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=stype_tumor, sname=tumorid) if tumorid else [],
        normaldedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=stype_normal, sname=normalid) if normalid else [],
        somalier_pairs = expand("{stype}/somalier/somalier.pairs.tsv", stype=stype_normal if normalid else stype_tumor),
        somalier_samples = expand("{stype}/somalier/somalier.samples.tsv", stype=stype_normal if normalid else stype_tumor),
        tmb = expand("{stype}/reports/{sname}_tmb.txt", stype=stype_tumor, sname=tumorid) if tumorid else [],
    output: 
        report(qcstats_wgs_admin_output)
    params:
        tumorid = tumorid if tumorid else '',
        normalid = normalid if normalid else '',
    run:
        # Parse somalier data
        somalier = SomalierParser(
            pairs_file=f"{input.somalier_pairs}",
            samples_file=f"{input.somalier_samples}",
            tumorstring=stype_tumor,
            normalstring=stype_normal,
            match_cutoff=filterconfig["somalier_filter"]["min_relatedness"]
        )
        create_qc_toaggregate(
            output=f"{output}",
            somalier_obj=somalier,
            tumorcov=f"{input.tumorcov}",
            normalcov=f"{input.normalcov}" ,
            tumordedup=f"{input.tumordedup}" ,
            normaldedup=f"{input.normaldedup}" ,
            tmb=f"{input.tmb}" ,
            tumorid = f"{params.tumorid}",
            normalid = f"{params.normalid}"
        )
