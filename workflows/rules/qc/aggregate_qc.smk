# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from workflows.scripts.create_qc_excel import create_excel_main

from pathlib import Path

def get_excel_qc_input(wcs):
    stype = "normal" if normalid else "tumor"
    sname = normalid or tumorid
    inputs = dict(
        ycov = f"{'normal' if normalid else 'tumor'}/reports/{normalid or tumorid}_Ycov.tsv",
        insilico_complete = expand(
            f"{'normal' if normalid else 'tumor'}/insilico/{{insiliconame}}/{normalid or tumorid}.complete",
            insiliconame=insilico_panels
        ),
    )
    if tumorid:
        inputs |= dict(
            tumorcov = f"tumor/reports/{tumorid}_WGScov.tsv",
            tumordedup = f"tumor/dedup/{tumorid}_DEDUP.txt",
            tumorvcf = f"tumor/dnascope/{tumorid}_germline_SNVsOnly.recode.vcf",
        )
    if normalid:
        inputs |= dict(
            normalcov = f"normal/reports/{normalid}_WGScov.tsv",
            normaldedup = f"normal/dedup/{normalid}_DEDUP.txt",
            normalvcf = f"normal/dnascope/{normalid}_germline_SNVsOnly.recode.vcf",
        )
    if tumorid and normalid:
        inputs |= dict(
            canvasvcf = f"tumor/canvas/{tumorid}_CNV_somatic.vcf",
        )
    
    return inputs

rule excel_qc:
    input:
        unpack(get_excel_qc_input)
    output:
        qc_stats = temp("qc_report/{tumorname}_qc_stats.xlsx") if tumorid else temp("qc_report/{normalname}_qc_stats.xlsx")
    params:
        insilicodir = lambda wildcards, input: str(Path(input.insilico_complete[0]).parent)
    run:
        create_excel_main(
            **{k: v for k, v in input.items() if k != "insilico_complete"},
            output = output.qc_stats,
            insilicodir = params.insilicodir
        ) 
