# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
from workflows.scripts.qc_eval import evaluate_qc

rule qc_pass:
    params:
        sname = "{sname}"
    input:
        coverage = "{stype}/reports/{sname}_WGScov.tsv",
        somalier = "{stype}/somalier/calculated_sex.txt"
    output:
        passfile = "{stype}/reports/{sname}/{sname}.QC_PASS"
    run:
        evaluate_qc(params.sname,
                    input.coverage,
                    input.somalier,
                    output.passfile)
