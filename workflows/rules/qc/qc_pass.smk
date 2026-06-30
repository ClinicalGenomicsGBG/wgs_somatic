from workflows.scripts.qc_eval import evaluate_qc

rule qc_pass:
    input:
        coverage = "{stype}/reports/{sname}_WGScov.tsv",
        somalier = "{stype}/somalier/calculated_sex.txt"
    output:
        passfile = "{stype}/reports/{sname}.QC_PASS"
    run:
        evaluate_qc(wildcards.stype,
                    wildcards.sname,
                    input.coverage,
                    input.somalier,
                    output.passfile,
                    gender)
