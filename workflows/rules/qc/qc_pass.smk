rule qc_pass:
    input:
        sname="{sname}"
	coverage="{stype}/reports/{sname}_WGScov.tsv",
        somalier=""
    output:
        passfile="{stype}/reports/{sname}/{sname}.QC_PASS"
    run:
        from workflows.scripts.qc_eval import evaluate_qc
        evaluate_qc(input.sname,
                    input.coverage,
                    input.somalier,
                    output.passfile)
