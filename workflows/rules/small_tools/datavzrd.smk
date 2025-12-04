rule datavzrd_report:
    input:
        admin_qc="qc_report/{sname}_qc_stats_wgsadmin.tsv",
        ascatstats="{stype}/ascat/{sname}_ascat_stats.tsv",
        tmb="{stype}/reports/{sname}_tmb.tsv",
        msi="{stype}/msi/{sname}_msi.tsv",
        config="configs/datavzrd_config.yaml",
    output:
        report(
            directory("{stype}/reports/{sname}_datavzrd_report"),
            htmlindex="index.html",
        ),
    singularity:
        pipeconfig["singularities"]["datavzrd"]["sing"]
    wrapper:
        "v8.0.3/utils/datavzrd"
