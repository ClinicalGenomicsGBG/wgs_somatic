# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from workflows.scripts.b_allele_igv_plot.plot_b_allele_freq import plot_freq

rule ballele_plot:
    input:
        "{stype}/dnascope/{sname}_germline.vcf"
    params:
        dbsnp = pipeconfig["rules"]["ballele_plot"]["dbsnp"],
        hg38ref = pipeconfig["rules"]["ballele_plot"]["hg38ref"]
    shadow:
        pipeconfig["rules"].get("ballele_plot", {}).get("shadow", pipeconfig.get("shadow", False))
    output:
        temp("{stype}/reports/{sname}_baf.igv")
    run:
        plot_freq(f"{input}", f"{output}", f"{params.dbsnp}", f"{params.hg38ref}")
