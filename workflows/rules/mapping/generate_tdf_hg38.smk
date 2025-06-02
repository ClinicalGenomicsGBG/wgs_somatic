# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
from tools.helpers import conditional_temp

rule generate_tdf_hg38:
    input:
        bam = "{stype}/realign/{sname}_REALIGNED.bam",
        bai = "{stype}/realign/{sname}_REALIGNED.bam.bai"
    params:
        igvtools_jar_path = pipeconfig["rules"]["generate_tdf"]["igvtools_jar_path"],
        igvtools_memory_limit = pipeconfig["rules"]["generate_tdf"]["igvtools_memory_limit"]
    shadow:
        pipeconfig["rules"].get("generate_tdf", {}).get("shadow", pipeconfig.get("shadow", False))
    output:
        conditional_temp("{stype}/reports/{sname}_REALIGNED.bam.tdf", keepfiles)
    run:
        shell("nohup java -Xmx{params.igvtools_memory_limit} -jar {params.igvtools_jar_path} count {input.bam} {output} hg38")
