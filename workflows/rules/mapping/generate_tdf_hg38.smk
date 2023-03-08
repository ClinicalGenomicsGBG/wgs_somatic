# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule generate_tdf_hg38:
    input:
        "{stype}/realign/{sname}_REALIGNED.bam"
    params:
        igvtools_jar_path = pipeconfig["rules"]["generate_tdf"]["igvtools_jar_path"],
        igvtools_memory_limit = pipeconfig["rules"]["generate_tdf"]["igvtools_memory_limit"]
    output:
        "{stype}/reports/{sname}_REALIGNED.bam.tdf"
    run:
        shell("nohup java -Xmx{params.igvtools_memory_limit} -jar {params.igvtools_jar_path} count {input} {output} hg38")
