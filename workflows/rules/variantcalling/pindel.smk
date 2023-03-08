# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

if normalid:
    rule pindelConfig:
        input:
            tumorbam = expand("{stype}/dedup/{sname}_DEDUP.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            normalbam = expand("{stype}/dedup/{sname}_DEDUP.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        output:
            pindelConfig = "{stype}/pindel/{sname}_pindelConfig.txt"
        shell:
            "echo '{input.tumorbam}\t300\t{tumorname}\n{input.normalbam}\t300\t{normalname}'>{output.pindelConfig}"
else:
    rule pindelConfig:
        input:
            tumorbam = expand("{stype}/dedup/{sname}_DEDUP.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        output:
            pindelConfig = "{stype}/pindel/{sname}_pindelConfig.txt"
        shell:
            "echo '{input.tumorbam}\t300\t{tumorname}'>{output.pindelConfig}"
if normalid:
    rule pindel:
        input:
            tumorbam = expand("{stype}/dedup/{sname}_DEDUP.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            normalbam = expand("{stype}/dedup/{sname}_DEDUP.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            pindelConfig = expand("{stype}/pindel/{sname}_pindelConfig.txt", sname=tumorid, stype=sampleconfig[tumorname]["stype"])
        params:
            bed = pipeconfig["rules"]["pindel"]["bed"],
            reference = pipeconfig["referencegenome"],
            threads = clusterconf["pindel"]["threads"],
            x = 2,
            B = 60
        singularity:
            pipeconfig["singularities"]["pindel"]["sing"]
        output:
            "{stype}/pindel/{sname}_BP",
            "{stype}/pindel/{sname}_CloseEndMapped",
            "{stype}/pindel/{sname}_D",
            "{stype}/pindel/{sname}_INT_final",
            "{stype}/pindel/{sname}_INV",
            "{stype}/pindel/{sname}_LI",
            "{stype}/pindel/{sname}_RP",
            "{stype}/pindel/{sname}_SI",
            "{stype}/pindel/{sname}_TD"
        shell:
            "echo $HOSTNAME;"
            " (pindel -f {params.reference} -i {input.pindelConfig} -T {params.threads} -x {params.x} -B {params.B} -j {params.bed} -o {wildcards.stype}/pindel/{wildcards.sname} ) "

else:
    rule pindel:
        input:
            tumorbam = expand("{stype}/dedup/{sname}_DEDUP.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            pindelConfig = expand("{stype}/pindel/{sname}_pindelConfig.txt", sname=tumorid, stype=sampleconfig[tumorname]["stype"])
        params:
            bed = pipeconfig["rules"]["pindel"]["bed"],
            reference = pipeconfig["referencegenome"],
            threads = clusterconf["pindel"]["threads"],
            x = 2,
            B = 60
        singularity:
            pipeconfig["singularities"]["pindel"]["sing"]
        output:
            "{stype}/pindel/{sname}_BP",
            "{stype}/pindel/{sname}_CloseEndMapped",
            "{stype}/pindel/{sname}_D",
            "{stype}/pindel/{sname}_INT_final",
            "{stype}/pindel/{sname}_INV",
            "{stype}/pindel/{sname}_LI",
            "{stype}/pindel/{sname}_RP",
            "{stype}/pindel/{sname}_SI",
            "{stype}/pindel/{sname}_TD"
        shell:
            "echo $HOSTNAME;"
            " (pindel -f {params.reference} -i {input.pindelConfig} -T {params.threads} -x {params.x} -B {params.B} -j {params.bed} -o {wildcards.stype}/pindel/{wildcards.sname} ) "

rule pindel2vcf:
    input:
        expand("{stype}/pindel/{sname}_BP", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/pindel/{sname}_CloseEndMapped", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/pindel/{sname}_D", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/pindel/{sname}_INT_final", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/pindel/{sname}_INV", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/pindel/{sname}_LI", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/pindel/{sname}_RP", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/pindel/{sname}_SI", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/pindel/{sname}_TD", sname=tumorid, stype=sampleconfig[tumorname]["stype"])
    params:
        threads = clusterconf["pindel"]["threads"],
        reference = pipeconfig["referencegenome"],
        refname = "GRCh38",
        refdate = 000000,
        e = 3, #e = 10,
        mc = 10,
        minsize = 5
    singularity:
        pipeconfig["singularities"]["pindel"]["sing"]
    output:
        "{stype}/pindel/{sname}_pindel_noDP_noContig.vcf"
    shell:
        "echo $HOSTNAME;"
        " (pindel2vcf -P {wildcards.stype}/pindel/{wildcards.sname} -r {params.reference} -R {params.refname} -d {params.refdate} -v {output} -e {params.e} -mc {params.mc} -G -is {params.minsize}) "

rule fixContigPindel:
    input:
        "{stype}/pindel/{sname}_pindel_noDP_noContig.vcf"
    output:
        "{stype}/pindel/{sname}_pindel_noDP.vcf"
    params: 
        referencefai = pipeconfig["referencefai"]
    shell:
        """/apps/bio/software/anaconda2/envs/wgs_somatic/bin/bcftools reheader --fai {params.referencefai} {input} > {output}"""

rule fixPindelDPoAF:
    input:
        "{stype}/pindel/{sname}_pindel_noDP.vcf"
    output:
        "{stype}/pindel/{sname}_pindel.vcf"
    params:
        python = pipeconfig["rules"]["pindel"]["python"],
        fix_DPoAF = pipeconfig["rules"]["pindel"].get("fix_DPoAF", f"{ROOT_DIR}/workflows/scripts/fix_pindelDPoAF.py") 
    run:
        shell(f"{params.python} {params.fix_DPoAF} {input} {output}")

rule pindel_xlsx:
    input:
        "{stype}/pindel/{sname}_pindel.vcf"
    output:
        "{stype}/pindel/{sname}_pindel.xlsx"
    params:
        python = pipeconfig["rules"]["pindel"]["python"],
        bed = pipeconfig["rules"]["pindel"]["bed"],
        pindel_excel = pipeconfig["rules"]["pindel"].get("pindel_excel", f"{ROOT_DIR}/workflows/scripts/pindel_excel.py")
    run:
        shell(f"{params.python} {params.pindel_excel} {input} {output} {params.bed}")
