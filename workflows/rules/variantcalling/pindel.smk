# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

if normalid:
    rule pindel:
        input:
            tumorbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            tumorbai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            normalbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            normalbai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        params:
            bed = pipeconfig["rules"]["pindel"]["bed"],
            reference = pipeconfig["referencegenome"],
            threads = clusterconf["pindel"]["threads"],
            x = 2,
            B = 60
        singularity:
            pipeconfig["singularities"]["pindel"]["sing"]
        output:
            temp("{stype}/pindel/{sname}_BP"),
            temp("{stype}/pindel/{sname}_CloseEndMapped"),
            temp("{stype}/pindel/{sname}_D"),
            temp("{stype}/pindel/{sname}_INT_final"),
            temp("{stype}/pindel/{sname}_INV"),
            temp("{stype}/pindel/{sname}_LI"),
            temp("{stype}/pindel/{sname}_RP"),
            temp("{stype}/pindel/{sname}_SI"),
            temp("{stype}/pindel/{sname}_TD"),
            pindelConfig = temp("{stype}/pindel/{sname}_pindelConfig.txt")
        shadow:
            pipeconfig["rules"].get("pindel", {}).get("shadow", pipeconfig.get("shadow", False))
        shell:
            "echo $HOSTNAME; "
            "echo '{input.tumorbam}\t300\t{tumorname}\n{input.normalbam}\t300\t{normalname}'>{output.pindelConfig}; "
            "(pindel -f {params.reference} -i {output.pindelConfig} -T {params.threads} -x {params.x} -B {params.B} -j {params.bed} -o {wildcards.stype}/pindel/{wildcards.sname} ) "

else:
    rule pindel:
        input:
            tumorbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            tumorbai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        params:
            bed = pipeconfig["rules"]["pindel"]["bed"],
            reference = pipeconfig["referencegenome"],
            threads = clusterconf["pindel"]["threads"],
            x = 2,
            B = 60
        singularity:
            pipeconfig["singularities"]["pindel"]["sing"]
        output:
            temp("{stype}/pindel/{sname}_BP"),
            temp("{stype}/pindel/{sname}_CloseEndMapped"),
            temp("{stype}/pindel/{sname}_D"),
            temp("{stype}/pindel/{sname}_INT_final"),
            temp("{stype}/pindel/{sname}_INV"),
            temp("{stype}/pindel/{sname}_LI"),
            temp("{stype}/pindel/{sname}_RP"),
            temp("{stype}/pindel/{sname}_SI"),
            temp("{stype}/pindel/{sname}_TD"),
            pindelConfig = temp("{stype}/pindel/{sname}_pindelConfig.txt")
        shadow:
            pipeconfig["rules"].get("pindel", {}).get("shadow", pipeconfig.get("shadow", False))
        shell:
            "echo $HOSTNAME; "
            "echo '{input.tumorbam}\t300\t{tumorname}'>{output.pindelConfig}; "
            "(pindel -f {params.reference} -i {output.pindelConfig} -T {params.threads} -x {params.x} -B {params.B} -j {params.bed} -o {wildcards.stype}/pindel/{wildcards.sname} ) "

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
        e = 1, #3, #e = 10,
        mc = 10,
        minsize = 5
    singularity:
        pipeconfig["singularities"]["pindel"]["sing"]
    output:
        temp("{stype}/pindel/{sname}_pindel_noDP_noContig.vcf")
    shadow:
        pipeconfig["rules"].get("pindel", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "echo $HOSTNAME;"
        " (pindel2vcf -P {wildcards.stype}/pindel/{wildcards.sname} -r {params.reference} -R {params.refname} -d {params.refdate} -v {output} -e {params.e} -mc {params.mc} -G -is {params.minsize}) "

rule fixContigPindel:
    input:
        "{stype}/pindel/{sname}_pindel_noDP_noContig.vcf"
    output:
        temp("{stype}/pindel/{sname}_pindel_noDP.vcf")
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
    shadow:
        pipeconfig["rules"].get("pindel", {}).get("shadow", pipeconfig.get("shadow", False))
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
    shadow:
        pipeconfig["rules"].get("pindel", {}).get("shadow", pipeconfig.get("shadow", False))
    run:
        shell(f"{params.python} {params.pindel_excel} {input} {output} {params.bed}")
