# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

def get_pindel_input(wcs):
    inputs = dict(
        tumorbam = f"tumor/dedup/{tumorid}_DEDUP.bam",
        tumorbai = f"tumor/dedup/{tumorid}_DEDUP.bam.bai",
    )
    if normalid:
        inputs |= dict(
            normalbam = f"normal/dedup/{normalid}_DEDUP.bam",
            normalbai = f"normal/dedup/{normalid}_DEDUP.bam.bai",
        )
    return inputs

def get_pindel_config_fmt(wcs, input):
    if normalid:
        return f"{input.tumorbam}\\t300\\t{tumorname}\\n{input.normalbam}\\t300\\t{normalname}"
    else:
        return f"{input.tumorbam}\\t300\\t{tumorname}"


rule pindel:
    input:
        unpack(get_pindel_input)
    params:
        config_fmt = get_pindel_config_fmt,
        bed = pipeconfig["rules"]["pindel"]["bed"],
        reference = pipeconfig["referencegenome"],
        threads = clusterconf["pindel"]["threads"],
        x = 2,
        B = 60
    singularity:
        pipeconfig["singularities"]["pindel"]["simg"]
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
        config = temp("{stype}/pindel/{sname}_pindelConfig.txt")
    wildcard_constraints:
        stype = "tumor"
    shadow:
        pipeconfig["rules"].get("pindel", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "echo $HOSTNAME; "
        "echo -e '{params.config_fmt}' > {output.config}; "
        "pindel "
            "-f {params.reference} "
            "-i {output.config} "
            "-T {params.threads} "
            "-x {params.x} "
            "-B {params.B} "
            "-j {params.bed} "
            "-o {wildcards.stype}/pindel/{wildcards.sname}"


rule pindel2vcf:
    input:
        f"tumor/pindel/{tumorid}_BP",
        f"tumor/pindel/{tumorid}_CloseEndMapped",
        f"tumor/pindel/{tumorid}_D",
        f"tumor/pindel/{tumorid}_INT_final",
        f"tumor/pindel/{tumorid}_INV",
        f"tumor/pindel/{tumorid}_LI",
        f"tumor/pindel/{tumorid}_RP",
        f"tumor/pindel/{tumorid}_SI",
        f"tumor/pindel/{tumorid}_TD",
    params:
        threads = clusterconf["pindel"]["threads"],
        reference = pipeconfig["referencegenome"],
        refname = "GRCh38",
        refdate = 000000,
        e = 3, #e = 10,
        mc = 10,
        minsize = 5
    singularity:
        pipeconfig["singularities"]["pindel"]["simg"]
    output:
        temp("{stype}/pindel/{sname}_pindel_noDP_noContig.vcf")
    wildcard_constraints:
        stype = "tumor"
    shadow:
        pipeconfig["rules"].get("pindel", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "echo $HOSTNAME; "
        "pindel2vcf "
            "-P {wildcards.stype}/pindel/{wildcards.sname} "
            "-r {params.reference} "
            "-R {params.refname} "
            "-d {params.refdate} "
            "-v {output} "
            "-e {params.e} "
            "-mc {params.mc} "
            "-G "
            "-is {params.minsize}"

rule fixContigPindel:
    input:
        "tumor/pindel/{sname}_pindel_noDP_noContig.vcf"
    output:
        temp("{stype}/pindel/{sname}_pindel_noDP.vcf")
    params: 
        referencefai = pipeconfig["referencefai"],
        bcftools = pipeconfig["rules"]["pindel"]["bcftools"]
    wildcard_constraints:
        stype = "tumor"
    shadow:
        pipeconfig["rules"].get("pindel", {}).get("shadow", pipeconfig.get("shadow", False))    
    shell:
        "{params.bcftools} reheader "
            "--fai {params.referencefai} "
            "{input} > {output}"

rule fixPindelDPoAF:
    input:
        "tumor/pindel/{sname}_pindel_noDP.vcf"
    output:
        temp("{stype}/pindel/{sname}_pindel.vcf")
    params:
        python = pipeconfig["rules"]["pindel"]["python"],
        fix_DPoAF = pipeconfig["rules"]["pindel"].get("fix_DPoAF", f"{ROOT_DIR}/workflows/scripts/fix_pindelDPoAF.py") 
    wildcard_constraints:
        stype = "tumor"
    shadow:
        pipeconfig["rules"].get("pindel", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "{params.python} {params.fix_DPoAF} {input} {output}"

rule pindel_xlsx:
    input:
        "tumor/pindel/{sname}_pindel.vcf"
    output:
        temp("{stype}/pindel/{sname}_pindel.xlsx")
    params:
        python = pipeconfig["rules"]["pindel"]["python"],
        bed = pipeconfig["rules"]["pindel"]["bed"],
        pindel_excel = pipeconfig["rules"]["pindel"].get("pindel_excel", f"{ROOT_DIR}/workflows/scripts/pindel_excel.py")
    wildcard_constraints:
        stype = "tumor"
    shadow:
        pipeconfig["rules"].get("pindel", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "{params.python} {params.pindel_excel} {input} {output} {params.bed}"
