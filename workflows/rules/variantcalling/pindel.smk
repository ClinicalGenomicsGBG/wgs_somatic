# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from workflows.scripts.pindel_util import pindel_fix_dpaf, pindel_write_xlsx


if normalid:
    rule pindel:
        input:
            tumorbam = expand("{stype}/dedup/{sname}_DEDUP.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            tumorbai = expand("{stype}/dedup/{sname}_DEDUP.bam.bai", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            normalbam = expand("{stype}/dedup/{sname}_DEDUP.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            normalbai = expand("{stype}/dedup/{sname}_DEDUP.bam.bai", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        params:
            bed = pipeconfig["rules"]["pindel"]["bed"],
            reference = pipeconfig["referencegenome"],
            threads = clusterconf["pindel"]["threads"],
            x = 2,
            B = 60,
            vstamp = f"{VDIR}/pindel.txt"
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
            """
            # Version info
            (pindel | sed -n 2p > {params.vstamp}) || true

            # Configure and run pindel
            echo $HOSTNAME
            echo '{input.tumorbam}\t300\t{tumorname}\n{input.normalbam}\t300\t{normalname}'>{output.pindelConfig}
            (pindel -f {params.reference} -i {output.pindelConfig} -T {params.threads} -x {params.x} -B {params.B} -j {params.bed} -o {wildcards.stype}/pindel/{wildcards.sname} )
            """

else:
    rule pindel:
        input:
            tumorbam = expand("{stype}/dedup/{sname}_DEDUP.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            tumorbai = expand("{stype}/dedup/{sname}_DEDUP.bam.bai", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        params:
            bed = pipeconfig["rules"]["pindel"]["bed"],
            reference = pipeconfig["referencegenome"],
            threads = clusterconf["pindel"]["threads"],
            x = 2,
            B = 60,
            vstamp = f"{VDIR}/pindel.txt"
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
            """
            # Version info
            (pindel | sed -n 2p > {params.vstamp}) || true

            # Configure and run pindel
            echo $HOSTNAME
            echo '{input.tumorbam}\t300\t{tumorname}'>{output.pindelConfig}
            (pindel -f {params.reference} -i {output.pindelConfig} -T {params.threads} -x {params.x} -B {params.B} -j {params.bed} -o {wildcards.stype}/pindel/{wildcards.sname} )
            """

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
        minsize = 5,
        vstamp = f"{VDIR}/pindel2vcf.txt"
    singularity:
        pipeconfig["singularities"]["pindel"]["sing"]
    output:
        temp("{stype}/pindel/{sname}_pindel_noDP_noContig.vcf")
    shadow:
        pipeconfig["rules"].get("pindel", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        """
        echo pindel2vcf $(pindel2vcf --help | grep Version) > {params.vstamp}
        echo $HOSTNAME
        (pindel2vcf -P {wildcards.stype}/pindel/{wildcards.sname} -r {params.reference} -R {params.refname} -d {params.refdate} -v {output} -e {params.e} -mc {params.mc} -G -is {params.minsize})
        """

rule pindel_fix_contig:
    input:
        "{stype}/pindel/{sname}_pindel_noDP_noContig.vcf"
    output:
        temp("{stype}/pindel/{sname}_pindel_noDP.vcf")
    params: 
        referencefai = pipeconfig["referencefai"],
        vstamp = f"{VDIR}/pindel_fix_contig.txt",
        bcftools = pipeconfig["rules"]["pindel_fix_contig"]["bcftools"]
    shell:
        """
        {params.bcftools} -v | head -n 2 > {params.vstamp}
        {params.bcftools} reheader --fai {params.referencefai} {input} > {output}
        """

rule pindel_fix_dpaf:
    input:
        "{stype}/pindel/{sname}_pindel_noDP.vcf"
    output:
        "{stype}/pindel/{sname}_pindel.vcf"
    shadow:
        pipeconfig["rules"].get("pindel", {}).get("shadow", pipeconfig.get("shadow", False))
    run:
        pindel_fix_dpaf(
            vcf_input = input[0],
            vcf_output = output[0]
        )

rule pindel_xlsx:
    input:
        "{stype}/pindel/{sname}_pindel.vcf"
    output:
        "{stype}/pindel/{sname}_pindel.xlsx"
    params:
        bed = pipeconfig["rules"]["pindel"]["bed"],
        tumorname = tumorname,
        normalname = normalname
    shadow:
        pipeconfig["rules"].get("pindel", {}).get("shadow", pipeconfig.get("shadow", False))
    run:
        pindel_write_xlsx(
            vcf_input = input[0],
            xlsx_output = output[0],
            bedfile = params.bed,
            tumor = params.tumorname,
            normal = params.normalname if params.normalname else None
        )
