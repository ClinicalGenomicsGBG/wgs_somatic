from workflows.scripts.parse_somalier import SomalierParser

rule somalier_extract:
    input:
        "{stype}/realign/{sname}_REALIGNED.bam"
    params:
        sites = pipeconfig["rules"]["somalier"]["sites"],
        reference = pipeconfig["referencegenome"],
        outdir = lambda wildcards: f"{wildcards.stype}/somalier",
        vstamp = f"{VDIR}/somalier_extract.txt"
    singularity:
        pipeconfig["singularities"]["somalier"]["sing"]
    threads:
        clusterconf["somalier_extract"]["threads"]
    output:
        temp("{stype}/somalier/{sname}.somalier")
    shell:
        """
        # Version info
        somalier 2>&1 | head -n 1 > {params.vstamp}

        # Run somalier extract
        somalier extract \
            -d {params.outdir} \
            --sites {params.sites} \
            -f {params.reference} \
            --sample-prefix={wildcards.stype}_ \
            {input}
        mv {params.outdir}/*.somalier {output}
        """

rule somalier_relate:
    input:
        normal_somalier = expand("{stype}/somalier/{sname}.somalier", sname=normalid, stype=stype_normal) if normalid else [],
        tumor_somalier = expand("{stype}/somalier/{sname}.somalier", sname=tumorid, stype=stype_tumor) if tumorid else [],
    params:
        outpath = lambda wildcards: f"{wildcards.stype}/somalier/somalier",
        vstamp = f"{VDIR}/somalier_relate.txt"
    singularity:
        pipeconfig["singularities"]["somalier"]["sing"]
    shadow:
        pipeconfig["rules"].get("somalier_relate", {}).get("shadow", pipeconfig.get("shadow", False))
    output:
        pairs = temp("{stype}/somalier/somalier.pairs.tsv"),
        samples = temp("{stype}/somalier/somalier.samples.tsv")
    shell:
        """
        # Version info
        somalier 2>&1 | head -n 1 > {params.vstamp}

        # Run somalier relate
        somalier relate \
            -i \
            -o {params.outpath} \
            {input.normal_somalier} {input.tumor_somalier}
        """

rule somalier_parse_sex:
    input:
        pairs = "{stype}/somalier/somalier.pairs.tsv",
        samples = "{stype}/somalier/somalier.samples.tsv"
    output:
        sex = temp("{stype}/somalier/calculated_sex.txt")
    run:
        parser = SomalierParser(
            pairs_file=f"{input.pairs}",
            samples_file=f"{input.samples}",
            tumorstring=stype_tumor,
            normalstring=stype_normal,
            match_cutoff=filterconfig["somalier_filter"]["min_relatedness"]
        )
        
        with open(output.sex, "w") as f:
            f.write(parser.sex + "\n")
