rule somalier_extract:
  input:
    "{stype}/realign/{sname}_REALIGNED.bam"
  params:
    sites = pipeconfig["rules"]["somalier"]["sites"],
    reference = pipeconfig["referencegenome"],
    outdir = lambda wildcards: f"{wildcards.stype}/somalier"
  singularity:
    pipeconfig["singularities"]["somalier"]["sing"]
  threads:
    clusterconf["somalier_extract"]["threads"]
  output:
    temp("{stype}/somalier/{sname}.somalier")
  shell:
    """
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
    normal_somalier = expand("{stype}/somalier/{sname}.somalier", sname=normalid, stype=normaltype) if normalid else [],
    tumor_somalier = expand("{stype}/somalier/{sname}.somalier", sname=tumorid, stype=tumortype) if tumorid else [],
  params:
    outpath = lambda wildcards: f"{wildcards.stype}/somalier/somalier"
  singularity:
    pipeconfig["singularities"]["somalier"]["sing"]
  shadow:
    pipeconfig["rules"].get("somalier_relate", {}).get("shadow", pipeconfig.get("shadow", False))
  output:
    pairs = temp("{stype}/somalier/somalier.pairs.tsv"),
    samples = temp("{stype}/somalier/somalier.samples.tsv")
  shell:
    """
    somalier relate \
      -i \
      -o {params.outpath} \
      {input.normal_somalier} {input.tumor_somalier}
    """

