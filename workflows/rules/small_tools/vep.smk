# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule vep_annotate:
  input:
    "{path}/{file}.vcf"
  params:
    reference = pipeconfig["referencegenome"],
    cache = pipeconfig["rules"]["vep_annotate"]["cache"],
    vstamp = f"{VDIR}/vep_annotate.txt"
  output:
    vcf = temp("{path}/{file}_vep.vcf"),
    warnings = temp("{path}/{file}_vep.vcf_warnings.txt")
  singularity:
    pipeconfig["singularities"]["vep"]["sing"]
  threads:
    clusterconf["vep_annotate"].get("threads", 4)
  shell:
    """
    vep --help | grep "ensembl-vep" | tr -s '[:space:]' > {params.vstamp}
    echo "cache version: $(ls {params.cache}/homo_sapiens)" >> {params.vstamp}
    vep --cache --offline --dir_cache {params.cache} --show_cache_info >> {params.vstamp}
    vep \
      --input_file {input} \
      --output_file {output.vcf} \
      --vcf \
      --cache \
      --dir_cache {params.cache} \
      --fasta {params.reference} \
      --offline \
      --pick \
      --no_stats \
      --fork {threads} \
      --af_gnomade \
      --af_gnomadg
    """

rule vep_gnomad_filter:
  input:
    "{path}/{file}_vep.vcf"
  params:
    max_af = filterconfig["gnomad_filter"]["max_allele_frequency"],
    vstamp = f"{VDIR}/vep_gnomad_filter.txt"
  output:
    temp("{path}/{file}_vep_gnomadfilt.vcf")
  singularity:
    pipeconfig["singularities"]["vep"]["sing"]
  threads:
    clusterconf["vep_gnomad_filter"].get("threads", 1)
  shell:
    """
    bcftools --version | head -n 2 > {params.vstamp}
    # Parse the gnomAD annotations and filter variants with AF > 1% (genome or exome)
    bcftools \
      +split-vep \
      -c gnomADe_AF:Float,gnomADg_AF:Float \
      -e 'gnomADe_AF>{params.max_af} || gnomADg_AF>{params.max_af}' \
      -o {output} {input}
    """

