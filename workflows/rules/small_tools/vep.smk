# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule vep_annotate:
  input:
    "{path}/{file}.vcf"
  params:
    reference = pipeconfig["referencegenome"],
    cache = pipeconfig["rules"]["vep_annotate"]["cache"]
  output:
    temp("{path}/{file}_vep.vcf"),
    temp("{path}/{file}_vep.vcf_warnings.txt")
  singularity:
    pipeconfig["singularities"]["vep"]["sing"]
  threads:
    clusterconf["vep_annotate"].get("threads", 4)
  shell:
    """
    vep \
      --input_file {input} \
      --output_file {output} \
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
  output:
    temp("{path}/{file}_vep_gnomadfilt.vcf")
  singularity:
    pipeconfig["singularities"]["vep"]["sing"]
  threads:
    clusterconf["vep_gnomad_filter"].get("threads", 1)
  shell:
    """
    # Parse the gnomAD annotations and filter variants with AF > 1% (genome or exome)
    bcftools \
      +split-vep \
      -c gnomADe_AF:Float,gnomADg_AF:Float \
      -e 'gnomADe_AF>0.01 || gnomADg_AF>0.01' \
      -o {output} {input}
    """

