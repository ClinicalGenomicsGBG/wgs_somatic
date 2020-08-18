# vim: syntax=python tabstop=4 expandtab
# # coding: utf-8
import os

rule manta_germline:
    input:
        expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=sampleconfig[normalname]["stype"])
    params:
        reference = pipeconfig["referencegenome"],
        svdb = pipeconfig["rules"]["manta"]["svdb"],
        mantaconf = pipeconfig["rules"]["manta"]["mantaconf"], 
        annotate = pipeconfig["rules"]["manta"]["annotate"],
        annotate_ref = pipeconfig["rules"]["manta"]["annotate_ref"],
        swegendb = pipeconfig["rules"]["manta"]["swegendb"],
        gnomaddb = pipeconfig["rules"]["manta"]["gnomaddb"],
        localdb = pipeconfig["rules"]["manta"]["localdb"],
        bcftools = pipeconfig["rules"]["manta"]["bcftools"] 
    output:
        "{stype}/manta/{sname}_germline_mantaSV.vcf",
        "{stype}/manta/{sname}_germline_mantaSV.vcf.xlsx",
        "{stype}/manta/{sname}_germline_MantaBNDs.vcf",
        "{stype}/manta/{sname}_germline_MantaNOBNDs.vcf"
    run:
        if not os.path.isfile(f"{wildcards.stype}/manta/runWorkflow.py"):
            shell("{params.mantaconf} --normalBam={input} --referenceFasta {params.reference} --runDir {wildcards.stype}/manta/") # Preparing Manta
        if not os.path.isfile(f"{wildcards.stype}/manta/results/variants/diploidSV.vcf.gz"):
            shell("{wildcards.stype}/manta/runWorkflow.py -m local") #Running Manta
        if not os.path.isfile(f"{wildcards.stype}/manta/results/variants/diploidSV.vcf"):
            shell("gunzip {wildcards.stype}/manta/results/variants/diploidSV.vcf.gz")
        shell("grep -e $'\t'PASS$'\t' -e '^#' {wildcards.stype}/manta/results/variants/diploidSV.vcf > {wildcards.stype}/manta/results/variants/diploidSV_PASS.vcf")
        shell("{params.svdb} --query --query_vcf {wildcards.stype}/manta/results/variants/diploidSV_PASS.vcf --db {params.swegendb} --out_frq SWEFRQ --in_frq FRQ --out_occ SWEOCC --in_occ OCC > {wildcards.stype}/manta/results/variants/diploidSV_PASS_swegen.vcf")
        shell("{params.svdb} --query --query_vcf {wildcards.stype}/manta/results/variants/diploidSV_PASS_swegen.vcf --db {params.gnomaddb} --out_frq GNOMAD_AF --in_frq AF --out_occ GNOMAD_AC --in_occ AC > {wildcards.stype}/manta/results/variants/diploidSV_PASS_swegen_gnomad.vcf")
        shell("{params.svdb} --query --query_vcf {wildcards.stype}/manta/results/variants/diploidSV_PASS_swegen_gnomad.vcf --out_frq KGGFRQ --out_occ KGGOCC --sqdb {params.localdb} > {wildcards.stype}/manta/results/variants/diploidSV_PASS_swegen_gnomad_kgg.vcf")
        shell("{params.bcftools} filter -e 'INFO/GNOMAD_AF >= 0.03 | INFO/SWEFRQ >= 0.03 | INFO/KGGFRQ >= 0.05' {wildcards.stype}/manta/results/variants/diploidSV_PASS_swegen_gnomad_kgg.vcf > {wildcards.stype}/manta/results/variants/diploidSV_PASS_swegen_gnomad_kgg_freqfiltered.vcf")
        shell("grep -Ev 'GL000|hs37d5' {wildcards.stype}/manta/results/variants/diploidSV_PASS_swegen_gnomad_kgg_freqfiltered.vcf > {wildcards.stype}/manta/{wildcards.sname}_germline_mantaSV.vcf")
        shell("grep -e '^#' -e 'MantaBND:' {wildcards.stype}/manta/{wildcards.sname}_germline_mantaSV.vcf > {wildcards.stype}/manta/{wildcards.sname}_germline_MantaBNDs.vcf")
        shell("grep -v 'MantaBND:' {wildcards.stype}/manta/{wildcards.sname}_germline_mantaSV.vcf > {wildcards.stype}/manta/{wildcards.sname}_germline_MantaNOBNDs.vcf")
        shell("{params.annotate} -v {wildcards.stype}/manta/{wildcards.sname}_germline_mantaSV.vcf -g {params.annotate_ref} -o {wildcards.stype}/manta")

rule manta_somatic:
    input:
        tumorbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        normalbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=sampleconfig[normalname]["stype"])
    params:
        reference = pipeconfig["referencegenome"],
        svdb = pipeconfig["rules"]["manta"]["svdb"],
        mantaconf = pipeconfig["rules"]["manta"]["mantaconf"],
        annotate = pipeconfig["rules"]["manta"]["annotate"],
        annotate_ref = pipeconfig["rules"]["manta"]["annotate_ref"],
        swegendb = pipeconfig["rules"]["manta"]["swegendb"],
        gnomaddb = pipeconfig["rules"]["manta"]["gnomaddb"],
        localdb = pipeconfig["rules"]["manta"]["localdb"],
        bcftools = pipeconfig["rules"]["manta"]["bcftools"]
    output:
        "{stype}/manta/{sname}_somatic_mantaSV.vcf",
        "{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx",
        "{stype}/manta/{sname}_somatic_MantaBNDs.vcf",
        "{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf"
    run:
        if not os.path.isfile(f"{wildcards.stype}/manta/runWorkflow.py"):
            shell("{params.mantaconf} --tumorBam={input.tumorbam} --normalBam={input.normalbam} --referenceFasta {params.reference} --runDir {wildcards.stype}/manta/") # Preparing Manta
        if not os.path.isfile(f"{wildcards.stype}/manta/results/variants/somaticSV.vcf.gz"):
            shell("{wildcards.stype}/manta/runWorkflow.py -m local") #Running Manta
        if not os.path.isfile(f"{wildcards.stype}/manta/results/variants/somaticSV.vcf"):
            shell("gunzip {wildcards.stype}/manta/results/variants/somaticSV.vcf.gz")
        shell("grep -e $'\t'PASS$'\t' -e '^#' {wildcards.stype}/manta/results/variants/somaticSV.vcf > {wildcards.stype}/manta/results/variants/somaticSV_PASS.vcf")
        shell("{params.svdb} --query --query_vcf {wildcards.stype}/manta/results/variants/somaticSV_PASS.vcf --db {params.swegendb} --out_frq SWEFRQ --in_frq FRQ --out_occ SWEOCC --in_occ OCC > {wildcards.stype}/manta/results/variants/somaticSV_PASS_swegen.vcf")
        shell("{params.svdb} --query --query_vcf {wildcards.stype}/manta/results/variants/somaticSV_PASS_swegen.vcf --db {params.gnomaddb} --out_frq GNOMAD_AF --in_frq AF --out_occ GNOMAD_AC --in_occ AC > {wildcards.stype}/manta/results/variants/somaticSV_PASS_swegen_gnomad.vcf")
        shell("{params.svdb} --query --query_vcf {wildcards.stype}/manta/results/variants/somaticSV_PASS_swegen_gnomad.vcf --out_frq KGGFRQ --out_occ KGGOCC --sqdb {params.localdb} > {wildcards.stype}/manta/results/variants/somaticSV_PASS_swegen_gnomad_kgg.vcf")
        shell("{params.bcftools} filter -e 'INFO/GNOMAD_AF >= 0.03 | INFO/SWEFRQ >= 0.03 | INFO/KGGFRQ >= 0.05' {wildcards.stype}/manta/results/variants/somaticSV_PASS_swegen_gnomad_kgg.vcf > {wildcards.stype}/manta/results/variants/somaticSV_PASS_swegen_gnomad_kgg_freqfiltered.vcf")
        shell("grep -Ev 'GL000|hs37d5' {wildcards.stype}/manta/results/variants/somaticSV_PASS_swegen_gnomad_kgg_freqfiltered.vcf > {wildcards.stype}/manta/{wildcards.sname}_somatic_mantaSV.vcf")
        shell("grep -e '^#' -e 'MantaBND:' {wildcards.stype}/manta/{wildcards.sname}_somatic_mantaSV.vcf > {wildcards.stype}/manta/{wildcards.sname}_somatic_MantaBNDs.vcf")
        shell("grep -v 'MantaBND:' {wildcards.stype}/manta/{wildcards.sname}_somatic_mantaSV.vcf > {wildcards.stype}/manta/{wildcards.sname}_somatic_MantaNOBNDs.vcf")
        shell("{params.annotate} -v {wildcards.stype}/manta/{wildcards.sname}_somatic_mantaSV.vcf -g {params.annotate_ref} -o {wildcards.stype}/manta")




