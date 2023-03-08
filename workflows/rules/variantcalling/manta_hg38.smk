# vim: syntax=python tabstop=4 expandtab
# # coding: utf-8
import os
from workflows.scripts.annotate_manta.manta_summary import manta_summary

rule manta_germline:
    input:
        expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=sampleconfig[normalname]["stype"])
    params:
        reference = pipeconfig["referencegenome"],
        svdb = pipeconfig["rules"]["manta"]["svdb"],
        mantaconf = pipeconfig["rules"]["manta"]["mantaconf"], 
        annotate = pipeconfig["rules"]["manta"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),
        annotate_ref = pipeconfig["rules"]["manta"]["annotate_ref"],
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
        shell("grep -Ev 'GL000|hs37d5' {wildcards.stype}/manta/results/variants/diploidSV_PASS.vcf > {wildcards.stype}/manta/{wildcards.sname}_germline_mantaSV.vcf")
        shell("grep -e '^#' -e 'MantaBND:' {wildcards.stype}/manta/{wildcards.sname}_germline_mantaSV.vcf > {wildcards.stype}/manta/{wildcards.sname}_germline_MantaBNDs.vcf")
        shell("grep -v 'MantaBND:' {wildcards.stype}/manta/{wildcards.sname}_germline_mantaSV.vcf > {wildcards.stype}/manta/{wildcards.sname}_germline_MantaNOBNDs.vcf")
        shell("{params.annotate} -v {wildcards.stype}/manta/{wildcards.sname}_germline_mantaSV.vcf -g {params.annotate_ref} -o {wildcards.stype}/manta")

if normalid:
    rule manta_somatic:
        input:
            tumorbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            normalbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=sampleconfig[normalname]["stype"])
        params:
            reference = pipeconfig["referencegenome"],
            svdb = pipeconfig["rules"]["manta"]["svdb"],
            mantaconf = pipeconfig["rules"]["manta"]["mantaconf"],
            annotate = pipeconfig["rules"]["manta"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),
            annotate_ref = pipeconfig["rules"]["manta"]["annotate_ref"],
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
            shell("grep -Ev 'GL000|hs37d5' {wildcards.stype}/manta/results/variants/somaticSV_PASS.vcf > {wildcards.stype}/manta/{wildcards.sname}_somatic_mantaSV.vcf")
            shell("grep -e '^#' -e 'MantaBND:' {wildcards.stype}/manta/{wildcards.sname}_somatic_mantaSV.vcf > {wildcards.stype}/manta/{wildcards.sname}_somatic_MantaBNDs.vcf")
            shell("grep -v 'MantaBND:' {wildcards.stype}/manta/{wildcards.sname}_somatic_mantaSV.vcf > {wildcards.stype}/manta/{wildcards.sname}_somatic_MantaNOBNDs.vcf")
            shell("{params.annotate} -v {wildcards.stype}/manta/{wildcards.sname}_somatic_mantaSV.vcf -g {params.annotate_ref} -o {wildcards.stype}/manta")
else:
    rule manta_somatic:
        input:
            tumorbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        params:
            reference = pipeconfig["referencegenome"],
            svdb = pipeconfig["rules"]["manta"]["svdb"],
            mantaconf = pipeconfig["rules"]["manta"]["mantaconf"],
            annotate = pipeconfig["rules"]["manta"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),
            annotate_ref = pipeconfig["rules"]["manta"]["annotate_ref"],
            bcftools = pipeconfig["rules"]["manta"]["bcftools"]
        output:
            "{stype}/manta/{sname}_somatic_mantaSV.vcf",
            "{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx",
            "{stype}/manta/{sname}_somatic_MantaBNDs.vcf",
            "{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf"
        run:
            if not os.path.isfile(f"{wildcards.stype}/manta/runWorkflow.py"):
                shell("{params.mantaconf} --tumorBam={input.tumorbam} --referenceFasta {params.reference} --runDir {wildcards.stype}/manta/") # Preparing Manta
            if not os.path.isfile(f"{wildcards.stype}/manta/results/variants/tumorSV.vcf.gz"):
                shell("{wildcards.stype}/manta/runWorkflow.py -m local") #Running Manta
            if not os.path.isfile(f"{wildcards.stype}/manta/results/variants/tumorSV.vcf"):
                shell("gunzip {wildcards.stype}/manta/results/variants/tumorSV.vcf.gz")
            shell("grep -e $'\t'PASS$'\t' -e '^#' {wildcards.stype}/manta/results/variants/tumorSV.vcf > {wildcards.stype}/manta/results/variants/tumorSV_PASS.vcf")
            shell("grep -Ev 'GL000|hs37d5' {wildcards.stype}/manta/results/variants/tumorSV_PASS.vcf > {wildcards.stype}/manta/{wildcards.sname}_somatic_mantaSV.vcf")
            shell("grep -e '^#' -e 'MantaBND:' {wildcards.stype}/manta/{wildcards.sname}_somatic_mantaSV.vcf > {wildcards.stype}/manta/{wildcards.sname}_somatic_MantaBNDs.vcf")
            shell("grep -v 'MantaBND:' {wildcards.stype}/manta/{wildcards.sname}_somatic_mantaSV.vcf > {wildcards.stype}/manta/{wildcards.sname}_somatic_MantaNOBNDs.vcf")
            shell("{params.annotate} -v {wildcards.stype}/manta/{wildcards.sname}_somatic_mantaSV.vcf -g {params.annotate_ref} -o {wildcards.stype}/manta")

if normalid:
    rule manta_summary:
        input:
            "{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx"
        params:
            genelist = pipeconfig["rules"]["manta_summary"]["genelist"]
        output:
            "{stype}/manta/{sname}_somatic_mantaSV_Summary.xlsx"
        run:
            manta_summary(input, output, tumorname, f"{params.genelist}", normalname)

else:
    rule manta_summary:
        input:
            "{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx"
        params:
            genelist = pipeconfig["rules"]["manta_summary"]["genelist_tumoronly"]
        output:
            "{stype}/manta/{sname}_somatic_mantaSV_Summary.xlsx"
        run:
            manta_summary(input, output, tumorname, f"{params.genelist}", normalname='')
