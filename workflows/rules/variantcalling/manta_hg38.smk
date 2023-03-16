# vim: syntax=python tabstop=4 expandtab
# # coding: utf-8
import os
from workflows.scripts.annotate_manta.manta_summary import manta_summary

rule manta_germline:
    input:
        normalbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        normalbai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=normalid, stype=sampleconfig[normalname]["stype"])
    params:
        reference = pipeconfig["referencegenome"],
        svdb = pipeconfig["rules"]["manta"]["svdb"],
        mantaconf = pipeconfig["rules"]["manta"]["mantaconf"], 
        annotate = pipeconfig["rules"]["manta"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),
        annotate_ref = pipeconfig["rules"]["manta"]["annotate_ref"],
        bcftools = pipeconfig["rules"]["manta"]["bcftools"] 
    output:
        sv_vcf = temp("{stype}/manta/{sname}_germline_mantaSV.vcf"),
        sv_xlsx = temp("{stype}/manta/{sname}_germline_mantaSV.vcf.xlsx"),
        bnd_vcf = temp("{stype}/manta/{sname}_germline_MantaBNDs.vcf"),
        nobnd_vcf = temp("{stype}/manta/{sname}_germline_MantaNOBNDs.vcf")
    shadow:
        pipeconfig["rules"].get("manta", {}).get("shadow", pipeconfig.get("shadow", False))
    run:
        if not os.path.isfile(f"{wildcards.stype}/manta/runWorkflow.py"):
            shell("{params.mantaconf} --normalBam={input.normalbam} --referenceFasta {params.reference} --runDir {wildcards.stype}/manta/") # Preparing Manta
        if not os.path.isfile(f"{wildcards.stype}/manta/results/variants/diploidSV.vcf.gz"):
            shell("{wildcards.stype}/manta/runWorkflow.py -m local") #Running Manta
        if not os.path.isfile(f"{wildcards.stype}/manta/results/variants/diploidSV.vcf"):
            shell("gunzip {wildcards.stype}/manta/results/variants/diploidSV.vcf.gz")
        shell("grep -e $'\t'PASS$'\t' -e '^#' {wildcards.stype}/manta/results/variants/diploidSV.vcf > {wildcards.stype}/manta/results/variants/diploidSV_PASS.vcf")
        shell("grep -Ev 'GL000|hs37d5' {wildcards.stype}/manta/results/variants/diploidSV_PASS.vcf > {output.sv_vcf}")
        shell("grep -e '^#' -e 'MantaBND:' {output.sv_vcf} > {output.bnd_vcf}")
        shell("grep -v 'MantaBND:' {output.sv_vcf} > {output.nobnd_vcf}")
        shell("{params.annotate} -v {output.sv_vcf} -g {params.annotate_ref} -o {wildcards.stype}/manta")

if normalid:
    rule manta_somatic:
        input:
            tumorbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            tumorbai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            normalbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            normalbai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=normalid, stype=sampleconfig[normalname]["stype"])
        params:
            reference = pipeconfig["referencegenome"],
            svdb = pipeconfig["rules"]["manta"]["svdb"],
            mantaconf = pipeconfig["rules"]["manta"]["mantaconf"],
            annotate = pipeconfig["rules"]["manta"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),
            annotate_ref = pipeconfig["rules"]["manta"]["annotate_ref"],
            bcftools = pipeconfig["rules"]["manta"]["bcftools"]
        output:
            sv_vcf = temp("{stype}/manta/{sname}_somatic_mantaSV.vcf"),
            sv_xlsx = temp("{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx"),
            bnd_vcf = temp("{stype}/manta/{sname}_somatic_MantaBNDs.vcf"),
            nobnd_vcf = temp("{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf")
        shadow:
            pipeconfig["rules"].get("manta", {}).get("shadow", pipeconfig.get("shadow", False))
        run:
            if not os.path.isfile(f"{wildcards.stype}/manta/runWorkflow.py"):
                shell("{params.mantaconf} --tumorBam={input.tumorbam} --normalBam={input.normalbam} --referenceFasta {params.reference} --runDir {wildcards.stype}/manta/") # Preparing Manta
            if not os.path.isfile(f"{wildcards.stype}/manta/results/variants/somaticSV.vcf.gz"):
                shell("{wildcards.stype}/manta/runWorkflow.py -m local") #Running Manta
            if not os.path.isfile(f"{wildcards.stype}/manta/results/variants/somaticSV.vcf"):
                shell("gunzip {wildcards.stype}/manta/results/variants/somaticSV.vcf.gz")
            shell("grep -e $'\t'PASS$'\t' -e '^#' {wildcards.stype}/manta/results/variants/somaticSV.vcf > {wildcards.stype}/manta/results/variants/somaticSV_PASS.vcf")
            shell("grep -Ev 'GL000|hs37d5' {wildcards.stype}/manta/results/variants/somaticSV_PASS.vcf > {output.sv_vcf}")
            shell("grep -e '^#' -e 'MantaBND:' {output.sv_vcf} > {output.bnd_vcf}")
            shell("grep -v 'MantaBND:' {output.sv_vcf} > {output.nobnd_vcf}")
            shell("{params.annotate} -v {output.sv_vcf} -g {params.annotate_ref} -o {wildcards.stype}/manta")
else:
    rule manta_somatic:
        input:
            tumorbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            tumorbai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        params:
            reference = pipeconfig["referencegenome"],
            svdb = pipeconfig["rules"]["manta"]["svdb"],
            mantaconf = pipeconfig["rules"]["manta"]["mantaconf"],
            annotate = pipeconfig["rules"]["manta"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),
            annotate_ref = pipeconfig["rules"]["manta"]["annotate_ref"],
            bcftools = pipeconfig["rules"]["manta"]["bcftools"]
        output:
            sv_vcf = temp("{stype}/manta/{sname}_somatic_mantaSV.vcf"),
            sv_xlsx = temp("{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx"),
            bnd_vcf = temp("{stype}/manta/{sname}_somatic_MantaBNDs.vcf"),
            nobnd_vcf = temp("{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf")
        shadow:
            pipeconfig["rules"].get("manta", {}).get("shadow", pipeconfig.get("shadow", False))
        run:
            if not os.path.isfile(f"{wildcards.stype}/manta/runWorkflow.py"):
                shell("{params.mantaconf} --tumorBam={input.tumorbam} --referenceFasta {params.reference} --runDir {wildcards.stype}/manta/") # Preparing Manta
            if not os.path.isfile(f"{wildcards.stype}/manta/results/variants/tumorSV.vcf.gz"):
                shell("{wildcards.stype}/manta/runWorkflow.py -m local") #Running Manta
            if not os.path.isfile(f"{wildcards.stype}/manta/results/variants/tumorSV.vcf"):
                shell("gunzip {wildcards.stype}/manta/results/variants/tumorSV.vcf.gz")
            shell("grep -e $'\t'PASS$'\t' -e '^#' {wildcards.stype}/manta/results/variants/tumorSV.vcf > {wildcards.stype}/manta/results/variants/tumorSV_PASS.vcf")
            shell("grep -Ev 'GL000|hs37d5' {wildcards.stype}/manta/results/variants/tumorSV_PASS.vcf > {output.sv_vcf}")
            shell("grep -e '^#' -e 'MantaBND:' {output.sv_vcf} > {output.bnd_vcf}")
            shell("grep -v 'MantaBND:' {output.sv_vcf} > {output.nobnd_vcf}")
            shell("{params.annotate} -v {output.sv_vcf} -g {params.annotate_ref} -o {wildcards.stype}/manta")

if normalid:
    rule manta_summary:
        input:
            "{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx"
        params:
            genelist = pipeconfig["rules"]["manta_summary"]["genelist"]
        output:
            temp("{stype}/manta/{sname}_somatic_mantaSV_Summary.xlsx")
        shadow:
            pipeconfig["rules"].get("manta_summary", {}).get("shadow", pipeconfig.get("shadow", False))
        run:
            manta_summary(input, output, tumorname, f"{params.genelist}", normalname)

else:
    rule manta_summary:
        input:
            "{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx"
        params:
            genelist = pipeconfig["rules"]["manta_summary"]["genelist_tumoronly"]
        output:
            temp("{stype}/manta/{sname}_somatic_mantaSV_Summary.xlsx")
        shadow:
            pipeconfig["rules"].get("manta_summary", {}).get("shadow", pipeconfig.get("shadow", False))
        run:
            manta_summary(input, output, tumorname, f"{params.genelist}", normalname='')
