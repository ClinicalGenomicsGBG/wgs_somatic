# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
from workflows.scripts.sex import calc_sex
from workflows.scripts.create_segfile import create_seg
from workflows.scripts.fix_sexploidyfile import mod_sex_vcf


if tumorid:
    if normalid:
        rule filter_canvas_somatic:
            input:
                expand("{stype}/canvas/{sname}_somatic_CNV.vcf.gz", sname=tumorid, stype=sampleconfig[tumorname]["stype"])
            params:
                annotate = pipeconfig["rules"]["canvas"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),
                annotate_ref = pipeconfig["rules"]["canvas"]["annotate_ref"]
            output:
                "{stype}/canvas/{sname}_CNV_somatic.vcf.xlsx",
                "{stype}/canvas/{sname}_CNV_somatic.vcf"
            run:
                shell("gunzip {input}")
                shell("grep -v 'Canvas:REF' {wildcards.stype}/canvas/{wildcards.sname}_somatic_CNV.vcf > {wildcards.stype}/canvas/{wildcards.sname}_CNV_somatic_noref.vcf")
                shell("{params.annotate} -v {wildcards.stype}/canvas/{wildcards.sname}_CNV_somatic_noref.vcf -g {params.annotate_ref} -o {wildcards.stype}/canvas/")
                os.rename(f"{wildcards.stype}/canvas/{wildcards.sname}_CNV_somatic_noref.vcf.xlsx", f"{wildcards.stype}/canvas/{wildcards.sname}_CNV_somatic.vcf.xlsx")
                os.rename(f"{wildcards.stype}/canvas/{wildcards.sname}_CNV_somatic_noref.vcf", f"{wildcards.stype}/canvas/{wildcards.sname}_CNV_somatic.vcf")

rule filter_canvas_germline:
    input:
        "{stype}/canvas/{sname}_germline_CNV.vcf.gz"
    params:
        annotate = pipeconfig["rules"]["canvas"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),
        annotate_ref = pipeconfig["rules"]["canvas"]["annotate_ref"]
    output:
        "{stype}/canvas/{sname}_CNV_germline.vcf.xlsx",
        "{stype}/canvas/{sname}_CNV_germline.vcf"
    run:
        shell("gunzip {input}")
        shell("grep -v 'Canvas:REF' {wildcards.stype}/canvas/{wildcards.sname}_germline_CNV.vcf > {wildcards.stype}/canvas/{wildcards.sname}_CNV_germline_noref.vcf")
        shell("{params.annotate} -v {wildcards.stype}/canvas/{wildcards.sname}_CNV_germline_noref.vcf -g {params.annotate_ref} -o {wildcards.stype}/canvas/")
        os.rename(f"{wildcards.stype}/canvas/{wildcards.sname}_CNV_germline_noref.vcf.xlsx", f"{wildcards.stype}/canvas/{wildcards.sname}_CNV_germline.vcf.xlsx")
        os.rename(f"{wildcards.stype}/canvas/{wildcards.sname}_CNV_germline_noref.vcf", f"{wildcards.stype}/canvas/{wildcards.sname}_CNV_germline.vcf")


if tumorid:
    if normalid:
        rule canvas_somatic:
            input:
                germline_snv_vcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                somatic_vcf = expand("{stype}/tnscope/{sname}_somatic.vcf", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                bamfile = "{stype}/realign/{sname}_REALIGNED.bam",
                normal_wgscov = expand("{stype}/reports/{sname}_WGScov.tsv", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                normal_ycov = expand("{stype}/reports/{sname}_Ycov.tsv", sname=normalid, stype=sampleconfig[normalname]["stype"])
            params:
                genomeversion = config["reference"],
                dll = pipeconfig["singularities"]["canvas"]["dll"],
                annotate = pipeconfig["rules"]["canvas"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),
                annotate_ref = pipeconfig["rules"]["canvas"]["annotate_ref"],
                genomedir = pipeconfig["singularities"]["canvas"]["reference"],
                kmerfile = pipeconfig["singularities"]["canvas"]["kmerfile"],
                run_py = pipeconfig["singularities"]["canvas"]["tool_path"],
                filter13 = pipeconfig["singularities"]["canvas"]["filter13"],
                samplename = sampleconfig["tumorname"]
            singularity:
                pipeconfig["singularities"]["canvas"]["sing"]
            output:
                "{stype}/canvas/{sname}_somatic_CNV.vcf.gz",
                "{stype}/canvas/{sname}_somatic_CNV_observed.seg",
                "{stype}/canvas/{sname}_somatic_CNV_called.seg"
            shell:
                "echo $HOSTNAME;"
                "{params.run_py} --genomeversion {params.genomeversion} --bam {input.bamfile} --normal_vcf {input.germline_snv_vcf} --o {wildcards.stype}/canvas/ -t TN --samplename {wildcards.sname} --wgscovfile {input.normal_wgscov} --ycovfile {input.normal_ycov} --somatic_vcf {input.somatic_vcf} --referencedir {params.genomedir} --kmerfile {params.kmerfile} --canvasdll {params.dll} --filterfile {params.filter13}"
 
if normalid:
    rule canvas_germline:
        input:
            germline_snv_vcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            bamfile = "{stype}/realign/{sname}_REALIGNED.bam",
            normal_wgscov = expand("{stype}/reports/{sname}_WGScov.tsv", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            normal_ycov = expand("{stype}/reports/{sname}_Ycov.tsv", sname=normalid, stype=sampleconfig[normalname]["stype"])
        params:
            genomeversion = config["reference"],
            dll = pipeconfig["singularities"]["canvas"]["dll"],
            annotate = pipeconfig["rules"]["canvas"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),
            annotate_ref = pipeconfig["rules"]["canvas"]["annotate_ref"],
            genomedir = pipeconfig["singularities"]["canvas"]["reference"],
            kmerfile = pipeconfig["singularities"]["canvas"]["kmerfile"],
            run_py = pipeconfig["singularities"]["canvas"]["tool_path"],
            filter13 = pipeconfig["singularities"]["canvas"]["filter13"],
            samplename = sampleconfig["normalname"]
        singularity:
            pipeconfig["singularities"]["canvas"]["sing"]
        output:
            "{stype}/canvas/{sname}_germline_CNV.vcf.gz",
            "{stype}/canvas/{sname}_germline_CNV_observed.seg",
            "{stype}/canvas/{sname}_germline_CNV_called.seg"
        shell:
            "echo $HOSTNAME;"
            "{params.run_py} --genomeversion {params.genomeversion} --bam {input.bamfile} --normal_vcf {input.germline_snv_vcf} --o {wildcards.stype}/canvas/ -t germline --samplename {wildcards.sname} --wgscovfile {input.normal_wgscov} --ycovfile {input.normal_ycov} --referencedir {params.genomedir} --kmerfile {params.kmerfile} --canvasdll {params.dll} --filterfile {params.filter13}"
else:
    rule canvas_germline:
        input:
            germline_snv_vcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            bamfile = "{stype}/realign/{sname}_REALIGNED.bam",
            tumor_wgscov = expand("{stype}/reports/{sname}_WGScov.tsv", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            tumor_ycov = expand("{stype}/reports/{sname}_Ycov.tsv", sname=tumorid, stype=sampleconfig[tumorname]["stype"])
        params:
            genomeversion = config["reference"],
            dll = pipeconfig["singularities"]["canvas"]["dll"],
            annotate = pipeconfig["rules"]["canvas"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),
            annotate_ref = pipeconfig["rules"]["canvas"]["annotate_ref"],
            genomedir = pipeconfig["singularities"]["canvas"]["reference"],
            kmerfile = pipeconfig["singularities"]["canvas"]["kmerfile"],
            run_py = pipeconfig["singularities"]["canvas"]["tool_path"],
            filter13 = pipeconfig["singularities"]["canvas"]["filter13"],
            samplename = sampleconfig["tumorname"]
        singularity:
            pipeconfig["singularities"]["canvas"]["sing"]
        output:
            "{stype}/canvas/{sname}_germline_CNV.vcf.gz",
            "{stype}/canvas/{sname}_germline_CNV_observed.seg",
            "{stype}/canvas/{sname}_germline_CNV_called.seg"
        shell:
            "echo $HOSTNAME;"
            "{params.run_py} --genomeversion {params.genomeversion} --bam {input.bamfile} --normal_vcf {input.germline_snv_vcf} --o {wildcards.stype}/canvas/ -t germline --samplename {wildcards.sname} --wgscovfile {input.tumor_wgscov} --ycovfile {input.tumor_ycov} --referencedir {params.genomedir} --kmerfile {params.kmerfile} --canvasdll {params.dll} --filterfile {params.filter13}"

rule convert_to_alissaformat:
    input:
        germline_cnv_vcf = expand("{stype}/canvas/{sname}_CNV_germline.vcf", sname=normalid, stype=sampleconfig[normalname]["stype"])
    params:
        python = pipeconfig["rules"]["convert_to_alissaformat"]["python"],
        converter = pipeconfig["rules"]["convert_to_alissaformat"].get("converter", f"{ROOT_DIR}/workflows/scripts/canvas_to_interpreter/canvasvcf_to_interpreter.py"),
        referencegenome = pipeconfig["referencegenome"]
    output:
        "{stype}/canvas/{sname}_CNV_germline_alissaformat.vcf"
    run:
        shell(f"{params.python} {params.converter} -l -q {input} {output} {params.referencegenome}")


rule merge_snvs_cnvs:
    input:
        snvs = expand("{stype}/dnascope/{sname}_germline_refseq3kfilt.vcf", stype=sampleconfig[normalname]["stype"], sname=normalid, hgX=reference), 
        cnvs = expand("{stype}/canvas/{sname}_CNV_germline_alissaformat.vcf", stype=sampleconfig[normalname]["stype"], sname=normalid, hgX=reference)
    params:
        bgzip = pipeconfig["rules"]["merge_snvs_cnvs"]["bgzip"],
        bcftools = pipeconfig["rules"]["merge_snvs_cnvs"]["bcftools"]
    output:
        "{stype}/dnascope/{sname}_{hgX}_SNV_CNV_germline.vcf.gz"
    run:
        # bgzip CNVcalls
        if not os.path.isfile(f"{input.cnvs}.gz"):
            shell(f"{params.bgzip} --stdout {input.cnvs} > {input.cnvs}.gz")
        # bgzip SNVcalls
        if not os.path.isfile(f"{input.snvs}.gz"):
            shell(f"{params.bgzip} --stdout {input.snvs} > {input.snvs}.gz")
        # index 
        shell(f"{params.bcftools} index -f {input.snvs}.gz")
        shell(f"{params.bcftools} index -f {input.cnvs}.gz")
        # merge calls
        shell(f"{params.bcftools} concat --allow-overlaps {input.snvs}.gz {input.cnvs}.gz -Oz -o {output}")
