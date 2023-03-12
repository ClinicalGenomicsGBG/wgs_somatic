# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
from workflows.scripts.sex import calc_sex
from workflows.scripts.create_segfile import create_seg
from workflows.scripts.fix_sexploidyfile import mod_sex_vcf

def get_canvas_input(wcs):
    inputs = dict(
        bam = "{stype}/realign/{sname}_REALIGNED.bam",
        bai = "{stype}/realign/{sname}_REALIGNED.bam.bai",
        germline_snv_vcf = f"{'normal' if normalid else 'tumor'}/dnascope/{normalid or tumorid}_germline_SNVsOnly.recode.vcf",
        wgscov = f"{'normal' if normalid else 'tumor'}/reports/{normalid or tumorid}_WGScov.tsv",
        ycov = f"{'normal' if normalid else 'tumor'}/reports/{normalid or tumorid}_Ycov.tsv",
    )

    if wcs.vartype == "somatic":
        inputs |= dict(somatic_vcf = f"tumor/tnscope/{tumorid}_somatic.vcf")
    
    return inputs


def get_canvas_extargs(wcs, input):
    if wcs.vartype == "somatic":
        return f"-t TN --somatic_vcf {input.somatic_vcf} "
    else:
        return "-t germline "


rule canvas:
    input:
        unpack(get_canvas_input)
    params:
        genomeversion = config["reference"],
        dll = pipeconfig["singularities"]["canvas"]["dll"],
        genomedir = pipeconfig["singularities"]["canvas"]["reference"],
        kmerfile = pipeconfig["singularities"]["canvas"]["kmerfile"],
        run_py = pipeconfig["singularities"]["canvas"]["tool_path"],
        filter13 = pipeconfig["singularities"]["canvas"]["filter13"],
        extargs = get_canvas_extargs
    singularity:
        pipeconfig["singularities"]["canvas"]["simg"]
    output:
        temp("{stype}/canvas/{sname}_{vartype}_CNV.vcf.gz"),
        temp("{stype}/canvas/{sname}_{vartype}_CNV_observed.seg"),
        temp("{stype}/canvas/{sname}_{vartype}_CNV_called.seg")
    shadow:
        pipeconfig["rules"].get("canvas", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "echo $HOSTNAME;"
        "{params.run_py} "
            "--genomeversion {params.genomeversion} "
            "--bam {input.bam} "
            "--normal_vcf {input.germline_snv_vcf} "
            "--o {wildcards.stype}/canvas/ "
            "--samplename {wildcards.sname} "
            "--wgscovfile {input.wgscov} "
            "--ycovfile {input.ycov} "
            "--referencedir {params.genomedir} "
            "--kmerfile {params.kmerfile} "
            "--canvasdll {params.dll} "
            "--filterfile {params.filter13} "
            "{params.extargs} "


rule filter_canvas:
    input:
        "{stype}/canvas/{sname}_{vartype}_CNV.vcf.gz"
    params:
        annotate = pipeconfig["rules"]["canvas"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),        
        annotate_ref = pipeconfig["rules"]["canvas"]["annotate_ref"]
    output:
        xlsx = temp("{stype}/canvas/{sname}_CNV_{vartype}.vcf.xlsx"),
        vcf = temp("{stype}/canvas/{sname}_CNV_{vartype}.vcf")
    shadow:
        pipeconfig["rules"].get("canvas", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "gunzip {input}; "
        "grep "
            "-v 'Canvas:REF' "
            "{wildcards.stype}/canvas/{wildcards.sname}_{wildcards.vartype}_CNV.vcf "
            "> {output.vcf}; "
        "{params.annotate} "
            "-v {output.vcf} "
            "-g {params.annotate_ref} "
            "-o {wildcards.stype}/canvas/"

rule convert_to_alissaformat:
    input:
        "{stype}/canvas/{sname}_CNV_germline.vcf"
    params:
        python = pipeconfig["rules"]["convert_to_alissaformat"]["python"],
        converter = pipeconfig["rules"]["convert_to_alissaformat"].get("converter", f"{ROOT_DIR}/workflows/scripts/canvas_to_interpreter/canvasvcf_to_interpreter.py"),
        referencegenome = pipeconfig["referencegenome"]
    output:
        temp("{stype}/canvas/{sname}_CNV_germline_alissaformat.vcf")
    shadow:
        pipeconfig["rules"].get("convert_to_alissaformat", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "{params.python} {params.converter} -l -q {input} {output} {params.referencegenome}"


rule merge_snvs_cnvs:
    input:
        snvs = "{stype}/dnascope/{sname}_germline_refseq3kfilt.vcf.gz",
        snvs_csi = "{stype}/dnascope/{sname}_germline_refseq3kfilt.vcf.gz.csi",
        cnvs = "{stype}/canvas/{sname}_CNV_germline_alissaformat.vcf.gz",
        cnvs_csi = "{stype}/canvas/{sname}_CNV_germline_alissaformat.vcf.gz.csi",
    output:
        temp("{stype}/dnascope/{sname}_{hgX}_SNV_CNV_germline.vcf")
    shadow:
        pipeconfig["rules"].get("merge_snvs_cnvs", {}).get("shadow", pipeconfig.get("shadow", False))
    params:
        bcftools = pipeconfig["rules"]["merge_snvs_cnvs"]["bcftools"]
    shell:
        "{params.bcftools} concat --allow-overlaps {input.snvs} {input.cnvs} -Ov -o {output}"
