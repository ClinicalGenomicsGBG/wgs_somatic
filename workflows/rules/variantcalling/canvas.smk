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
                temp("{stype}/canvas/{sname}_CNV_somatic.vcf.xlsx"),
                temp("{stype}/canvas/{sname}_CNV_somatic.vcf")
            shadow:
                pipeconfig["rules"].get("canvas", {}).get("shadow", pipeconfig.get("shadow", False))
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
        temp("{stype}/canvas/{sname}_CNV_germline.vcf.xlsx"),
        temp("{stype}/canvas/{sname}_CNV_germline.vcf")
    shadow:
        pipeconfig["rules"].get("canvas", {}).get("shadow", pipeconfig.get("shadow", False))
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
                bam = "{stype}/realign/{sname}_REALIGNED.bam",
                bai = "{stype}/realign/{sname}_REALIGNED.bam.bai",
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
                temp("{stype}/canvas/{sname}_somatic_CNV.vcf.gz"),
                temp("{stype}/canvas/{sname}_somatic_CNV_observed.seg"),
                temp("{stype}/canvas/{sname}_somatic_CNV_called.seg")
            shadow:
                pipeconfig["rules"].get("canvas", {}).get("shadow", pipeconfig.get("shadow", False))
            shell:
                "echo $HOSTNAME;"
                "{params.run_py} --genomeversion {params.genomeversion} --bam {input.bam} --normal_vcf {input.germline_snv_vcf} --o {wildcards.stype}/canvas/ -t TN --samplename {wildcards.sname} --wgscovfile {input.normal_wgscov} --ycovfile {input.normal_ycov} --somatic_vcf {input.somatic_vcf} --referencedir {params.genomedir} --kmerfile {params.kmerfile} --canvasdll {params.dll} --filterfile {params.filter13}"
 
if normalid:
    rule canvas_germline:
        input:
            germline_snv_vcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            bam = "{stype}/realign/{sname}_REALIGNED.bam",
            bai = "{stype}/realign/{sname}_REALIGNED.bam.bai",
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
            temp("{stype}/canvas/{sname}_germline_CNV.vcf.gz"),
            temp("{stype}/canvas/{sname}_germline_CNV_observed.seg"),
            temp("{stype}/canvas/{sname}_germline_CNV_called.seg")
        shadow:
            pipeconfig["rules"].get("canvas", {}).get("shadow", pipeconfig.get("shadow", False))
        shell:
            "echo $HOSTNAME;"
            "{params.run_py} --genomeversion {params.genomeversion} --bam {input.bam} --normal_vcf {input.germline_snv_vcf} --o {wildcards.stype}/canvas/ -t germline --samplename {wildcards.sname} --wgscovfile {input.normal_wgscov} --ycovfile {input.normal_ycov} --referencedir {params.genomedir} --kmerfile {params.kmerfile} --canvasdll {params.dll} --filterfile {params.filter13}"
else:
    rule canvas_germline:
        input:
            germline_snv_vcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            bam = "{stype}/realign/{sname}_REALIGNED.bam",
            bai = "{stype}/realign/{sname}_REALIGNED.bam.bai",
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
            temp("{stype}/canvas/{sname}_germline_CNV.vcf.gz"),
            temp("{stype}/canvas/{sname}_germline_CNV_observed.seg"),
            temp("{stype}/canvas/{sname}_germline_CNV_called.seg")
        shadow:
            pipeconfig["rules"].get("canvas", {}).get("shadow", pipeconfig.get("shadow", False))
        shell:
            "echo $HOSTNAME;"
            "{params.run_py} --genomeversion {params.genomeversion} --bam {input.bam} --normal_vcf {input.germline_snv_vcf} --o {wildcards.stype}/canvas/ -t germline --samplename {wildcards.sname} --wgscovfile {input.tumor_wgscov} --ycovfile {input.tumor_ycov} --referencedir {params.genomedir} --kmerfile {params.kmerfile} --canvasdll {params.dll} --filterfile {params.filter13}"

rule convert_to_alissaformat:
    input:
        germline_cnv_vcf = expand("{stype}/canvas/{sname}_CNV_germline.vcf", sname=normalid, stype=sampleconfig[normalname]["stype"])
    params:
        python = pipeconfig["rules"]["convert_to_alissaformat"]["python"],
        converter = pipeconfig["rules"]["convert_to_alissaformat"].get("converter", f"{ROOT_DIR}/workflows/scripts/canvas_to_interpreter/canvasvcf_to_interpreter.py"),
        referencegenome = pipeconfig["referencegenome"]
    output:
        temp("{stype}/canvas/{sname}_CNV_germline_alissaformat.vcf")
    shadow:
        pipeconfig["rules"].get("convert_to_alissaformat", {}).get("shadow", pipeconfig.get("shadow", False))
    run:
        shell(f"{params.python} {params.converter} -l -q {input} {output} {params.referencegenome}")


rule merge_snvs_cnvs:
    input:
        snvs = expand("{stype}/dnascope/{sname}_germline_refseq3kfilt.vcf.gz", stype=sampleconfig[normalname]["stype"], sname=normalid, hgX=reference), 
        snvs_csi = expand("{stype}/dnascope/{sname}_germline_refseq3kfilt.vcf.gz.csi", stype=sampleconfig[normalname]["stype"], sname=normalid, hgX=reference), 
        cnvs = expand("{stype}/canvas/{sname}_CNV_germline_alissaformat.vcf.gz", stype=sampleconfig[normalname]["stype"], sname=normalid, hgX=reference),
        cnvs_csi = expand("{stype}/canvas/{sname}_CNV_germline_alissaformat.vcf.gz.csi", stype=sampleconfig[normalname]["stype"], sname=normalid, hgX=reference),
    output:
        temp("{stype}/dnascope/{sname}_{hgX}_SNV_CNV_germline.vcf")
    shadow:
        pipeconfig["rules"].get("merge_snvs_cnvs", {}).get("shadow", pipeconfig.get("shadow", False))
    params:
        bcftools = pipeconfig["rules"]["merge_snvs_cnvs"]["bcftools"]
    run:
        shell(f"{params.bcftools} concat --allow-overlaps {input.snvs} {input.cnvs} -Ov -o {output}")
