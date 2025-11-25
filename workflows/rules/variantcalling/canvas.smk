# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
from workflows.scripts.create_segfile import create_seg
from workflows.scripts.fix_sexploidyfile import mod_sex_vcf


rule filter_canvas:
    input:
        vcf = "{stype}/canvas/{sname}_canvas_{mode}.vcf.gz"
    params:
        annotate = pipeconfig["rules"]["canvas"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),
        annotate_ref = pipeconfig["rules"]["canvas"]["annotate_ref"]
    output:
        xlsx = "{stype}/canvas/{sname}_canvas_{mode}_filt.vcf.xlsx",
        vcf_out = "{stype}/canvas/{sname}_canvas_{mode}_filt.vcf"
    shadow:
        pipeconfig["rules"].get("canvas", {}).get("shadow", pipeconfig.get("shadow", False))
    run:
        vcf_unzipped = input.vcf[:-3]
        shell("gunzip {input.vcf}")
        shell("grep -v 'Canvas:REF' {vcf_unzipped} > {output.vcf_out}")
        shell("{params.annotate} -v {output.vcf_out} -g {params.annotate_ref} -o {wildcards.stype}/canvas/")

if tumorid:
    if normalid:
        rule canvas_somatic:
            input:
                germline_snv_vcf = expand("{stype}/dnascope/{sname}_DNAscope_germline_SNVsOnly.recode.vcf", sname=normalid, stype=stype_normal),
                somatic_vcf = expand("{stype}/tnscope/{sname}_TNscope_somatic.vcf", sname=tumorid, stype=stype_tumor),
                bam = "{stype}/realign/{sname}_REALIGNED.bam",
                bai = "{stype}/realign/{sname}_REALIGNED.bam.bai",
                somalier_sex = "{stype}/somalier/calculated_sex.txt",
            params:
                genomeversion = config["reference"],
                dll = pipeconfig["singularities"]["canvas"]["dll"],
                annotate = pipeconfig["rules"]["canvas"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),
                annotate_ref = pipeconfig["rules"]["canvas"]["annotate_ref"],
                genomedir = pipeconfig["singularities"]["canvas"]["reference"],
                kmerfile = pipeconfig["singularities"]["canvas"]["kmerfile"],
                run_py = pipeconfig["singularities"]["canvas"]["tool_path"],
                filter13 = pipeconfig["singularities"]["canvas"]["filter13"],
                samplename = sampleconfig["tumorname"],
                intermediate_vcf = "{stype}/canvas/{sname}_somatic_CNV.vcf.gz",
                intermediate_observed = "{stype}/canvas/{sname}_somatic_CNV_observed.seg",
                intermediate_called = "{stype}/canvas/{sname}_somatic_CNV_called.seg"
            singularity:
                pipeconfig["singularities"]["canvas"]["sing"]
            output:
                out_vcf = "{stype}/canvas/{sname}_canvas_somatic.vcf.gz",
                out_observed = "{stype}/canvas/{sname}_canvas_somatic_observed.seg",
                out_called = "{stype}/canvas/{sname}_canvas_somatic_called.seg"
            shadow:
                pipeconfig["rules"].get("canvas", {}).get("shadow", pipeconfig.get("shadow", False))
            shell:
                """
                echo $HOSTNAME;
                SEX=$(cat {input.somalier_sex})
                {params.run_py} --genomeversion {params.genomeversion} --bam {input.bam} --normal_vcf {input.germline_snv_vcf} --o {wildcards.stype}/canvas/ -t TN --samplename {wildcards.sname} --somatic_vcf {input.somatic_vcf} --sex "$SEX" --referencedir {params.genomedir} --kmerfile {params.kmerfile} --canvasdll {params.dll} --filterfile {params.filter13}
                mv {params.intermediate_vcf} {output.out_vcf}
                mv {params.intermediate_observed} {output.out_observed}
                mv {params.intermediate_called} {output.out_called}
                """
    else:
        rule canvas_tumoronly:
            # Note: for tumor-only we use "-t germline" to call variants in the tumor
            input:
                germline_snv_vcf = expand("{stype}/dnascope/{sname}_DNAscope_germline_SNVsOnly.recode.vcf", sname=tumorid, stype=stype_tumor),
                bam = "{stype}/realign/{sname}_REALIGNED.bam",
                bai = "{stype}/realign/{sname}_REALIGNED.bam.bai",
                somalier_sex = "{stype}/somalier/calculated_sex.txt",
            params:
                genomeversion = config["reference"],
                dll = pipeconfig["singularities"]["canvas"]["dll"],
                annotate = pipeconfig["rules"]["canvas"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),
                annotate_ref = pipeconfig["rules"]["canvas"]["annotate_ref"],
                genomedir = pipeconfig["singularities"]["canvas"]["reference"],
                kmerfile = pipeconfig["singularities"]["canvas"]["kmerfile"],
                run_py = pipeconfig["singularities"]["canvas"]["tool_path"],
                filter13 = pipeconfig["singularities"]["canvas"]["filter13"],
                samplename = sampleconfig["tumorname"],
                intermediate_vcf = "{stype}/canvas/{sname}_germline_CNV.vcf.gz",
                intermediate_observed = "{stype}/canvas/{sname}_germline_CNV_observed.seg",
                intermediate_called = "{stype}/canvas/{sname}_germline_CNV_called.seg"
            singularity:
                pipeconfig["singularities"]["canvas"]["sing"]
            output:
                out_vcf = "{stype}/canvas/{sname}_canvas_tumoronly.vcf.gz",
                out_observed = "{stype}/canvas/{sname}_canvas_tumoronly_observed.seg",
                out_called = "{stype}/canvas/{sname}_canvas_tumoronly_called.seg"
            shadow:
                pipeconfig["rules"].get("canvas", {}).get("shadow", pipeconfig.get("shadow", False))
            shell:
                """
                echo $HOSTNAME;
                SEX=$(cat {input.somalier_sex})
                {params.run_py} --genomeversion {params.genomeversion} --bam {input.bam} --normal_vcf {input.germline_snv_vcf} --o {wildcards.stype}/canvas/ -t germline --samplename {wildcards.sname} --sex "$SEX" --referencedir {params.genomedir} --kmerfile {params.kmerfile} --canvasdll {params.dll} --filterfile {params.filter13}
                mv {params.intermediate_vcf} {output.out_vcf}
                mv {params.intermediate_observed} {output.out_observed}
                mv {params.intermediate_called} {output.out_called}
                """

if normalid:
    rule canvas_germline:
        input:
            germline_snv_vcf = expand("{stype}/dnascope/{sname}_DNAscope_germline_SNVsOnly.recode.vcf", sname=normalid, stype=stype_normal),
            bam = "{stype}/realign/{sname}_REALIGNED.bam",
            bai = "{stype}/realign/{sname}_REALIGNED.bam.bai",
            somalier_sex = "{stype}/somalier/calculated_sex.txt",
        params:
            genomeversion = config["reference"],
            dll = pipeconfig["singularities"]["canvas"]["dll"],
            annotate = pipeconfig["rules"]["canvas"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),
            annotate_ref = pipeconfig["rules"]["canvas"]["annotate_ref"],
            genomedir = pipeconfig["singularities"]["canvas"]["reference"],
            kmerfile = pipeconfig["singularities"]["canvas"]["kmerfile"],
            run_py = pipeconfig["singularities"]["canvas"]["tool_path"],
            filter13 = pipeconfig["singularities"]["canvas"]["filter13"],
            samplename = sampleconfig["normalname"],
            intermediate_vcf = "{stype}/canvas/{sname}_germline_CNV.vcf.gz",
            intermediate_observed = "{stype}/canvas/{sname}_germline_CNV_observed.seg",
            intermediate_called = "{stype}/canvas/{sname}_germline_CNV_called.seg"
        singularity:
            pipeconfig["singularities"]["canvas"]["sing"]
        output:
            out_vcf = "{stype}/canvas/{sname}_canvas_germline.vcf.gz",
            out_observed = "{stype}/canvas/{sname}_canvas_germline_observed.seg",
            out_called = "{stype}/canvas/{sname}_canvas_germline_called.seg"
        shadow:
            pipeconfig["rules"].get("canvas", {}).get("shadow", pipeconfig.get("shadow", False))
        shell:
            """
            echo $HOSTNAME;
            SEX=$(cat {input.somalier_sex})
            {params.run_py} --genomeversion {params.genomeversion} --bam {input.bam} --normal_vcf {input.germline_snv_vcf} --o {wildcards.stype}/canvas/ -t germline --samplename {wildcards.sname} --sex "$SEX" --referencedir {params.genomedir} --kmerfile {params.kmerfile} --canvasdll {params.dll} --filterfile {params.filter13}
            mv {params.intermediate_vcf} {output.out_vcf}
            mv {params.intermediate_observed} {output.out_observed}
            mv {params.intermediate_called} {output.out_called}
            """
