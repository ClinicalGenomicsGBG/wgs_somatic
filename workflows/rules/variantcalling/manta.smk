# vim: syntax=python tabstop=4 expandtab
# # coding: utf-8
import os
from workflows.scripts.annotate_manta.manta_summary import manta_summary
from workflows.scripts.filter_manta import filter_vcf

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
            min_tumor_support = filterconfig["sv_filter"]["min_tumor_support"],
            max_normal_support = filterconfig["sv_filter"]["max_normal_support"],
        output:
            sv_vcf = "{stype}/manta/{sname}_manta_somatic.vcf",
            sv_xlsx = "{stype}/manta/{sname}_manta_somatic.vcf.xlsx",
        shadow:
            pipeconfig["rules"].get("manta", {}).get("shadow", pipeconfig.get("shadow", False))
        run:
            if not os.path.isfile(f"{wildcards.stype}/manta/runWorkflow.py"):
                shell("{params.mantaconf} --tumorBam={input.tumorbam} --normalBam={input.normalbam} --referenceFasta {params.reference} --runDir {wildcards.stype}/manta/") # Preparing Manta
            if not os.path.isfile(f"{wildcards.stype}/manta/results/variants/somaticSV.vcf.gz"):
                shell("{wildcards.stype}/manta/runWorkflow.py -m local") #Running Manta
            if not os.path.isfile(f"{wildcards.stype}/manta/results/variants/somaticSV.vcf"):
                shell("gunzip {wildcards.stype}/manta/results/variants/somaticSV.vcf.gz")
            filter_vcf(
                f"{wildcards.stype}/manta/results/variants/somaticSV.vcf",
                f"{output.sv_vcf}",
                tumor_name=tumorname,
                normal_name=normalname,
                min_tumor_support=params.min_tumor_support,
                max_normal_support=params.max_normal_support
                )
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
            min_tumor_support = filterconfig["sv_filter"]["min_tumor_support"],
        output:
            sv_vcf = "{stype}/manta/{sname}_manta_somatic.vcf",
            sv_xlsx = "{stype}/manta/{sname}_manta_somatic.vcf.xlsx",
        shadow:
            pipeconfig["rules"].get("manta", {}).get("shadow", pipeconfig.get("shadow", False))
        run:
            if not os.path.isfile(f"{wildcards.stype}/manta/runWorkflow.py"):
                shell("{params.mantaconf} --tumorBam={input.tumorbam} --referenceFasta {params.reference} --runDir {wildcards.stype}/manta/") # Preparing Manta
            if not os.path.isfile(f"{wildcards.stype}/manta/results/variants/tumorSV.vcf.gz"):
                shell("{wildcards.stype}/manta/runWorkflow.py -m local") #Running Manta
            if not os.path.isfile(f"{wildcards.stype}/manta/results/variants/tumorSV.vcf"):
                shell("gunzip {wildcards.stype}/manta/results/variants/tumorSV.vcf.gz")
            filter_vcf(
                f"{wildcards.stype}/manta/results/variants/tumorSV.vcf",
                f"{output.sv_vcf}",
                tumor_name=tumorname,
                min_tumor_support=params.min_tumor_support,
                )
            shell("{params.annotate} -v {output.sv_vcf} -g {params.annotate_ref} -o {wildcards.stype}/manta")

if normalid:
    rule manta_summary:
        input:
            "{stype}/manta/{sname}_manta_somatic.vcf.xlsx"
        params:
            genelist = pipeconfig["rules"]["manta_summary"]["genelist"]
        output:
            "{stype}/manta/{sname}_manta_somatic_Summary.xlsx"
        shadow:
            pipeconfig["rules"].get("manta_summary", {}).get("shadow", pipeconfig.get("shadow", False))
        run:
            manta_summary(input, output, tumorname, f"{params.genelist}", normalname)

else:
    rule manta_summary:
        input:
            "{stype}/manta/{sname}_manta_somatic.vcf.xlsx"
        params:
            genelist = pipeconfig["rules"]["manta_summary"]["genelist_tumoronly"]
        output:
            "{stype}/manta/{sname}_manta_somatic_Summary.xlsx"
        shadow:
            pipeconfig["rules"].get("manta_summary", {}).get("shadow", pipeconfig.get("shadow", False))
        run:
            manta_summary(input, output, tumorname, f"{params.genelist}", normalname='')
