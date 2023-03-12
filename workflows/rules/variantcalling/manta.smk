# vim: syntax=python tabstop=4 expandtab
# # coding: utf-8
import os
from workflows.scripts.annotate_manta.manta_summary import manta_summary

def get_manta_input(wcs):
    inputs = dict()
    if tumorid and wcs.vartype != "germline":
        inputs |= dict(
            tumorbam = f"tumor/realign/{tumorid}_REALIGNED.bam",
            tumorbai = f"tumor/realign/{tumorid}_REALIGNED.bam.bai",
        )
    if normalid:
        inputs |= dict(
            normalbam = f"normal/realign/{normalid}_REALIGNED.bam",
            normalbai = f"normal/realign/{normalid}_REALIGNED.bam.bai",
        )

    return inputs

def get_manta_varfile(wcs):
    if wcs.vartype == "germline":
        return "diploidSV.vcf"
    elif normalid and tumorid:
        return "somaticSV.vcf"
    else:
        return "tumorSV.vcf"


def get_manta_extargs(wcs, input):
    extargs = ""
    if tumorid and wcs.vartype == "somatic":
        extargs += f"--tumorBam={input.tumorbam} "
    if normalid and tumorid or wcs.vartype == "germline":
        extargs += f"--normalBam={input.normalbam} "
    return extargs
  

rule manta:
    input:
        unpack(get_manta_input)
    params:
        annotate = pipeconfig["rules"]["manta"].get("annotate", f"{ROOT_DIR}/workflows/scripts/annotate_manta_canvas/annotate_manta_canvas.py"),
        reference = pipeconfig["referencegenome"],
        svdb = pipeconfig["rules"]["manta"]["svdb"],
        mantaconf = pipeconfig["rules"]["manta"]["mantaconf"],
        varfile = get_manta_varfile,
        annotate_ref = pipeconfig["rules"]["manta"]["annotate_ref"],
        swegendb = pipeconfig["rules"]["manta"]["swegendb"] if reference == "hg19" else None,
        gnomaddb = pipeconfig["rules"]["manta"]["gnomaddb"] if reference == "hg19" else None,
        localdb = pipeconfig["rules"]["manta"]["localdb"] if reference == "hg19" else None,
        bcftools = pipeconfig["rules"]["manta"]["bcftools"],
        extargs = get_manta_extargs
    output:
        sv_vcf = temp("{stype}/manta/{sname}_{vartype}_mantaSV.vcf"),
        sv_xlsx = temp("{stype}/manta/{sname}_{vartype}_mantaSV.vcf.xlsx"),
        bnd_vcf = temp("{stype}/manta/{sname}_{vartype}_MantaBNDs.vcf"),
        nobnd_vcf = temp("{stype}/manta/{sname}_{vartype}_MantaNOBNDs.vcf")
    shadow:
        pipeconfig["rules"].get("manta", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        (
            "{params.mantaconf} {params.extargs} --referenceFasta {params.reference} --runDir {wildcards.stype}/manta/; " # Preparing Manta
            "{wildcards.stype}/manta/runWorkflow.py -m local; " #Running Manta
            "gunzip {wildcards.stype}/manta/results/variants/{params.varfile}.gz -c "
                "> {wildcards.stype}/manta/results/variants/SV.vcf; "
            "grep -e $'\t'PASS$'\t' -e '^#' {wildcards.stype}/manta/results/variants/SV.vcf "
                "> {wildcards.stype}/manta/results/variants/SV_PASS.vcf; "
            + (
                "{params.svdb} --query --query_vcf {wildcards.stype}/manta/results/variants/SV_PASS.vcf "
                    "--db {params.swegendb} --out_frq SWEFRQ --in_frq FRQ --out_occ SWEOCC --in_occ OCC "
                    "> {wildcards.stype}/manta/results/variants/SV_PASS_swegen.vcf; "
                "{params.svdb} --query --query_vcf {wildcards.stype}/manta/results/variants/SV_PASS_swegen.vcf "
                    "--db {params.gnomaddb} --out_frq GNOMAD_AF --in_frq AF --out_occ GNOMAD_AC --in_occ AC "
                    "> {wildcards.stype}/manta/results/variants/SV_PASS_swegen_gnomad.vcf; "
                "{params.svdb} --query --query_vcf {wildcards.stype}/manta/results/variants/SV_PASS_swegen_gnomad.vcf "
                    "--sqdb {params.localdb} --out_frq KGGFRQ --out_occ KGGOCC "
                    "> {wildcards.stype}/manta/results/variants/SV_PASS_swegen_gnomad_kgg.vcf; "
                "{params.bcftools} filter "
                    "-e 'INFO/GNOMAD_AF >= 0.03 | INFO/SWEFRQ >= 0.03 | INFO/KGGFRQ >= 0.05' "
                    "{wildcards.stype}/manta/results/variants/SV_PASS_swegen_gnomad_kgg.vcf "
                    "> {wildcards.stype}/manta/results/variants/SV_PASS_swegen_gnomad_kgg_freqfiltered.vcf; "
                "grep -Ev 'GL000|hs37d5' {wildcards.stype}/manta/results/variants/SV_PASS_swegen_gnomad_kgg_freqfiltered.vcf "
                    "> {output.sv_vcf}; "
                if reference == "hg19" else 
                "grep -Ev 'GL000|hs37d5' {wildcards.stype}/manta/results/variants/SV_PASS.vcf "
                    "> {output.sv_vcf}; "
            ) +
            "grep -e '^#' -e 'MantaBND:' {output.sv_vcf} > {output.bnd_vcf}; "
            "grep -v 'MantaBND:' {output.sv_vcf} > {output.nobnd_vcf}; "
            "{params.annotate} "
                "-v {output.sv_vcf} "
                "-g {params.annotate_ref} "
                "-o {wildcards.stype}/manta"
        )

rule manta_summary:
    input:
        "{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx"
    params:
        genelist = pipeconfig["rules"]["manta_summary"]["genelist" if normalid else "genelist_tumoronly"]
    output:
        temp("{stype}/manta/{sname}_somatic_mantaSV_Summary.xlsx")
    shadow:
        pipeconfig["rules"].get("manta_summary", {}).get("shadow", pipeconfig.get("shadow", False))
    run:
        manta_summary(input, output, tumorname, f"{params.genelist}", normalname)