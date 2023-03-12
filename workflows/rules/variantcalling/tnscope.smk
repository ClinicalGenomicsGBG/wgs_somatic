# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

def get_tnscope_input(wcs):
    inputs = dict(
        tumorbam = f"tumor/realign/{tumorid}_REALIGNED.bam",
        tumorbai = f"tumor/realign/{tumorid}_REALIGNED.bam.bai",
        tumortable = f"tumor/recal/{tumorid}_RECAL_DATA.TABLE",
        
    )
    if normalid:
        inputs |= dict(
            normalbam = f"normal/realign/{normalid}_REALIGNED.bam",
            normalbai = f"normal/realign/{normalid}_REALIGNED.bam.bai",
            normaltable = f"normal/recal/{normalid}_RECAL_DATA.TABLE",
        )
    return inputs


def get_tnscopevcf_ml(wcs):
    if normalid:
        return "tumor/tnscope/{sname}_TNscope_tn_ML.vcf"
    else:
        return "tumor/tnscope/{sname}_TNscope_tn.vcf"


rule tnscope:
    input:
        unpack(get_tnscope_input)
    params:
        threads = clusterconf["tnscope"]["threads"],
        tumorname = tumorname,
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        reference = pipeconfig["singularities"]["sentieon"]["reference"],
        dbsnp = pipeconfig["singularities"]["sentieon"]["dbsnp"],
        callsettings = pipeconfig["rules"]["tnscope"]["settings"],
        normalname = normalname if normalid else None,
        pon = pipeconfig["rules"]["tnscope"]["pon"] if not normalid else None,
    singularity:
        pipeconfig["singularities"]["sentieon"]["simg"]
    output:
        tnscope_vcf = temp("{stype}/tnscope/{sname}_TNscope_tn.vcf"),
        tnscope_idx = temp("{stype}/tnscope/{sname}_TNscope_tn.vcf.idx"),
        tnscope_bam = temp("{stype}/tnscope/{sname}_REALIGNED_realignedTNscope.bam"),
        tnscope_bai = temp("{stype}/tnscope/{sname}_REALIGNED_realignedTNscope.bam.bai")
    wildcard_constraints:
        stype = "tumor"
    shadow:
        pipeconfig["rules"].get("tnscope", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        (
            "echo $HOSTNAME; "
            "{params.sentieon} driver "
                "-t {params.threads} "
                "-r {params.reference} "
                "-i {input.tumorbam} "
                "-q {input.tumortable} "
                + (
                    "-i {input.normalbam} "
                    "-q {input.normaltable} "
                    if normalid else
                    ""
                ) +
                "--algo TNscope "
                "--tumor_sample {params.tumorname} "
                + (
                    "--normal_sample {params.normalname} "
                    if normalid else
                    "--pon {params.pon}"
                ) +
                "--bam_output {output.tnscope_bam} "
                "{params.callsettings} "
                "{output.tnscope_vcf}"
        )

rule tnscope_modelfilter:
    input:
        tnscopevcf = "tumor/tnscope/{sname}_TNscope_tn.vcf",
        tnscopeidx = "tumor/tnscope/{sname}_TNscope_tn.vcf.idx"
    params:
        threads = clusterconf["tnscope_modelfilter"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        reference = pipeconfig["singularities"]["sentieon"]["reference"],
        modelpath = pipeconfig["singularities"]["sentieon"]["tnscope_m"],
    singularity:
        pipeconfig["singularities"]["sentieon"]["simg"]
    output:
        vcf = temp("{stype}/tnscope/{sname}_TNscope_tn_ML.vcf"),
        idx = temp("{stype}/tnscope/{sname}_TNscope_tn_ML.vcf.idx")
    wildcard_constraints:
        stype = "tumor"
    shadow:
        pipeconfig["rules"].get("tnscope_modelfilter", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        "echo $HOSTNAME; "
        "{params.sentieon} driver "
            "-t {params.threads} "
            "-r {params.reference} "
            "--algo TNModelApply "
            "-m {params.modelpath} "
            "-v {input.tnscopevcf} "
            "{output.vcf}"


rule tnscope_vcffilter:
    input:
        tnscopevcf_ml = get_tnscopevcf_ml
    params:
        threads = clusterconf["tnscope_vcffilter"]["threads"],
        outputdir = pipeconfig["rules"]["tnscope_vcffilter"]["outputdir"],
        bcftools = pipeconfig["rules"]["tnscope_vcffilter"]["bcftools"],
        vcfname = lambda wildcards: f"{wildcards.sname}_TNscope_tn_ML"
    output:
        somatic_n = temp("{stype}/tnscope/{sname}_somatic_w_normal.vcf"),
        somatic = temp("{stype}/tnscope/{sname}_somatic.vcf")
    wildcard_constraints:
        stype = "tumor"
    shadow:
        pipeconfig["rules"].get("tnscope_vcffilter", {}).get("shadow", pipeconfig.get("shadow", False))
    shell:
        (
            "{params.bcftools} filter -s LowQual -e 'QUAL < 1' -m + {input.tnscopevcf_ml} "
                "> {params.outputdir}/{params.vcfname}_lowqual.vcf; "
            "{params.bcftools} annotate -x FILTER/triallelic_site {params.outputdir}/{params.vcfname}_lowqual.vcf "
                "> {params.outputdir}/{params.vcfname}_triallelic.vcf; "
            "{params.bcftools} annotate -x FILTER/alt_allele_in_normal {params.outputdir}/{params.vcfname}_triallelic.vcf "
                "> {params.outputdir}/{params.vcfname}_altalleleinnormal.vcf; "
            "{params.bcftools} filter -s uncertainAF -e 'FORMAT/AF[0]<0.045 && FORMAT/AD[0:1]<4' -m + {params.outputdir}/{params.vcfname}_altalleleinnormal.vcf "
                "> {params.outputdir}/{params.vcfname}_uncertainaf.vcf; "
            + (
                "{params.bcftools} filter -s likely_artifact -e 'FORMAT/AF[0]<0.1 && FORMAT/AF[1]>0.06' -m + {params.outputdir}/{params.vcfname}_uncertainaf.vcf "
                    "> {params.outputdir}/{params.vcfname}_likelyartifact.vcf; "
                "{params.bcftools} filter -s lowAD -e 'FORMAT/AD[0:1] < 3' {params.outputdir}/{params.vcfname}_likelyartifact.vcf "
                    "> {params.outputdir}/{params.vcfname}_lowad.vcf; "
                "{params.bcftools} filter -s MLrejected -e 'INFO/ML_PROB<0.37' -m + {params.outputdir}/{params.vcfname}_lowad.vcf "
                    "> {params.outputdir}/{params.vcfname}_mladjusted.vcf; "
                "{params.bcftools} filter -s 'orientation_bias' -e 'FMT/FOXOG[0] == 1' -m + {params.outputdir}/{params.vcfname}_mladjusted.vcf "
                    "> {params.outputdir}/{params.vcfname}_foxogadj.vcf; "
                if normalid else
                "{params.bcftools} filter -s 'orientation_bias' -e 'FMT/FOXOG[0] == 1' -m + {params.outputdir}/{params.vcfname}_mladjusted.vcf "
                    "> {params.outputdir}/{params.vcfname}_foxogadj.vcf; "
            ) +
            "{params.bcftools} filter -s 'strand_bias' -e 'SOR > 3' -m + {params.outputdir}/{params.vcfname}_foxogadj.vcf "
                "> {params.outputdir}/{params.vcfname}_strandbiasadj.vcf; "
            "{params.bcftools} filter -i 'FILTER=\"PASS\"' {params.outputdir}/{params.vcfname}_strandbiasadj.vcf "
                "> {output.somatic_n}; "
            "{params.bcftools} view -s {tumorname} {output.somatic_n} "
                "> {output.somatic}; "
        )
