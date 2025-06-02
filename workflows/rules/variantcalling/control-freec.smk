from tools.helpers import conditional_temp

if normalid:
    rule control_freec_run:
        input:
            tumor_bam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            tumor_bai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            normal_bam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            normal_bai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            wgscovfile = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[normalname]["stype"], sname=normalid),
            ycov = expand("{stype}/reports/{sname}_Ycov.tsv", stype=sampleconfig[normalname]["stype"], sname=normalid),
        params:
            config_template = pipeconfig["rules"]["control-freec"].get("config_template", f"{ROOT_DIR}/workflows/scripts/control_freec_config.txt"),
            edit_config = pipeconfig["rules"]["control-freec"].get("edit_config", f"{ROOT_DIR}/workflows/scripts/control_freec_edit_config.py"),
            chrLenFile = pipeconfig["rules"]["control-freec"]["chrLenFile"],
            chrFiles = pipeconfig["rules"]["control-freec"]["chrFiles"],
            mappability = pipeconfig["rules"]["control-freec"]["mappability"],
            normalid = normalid,
        singularity:
            pipeconfig["singularities"]["control-freec"]["sing"]
        threads:
            clusterconf["control_freec_run"]["threads"]
        output:
            config = conditional_temp("{stype}/control-freec_{ploidy}/{sname}.config", keepfiles),
            tumor_ratio = conditional_temp("{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_ratio.txt", keepfiles),
            info = conditional_temp("{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_info.txt", keepfiles),
        shadow:
            pipeconfig["rules"].get("control-freec", {}).get("shadow", pipeconfig.get("shadow", False))
        shell:
            # The output files are created by freec, but we also create them if freec fails
            """
            set -e
            (python {params.edit_config} {params.config_template} {wildcards.ploidy} {input.tumor_bam} {input.normal_bam} {params.chrLenFile} {params.chrFiles} {params.mappability} {threads} {output.config} {input.wgscovfile} {input.ycov} &&
            freec -conf {output.config}) || (touch {output.config} {output.tumor_ratio} {output.info}; echo "[control-freec] WARNING: FREEC failed with ploidy {wildcards.ploidy}, creating empty outputs and continuing." >&2)
            """

else:
    rule control_freec_run:
        input:
            tumor_bam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            tumor_bai = expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            wgscovfile = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[tumorname]["stype"], sname=tumorid),  # For tumor only, ycov and wgscov are from tumor
            ycov = expand("{stype}/reports/{sname}_Ycov.tsv", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
        params:
            config_template = pipeconfig["rules"]["control-freec"].get("config_template", f"{ROOT_DIR}/workflows/scripts/control_freec_config.txt"),
            edit_config = pipeconfig["rules"]["control-freec"].get("edit_config", f"{ROOT_DIR}/workflows/scripts/control_freec_edit_config.py"),
            chrLenFile = pipeconfig["rules"]["control-freec"]["chrLenFile"],
            chrFiles = pipeconfig["rules"]["control-freec"]["chrFiles"],
            mappability = pipeconfig["rules"]["control-freec"]["mappability"],
            normal_bam = "None",
        singularity:
            pipeconfig["singularities"]["control-freec"]["sing"]
        threads:
            clusterconf["control_freec_run"]["threads"]
        output:
            config = conditional_temp("{stype}/control-freec_{ploidy}/{sname}.config", keepfiles),
            tumor_ratio = conditional_temp("{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_ratio.txt", keepfiles),
            info = conditional_temp("{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_info.txt", keepfiles),
        shadow:
            pipeconfig["rules"].get("control-freec", {}).get("shadow", pipeconfig.get("shadow", False))
        shell:
            """
            set -e
            (python {params.edit_config} {params.config_template} {wildcards.ploidy} {input.tumor_bam} {params.normal_bam} {params.chrLenFile} {params.chrFiles} {params.mappability} {threads} {output.config} {input.wgscovfile} {input.ycov} &&
            freec -conf {output.config}) || (touch {output.config} {output.tumor_ratio} {output.info}; echo "[control-freec] WARNING: FREEC failed with ploidy {wildcards.ploidy}, creating empty outputs and continuing." >&2)
            """

rule control_freec_plot:
    input:
        ratio = "{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_ratio.txt",
        BAF = expand("{stype}/reports/{sname}_baf.igv", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
    params:
        plot_script = pipeconfig["rules"]["control-freec"].get("plot_script", f"{ROOT_DIR}/workflows/scripts/control_freec_makeGraph3.0.R"),
        fai =  pipeconfig["referencefai"],
        cytoBandIdeo = pipeconfig["rules"]["control-freec"]["cytoBandIdeo"],
    singularity:
        pipeconfig["singularities"]["control-freec"]["sing"]
    output:
        ratio_plot = conditional_temp("{stype}/control-freec_{ploidy}/{sname}_controlfreec_ploidy{ploidy}.png", keepfiles),
        ratio_seg = conditional_temp("{stype}/control-freec_{ploidy}/{sname}_controlfreec_ploidy{ploidy}.seg", keepfiles),
    shell:
        # The Rscript accepts empty input files, so we can run it even if the previous step failed
        """
        Rscript {params.plot_script} {wildcards.sname} {input.ratio} {input.BAF} {params.fai} {params.cytoBandIdeo} {output.ratio_plot} {output.ratio_seg}
        """