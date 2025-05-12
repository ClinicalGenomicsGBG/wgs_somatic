if normalid:
    rule control_freec_run:
        input:
            tumor_pileup = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            normal_pileup = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
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
            config = temp("{stype}/control-freec_{ploidy}/{sname}.config"),
            tumor_ratio = temp("{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_ratio.txt"),
            #normal_ratio = temp("{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_normal_ratio.txt"),
            info = temp("{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_info.txt"),
        #shadow:
        #    pipeconfig["rules"].get("control-freec", {}).get("shadow", pipeconfig.get("shadow", False))
        shell:
            """
            python {params.edit_config} {params.config_template} {wildcards.ploidy} {input.tumor_pileup} {input.normal_pileup} {params.chrLenFile} {params.chrFiles} {params.mappability} {threads} {output.config}
            freec -conf {output.config}
            """

    rule control_freec_plot:
        input:
            ratio = "{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_ratio.txt",
            #normal_ratio = "{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_normal_ratio.txt",
            BAF = expand("{stype}/reports/{sname}_baf.igv", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            normal_BAF = expand("{stype}/reports/{sname}_baf.igv", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        params:
            plot_script = pipeconfig["rules"]["control-freec"].get("plot_script", f"{ROOT_DIR}/workflows/scripts/control_freec_makeGraph3.0.R"),
            fai =  pipeconfig["referencefai"],
            cytoBandIdeo = pipeconfig["rules"]["control-freec"]["cytoBandIdeo"],
        singularity:
            pipeconfig["singularities"]["control-freec"]["sing"]
        output:
            ratio_plot = temp("{stype}/control-freec_{ploidy}/{sname}_controlfreec_ploidy{ploidy}.png"),
            ratio_seg = temp("{stype}/control-freec_{ploidy}/{sname}_controlfreec_ploidy{ploidy}.seg"),
            #normal_ratio_plot = temp("{stype}/control-freec_{ploidy}/{sname}_ploidy{ploidy}_normal_ratio.png"),
            #normal_ratio_seg = temp("{stype}/control-freec_{ploidy}/{sname}_ploidy{ploidy}_normal_ratio.seg"),
        shadow:
            pipeconfig["rules"].get("control-freec", {}).get("shadow", pipeconfig.get("shadow", False))
        shell:
            """
            Rscript {params.plot_script} {wildcards.sname} {input.ratio} {input.BAF} {params.fai} {params.cytoBandIdeo} {output.ratio_plot} {output.ratio_seg}
            """
            #Rscript {params.plot_script} {wildcards.sname} {input.normal_ratio} {input.normal_BAF} {params.fai} {params.cytoBandIdeo} {output.normal_ratio_plot} {output.normal_ratio_seg}
            


else:
    rule control_freec_run:
        input:
            tumor_pileup = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        params:
            config_template = pipeconfig["rules"]["control-freec"].get("config_template", f"{ROOT_DIR}/workflows/scripts/control_freec_config.txt"),
            edit_config = pipeconfig["rules"]["control-freec"].get("edit_config", f"{ROOT_DIR}/workflows/scripts/control_freec_edit_config.py"),
            chrLenFile = pipeconfig["rules"]["control-freec"]["chrLenFile"],
            chrFiles = pipeconfig["rules"]["control-freec"]["chrFiles"],
            mappability = pipeconfig["rules"]["control-freec"]["mappability"],
            normal_pileup = "None",
        singularity:
            pipeconfig["singularities"]["control-freec"]["sing"]
        threads:
            clusterconf["control_freec_run"]["threads"]
        output:
            config = temp("{stype}/control-freec_{ploidy}/{sname}.config"),
            tumor_ratio = temp("{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_ratio.txt"),
            tumor_BAF = temp("{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_BAF.txt"),
            info = temp("{stype}/control-freec_{ploidy}/{sname}_REALIGNED.bam_info.txt"),
        shadow:
            pipeconfig["rules"].get("control-freec", {}).get("shadow", pipeconfig.get("shadow", False))
        shell:
            """
            python {params.edit_config} {params.config_template} {wildcards.ploidy} {input.tumor_pileup} {params.normal_pileup} {params.chrLenFile} {params.chrFiles} {params.mappability} {threads} {output.config}
            freec -conf {output.config}
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
            ratio_plot = temp("{stype}/control-freec_{ploidy}/{sname}_controlfreec_ploidy{ploidy}.png"),
            ratio_seg = temp("{stype}/control-freec_{ploidy}/{sname}_controlfreec_ploidy{ploidy}.seg"),
        shadow:
            pipeconfig["rules"].get("control-freec", {}).get("shadow", pipeconfig.get("shadow", False))
        shell:
            """
            Rscript {params.plot_script} {wildcards.sname} {input.ratio} {input.BAF} {params.fai} {params.cytoBandIdeo} {output.ratio_plot} {output.ratio_seg}
            """