rule control_freec_pileup:
    input:
        bam = "{stype}/realign/{sname}_REALIGNED.bam"
    params:
        ref = pipeconfig["referencegenome"],
    singularity:
        pipeconfig["singularities"]["control-freec"]["sing"]
    threads:
        clusterconf["control_freec_pileup"]["threads"]
    output:
        pileup = temp("{stype}/realign/{sname}.pileup")
    shell:
        """
        mkdir -p {wildcards.stype}/control-freec/
        sambamba mpileup -t {threads} -o {output.pileup} {input.bam} --samtools -f {params.ref} -d 8000 -Q 0
        """

if normalid:
    rule control_freec_run:
        input:
            tumor_pileup = expand("{stype}/realign/{sname}.pileup", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            normal_pileup = expand("{stype}/realign/{sname}.pileup", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        params:
            config_template = pipeconfig["rules"]["control-freec"].get("config_template", f"{ROOT_DIR}/workflows/scripts/control_freec_config.txt"),
            edit_config = pipeconfig["rules"]["control-freec"].get("edit_config", f"{ROOT_DIR}/workflows/scripts/control_freec_edit_config.py"),
            chrLenFile = pipeconfig["rules"]["control-freec"]["chrLenFile"],
            chrFiles = pipeconfig["rules"]["control-freec"]["chrFiles"],
            SNPfile = pipeconfig["rules"]["control-freec"]["SNPfile"],
        singularity:
            pipeconfig["singularities"]["control-freec"]["sing"]
        threads:
            clusterconf["control_freec_run"]["threads"]
        output:
            config = temp("{stype}/control-freec_{ploidy}/{sname}.config"),
            tumor_ratio = temp("{stype}/control-freec_{ploidy}/{sname}.pileup_ratio.txt"),
            tumor_BAF = temp("{stype}/control-freec_{ploidy}/{sname}.pileup_BAF.txt"),
        shell:
            """
            python {params.edit_config} {params.config_template} {wildcards.ploidy} {input.tumor_pileup} {input.normal_pileup} {params.chrLenFile} {params.chrFiles} {params.SNPfile} {threads} {output.config}
            freec -conf {output.config}
            """
else:
    rule control_freec_run:
        input:
            tumor_pileup = expand("{stype}/realign/{sname}.pileup", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        params:
            config_template = pipeconfig["rules"]["control-freec"].get("config_template", f"{ROOT_DIR}/workflows/scripts/control_freec_config.txt"),
            edit_config = pipeconfig["rules"]["control-freec"].get("edit_config", f"{ROOT_DIR}/workflows/scripts/control_freec_edit_config.py"),
            chrLenFile = pipeconfig["rules"]["control-freec"]["chrLenFile"],
            chrFiles = pipeconfig["rules"]["control-freec"]["chrFiles"],
            SNPfile = pipeconfig["rules"]["control-freec"]["SNPfile"],
            normal_pileup = "None",
        singularity:
            pipeconfig["singularities"]["control-freec"]["sing"]
        threads:
            clusterconf["control_freec_run"]["threads"]
        output:
            config = temp("{stype}/control-freec_{ploidy}/{sname}.config"),
            tumor_ratio = temp("{stype}/control-freec_{ploidy}/{sname}.pileup_ratio.txt"),
            tumor_BAF = temp("{stype}/control-freec_{ploidy}/{sname}.pileup_BAF.txt"),
        shell:
            """
            python {params.edit_config} {params.config_template} {wildcards.ploidy} {input.tumor_pileup} {params.normal_pileup} {params.chrLenFile} {params.chrFiles} {params.SNPfile} {threads} {output.config}
            freec -conf {output.config}
            """

rule control_freec_plot:
    input:
        ratio = "{stype}/control-freec_{ploidy}/{sname}.pileup_ratio.txt",
        BAF = "{stype}/control-freec_{ploidy}/{sname}.pileup_BAF.txt"
    params:
        plot_script = pipeconfig["rules"]["control-freec"].get("plot_script", f"{ROOT_DIR}/workflows/scripts/control_freec_makeGraph3.0.R"),
        fai =  pipeconfig["referencefai"],
    singularity:
        pipeconfig["singularities"]["control-freec"]["sing"]
    output:
        ratio_plot = temp("{stype}/control-freec_{ploidy}/{sname}_ploidy{ploidy}_ratio.png"),
        ratio_seg = temp("{stype}/control-freec_{ploidy}/{sname}_ploidy{ploidy}_ratio.seg"),
        BAF_igv = temp("{stype}/control-freec_{ploidy}/{sname}_ploidy{ploidy}_BAF.igv"),
    shell:
        """
        Rscript {params.plot_script} {wildcards.sname} {input.ratio} {input.BAF} {params.fai} {output.ratio_plot} {output.ratio_seg} {output.BAF_igv}
        """