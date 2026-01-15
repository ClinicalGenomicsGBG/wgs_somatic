rule ascat_run:
    input:
        tumor_bam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=stype_tumor),
        normal_bam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=stype_normal) if normalid else [],
        somalier_sex = "{stype}/somalier/calculated_sex.txt",
        # SomalierParser outputs male / female, ascat Rscript accepts male/XY or female/XX
        # ascat will currently not calculate on chrY but will add it to the plot depending on sex
        # see also: https://github.com/VanLoo-lab/ascat/issues/125
    params:
        # alleleCounter executable is in conda bin directory in the ascat container
        allelecounter_exe = pipeconfig["rules"]["ascat_run"]["allelecounter_exe"],
        alleles_prefix = pipeconfig["rules"]["ascat_run"]["alleles_prefix"],
        loci_prefix = pipeconfig["rules"]["ascat_run"]["loci_prefix"],
        genome_version = pipeconfig["rules"]["ascat_run"]["genome_version"],
        gc_content_file = pipeconfig["rules"]["ascat_run"]["gc_content_file"],
        replic_timing_file = pipeconfig["rules"]["ascat_run"]["replic_timing_file"],
        ascat_run_script = f"{ROOT_DIR}/workflows/scripts/ascat_run.R",
        tumoronly = "TRUE" if not normalid else "FALSE",
        normal_bam_arg = lambda wildcards, input: f"--normal-bam {input.normal_bam}" if input.normal_bam else "",
        normal_name_arg = lambda wildcards, input: f"--normal-name {normalid}" if normalid else "",
        vstamp = f"{VDIR}/ascat_run.txt"
    output:
        # All output will be stored in the temporary output_directory. 
        # The Rdata and segments files are moved to the below output locations for further processing.
        output_dir = temp(directory("{stype}/ascat/{sname}_run_output")),
        rdata_file = temp("{stype}/ascat/{sname}_ascat_bc.Rdata"),
        stats = temp("{stype}/ascat/{sname}_ascat_stats.tsv"),
    singularity:
        pipeconfig["singularities"]["ascat"]["sing"]
    threads:
        clusterconf["ascat_run"]["threads"]
    shell:
        """
        # Set cache to avoid protected files
        export XDG_CACHE_HOME="${{TMPDIR:-/tmp}}";

        # Version info
        Rscript --version 2>&1 | head -n 1 > {params.vstamp}
        Rscript -e "library(ASCAT); cat(paste('ASCAT version:', packageVersion('ASCAT')),'\n')" >> {params.vstamp}

        # Run ascat
        SEX=$(cat {input.somalier_sex})
        Rscript {params.ascat_run_script} \
            --tumor-bam {input.tumor_bam} \
            --tumor-name {wildcards.sname} \
            {params.normal_bam_arg} \
            {params.normal_name_arg} \
            --allelecounter-exe {params.allelecounter_exe} \
            --alleles-prefix {params.alleles_prefix} \
            --loci-prefix {params.loci_prefix} \
            --gender "$SEX" \
            --genome-version {params.genome_version} \
            --nthreads {threads} \
            --gc-content-file {params.gc_content_file} \
            --replic-timing-file {params.replic_timing_file} \
            --output-dir {output.output_dir} \
            --tumoronly {params.tumoronly}
        mv {output.output_dir}/{wildcards.sname}_ascat_bc.Rdata {output.rdata_file}
        mv {output.output_dir}/{wildcards.sname}_ascat_stats.tsv {output.stats}
        """

rule ascat_plot:
    input:
        rdata_file = "{stype}/ascat/{sname}_ascat_bc.Rdata",
        somalier_sex = "{stype}/somalier/calculated_sex.txt"
    params:
        tumorname = tumorname,
        genome_fai = pipeconfig["referencefai"],
        ascat_plot_script = f"{ROOT_DIR}/workflows/scripts/ascat_custom_plot.R",
        cytoBandIdeo = pipeconfig["rules"]["ascat_run"]["cytoBandIdeo"],
        vstamp = f"{VDIR}/ascat_plot.txt"
    output:
        plot = report("{stype}/ascat/{sname}_ascat_plot.pdf", labels={"result": "ascat plot"}),
        seg_smooth = "{stype}/ascat/{sname}_ascat_CN_smooth_IGV.seg",
        seg_call = "{stype}/ascat/{sname}_ascat_CN_call_IGV.seg",
        BAF = "{stype}/ascat/{sname}_ascat_BAF_IGV.seg",
    singularity:
        pipeconfig["singularities"]["ascat"]["sing"]
    shell:
        """
        # Set cache to avoid protected files
        export XDG_CACHE_HOME="${{TMPDIR:-/tmp}}";

        # Version info
        Rscript --version 2>&1 | head -n 1 > {params.vstamp}
        Rscript -e "library(ggplot2); cat(paste('ggplot2 version:', packageVersion('ggplot2')),'\n')" >> {params.vstamp}
            
        # Plot ascat
        SEX=$(cat {input.somalier_sex})
        Rscript {params.ascat_plot_script} \
            --tumorname {params.tumorname} \
            --gender "$SEX" \
            --genome-fai {params.genome_fai} \
            --Rdata-file {input.rdata_file} \
            --cytoband {params.cytoBandIdeo} \
            --output-plot {output.plot} \
            --output-seg-smooth {output.seg_smooth} \
            --output-seg-call {output.seg_call} \
            --output-baf {output.BAF}
        """
