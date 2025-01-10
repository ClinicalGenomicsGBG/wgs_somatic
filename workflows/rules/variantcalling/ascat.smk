rule ascat_allele_count:
    input:
        bam = "{stype}/realign/{sname}_REALIGNED.bam"
    singularity:
        pipeconfig['rules']['ascat']['singularity']
    params:
        allele_count_script = f"{ROOT_DIR}/workflows/scripts/ascat/allele_count.R",
        loci_prefix = f"{pipeconfig['rules']['ascat']['resources']}/G1000_loci_hg38_chr",
        allelecounter_exe = "/opt/conda/bin/alleleCounter"
    output:
        allele_freq = temp("{stype}/ascat/{sname}_alleleFrequencies_chr{chr_num}.txt")
    shell:
        "Rscript {params.allele_count_script} {input.bam} {wildcards.sname} {wildcards.chr_num} {output.allele_freq} {params.loci_prefix} {params.allelecounter_exe}"

if tumorid and normalid:
    rule ascat_get_BAF_LogR:
        input:
            in_tumor = expand("{stype}/ascat/{sname}_alleleFrequencies_chr{chr_num}.txt", chr_num=[str(i) for i in range(1, 23)] + ["X"], sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            in_normal = expand("{stype}/ascat/{sname}_alleleFrequencies_chr{chr_num}.txt", chr_num=[str(i) for i in range(1, 23)] + ["X"], sname=normalid, stype=sampleconfig[normalname]["stype"]),
        singularity:
            pipeconfig['rules']['ascat']['singularity']
        params:
            get_BAF_LogR_script = f"{ROOT_DIR}/workflows/scripts/ascat/get_BAF_LogR.R",
            tumorid = tumorid,
            allele_prefix = f"{pipeconfig['rules']['ascat']['resources']}/G1000_alleles_hg38_chr",
            genome_id = f"{pipeconfig['rules']['ascat']['genome_id']}",
        output:
            raw_BAF_file = temp("{stype}/ascat/{sname}_BAF_rawBAF.txt"),
            BAF_file = temp("{stype}/ascat/{sname}_BAF.txt"),
            LogR_file = temp("{stype}/ascat/{sname}_LogR.txt"),
        shell:
            "Rscript {params.get_BAF_LogR_script} {params.tumorid} {input.in_tumor[0]} {input.in_normal[0]} {output.BAF_file} {output.LogR_file} {params.allele_prefix} {wildcards.stype} {params.genome_id}"

if tumorid and not normalid:
    rule ascat_get_BAF_LogR:
        input:
            in_tumor = expand("{stype}/ascat/{sname}_alleleFrequencies_chr{chr_num}.txt", chr_num=[str(i) for i in range(1, 23)] + ["X"], sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        singularity:
            pipeconfig['rules']['ascat']['singularity']
        params:
            get_BAF_LogR_script = f"{ROOT_DIR}/workflows/scripts/ascat/get_BAF_LogR.R",
            tumorid = tumorid,
            allele_prefix = f"{pipeconfig['rules']['ascat']['resources']}/G1000_alleles_hg38_chr",
            tumoronly = "TRUE",
            genome_id = f"{pipeconfig['rules']['ascat']['genome_id']}",
        output:
            raw_BAF_file = temp("{stype}/ascat/{sname}_BAF_rawBAF.txt"),
            BAF_file = temp("{stype}/ascat/{sname}_BAF.txt"),
            LogR_file = temp("{stype}/ascat/{sname}_LogR.txt"),
        shell:
            "Rscript {params.get_BAF_LogR_script} {params.tumorid} {input.in_tumor[0]} {input.in_tumor[0]} {output.BAF_file} {output.LogR_file} {params.allele_prefix} {wildcards.stype} {params.genome_id} {params.tumoronly}"

if tumorid and normalid:
    rule ascat_run_main:
        input:
            BAF_tumor_file = expand("{stype}/ascat/{sname}_BAF.txt", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            LogR_tumor_file = expand("{stype}/ascat/{sname}_LogR.txt", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            BAF_normal_file = expand("{stype}/ascat/{sname}_BAF.txt", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            LogR_normal_file = expand("{stype}/ascat/{sname}_LogR.txt", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        singularity:
            pipeconfig['rules']['ascat']['singularity']
        params:
            tumorid = tumorid,
            run_ascat_script = f"{ROOT_DIR}/workflows/scripts/ascat/run_ascat.R",
            GCcontentfile = f"{pipeconfig['rules']['ascat']['resources']}/GC_G1000_hg38.txt",
            replictimingfile = f"{pipeconfig['rules']['ascat']['resources']}/RT_G1000_hg38.txt",
            genome_id = f"{pipeconfig['rules']['ascat']['genome_id']}",
            genome_fai = f"{pipeconfig['referencefai']}"
        output:
            stats_file = temp("{stype}/ascat/{sname}_ascat_stats.tsv"),
            sunrise_plot = temp("{stype}/ascat/{sname}.sunrise.png"),
            ASCATprofile_plot = temp("{stype}/ascat/{sname}.ASCATprofile.png"),
            ASPCF_plot = temp("{stype}/ascat/{sname}.ASPCF.png"),
            rawprofile_plot = temp("{stype}/ascat/{sname}.rawprofile.png"),
            segments_raw = temp("{stype}/ascat/{sname}.segments_raw.txt"),
            segments = temp("{stype}/ascat/{sname}.segments.txt"),
            after_corr_germline = temp("{stype}/ascat/After_correction_{sname}.germline.png"),
            after_corr_tumour = temp("{stype}/ascat/After_correction_{sname}.tumour.png"),
            before_corr_germline = temp("{stype}/ascat/Before_correction_{sname}.germline.png"),
            before_corr_tumour = temp("{stype}/ascat/Before_correction_{sname}.tumour.png"),
            ascat_plot = temp("{stype}/ascat/{sname}.ascat_out.png"),
            ascat_interactive = temp("{stype}/ascat/{sname}.ascat_interactive.html"),
        shell:
            "Rscript {params.run_ascat_script} {params.tumorid} {input.LogR_tumor_file} {input.BAF_tumor_file} {input.LogR_normal_file} {input.BAF_normal_file} {params.GCcontentfile} {params.replictimingfile} {output.stats_file} {params.genome_id} {params.genome_fai}"

if tumorid and not normalid:
    rule ascat_run_main:
        input:
            BAF_tumor_file = expand("{stype}/ascat/{sname}_BAF.txt", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            LogR_tumor_file = expand("{stype}/ascat/{sname}_LogR.txt", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        singularity:
            pipeconfig['rules']['ascat']['singularity']
        params:
            tumorid = tumorid,
            run_ascat_script = f"{ROOT_DIR}/workflows/scripts/ascat/run_ascat.R",
            GCcontentfile = f"{pipeconfig['rules']['ascat']['resources']}/GC_G1000_hg38.txt",
            replictimingfile = f"{pipeconfig['rules']['ascat']['resources']}/RT_G1000_hg38.txt",
            genome_id = f"{pipeconfig['rules']['ascat']['genome_id']}",
            genome_fai = f"{pipeconfig['referencefai']}"
        output:
            stats_file = temp("{stype}/ascat/{sname}_ascat_stats.tsv"),
            sunrise_plot = temp("{stype}/ascat/{sname}.sunrise.png"),
            ASCATprofile_plot = temp("{stype}/ascat/{sname}.ASCATprofile.png"),
            ASPCF_plot = temp("{stype}/ascat/{sname}.ASPCF.png"),
            rawprofile_plot = temp("{stype}/ascat/{sname}.rawprofile.png"),
            segments_raw = temp("{stype}/ascat/{sname}.segments_raw.txt"),
            segments = temp("{stype}/ascat/{sname}.segments.txt"),
            after_corr_germline = temp("{stype}/ascat/After_correction_{sname}.germline.png"),
            after_corr_tumour = temp("{stype}/ascat/After_correction_{sname}.tumour.png"),
            before_corr_germline = temp("{stype}/ascat/Before_correction_{sname}.germline.png"),
            before_corr_tumour = temp("{stype}/ascat/Before_correction_{sname}.tumour.png"),
            ascat_plot = temp("{stype}/ascat/{sname}.ascat_out.png"),
            ascat_interactive = temp("{stype}/ascat/{sname}.ascat_interactive.html"),
        shell:
            "Rscript {params.run_ascat_script} {params.tumorid} {input.LogR_tumor_file} {input.BAF_tumor_file} {input.LogR_tumor_file} {input.BAF_tumor_file} {params.GCcontentfile} {params.replictimingfile} {output.stats_file} {params.genome_id} {params.genome_fai}"
