import os
from workflows.scripts.sex import calc_sex

# Helper function to resolve file paths dynamically
def get_cov_files(wildcards):
    if normalid:
        # If normal is included in the run, the Ycov and WGScov files are taken from the normal sample
        return {
            "wgscov": "{stype}/reports/{sname}_WGScov.tsv".format(
                stype=sampleconfig[normalname]["stype"], sname=normalid),
            "ycov": "{stype}/reports/{sname}_Ycov.tsv".format(
                stype=sampleconfig[normalname]["stype"], sname=normalid)
        }
    else:
        # If no normal is included, the Ycov and WGScov files are taken from the tumor sample
        return {
            "wgscov": "{stype}/reports/{sname}_WGScov.tsv".format(
                stype=sampleconfig[tumorname]["stype"], sname=tumorid),
            "ycov": "{stype}/reports/{sname}_Ycov.tsv".format(
                stype=sampleconfig[tumorname]["stype"], sname=tumorid)
        }

if normalid:
    rule ascat_run:
        input:
            tumor_bam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            normal_bam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            wgscov = lambda wildcards: get_cov_files(wildcards)["wgscov"],
            ycov = lambda wildcards: get_cov_files(wildcards)["ycov"]
        params:
            allelecounter_exe = pipeconfig["rules"]["ascat_run"]["allelecounter_exe"],
            alleles_prefix = pipeconfig["rules"]["ascat_run"]["alleles_prefix"],
            loci_prefix = pipeconfig["rules"]["ascat_run"]["loci_prefix"],
            # calc_sex outputs male / female, ascat Rscript accepts male/XY or female/XX
            # ascat will currently not calculate on chrY but will add it to the plot depending on sex
            # see also: https://github.com/VanLoo-lab/ascat/issues/125
            gender = lambda wildcards: calc_sex(
                get_cov_files(wildcards)["wgscov"],
                get_cov_files(wildcards)["ycov"]
            ),            
            genome_version = pipeconfig["rules"]["ascat_run"]["genome_version"],
            gc_content_file = pipeconfig["rules"]["ascat_run"]["gc_content_file"],
            replic_timing_file = pipeconfig["rules"]["ascat_run"]["replic_timing_file"],
            ascat_run_script = f"{ROOT_DIR}/workflows/scripts/ascat_run.R",
        output:
            # All output will be stored in the temporary output_directory except the Rdata and segments files
            output_dir = temp(directory("{stype}/ascat/{sname}_run_output")),
            rdata_file = temp("{stype}/ascat/{sname}_ascat_bc.Rdata"),
            segments = temp("{stype}/ascat/{sname}.segments.txt"),
        singularity:
            pipeconfig["singularities"]["ascat"]["sing"]
        threads:
            clusterconf["ascat_run"]["threads"]
        shell:
            """
            Rscript {params.ascat_run_script} \
                --tumor-bam {input.tumor_bam} \
                --tumor-name {wildcards.sname} \
                --normal-bam {input.normal_bam} \
                --normal-name {wildcards.sname}_normal \
                --allelecounter-exe {params.allelecounter_exe} \
                --alleles-prefix {params.alleles_prefix} \
                --loci-prefix {params.loci_prefix} \
                --gender {params.gender} \
                --genome-version {params.genome_version} \
                --nthreads {threads} \
                --gc-content-file {params.gc_content_file} \
                --replic-timing-file {params.replic_timing_file} \
                --output-dir {output.output_dir} \
                --tumoronly FALSE
            mv {output.output_dir}/{wildcards.sname}_ascat_bc.Rdata {output.rdata_file}
            mv {output.output_dir}/{wildcards.sname}.segments.txt {output.segments}
            """

else:
    rule ascat_run:
        input:
            tumor_bam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            wgscov = lambda wildcards: get_cov_files(wildcards)["wgscov"],
            ycov = lambda wildcards: get_cov_files(wildcards)["ycov"]
        params:
            allelecounter_exe = pipeconfig["rules"]["ascat_run"]["allelecounter_exe"],
            alleles_prefix = pipeconfig["rules"]["ascat_run"]["alleles_prefix"],
            loci_prefix = pipeconfig["rules"]["ascat_run"]["loci_prefix"],
            gender = lambda wildcards: calc_sex(
                get_cov_files(wildcards)["wgscov"],
                get_cov_files(wildcards)["ycov"]
            ),            
            genome_version = pipeconfig["rules"]["ascat_run"]["genome_version"],
            gc_content_file = pipeconfig["rules"]["ascat_run"]["gc_content_file"],
            replic_timing_file = pipeconfig["rules"]["ascat_run"]["replic_timing_file"],
            ascat_run_script = f"{ROOT_DIR}/workflows/scripts/ascat_run.R",
        output:
            output_dir = temp(directory("{stype}/ascat/{sname}_run_output")),
            rdata_file = temp("{stype}/ascat/{sname}_ascat_bc.Rdata"),
            segments = temp("{stype}/ascat/{sname}.segments.txt"),
        singularity:
            pipeconfig["singularities"]["ascat"]["sing"]
        threads:
            clusterconf["ascat_run"]["threads"]
        shell:
            """
            Rscript {params.ascat_run_script} \
                --tumor-bam {input.tumor_bam} \
                --tumor-name {wildcards.sname} \
                --allelecounter-exe {params.allelecounter_exe} \
                --alleles-prefix {params.alleles_prefix} \
                --loci-prefix {params.loci_prefix} \
                --gender {params.gender} \
                --genome-version {params.genome_version} \
                --nthreads {threads} \
                --gc-content-file {params.gc_content_file} \
                --replic-timing-file {params.replic_timing_file} \
                --output-dir {output.output_dir} \
                --tumoronly TRUE
            mv {output.output_dir}/{wildcards.sname}_ascat_bc.Rdata {output.rdata_file}
            mv {output.output_dir}/{wildcards.sname}.segments.txt {output.segments}
            """

rule ascat_plot:
    input:
        rdata_file = "{stype}/ascat/{sname}_ascat_bc.Rdata",
        segments = "{stype}/ascat/{sname}.segments.txt",
        wgscov = lambda wildcards: get_cov_files(wildcards)["wgscov"],
        ycov = lambda wildcards: get_cov_files(wildcards)["ycov"]
    params:
        tumorname = tumorname,
        gender = lambda wildcards: calc_sex(
                get_cov_files(wildcards)["wgscov"],
                get_cov_files(wildcards)["ycov"]
        ),        
        genome_fai = pipeconfig["referencefai"],
        ascat_plot_script = f"{ROOT_DIR}/workflows/scripts/ascat_custom_plot.R",
    output:
        plot = temp("{stype}/ascat/{sname}_ascat_plot.png"),
        seg = temp("{stype}/ascat/{sname}_ascat_copynumber_IGV.seg"),
        BAF = temp("{stype}/ascat/{sname}_ascat_BAF_IGV.seg")
    singularity:
        pipeconfig["singularities"]["ascat"]["sing"]
    shell:
        """
        Rscript {params.ascat_plot_script} \
            --tumorname {params.tumorname} \
            --gender {params.gender} \
            --genome-fai {params.genome_fai} \
            --segments {input.segments} \
            --Rdata-file {input.rdata_file} \
            --output-plot {output.plot} \
            --output-seg {output.seg} \
            --output-baf {output.BAF}
        """
