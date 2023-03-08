# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
from shutil import copyfile


if tumorid:
    if normalid:
        rule share_to_resultdir:
            input:
                expand("{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.vcf", stype=sampleconfig[tumorname]["stype"], caller="tnscope", sname=tumorid, vcftype="somatic"),
                expand("{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.vcf", stype=sampleconfig[normalname]["stype"], caller="dnascope", sname=normalid, vcftype="germline"),
                expand("qc_report/{tumorname}_qc_stats.xlsx", tumorname=tumorname),
                expand("{stype}/canvas/{sname}_CNV_somatic.vcf.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/canvas/{sname}_CNV_germline.vcf.xlsx", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/manta/{sname}_somatic_mantaSV_Summary.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/dnascope/{sname}_{hgX}_SNV_CNV_germline.vcf.gz", sname=normalid, stype=sampleconfig[normalname]["stype"], hgX=reference),            
                expand("{stype}/pindel/{sname}_pindel.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/canvas/{sname}_CNV_germline.vcf", sname=normalid, vartype="germline", stype=sampleconfig[normalname]["stype"]),
                expand("{stype}/canvas/{sname}_CNV_somatic.vcf", sname=tumorid, vartype="somatic", stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/canvas/{sname}_somatic_CNV_observed.seg", vartype="somatic", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/canvas/{sname}_somatic_CNV_called.seg", vartype="somatic", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/canvas/{sname}_germline_CNV_observed.seg", vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{stype}/canvas/{sname}_germline_CNV_called.seg", vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                #expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                #expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{stype}/reports/{sname}_REALIGNED.bam.tdf", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/reports/{sname}_REALIGNED.bam.tdf",  sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{stype}/reports/{sname}_baf.igv", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/manta/{sname}_somatic_MantaBNDs.vcf", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                #expand("{stype}/manta/{sname}_somatic_MantaBNDs.vcf.gzi.csi", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                #expand("{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf.gzi.csi", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/manta/{sname}_germline_MantaBNDs.vcf", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                #expand("{stype}/manta/{sname}_germline_MantaBNDs.vcf.gzi.csi", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{stype}/manta/{sname}_germline_MantaNOBNDs.vcf", sname=normalid, stype=sampleconfig[normalname]["stype"])
                #expand("{stype}/manta/{sname}_germline_MantaNOBNDs.vcf.gz.csi", sname=normalid, stype=sampleconfig[normalname]["stype"])
            params:
                bgzip = pipeconfig["rules"]["share_to_resultdir"]["bgzip"],
                bcftools = pipeconfig["rules"]["share_to_resultdir"]["bcftools"]
            output:
                "reporting/shared_result_files.txt"
            run:
                for resultfile in input:
                    if resultfile.endswith(".vcf"):
                        if not os.path.isfile(f"{resultfile}.gz"):
                            shell(f"{params.bgzip} --stdout {resultfile} > {resultfile}.gz")
                        shell(f"{params.bcftools} index -f {resultfile}.gz")
                        resultfile = os.path.abspath(f"{resultfile}.gz")
                    filebase = os.path.basename(f"{resultfile}")
                    if resultfile.endswith("REALIGNED.bam"):
                        copyfile(f"{resultfile}", f"{filebase}")
                        copyfile(f"{resultfile}.bai", f"{filebase}.bai")
                        shell("echo {resultfile}.bai >> {output}")
                    if resultfile.endswith(".vcf.gz"):
                        if not os.path.isfile(f"{resultfile}.csi"):
                            shell(f"{params.bcftools} index -f {resultfile}")
                        copyfile(f"{resultfile}.csi", f"{filebase}.csi")
                        shell("echo {resultfile}.csi >> {output}")
                    copyfile(f"{resultfile}", f"{filebase}")
                    shell("echo {resultfile} >> {output}")
    else:
        rule share_to_resultdir:
            input:
                expand("{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.vcf", stype=sampleconfig[tumorname]["stype"], caller="tnscope", sname=tumorid, vcftype="somatic"),
                expand("qc_report/{tumorname}_qc_stats.xlsx", tumorname=tumorname),
                expand("{stype}/canvas/{sname}_CNV_germline.vcf.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/manta/{sname}_somatic_mantaSV_Summary.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/pindel/{sname}_pindel.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/canvas/{sname}_CNV_germline.vcf", sname=tumorid, vartype="germline", stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/canvas/{sname}_germline_CNV_observed.seg", vartype="germline", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/canvas/{sname}_germline_CNV_called.seg", vartype="germline", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                #expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/reports/{sname}_REALIGNED.bam.tdf", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/reports/{sname}_baf.igv", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/manta/{sname}_somatic_MantaBNDs.vcf", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                #expand("{stype}/manta/{sname}_somatic_MantaBNDs.vcf.gz.csi", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf", sname=tumorid, stype=sampleconfig[tumorname]["stype"])
                #expand("{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf.gzi.csi", sname=tumorid, stype=sampleconfig[tumorname]["stype"])
            params:
                bgzip = pipeconfig["rules"]["share_to_resultdir"]["bgzip"],
                bcftools = pipeconfig["rules"]["share_to_resultdir"]["bcftools"]
            output:
                "reporting/shared_result_files.txt"
            run:
                for resultfile in input:
                    if resultfile.endswith(".vcf"):
                        if not os.path.isfile(f"{resultfile}.gz"):
                            shell(f"{params.bgzip} --stdout {resultfile} > {resultfile}.gz")
                        shell(f"{params.bcftools} index -f {resultfile}.gz")
                        resultfile = os.path.abspath(f"{resultfile}.gz")
                    filebase = os.path.basename(f"{resultfile}")
                    if resultfile.endswith("REALIGNED.bam"):
                        copyfile(f"{resultfile}", f"{filebase}")
                        copyfile(f"{resultfile}.bai", f"{filebase}.bai")
                        shell("echo {resultfile}.bai >> {output}")
                    if resultfile.endswith(".vcf.gz"):
                        if not os.path.isfile(f"{resultfile}.csi"):
                            shell(f"{params.bcftools} index -f {resultfile}")
                        copyfile(f"{resultfile}.csi", f"{filebase}.csi")
                        shell("echo {resultfile}.csi >> {output}")
                    copyfile(f"{resultfile}", f"{filebase}")
                    shell("echo {resultfile} >> {output}")

else:
    rule share_to_resultdir:
        input:
            expand("{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.vcf", stype=sampleconfig[normalname]["stype"], caller="dnascope", sname=normalid, vcftype="germline"),
            expand("qc_report/{normalname}_qc_stats.xlsx", normalname=normalname),
            expand("{stype}/canvas/{sname}_CNV_germline.vcf.xlsx", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/dnascope/{sname}_{hgX}_SNV_CNV_germline.vcf.gz", sname=normalid, stype=sampleconfig[normalname]["stype"], hgX=reference),
            expand("{stype}/canvas/{sname}_CNV_germline.vcf", sname=normalid, vartype="germline", stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/canvas/{sname}_germline_CNV_observed.seg", vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/canvas/{sname}_germline_CNV_called.seg", vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            #expand("{stype}/realign/{sname}_REALIGNED.bam.bai", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/reports/{sname}_REALIGNED.bam.tdf", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/manta/{sname}_germline_MantaBNDs.vcf", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            #expand("{stype}/manta/{sname}_germline_MantaBNDs.vcf.gz.csi", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/manta/{sname}_germline_MantaNOBNDs.vcf", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            #expand("{stype}/manta/{sname}_germline_MantaNOBNDs.vcf.csi", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/reports/{sname}_baf.igv", sname=normalid, stype=sampleconfig[normalname]["stype"])
        params:
            bgzip = pipeconfig["rules"]["share_to_resultdir"]["bgzip"],
            bcftools = pipeconfig["rules"]["share_to_resultdir"]["bcftools"]
        output:
            "reporting/shared_result_files.txt"
        run:
                for resultfile in input:
                    if resultfile.endswith(".vcf"):
                        if not os.path.isfile(f"{resultfile}.gz"):
                            shell(f"{params.bgzip} --stdout {resultfile} > {resultfile}.gz")
                        shell(f"{params.bcftools} index -f {resultfile}.gz")
                        resultfile = os.path.abspath(f"{resultfile}.gz")
                    filebase = os.path.basename(f"{resultfile}")
                    if resultfile.endswith("REALIGNED.bam"):
                        copyfile(f"{resultfile}", f"{filebase}")
                        copyfile(f"{resultfile}.bai", f"{filebase}.bai")
                        shell("echo {resultfile}.bai >> {output}")
                    if resultfile.endswith(".vcf.gz"):
                        if not os.path.isfile(f"{resultfile}.csi"):
                            shell(f"{params.bcftools} index -f {resultfile}")
                        copyfile(f"{resultfile}.csi", f"{filebase}.csi")
                        shell("echo {resultfile}.csi >> {output}")
                    copyfile(f"{resultfile}", f"{filebase}")
                    shell("echo {resultfile} >> {output}")
