# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
from shutil import copyfile


if tumorid:
    if normalid:
        rule share_to_resultdir:
            input:
                #expand("{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], stype=sampleconfig[tumorname]["stype"], caller="tnscope", sname=tumorid, vcftype="somatic"),
                expand("{stype}/{caller}/{sname}_{vcftype}.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], stype=sampleconfig[tumorname]["stype"], caller="tnscope", sname=tumorid, vcftype="somatic"),
                expand("{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], stype=sampleconfig[normalname]["stype"], caller="dnascope", sname=normalid, vcftype="germline"),
                expand("qc_report/{tumorname}_qc_stats.xlsx", tumorname=tumorname),
                expand("{stype}/canvas/{sname}_CNV_somatic.vcf.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/canvas/{sname}_CNV_germline.vcf.xlsx", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/manta/{sname}_somatic_mantaSV_Summary.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/dnascope/{sname}_{hgX}_SNV_CNV_germline.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], sname=normalid, stype=sampleconfig[normalname]["stype"], hgX=reference),            
                expand("{stype}/pindel/{sname}_pindel.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/canvas/{sname}_CNV_germline.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], sname=normalid, vartype="germline", stype=sampleconfig[normalname]["stype"]),
                expand("{stype}/canvas/{sname}_CNV_somatic.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], sname=tumorid, vartype="somatic", stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/canvas/{sname}_somatic_CNV_observed.seg", vartype="somatic", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/canvas/{sname}_somatic_CNV_called.seg", vartype="somatic", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/canvas/{sname}_germline_CNV_observed.seg", vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{stype}/canvas/{sname}_germline_CNV_called.seg", vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{stype}/realign/{sname}_REALIGNED.{fmt}", fmt=["bam", "bam.bai", "cram", "cram.crai"], sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/reports/{sname}_REALIGNED.bam.tdf", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/realign/{sname}_REALIGNED.{fmt}", fmt=["bam", "bam.bai", "cram", "cram.crai"], sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{stype}/reports/{sname}_REALIGNED.bam.tdf", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{stype}/reports/{sname}_baf.igv", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/reports/{sname}_baf.igv", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{stype}/manta/{sname}_somatic_MantaBNDs.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                #expand("{stype}/manta/{sname}_somatic_MantaBNDs.vcf.gzi.csi", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/manta/{sname}_somatic_MantaNOBNDs.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                #expand("{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf.gzi.csi", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/manta/{sname}_germline_MantaBNDs.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], sname=normalid, stype=sampleconfig[normalname]["stype"]),
                #expand("{stype}/manta/{sname}_germline_MantaBNDs.vcf.gzi.csi", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{stype}/manta/{sname}_germline_MantaNOBNDs.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], sname=normalid, stype=sampleconfig[normalname]["stype"])
                #expand("{stype}/manta/{sname}_germline_MantaNOBNDs.vcf.gz.csi", sname=normalid, stype=sampleconfig[normalname]["stype"])
            output:
                "reporting/shared_result_files.txt"
            run:
                for resultfile in input:
                    filebase = os.path.basename(f"{resultfile}")
                    copyfile(f"{resultfile}", f"{filebase}")
                    shell("echo {resultfile} >> {output}")
    else:
        rule share_to_resultdir:
            input:
                expand("{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], stype=sampleconfig[tumorname]["stype"], caller="tnscope", sname=tumorid, vcftype="somatic"),
                expand("qc_report/{tumorname}_qc_stats.xlsx", tumorname=tumorname),
                expand("{stype}/canvas/{sname}_CNV_germline.vcf.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/manta/{sname}_somatic_mantaSV_Summary.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/pindel/{sname}_pindel.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/canvas/{sname}_CNV_germline.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], sname=tumorid, vartype="germline", stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/canvas/{sname}_germline_CNV_observed.seg", vartype="germline", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/canvas/{sname}_germline_CNV_called.seg", vartype="germline", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/realign/{sname}_REALIGNED.{fmt}", fmt=["bam", "bam.bai", "cram", "cram.crai"], sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/reports/{sname}_REALIGNED.bam.tdf", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/reports/{sname}_baf.igv", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/manta/{sname}_somatic_MantaBNDs.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                #expand("{stype}/manta/{sname}_somatic_MantaBNDs.vcf.gz.csi", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{stype}/manta/{sname}_somatic_MantaNOBNDs.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], sname=tumorid, stype=sampleconfig[tumorname]["stype"])
                #expand("{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf.gzi.csi", sname=tumorid, stype=sampleconfig[tumorname]["stype"])
            output:
                "reporting/shared_result_files.txt"
            run:
                for resultfile in input:
                    filebase = os.path.basename(f"{resultfile}")
                    copyfile(f"{resultfile}", f"{filebase}")
                    shell("echo {resultfile} >> {output}")

else:
    rule share_to_resultdir:
        input:
            expand("{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], stype=sampleconfig[normalname]["stype"], caller="dnascope", sname=normalid, vcftype="germline"),
            expand("qc_report/{normalname}_qc_stats.xlsx", normalname=normalname),
            expand("{stype}/canvas/{sname}_CNV_germline.vcf.xlsx", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/dnascope/{sname}_{hgX}_SNV_CNV_germline.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], sname=normalid, stype=sampleconfig[normalname]["stype"], hgX=reference),
            expand("{stype}/canvas/{sname}_CNV_germline.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], sname=normalid, vartype="germline", stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/canvas/{sname}_germline_CNV_observed.seg", vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/canvas/{sname}_germline_CNV_called.seg", vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/realign/{sname}_REALIGNED.{fmt}", fmt=["bam", "bam.bai", "cram", "cram.crai"], sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/reports/{sname}_REALIGNED.bam.tdf", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/manta/{sname}_germline_MantaBNDs.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], sname=normalid, stype=sampleconfig[normalname]["stype"]),
            #expand("{stype}/manta/{sname}_germline_MantaBNDs.vcf.gz.csi", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/manta/{sname}_germline_MantaNOBNDs.{fmt}", fmt=["vcf.gz", "vcf.gz.csi"], sname=normalid, stype=sampleconfig[normalname]["stype"]),
            #expand("{stype}/manta/{sname}_germline_MantaNOBNDs.vcf.csi", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/reports/{sname}_baf.igv", sname=normalid, stype=sampleconfig[normalname]["stype"])
        output:
            "reporting/shared_result_files.txt"
        run:
            for resultfile in input:
                filebase = os.path.basename(f"{resultfile}")
                copyfile(f"{resultfile}", f"{filebase}")
                shell("echo {resultfile} >> {output}")
