# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
from workflows.scripts.vcf_conversion import vcf_conversion

rule alissa_vcf:
    input:
        "{workingdir}/tumor/tnscope/{sname}_somatic_refseq3kfilt.vcf"
    output:
        "{workingdir}/{sname}_somatic_refseq3kfilt_Alissa.vcf"
    run:
        shell("cp {input} {output}")
        vcf_conversion(output)
