# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
from shutil import copyfile

rule share_to_resultdir:
    input:
        get_result_input
    output:
        "reporting/shared_result_files.txt"
    run:
        for resultfile in input:
            filebase = os.path.basename(f"{resultfile}")
            try:
                os.link(f"{resultfile}", f"{filebase}")
            except:
                shutil.copyfile(f"{resultfile}", f"{filebase}")
            shell("echo {resultfile} >> {output}")
