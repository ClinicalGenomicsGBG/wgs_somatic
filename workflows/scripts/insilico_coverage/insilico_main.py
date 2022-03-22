#!/usr/bin/env python
import argparse
import subprocess
import os

from insilico_annotate_mpileup import annotate_coverage
from insilico_per_region_stats import coverage_stats
from excel_summary import write_excel
from csv_to_excel import csv_to_excel
from anybelow_sample import alt_main


def run_insilico_main(bam, bed, version, outputdir, level, sample, samtools='samtools'):
    # handle outputdir input
    if outputdir.endswith("/"):
        outputdir = outputdir[:-1]
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    # bambase = os.path.basename(bam)
    bedbase = os.path.basename(bed)
    bedbase = bedbase.replace(".bed", "")
    bedbase = f'{bedbase}_v{version}'

    # run samtools (should be installed through snakemake envs)
    coverage_output = f"{outputdir}/{sample}_{bedbase}_cov.tsv"
    samtools_args = [samtools, "mpileup", "-a", "-l", bed, "-Q", "1", bam, "|", "cut", "-f1,2,4", ">", coverage_output]
    samtools_args = " ".join(samtools_args)
    subprocess.call(samtools_args, shell=True)

    # annotate mpileup
    annotated_coverage = f"{coverage_output}_annotated.tsv"
    if os.path.isfile(annotated_coverage):
        os.remove(annotated_coverage)
    annotate_coverage(coverage_output, bed, annotated_coverage)

    # insilico_per_region_stats
    per_regionstats_file = f"{outputdir}/{sample}_{bedbase}.csv"
    if os.path.isfile(per_regionstats_file):
        os.remove(per_regionstats_file)
    coverage_stats(level, annotated_coverage, "10,20", per_regionstats_file)

    # write_excel
    excel_out = per_regionstats_file.replace(".csv", "_genes_below10x.xlsx")
    write_excel(per_regionstats_file, excel_out)

    # any below regions
    below10x = f"{outputdir}/{sample}_{bedbase}_10x.csv"
    below20x = f"{outputdir}/{sample}_{bedbase}_20x.csv"
    alt_main(annotated_coverage, below10x, below20x)

    # convert to excel
    below10x_xlsx = below10x.replace(".csv", ".xlsx")
    below20x_xlsx = below20x.replace(".csv", ".xlsx")
    csv_to_excel(below10x, below10x_xlsx)
    csv_to_excel(below20x, below20x_xlsx)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-bam', '--bamfile', nargs='?', help='bamfile input to perform coverage analysis on', required=True)
    parser.add_argument('-bed', '--bedfile', nargs='?', help='bedfile with regions to calculate coverage within', required=True)
    parser.add_argument('-v', '--version', nargs='?', help='Version of supplied bedfile', required=True)
    parser.add_argument('-output', '--outputdir', nargs='?', help='where to put output', required=True)
    parser.add_argument('-level', '--annotationlevel', nargs='?', help='annotationlevel of bedfile, example: SDHB_NM_003000.2_exon_1 = 2, ISCA2 = 1 (gene only)' , required=True)
    parser.add_argument('-sample', '--samplename', nargs='?', help='name of sample to be analysed (provide full sampleid)', required=True)
    parser.add_argument('-sam', '--samtools', nargs='?', help='if samtools is not set in environment, supply here', required=False)
    args = parser.parse_args()
    run_insilico_main(args.bamfile, args.bedfile, args.version, args.outputdir, args.annotationlevel, args.samplename, args.samtools)
