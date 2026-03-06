#!/bin/python3.6
#import sys
import argparse
import xlsxwriter

from pysam import VariantFile


def position_gene(chr, position_start, position_stop, bed):
    """Function to get gene name based on start/stop position"""
    gene = None
    with open(bed) as bed:
        for line in bed:
            chr_bed, start, end, gene = line.split()[:4]
            if chr_bed != chr:
                continue
            if int(start) <= position_start <= int(end):
                return gene
            elif int(start) <= position_stop <= int(end):
                return gene
    return gene


def write_pindel_xlsx(vcf_input, xlsx_output, bedfile, tumor=None, normal=None):
    # Create excel file
    workbook = xlsxwriter.Workbook(xlsx_output)

    # Define formats to be used.
    headingFormat = workbook.add_format({"bold": True, "font_size": 18})
    tableHeadFormat = workbook.add_format({"bold": True, "text_wrap": True})

    # Create worksheet
    worksheet = workbook.add_worksheet("Pindel")
    worksheet.write("A1", "Pindel results", headingFormat)

    vcf_input = VariantFile(vcf_input, "r")

    # Get reference used
    ref = None
    for x in vcf_input.header.records:
        if x.key == "reference":
            ref = x.value

    worksheet.write("A4", "Reference used: " + str(ref))

    # Get genes in bedfile
    with open(bedfile) as bed:
        genesDup = [line.split()[3] for line in bed]
        genes = set(genesDup)

    # Write genes included in excel file
    worksheet.write("A6", "Genes included: ")
    row = 7
    for gene in genes:
        worksheet.write("A" + str(row), gene)
        row += 1

    row += 1
    samples = list(vcf_input.header.samples)
    # Write info from vcf to excel file
    # If only tumor sample in vcf
    if len(samples) == 1:
        sample_tumor = tumor if tumor else samples[0]
        if sample_tumor not in samples:
            raise ValueError(f"Tumor sample {sample_tumor} is not in vcf file")
        if normal:
            raise ValueError("Normal sample specified but only one sample in vcf file")
        tableheading = [
            "Chr",
            "Gene",
            "Start",
            "Stop",
            "SV length",
            "Ref",
            "Alt",
            sample_tumor + "\nDP",
            sample_tumor + "\nAF",
            "Supporting reads",
        ]
        worksheet.set_column("H:I", 10)
        worksheet.write_row("A" + str(row), tableheading, tableHeadFormat)
        for indel in vcf_input.fetch():
            svlen = indel.info["SVLEN"]
            alt = indel.alts[0] if indel.alts and len(indel.alts) == 1 else None
            s_dp = indel.samples[sample_tumor].get("DP")
            s_af = indel.samples[sample_tumor].get("AF")
            s_reads = indel.samples[sample_tumor].get("AD")[1]
            worksheet.write_row(
                row,
                0,
                [
                    indel.contig,
                    position_gene(indel.contig, indel.pos, indel.stop, bedfile),
                    indel.pos,
                    indel.stop,
                    svlen,
                    indel.ref,
                    alt,
                    s_dp,
                    s_af,
                    s_reads,
                ],
            )
            row += 1
    # If tumor and normal sample in vcf
    if len(samples) == 2:
        if not tumor:
            raise ValueError("Tumor sample not specified but two samples in vcf file")
        sample_tumor = samples[samples.index(tumor)]
        if normal:
            sample_normal = samples[samples.index(normal)]
        else:
            sample_normal = samples[0] if sample_tumor == samples[1] else samples[1]

        tableheading = [
            "Chr",
            "Gene",
            "Start",
            "Stop",
            "SV length",
            "Ref",
            "Alt",
            sample_tumor + "\nDP",
            sample_tumor + "\nAF",
            sample_normal + "\nDP",
            sample_normal + "\nAF",
            sample_tumor + "\nSupporting reads",
            sample_normal + "\nSupporting reads",
        ]
        worksheet.set_column("H:K", 10)
        worksheet.write_row("A" + str(row), tableheading, tableHeadFormat)
        for indel in vcf_input.fetch():
            svlen = indel.info["SVLEN"]
            alt = indel.alts[0] if indel.alts and len(indel.alts) == 1 else None
            s1_dp = indel.samples[sample_tumor].get("DP")
            s1_af = indel.samples[sample_tumor].get("AF")
            s2_dp = indel.samples[sample_normal].get("DP")
            s2_af = indel.samples[sample_normal].get("AF")
            s1_reads = indel.samples[sample_tumor].get("AD")[1]
            s2_reads = indel.samples[sample_normal].get("AD")[1]

            # Only write variant to excel if it exists in tumor sample and in less than 3 reads in normal sample
            if s1_reads != 0 and s2_reads <= 3:
                worksheet.write_row(
                    row,
                    0,
                    [
                        indel.contig,
                        position_gene(indel.contig, indel.pos, indel.stop, bedfile),
                        indel.pos,
                        indel.stop,
                        svlen,
                        indel.ref,
                        alt,
                        s1_dp,
                        s1_af,
                        s2_dp,
                        s2_af,
                        s1_reads,
                        s2_reads,
                    ],
                )
                row += 1

    workbook.close()


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--vcf_input", nargs="?", type=str, help=':Full path to your vcf input file')
    parser.add_argument("--xlsx_output", nargs="?", type=str, help=':Specify full output path and filename of the vcf file you want to create')
    parser.add_argument("--bedfile", nargs="?", type=str, help=':Full path to bedfile used by pindel')
    parser.add_argument("--tumor", nargs="?", type=str, help=':Name of tumor sample in vcf file')
    parser.add_argument("--normal", nargs="?", type=str, help=':Name of normal sample in vcf file', required=False)
    args = parser.parse_args()
    write_pindel_xlsx(args.vcf_input, args.xlsx_output, args.bedfile, args.tumor, args.normal)


if __name__ == "__main__":
    main()
