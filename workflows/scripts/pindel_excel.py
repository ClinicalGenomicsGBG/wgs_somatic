#!/bin/python3.6
#import sys
import argparse
import xlsxwriter

from pysam import VariantFile

def position_gene(position_start, position_stop, bed):
    '''Function to get gene name based on start/stop position'''
    gene_start = None
    gene_stop = None
    with open(bed) as bed:
        for line in bed:
            if position_start >= int(line.split(" ")[1]) and position_start <= int(line.split(" ")[2]):
                gene_start = line.split(" ")[3]
            if position_stop >= int(line.split(" ")[1]) and position_stop <= int(line.split(" ")[2]):
                gene_stop = line.split(" ")[3]
    if gene_stop == None:
        gene = gene_start
    elif gene_start == None:
        gene = gene_stop
    elif gene_start == gene_stop:
        gene = gene_start
    elif gene_start != gene_stop:
        gene = "two different genes = weird"
    return(gene)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("vcf_input", nargs="?", type=str, help=':Full path to your vcf input file')
    parser.add_argument("xlsx_output", nargs="?", type=str, help=':Specify full output path and filename of the vcf file you want to create')
    parser.add_argument("bedfile", nargs="?", type=str, help=':Full path to bedfile used by pindel')

    args = parser.parse_args()

    # Create excel file
    workbook = xlsxwriter.Workbook(args.xlsx_output)

    # Define formats to be used.
    headingFormat = workbook.add_format({'bold': True, 'font_size': 18})
    tableHeadFormat = workbook.add_format({'bold': True, 'text_wrap': True})

    # Create worksheet
    worksheet = workbook.add_worksheet('Pindel')
    worksheet.write('A1', 'Pindel results', headingFormat)

    vcf_input = args.vcf_input
    vcf_input = VariantFile(vcf_input, 'r')

    # Get reference used
    for x in vcf_input.header.records:
        if x.key=='reference':
            ref = x.value

    worksheet.write('A4', 'Reference used: '+str(ref))


    # Get genes in bedfile
    with open(args.bedfile) as bed:
        genesDup = [line.split(" ")[3] for line in bed]
        genes = set(genesDup)


    # Write genes included in excel file
    worksheet.write('A6', 'Genes included: ')
    row = 7
    for gene in genes:
        worksheet.write('A'+str(row), gene)
        row += 1

    row += 1
    samples = list(vcf_input.header.samples)
    # Write info from vcf to excel file
    # If only tumor sample in vcf
    if len(samples) == 1:
        sample = samples[0]
        tableheading = ['Chr', 'Gene', 'Start', 'Stop', 'SV length', 'Ref', 'Alt', sample+'\nDP', sample+'\nAF', 'Supporting reads']
        worksheet.set_column('H:I', 10)
        worksheet.write_row('A'+str(row), tableheading, tableHeadFormat)
        for indel in vcf_input.fetch():
            svlen = indel.info["SVLEN"]
            if len(indel.alts) == 1:
                alt=indel.alts[0]
            s_dp = indel.samples[sample].get("DP")
            s_af = indel.samples[sample].get("AF")
            s_reads = indel.samples[sample].get("AD")[1]
            worksheet.write_row(row,0,[indel.contig, position_gene(indel.pos, indel.stop, args.bedfile),indel.pos, indel.stop, svlen, indel.ref, alt, s_dp, s_af, s_reads])
            row +=1
    # If tumor and normal sample in vcf
    if len(samples) == 2:
        sample1 = samples[0]
        sample2 = samples[1]
        tableheading = ['Chr', 'Gene', 'Start', 'Stop', 'SV length', 'Ref', 'Alt', sample1+'\nDP', sample1+'\nAF', sample2+'\nDP', sample2+'\nAF', sample1+'\nSupporting reads', sample2+'\nSupporting reads']
        worksheet.set_column('H:K', 10)
        worksheet.write_row('A'+str(row), tableheading, tableHeadFormat)
        for indel in vcf_input.fetch():
            svlen = indel.info["SVLEN"]
            if len(indel.alts) == 1:
                alt=indel.alts[0]
            s1_dp = indel.samples[sample1].get("DP")
            s1_af = indel.samples[sample1].get("AF")
            s2_dp = indel.samples[sample2].get("DP")
            s2_af = indel.samples[sample2].get("AF")
            s1_reads = indel.samples[sample1].get("AD")[1]
            s2_reads = indel.samples[sample2].get("AD")[1]

            # Only write variant to excel if it exists in tumor sample and in less than 3 reads in normal sample
            if s1_reads !=0 and s2_reads <=3:
                worksheet.write_row(row,0,[indel.contig, position_gene(indel.pos, indel.stop, args.bedfile),indel.pos, indel.stop, svlen, indel.ref, alt, s1_dp, s1_af, s2_dp, s2_af, s1_reads, s2_reads])
                row += 1


    workbook.close()

if __name__ == "__main__":
    main()
