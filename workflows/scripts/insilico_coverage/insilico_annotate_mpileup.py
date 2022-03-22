#!/usr/bin/env python
import os
import csv
import sys
import argparse

from collections import defaultdict


def prepare_insilico_list(bedfile):
    """Generate a list of annotation rows from the given bedfile."""
    with open(bedfile, "r") as bedfile:
        insilico_list = []
        for annotationrow in bedfile:
            annotation_dict = {}
            annotation_info = annotationrow.split("\t")[:-1]
            annotation_dict["chrom"] = annotation_info[0]
            annotation_dict["start"] = annotation_info[1]
            annotation_dict["stop"] = annotation_info[2]
            annotation_name = annotation_info[3].split("_")
            # Possibly annotation formats
            # GENE_NM_XXXXXXXX_Exon_X
            # GENE_Exon_X
            # NM_XXXXXXXX_Exon_X
            # NM_XXXXXXXX
            # GENE
            # Discover format
            annotation_len = len(annotation_name)
            if annotation_len == 1:
                # Only Genes
                annotation_dict["gene"] = annotation_name[0]
            elif annotation_len == 2:
                # Only Transcript
                transcript = '_'.join(annotation_name[0:1])
                annotation_dict["transcript"] = transcript
            elif annotation_len == 3:
                # Gene + Exon
                annotation_dict["gene"] = annotation_name[0]
                exon = '_'.join(annotation_name[1:2])
                annotation_dict["exon"] = exon
            elif annotation_len == 4:
                # Transcript + Exon
                transcript = '_'.join(annotation_name[0:1])
                exon = '_'.join(annotation_name[1:2])
                annotation_dict["transcript"] = transcript
                annotation_dict["exon"] = exon
            elif annotation_len == 5:
                # Gene + Transcript + Exon
                annotation_dict["gene"] = annotation_name[0]
                transcript = '_'.join(annotation_name[1:3])
                exon = '_'.join(annotation_name[3:5])
                annotation_dict["transcript"] = transcript
                annotation_dict["exon"] = exon
            else:
                raise Exception(f"unexpected annotation format: {annotation_name}")

            insilico_list.append(annotation_dict)
    return insilico_list


def prepare_insilico_dict(insilico_list):
    """Divide the list on chromosomes for quicker lookup."""
    insilico_dict = defaultdict(list)
    for element in insilico_list:
        insilico_dict[element['chrom']].append(element)

    return insilico_dict


def annotate_coverage(coverage, bedfile, output):
    """Annotate the coverage file with entries from the bedfile."""
    insilico_list = prepare_insilico_list(bedfile)
    insilico_dict = prepare_insilico_dict(insilico_list)
    if output.endswith('/'):
        output = output[:-1]

    with open(output, 'w') as csvout:
        csv_writer = csv.writer(csvout)
        with open(coverage, "r") as coverage:
            for covrow in coverage:
                cov_list = covrow.rstrip('\n').split("\t")
                chrom = cov_list[0]
                # Only look through entries on the specified chromosome
                for region in insilico_dict[chrom]:
                    if int(region['start']) <= int(cov_list[1]) <= int(region['stop']):
                        gene = region.get('gene', '')
                        transcript = region.get('transcript', '')
                        exon = region.get('exon', '')
                        # NOTE: This is the same logic as before. Can only 1 be missing?
                        cov_list += [x for x in [gene, transcript, exon] if x]
                csv_writer.writerow(cov_list)  # Write the entry directly to file


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--coverage', nargs='?', help='(from samtools mpileup -l $BED -Q 1 $BAM | cut -f1,2,4)', required=True)
    parser.add_argument('-b', '--bedfile', nargs='?', help='(annotate coverage-positions with Gene Transcript Exonnumber)', required=True)
    parser.add_argument('-o', '--output', nargs='?', help='(Output location)', required=True)
    args = parser.parse_args()
    coverage = args.coverage
    bedfile = args.bedfile
    output = args.output
    annotate_coverage(coverage, bedfile, output)

