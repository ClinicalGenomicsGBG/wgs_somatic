#!/apps/bio/software/anaconda2/envs/mathias_general/bin/python3.6
import argparse
import xlsxwriter
import os
import re
import json

from prepare_dicts.extract_dicts import extract_variantlist
from prepare_dicts.extract_dicts import prepare_variantdict
from prepare_dicts.extract_dicts import create_refseq_dict

from prepare_dicts.annotate_dicts import find_overlapping_gene
from prepare_dicts.annotate_dicts import find_nearby_gene
from prepare_dicts.annotate_dicts import find_crossing_genes

from prepare_dicts.prepare_excel_dict import prepare_excel_variantdict

from write_excelfile.write_to_excel import write_to_excel


def annotate_vcf(vcf, refseq, output):

    # Prepare OutputNames
    vcfname = os.path.basename(vcf)
    if output.endswith("/"):
        output = output[:-1]
    # Prepare Dict of Variants
    variantlist, vcf_header, canvas_info = extract_variantlist(vcf)
    variant_dict_list, unique_info_columns = prepare_variantdict(
        variantlist, vcf_header
    )

    # Prepare Header for ExcelFile
    variant_write_table_header = [
        "Varianttype",
        "Breakpoint 1",
        "GeneInfo 1",
        "Breakpoint 2",
        "GeneInfo 2",
    ]
    # If no variants, write empty header and return
    if not variant_dict_list:
        print("No variants found in the VCF file.")
        write_to_excel(output, vcfname, variant_write_table_header, [], canvas_info)
        return

    for columnname in variant_dict_list[0]:
        if columnname == "INFO":
            for info_columnname in unique_info_columns:
                variant_write_table_header.append(info_columnname)
        else:
            variant_write_table_header.append(columnname)

    # Prepare Dict of RefseqTranscripts
    refseq_dict = create_refseq_dict(refseq)
    # with open("refseqdict.json", "w") as jsondump:
    #    json.dump(refseq_dict, jsondump, ensure_ascii=False, indent=4)
    # Chrom List
    chrom_list = []
    for chrom in refseq_dict:
        chrom_list.append(chrom)

    #    count = 0
    #    for variant in variant_dict_list:
    #        count += 1
    #        variant_gene_dict = find_crossing_genes(refseq_dict, variant)
    #
    #        with open(f"variant_dicts_examples/variant_{count}.json", "w") as jsondump:
    #            json.dump(variant_gene_dict, jsondump, ensure_ascii=False, indent=4)
    #        if count > 100:
    #            break

    variant_dict_list_new = []
    for variant in variant_dict_list:
        variant_gene_dict = find_crossing_genes(refseq_dict, variant)
        variant_dict_list_new.append(variant_gene_dict)
    variant_dict_list = variant_dict_list_new
    variant_write_table_header.append("DEL/DUP Genecrossings")

    # Loop Through each Variant in dict and find overlapping / nearby genes
    variant_write_table = []
    for variant in variant_dict_list:
        variant_write = []
        variant_chrom = variant["#CHROM"]
        variant_pos = int(variant["POS"])
        variant_type = variant["ID"].split(":")[0]

        if variant_type == "MantaBND":
            # Handle breakends
            regex = re.compile(r"[\[\]]([\w\.]*:[0-9]*)[\[\]]")
            end_location = regex.findall(variant["ALT"])[0]
            end_chrom = end_location.split(":")[0]
            end_pos = int(end_location.split(":")[1])
        else:
            end_chrom = variant_chrom
            # Handle missing or invalid END field
            end_value = variant["INFO"].get("END", None)
            try:
                end_pos = int(end_value)
            except (ValueError, TypeError):
                # If END is missing or invalid, calculate it for deletions using SVLEN
                if variant_type == "MantaDEL" and "SVLEN" in variant["INFO"]:
                    svlen = int(variant["INFO"]["SVLEN"])
                    end_pos = variant_pos + abs(svlen)
                elif variant_type == "MantaDUP" and "SVLEN" in variant["INFO"]:
                    svlen = int(variant["INFO"]["SVLEN"])
                    end_pos = variant_pos + svlen
                elif variant_type == "MantaINS":
                    end_pos = variant_pos  # For insertions, END is the same as POS
                else:
                    print(
                        f"Warning: Missing or invalid END for variant {variant['ID']} at {variant_chrom}:{variant_pos}"
                    )
                    end_pos = variant_pos  # Default to POS if END is missing or invalid

        if variant_chrom not in chrom_list:
            pos_gene_info = ["Not available for chromosome"]
        else:
            pos_gene_info = find_overlapping_gene(
                refseq_dict, variant_chrom, variant_pos
            )
            if not pos_gene_info:
                pos_gene_info = find_nearby_gene(
                    refseq_dict, variant_chrom, variant_pos
                )

        if not variant_type == "MantaINS":
            if end_chrom in chrom_list:
                endpos_gene_info = find_overlapping_gene(
                    refseq_dict, end_chrom, end_pos
                )
                if not endpos_gene_info:
                    endpos_gene_info = find_nearby_gene(refseq_dict, end_chrom, end_pos)
            else:
                endpos_gene_info = ["Not available for chromosome"]
        else:
            endpos_gene_info = ["Not available for INS"]

        variant_write.extend(
            [
                variant_type,
                f"{variant_chrom}:{variant_pos}",
                pos_gene_info,
                f"{end_chrom}:{end_pos}",
                endpos_gene_info,
            ]
        )
        for stat in variant:
            if stat == "INFO":
                for unique_info_name in unique_info_columns:
                    variant_write.extend([variant[stat][unique_info_name]])
            else:
                variant_write.extend([variant[stat]])
        variant_write_table.append(variant_write)
    #######################################################
    # Write ExcelFile

    variant_dict = prepare_excel_variantdict(variant_write_table)

    #    with open("variantdict.json", "w") as jsondump:
    #        json.dump(variant_dict, jsondump, ensure_ascii=False, indent=4)
    write_to_excel(
        output, vcfname, variant_write_table_header, variant_dict, canvas_info
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v", "--vcf", nargs="?", help="Input MantaVCF to Annotate", required=True
    )
    parser.add_argument(
        "-g", "--refseqgtf", nargs="?", help="Input Refseq GTF file", required=True
    )
    parser.add_argument(
        "-o", "--output", nargs="?", help="location to output results", required=True
    )
    args = parser.parse_args()
    annotate_vcf(args.vcf, args.refseqgtf, args.output)
