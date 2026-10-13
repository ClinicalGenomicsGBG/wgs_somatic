#!/apps/bio/software/anaconda2/envs/mathias_general/bin/python3.6
import argparse
import os
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def extract_variantlist(vcf):
    with open(vcf, "r") as vcffile:
        variantlist = []
        for variant in vcffile:
            variant = variant.rstrip("\n")
            variant_info = variant.split("\t")
            if not variant_info[0].startswith("#"):
                variantlist.append(variant_info)
            else:
                if variant_info[0].startswith("#CHROM"):
                    vcf_header = variant_info
    return variantlist, vcf_header


def prepare_variantdict(variantlist, vcf_header):
    variant_dict_list = []
    all_info_columns = []
    for variant in variantlist:
        variant_dict = {}
        for column_index, column_name in enumerate(vcf_header):
            variant_dict[column_name] = variant[column_index]
            # Collect all info-columnnames
            if column_name == "INFO":
                info_columns = [
                    info_column.split("=")[0]
                    for info_column in variant_dict["INFO"].split(";")
                ]
                all_info_columns.extend(info_columns)
        variant_dict_list.append(variant_dict)

    # Remove all duplicate info-columnnames
    unique_info_columns = list(set(all_info_columns))

    # Replace variant info string with a variant info dict
    final_variant_dict_list = []
    for variant_dict in variant_dict_list:
        variant_info_dict = {}
        variant_info_list = [
            info_column.split("=") for info_column in variant_dict["INFO"].split(";")
        ]
        for info_type in variant_info_list:
            if len(info_type) < 2:
                variant_info_dict[info_type[0]] = "yes"
            else:
                variant_info_dict[info_type[0]] = info_type[1]
        for info_column in unique_info_columns:
            if info_column in variant_info_dict:
                continue
            else:
                variant_info_dict[info_column] = "N/A"
        variant_dict["INFO"] = variant_info_dict
        final_variant_dict_list.append(variant_dict)
    return final_variant_dict_list, unique_info_columns


def plot_freq(vcf, output, dbsnp, hg38ref=False):

    # Prepare OutputNames
    vcfname = os.path.basename(vcf)
    if output.endswith("/"):
        output = output[:-1]
    # Prepare Dict of Variants
    variantlist, vcf_header = extract_variantlist(vcf)
    variant_dict_list, unique_info_columns = prepare_variantdict(
        variantlist, vcf_header
    )
    samplename = vcf_header[-1]
    with open(f"{output}", "w") as igvallelefile:
        igvallelefile.write("#type=GENE_EXPRESSION\n")
        igvallelefile.write(
            f'#track graphtype=points name="{samplename}" color=0,0,255 altColor=255,0,0 maxHeightPixels=80:80:80 viewLimits=-1:1\n'
        )
        igvallelefile.write("#Chromosome\tStart\tEnd\tFeatures\tvalues\n")
        if dbsnp:
            for variant in variant_dict_list:
                if variant["ID"] != ".":
                    samplecolumn = variant[samplename]
                    try:
                        allele_count = samplecolumn.split(":")[1].split(",")[1]
                        coverage = samplecolumn.split(":")[2]
                        if int(coverage) == 0:
                            logging.info(
                                f"Skipping variant with 0 depth: {samplecolumn}"
                            )
                            continue
                        fraction = round(float(allele_count) / float(coverage), 3)
                        if not hg38ref == "yes":
                            igvallelefile.write(
                                f"chr{variant['#CHROM']}\t{variant['POS']}\t{variant['POS']}\t{variant['ALT']}\t{fraction}\n"
                            )
                        else:
                            igvallelefile.write(
                                f"{variant['#CHROM']}\t{variant['POS']}\t{variant['POS']}\t{variant['ALT']}\t{fraction}\n"
                            )
                    except Exception as e:
                        logging.error(f"Error processing line: {samplecolumn}")
                        logging.error(f"Exception: {e}")
        else:
            for variant in variant_dict_list:
                samplecolumn = variant[samplename]
                try:
                    allele_count = samplecolumn.split(":")[1].split(",")[1]
                    coverage = samplecolumn.split(":")[2]
                    if int(coverage) == 0:
                        logging.info(f"Skipping variant with 0 depth: {samplecolumn}")
                        continue
                    fraction = round(float(allele_count) / float(coverage), 3)
                    if not hg38ref == "yes":
                        igvallelefile.write(
                            f"chr{variant['#CHROM']}\t{variant['POS']}\t{variant['POS']}\t{variant['ALT']}\t{fraction}\n"
                        )
                    else:
                        igvallelefile.write(
                            f"{variant['#CHROM']}\t{variant['POS']}\t{variant['POS']}\t{variant['ALT']}\t{fraction}\n"
                        )
                except Exception as e:
                    logging.error(f"Error processing line: {samplecolumn}")
                    logging.error(f"Exception: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--dbsnp",
        nargs="?",
        help="Restrict variants to those present in dbSNP (must be annotated with dbsnp)",
    )
    parser.add_argument(
        "-v", "--vcf", nargs="?", help="Input MantaVCF to Annotate", required=True
    )
    parser.add_argument(
        "-o",
        "--output",
        nargs="?",
        help="location to output results (directory)",
        required=True,
    )
    parser.add_argument(
        "-hg38",
        "--hg38ref",
        nargs="?",
        help="not adding extra chr if using hg38 reference",
        required=False,
    )
    args = parser.parse_args()
    plot_freq(args.vcf, args.output, args.dbsnp, args.hg38ref)

# type=GENE_EXPRESSION
# track graphtype=points name="Sample1.igv" color=0,0,255 altColor=255,0,0 maxHeightPixels=80:80:80 viewLimits=-1:1
# Chromosome  Start   End Features    values
# chr1    2488157 2488157        variant1 0.4818
# chr1    2491207 2491207        variant2 0.537
