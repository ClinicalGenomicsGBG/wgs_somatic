import pysam
import argparse


def is_reciprocal_pair(bnd_pairs, key, mate_key):
    """
    Check if the given key:value pair already exists as a value:key in the dictionary.
    """
    return mate_key in bnd_pairs and bnd_pairs[mate_key] == key


def filter_vcf(input_vcf, output_vcf, tumor_name=None, normal_name=None, min_tumor_support=3, max_normal_support=2):
    """
    Filters a VCF file based on the following criteria:
    - FILTER column must be "PASS".
    - Both PR and SR fields must be present.
    - At least `min_tumor_support` reads supporting the variant in the tumor.
    - If a normal sample is provided:
        - No more than `max_normal_support` reads supporting the variant in the normal.
    - Removes duplicate breakends (reciprocal variants).
    """
    # Dictionary to track reciprocal breakends
    seen_bnd_ids = set()

    with pysam.VariantFile(input_vcf) as vcf_in, pysam.VariantFile(output_vcf, "w", header=vcf_in.header) as vcf_out:
        for record in vcf_in:
            # Check if the variant passes the FILTER column
            if "PASS" not in record.filter.keys():
                continue

            # Extract tumor sample data
            if tumor_name:
                try:
                    tumor_sample = record.samples[tumor_name]
                except KeyError:
                    continue  # Skip if tumor sample is missing

                # Ensure both PR and SR are present in the tumor sample
                if "PR" not in tumor_sample or "SR" not in tumor_sample:
                    continue

                # Extract PR and SR values for the tumor sample
                tumor_pr, tumor_sr = tumor_sample["PR"][1], tumor_sample["SR"][1]

                # Require both PR and SR supporting the variant
                if tumor_pr == 0 or tumor_sr == 0:
                    continue
                tumor_support = tumor_pr + tumor_sr

                # Apply the tumor support filter
                if tumor_support < min_tumor_support:
                    continue

            # If a normal sample is provided, apply normal-specific filters
            if normal_name:
                try:
                    normal_sample = record.samples[normal_name]
                except KeyError:
                    continue  # Skip if normal sample is missing

                # Combine PR and SR support from the normal sample
                normal_support = 0
                if "PR" in normal_sample:
                    normal_support += normal_sample["PR"][1]
                if "SR" in normal_sample:
                    normal_support += normal_sample["SR"][1]

                if tumor_name:
                    # Apply the normal support filter
                    if normal_support > max_normal_support:
                        continue

            # Handle reciprocal breakends (BNDs)
            if record.info.get("SVTYPE") == "BND" and "MATEID" in record.info:
                id1 = record.id
                id2 = record.info["MATEID"]
                if isinstance(id2, (tuple, list)):
                    id2 = id2[0]
                pair = tuple(sorted([id1, id2]))
                if pair in seen_bnd_ids:
                    continue  # Skip duplicate
                seen_bnd_ids.add(pair)

            # Write the record to the output VCF
            vcf_out.write(record)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter a VCF file based on tumor and normal sample criteria.")
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument("output_vcf", help="Output VCF file")
    parser.add_argument("--tumor_name", help="Tumor sample name")
    parser.add_argument("--normal_name", help="Normal sample name (optional)", default=None)
    parser.add_argument("--min_tumor_support", type=int, help="Minimum tumor support (default: 3)", default=3)
    parser.add_argument("--max_normal_support", type=int, help="Maximum normal support (default: 2)", default=2)

    args = parser.parse_args()

    filter_vcf(
        input_vcf=args.input_vcf,
        output_vcf=args.output_vcf,
        tumor_name=args.tumor_name,
        normal_name=args.normal_name,
        min_tumor_support=args.min_tumor_support,
        max_normal_support=args.max_normal_support,
    )
