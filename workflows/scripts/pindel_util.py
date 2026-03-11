import xlsxwriter

from pysam import VariantFile


def _position_gene(chrom, position_start, position_stop, bed):
    """Get gene name based on start/stop position."""

    # Normalize incoming positions.
    position_start = int(position_start) if position_start is not None else None
    position_stop = int(position_stop) if position_stop is not None else None

    with open(bed) as bed_handle:
        for line in bed_handle:
            chr_bed, start, end, gene = line.split()[:4]
            start = int(start)
            end = int(end)

            if chr_bed != chrom:
                continue
            if position_start is not None and start <= position_start <= end:
                return gene
            if position_stop is not None and start <= position_stop <= end:
                return gene

    return None

def pindel_write_xlsx(vcf_input, xlsx_output, bedfile, tumor=None, normal=None):
    """ 
    Write pindel results to excel file. 

    If tumor and normal sample in vcf, only write variants that 
    - exist in tumor sample
    - less than 3 reads in normal sample. 
    """
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
        if tumor:
            sample_tumor = tumor
        else:
            print(f"Warning: Tumor sample not specified but only one sample in vcf file. Using {samples[0]} as tumor sample.")
            sample_tumor = samples[0]

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
                    _position_gene(indel.contig, indel.pos, indel.stop, bedfile),
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
            st_dp = indel.samples[sample_tumor].get("DP")
            st_af = indel.samples[sample_tumor].get("AF")
            sn_dp = indel.samples[sample_normal].get("DP")
            sn_af = indel.samples[sample_normal].get("AF")
            st_reads = indel.samples[sample_tumor].get("AD")[1]
            sn_reads = indel.samples[sample_normal].get("AD")[1]

            # Only write variant to excel if it exists in tumor sample and in less than 3 reads in normal sample
            if st_reads != 0 and sn_reads <= 3:
                worksheet.write_row(
                    row,
                    0,
                    [
                        indel.contig,
                        _position_gene(indel.contig, indel.pos, indel.stop, bedfile),
                        indel.pos,
                        indel.stop,
                        svlen,
                        indel.ref,
                        alt,
                        st_dp,
                        st_af,
                        sn_dp,
                        sn_af,
                        st_reads,
                        sn_reads,
                    ],
                )
                row += 1

    workbook.close()


def pindel_fix_dpaf(vcf_input, vcf_output):
    """
    Add/repair per-sample DP and AF in a Pindel VCF.

    DP = sum(AD)
    AF = alt AD / sum(AD)

    This is applied to every sample in each record, regardless of whether
    the sample is tumor, normal, or something else.
    """
    vcf_in = VariantFile(vcf_input, "r")
    new_header = vcf_in.header

    # Add FORMAT fields if they do not already exist
    if "DP" not in new_header.formats:
        new_header.formats.add("DP", "1", "Integer", "Sum of AD fields")
    if "AF" not in new_header.formats:
        new_header.formats.add("AF", "1", "Float", "Alt AD / DP")

    vcf_out = VariantFile(vcf_output, "w", header=new_header)

    for record in vcf_in.fetch():
        for sample_name, sample_data in record.samples.items():
            ad = sample_data.get("AD")

            # Handle missing or malformed AD safely
            if not ad:
                print(f"Warning: {sample_name} has no AD at {record.contig}:{record.pos}")
                dp = 0
                af = 0.0
            else:
                dp = sum(reads for reads in ad if reads is not None)

                # Alt depth is the second element of AD if it exists and is not None
                if len(ad) > 1 and ad[1] is not None:
                    alt_depth = ad[1]
                else:
                    print(f"Warning: {sample_name} has no ALT depth in AD at {record.contig}:{record.pos}")
                    alt_depth = 0
                af = 0.0 if dp == 0 else alt_depth / dp

            sample_data["DP"] = dp
            sample_data["AF"] = af

        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()
