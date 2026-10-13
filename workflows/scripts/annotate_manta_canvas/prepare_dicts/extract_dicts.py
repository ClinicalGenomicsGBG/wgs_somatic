def extract_variantlist(vcf):
    canvas_infofields = [
        "##OverallPloidy",
        "##DiploidCoverage",
        "##EstimatedTumorPurity",
        "##PurityModelFit",
        "##InterModelDistance",
        "##LocalSDmetric",
        "##EvennessScore",
        "##HeterogeneityProportion",
        "##EstimatedChromosomeCount",
    ]
    canvas_info = []
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
                if variant_info[0].split("=")[0] in canvas_infofields:
                    canvas_info.append(variant_info[0])

    return variantlist, vcf_header, canvas_info


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


def create_refseq_dict(refseqgtf):
    refseqgene_dict = {}
    with open(refseqgtf, "r") as refseq:
        header = refseq.readline().split("\t")
        keep_cols = ["name", "chrom", "strand", "exonStarts", "exonEnds", "name2"]
        keep_cols_dict = {}
        for keepcol in keep_cols:
            for colnumber, column in enumerate(header):
                if column == keepcol:
                    keep_cols_dict.update({keepcol: colnumber})

        for transcript in refseq:
            transcript_dict = {}
            transcript_info_list = transcript.split("\t")
            chrom = transcript_info_list[keep_cols_dict["chrom"]]
            name = transcript_info_list[keep_cols_dict["name"]]
            strand = transcript_info_list[keep_cols_dict["strand"]]
            name2 = transcript_info_list[keep_cols_dict["name2"]]
            if chrom not in refseqgene_dict:
                refseqgene_dict[chrom] = {}
            if name2 not in refseqgene_dict[chrom]:
                refseqgene_dict[chrom][name2] = {}
            refseqgene_dict[chrom][name2][name] = {}
            refseqgene_dict[chrom][name2][name]["strand"] = strand
            refseqgene_dict[chrom][name2][name]["exon"] = {}

            exon_start_list = transcript_info_list[keep_cols_dict["exonStarts"]].split(
                ","
            )[:-1]
            exon_stop_list = transcript_info_list[keep_cols_dict["exonEnds"]].split(
                ","
            )[:-1]

            refseqgene_dict[chrom][name2][name]["smallest"] = exon_start_list[0]
            refseqgene_dict[chrom][name2][name]["largest"] = exon_stop_list[-1]

            for exonnum, exon_start in enumerate(exon_start_list):
                real_exon_num = exonnum + 1
                if strand == "-":
                    real_exon_num = len(exon_start_list) - exonnum
                refseqgene_dict[chrom][name2][name]["exon"][real_exon_num] = {}
                exon_stop = exon_stop_list[exonnum]
                refseqgene_dict[chrom][name2][name]["exon"][real_exon_num]["start"] = (
                    exon_start
                )
                refseqgene_dict[chrom][name2][name]["exon"][real_exon_num]["stop"] = (
                    exon_stop
                )
        return refseqgene_dict
