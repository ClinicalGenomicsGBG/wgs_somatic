import re


def prepare_excel_variantdict(variantinfo):
    variant_dict = {}
    for variantnum, variant in enumerate(variantinfo):
        variant_dict[variantnum] = {}
        variant_dict[variantnum]["row1"] = {}
        for columnpos, info in enumerate(variant):
            if columnpos == 2 or columnpos == 4:
                continue
            else:
                variant_dict[variantnum]["row1"][columnpos] = info
        variant_dict[variantnum]["break1"] = {}
        variant_dict[variantnum]["break2"] = {}
        variant_dict[variantnum]["countbreak1"] = len(variant[2])
        variant_dict[variantnum]["countbreak2"] = len(variant[4])
        # for annotation in breakpoint 1
        for transcriptnum, geneinfo in enumerate(variant[2]):
            variant_dict[variantnum]["break1"][geneinfo] = {}
            try:
                gene_info_break2 = variant[4][transcriptnum]
            except:
                gene_info_break2 = variant[4][0]
            variant_dict[variantnum]["break1"][geneinfo]["match"] = gene_info_break2
            gene_info_break1 = re.sub("(\([^\)]*\))", "", geneinfo)
            gene_info_break2 = re.sub("(\([^\)]*\))", "", gene_info_break2)
            if "Not available for chromosome" in gene_info_break2:
                if ":exon" in geneinfo and "-exon" not in geneinfo:
                    variant_dict[variantnum]["break1"][geneinfo]["status"] = "exonic"
                elif "Not available for chromosome" in geneinfo:
                    variant_dict[variantnum]["break1"][geneinfo]["status"] = "strange"
                else:
                    variant_dict[variantnum]["break1"][geneinfo]["status"] = (
                        "likelyartifact"
                    )
            elif "Not available for chromosome" in geneinfo:
                variant_dict[variantnum]["break1"][geneinfo]["status"] = "strange"
            elif gene_info_break1 == gene_info_break2:
                if ":exon" in geneinfo and "-exon" not in geneinfo:
                    variant_dict[variantnum]["break1"][geneinfo]["status"] = "exonic"
                else:
                    variant_dict[variantnum]["break1"][geneinfo]["status"] = "intronic"
            elif "MantaINS" in variant[0]:
                if ":exon" in geneinfo and "-exon" not in geneinfo:
                    variant_dict[variantnum]["break1"][geneinfo]["status"] = "exonic"
                else:
                    variant_dict[variantnum]["break1"][geneinfo]["status"] = "intronic"
            else:
                if ":exon" in geneinfo:
                    if "-exon" not in geneinfo:
                        variant_dict[variantnum]["break1"][geneinfo]["status"] = (
                            "exonic"
                        )
                    else:
                        variant_dict[variantnum]["break1"][geneinfo]["status"] = (
                            "crossintronic"
                        )
                else:
                    variant_dict[variantnum]["break1"][geneinfo]["status"] = (
                        "crossgenic"
                    )

        # for annotation in breakpoint 2
        for transcriptnum, geneinfo in enumerate(variant[4]):
            variant_dict[variantnum]["break2"][geneinfo] = {}
            try:
                gene_info_break1 = variant[2][transcriptnum]
            except:
                gene_info_break1 = variant[2][0]
            variant_dict[variantnum]["break2"][geneinfo]["match"] = gene_info_break1
            gene_info_break2 = re.sub("(\([^\)]*\))", "", geneinfo)
            gene_info_break1 = re.sub("(\([^\)]*\))", "", gene_info_break1)
            if "Not available for chromosome" in gene_info_break1:
                if ":exon" in geneinfo and "-exon" not in geneinfo:
                    variant_dict[variantnum]["break2"][geneinfo]["status"] = "exonic"
                elif "Not available for chromosome" in geneinfo:
                    variant_dict[variantnum]["break2"][geneinfo]["status"] = "strange"
                else:
                    variant_dict[variantnum]["break2"][geneinfo]["status"] = (
                        "likelyartifact"
                    )
            elif "Not available for chromosome" in geneinfo:
                variant_dict[variantnum]["break2"][geneinfo]["status"] = "strange"
            elif gene_info_break1 == gene_info_break2:
                if ":exon" in geneinfo and "-exon" not in geneinfo:
                    variant_dict[variantnum]["break2"][geneinfo]["status"] = "exonic"
                else:
                    variant_dict[variantnum]["break2"][geneinfo]["status"] = "intronic"
            elif "MantaINS" in variant[0]:
                if ":exon" in geneinfo and "-exon" not in geneinfo:
                    variant_dict[variantnum]["break2"][geneinfo]["status"] = "exonic"
                else:
                    variant_dict[variantnum]["break2"][geneinfo]["status"] = "intronic"
            else:
                if ":exon" in geneinfo:
                    if "-exon" not in geneinfo:
                        variant_dict[variantnum]["break2"][geneinfo]["status"] = (
                            "exonic"
                        )
                    else:
                        variant_dict[variantnum]["break2"][geneinfo]["status"] = (
                            "crossintronic"
                        )
                else:
                    variant_dict[variantnum]["break2"][geneinfo]["status"] = (
                        "crossgenic"
                    )
    return variant_dict
