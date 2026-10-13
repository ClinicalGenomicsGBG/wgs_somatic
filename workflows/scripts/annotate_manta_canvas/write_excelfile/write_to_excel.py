import xlsxwriter
import json


def divide_into_sheets(variantdict):
    variantdict_tier1 = {}
    variantdict_tier2 = {}
    tier2_status_list = ["strange", "intronic", "likelyartifact"]
    for variant in variantdict:
        break1_annos = variantdict[variant]["break1"]
        break2_annos = variantdict[variant]["break2"]
        for an_count, anno1 in enumerate(break1_annos):
            status = break1_annos[anno1]["status"]
            break2_match = break1_annos[anno1]["match"]
            break2_status = variantdict[variant]["break2"][break2_match]["status"]
            if break2_status != "exonic":
                break2_status = status
            if status in tier2_status_list:
                if variant not in variantdict_tier2:
                    variantdict_tier2[variant] = {}
                    variantdict_tier2[variant]["row1"] = variantdict[variant]["row1"]
                    variantdict_tier2[variant]["breakcolumns"] = {}
                variantdict_tier2[variant]["breakcolumns"][an_count] = {}
                variantdict_tier2[variant]["breakcolumns"][an_count]["break1"] = anno1
                variantdict_tier2[variant]["breakcolumns"][an_count][
                    "break1_status"
                ] = status
                variantdict_tier2[variant]["breakcolumns"][an_count]["break2"] = (
                    break2_match
                )
                variantdict_tier2[variant]["breakcolumns"][an_count][
                    "break2_status"
                ] = break2_status
            else:
                if variant not in variantdict_tier1:
                    variantdict_tier1[variant] = {}
                    variantdict_tier1[variant]["row1"] = variantdict[variant]["row1"]
                    variantdict_tier1[variant]["breakcolumns"] = {}
                variantdict_tier1[variant]["breakcolumns"][an_count] = {}
                variantdict_tier1[variant]["breakcolumns"][an_count]["break1"] = anno1
                variantdict_tier1[variant]["breakcolumns"][an_count][
                    "break1_status"
                ] = status
                variantdict_tier1[variant]["breakcolumns"][an_count]["break2"] = (
                    break2_match
                )
                variantdict_tier1[variant]["breakcolumns"][an_count][
                    "break2_status"
                ] = break2_status
    return variantdict_tier1, variantdict_tier2


def write_to_sheet(sheetname, variantdict, header, canvas_info, formatdict, excelfile):
    worksheet = excelfile.add_worksheet(sheetname)
    worksheet.set_column(2, 2, 70)
    worksheet.set_column(4, 4, 70)
    row = 0
    column = 0

    for column_name in header:
        worksheet.write(row, column, column_name, formatdict["header"])
        column += 1

    row += 1
    for variant in variantdict:
        for column in variantdict[variant]["row1"]:
            worksheet.write(row, column, variantdict[variant]["row1"][column])
        for break_annotation in variantdict[variant]["breakcolumns"]:
            break1 = variantdict[variant]["breakcolumns"][break_annotation]["break1"]
            break1_s = variantdict[variant]["breakcolumns"][break_annotation][
                "break1_status"
            ]
            break2 = variantdict[variant]["breakcolumns"][break_annotation]["break2"]
            break2_s = variantdict[variant]["breakcolumns"][break_annotation][
                "break2_status"
            ]
            worksheet.write(row, 2, break1, formatdict[break1_s])
            worksheet.write(row, 4, break2, formatdict[break2_s])
            row += 1
    return excelfile


def write_to_excel(output, vcfname, header, variantinfo, canvas_info):
    variantdict_tier1, variantdict_tier2 = divide_into_sheets(variantinfo)

    #    with open("variantdict_tier1.json", "w") as jsondump:
    #        json.dump(variantdict_tier1, jsondump, ensure_ascii=False, indent=4)
    #    with open("variantdict_tier2.json", "w") as jsondump:
    #        json.dump(variantdict_tier2, jsondump, ensure_ascii=False, indent=4)

    excelfile = xlsxwriter.Workbook(f"{output}/{vcfname}.xlsx")

    formatdict = {}
    formatdict["header"] = excelfile.add_format(
        {"bold": True, "font_color": "white", "bg_color": "black"}
    )
    formatdict["crossintronic"] = excelfile.add_format({"bg_color": "yellow"})
    formatdict["crossgenic"] = excelfile.add_format({"bg_color": "yellow"})
    formatdict["exonic"] = excelfile.add_format({"bg_color": "lime"})
    formatdict["intronic"] = excelfile.add_format({"bg_color": "silver"})
    formatdict["likelyartifact"] = excelfile.add_format({"bg_color": "gray"})
    formatdict["strange"] = excelfile.add_format({"bg_color": "5422AB"})

    excelfile = write_to_sheet(
        "tier1", variantdict_tier1, header, canvas_info, formatdict, excelfile
    )
    excelfile = write_to_sheet(
        "tier2", variantdict_tier2, header, canvas_info, formatdict, excelfile
    )

    excelfile.close()

    return

    excelfile = xlsxwriter.Workbook(f"{output}/{vcfname}.xlsx")
    worksheet = excelfile.add_worksheet("tier1")
    worksheet2 = excelfile.add_worksheet("tier2")
    worksheet.set_column(2, 2, 70)
    worksheet.set_column(4, 4, 70)
    worksheet2.set_column(2, 2, 70)
    worksheet2.set_column(4, 4, 70)
    row = 0
    column = 0
    row2 = 0
    column2 = 0

    header_format = excelfile.add_format(
        {"bold": True, "font_color": "white", "bg_color": "black"}
    )
    inside_gene_format = excelfile.add_format({"bg_color": "yellow"})
    inside_exon_format = excelfile.add_format({"bg_color": "lime"})
    intronic = excelfile.add_format({"bg_color": "silver"})
    likelyartifact = excelfile.add_format({"bg_color": "gray"})
    strange = excelfile.add_format({"bg_color": "5422AB"})

    for column_name in header:
        worksheet.write(row, column, column_name, header_format)
        worksheet.write(row2, column2, column_name, header_format)
        column += 1
        column2 += 1

    column = 0
    column2 = 0
    row += 1
    row2 += 1
    for variant in variantdict:
        annotations1 = variantdict[variant]["countbreak1"]
        annotations2 = variantdict[variant]["countbreak2"]
        writesheet = "Unknown"
        for annotation in variantdict[variant]["break1"]:
            if (
                "exonic" in variantdict[variant]["break1"][annotation]
                or "insidegene" in variantdict[variant]["break1"][annotation]
                or "crossgenic" in variantdict[variant]["break1"][annotation]
            ):
                writesheet = "tier1"
            else:
                writesheet = "tier2"
        if writesheet != "tier1":
            for annotation in variantdict[variant]["break2"]:
                if (
                    "exonic" in variantdict[variant]["break2"][annotation]
                    or "insidegene" in variantdict[variant]["break2"][annotation]
                    or "crossgenic" in variantdict[variant]["break2"][annotation]
                ):
                    writesheet = "tier1"
                else:
                    writesheet = "tier2"
        if writesheet == "tier1":
            for column in variantdict[variant]["row1"]:
                worksheet.write(row, column, variantdict[variant]["row1"][column])
            for annotation in variantdict[variant]["break1"]:
                if variantdict[variant]["break1"][annotation] == "exonic":
                    writestyle = inside_exon_format
                elif (
                    variantdict[variant]["break1"][annotation] == "crossgenic"
                    or variantdict[variant]["break1"][annotation] == "insidegene"
                ):
                    writestyle = inside_gene_format
                elif variantdict[variant]["break1"][annotation] == "intronic":
                    writestyle = intronic
                elif variantdict[variant]["break1"][annotation] == "likelyartifact":
                    writestyle = likelyartifact
                else:
                    writestyle = strange
                worksheet.write(row, 2, annotation, writestyle)
                row += 1
            row -= annotations1
            for annotation in variantdict[variant]["break2"]:
                if variantdict[variant]["break2"][annotation] == "exonic":
                    writestyle = inside_exon_format
                elif (
                    variantdict[variant]["break2"][annotation] == "crossgenic"
                    or variantdict[variant]["break2"][annotation] == "insidegene"
                ):
                    writestyle = inside_gene_format
                elif variantdict[variant]["break2"][annotation] == "intronic":
                    writestyle = intronic
                elif variantdict[variant]["break2"][annotation] == "likelyartifact":
                    writestyle = likelyartifact
                else:
                    writestyle = strange
                worksheet.write(row, 4, annotation, writestyle)
                row += 1
            row -= annotations2
            if annotations1 >= annotations2:
                row += annotations1
            else:
                row += annotations2
        else:
            for column in variantdict[variant]["row1"]:
                worksheet2.write(row, column, variantdict[variant]["row1"][column])

            for annotation in variantdict[variant]["break1"]:
                if variantdict[variant]["break1"][annotation] == "exonic":
                    writestyle = inside_exon_format
                elif (
                    variantdict[variant]["break1"][annotation] == "crossgenic"
                    or variantdict[variant]["break1"][annotation] == "insidegene"
                ):
                    writestyle = inside_gene_format
                elif variantdict[variant]["break1"][annotation] == "intronic":
                    writestyle = intronic
                elif variantdict[variant]["break1"][annotation] == "likelyartifact":
                    writestyle = likelyartifact
                else:
                    writestyle = strange
                worksheet2.write(row, 2, annotation, writestyle)
                row2 += 1
            row2 -= annotations1
            for annotation in variantdict[variant]["break2"]:
                if variantdict[variant]["break2"][annotation] == "exonic":
                    writestyle = inside_exon_format
                elif (
                    variantdict[variant]["break2"][annotation] == "crossgenic"
                    or variantdict[variant]["break2"][annotation] == "insidegene"
                ):
                    writestyle = inside_gene_format
                elif variantdict[variant]["break2"][annotation] == "intronic":
                    writestyle = intronic
                elif variantdict[variant]["break2"][annotation] == "likelyartifact":
                    writestyle = likelyartifact
                else:
                    writestyle = strange
                worksheet2.write(row, 4, annotation, writestyle)
                row2 += 1
            row2 -= annotations2
            if annotations1 >= annotations2:
                row2 += annotations1
            else:
                row2 += annotations2

    row += 3
    #    worksheet.write(row, 2, "breakpoint between exons", inside_gene_format)
    row += 1
    #    worksheet.write(row, 2, "breakpoint inside exon", inside_exon_format)
    row += 1
    #    worksheet.write(row, 2, "both breakpoints within same intron, or between same genes, or second breakpoint inside non-standard chromosomes", both_breaks_interexonic_intergenic)
    row += 1
    #    worksheet.write(row, 2, "insertion inside intron", insertion_intron)
    row += 1
    if canvas_info:
        for info in canvas_info:
            name, value = info.split("=")
            row += 1
            worksheet.write(row, 2, name, header_format)
            worksheet.write(row, 3, value, header_format)

    excelfile.close()
