#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python
# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import argparse
import xlsxwriter
import os
from workflows.scripts.determine_match import determine_match
from tools.git_versions import get_git_commit, get_git_tag, get_git_reponame
from workflows.scripts.sex import calc_sex
import time
import pandas as pd
import gzip


def extract_stats(statsfile, statstype, sampletype, statsdict):
    with open(statsfile, 'r') as statsfile:
        if statstype not in statsdict:
            statsdict[statstype] = {}
        statsdict[statstype][sampletype] = {}
        
        for row in statsfile:
            row = row.rstrip()
            row_list = row.split("\t")
            if 0 in statsdict[statstype][sampletype]:
                statsvalues = row_list
                for column, statsvalue in enumerate(statsvalues):
                    statsvalue = statsvalue.replace(".", ",")
                    statsdict[statstype][sampletype][column]["colvalue"] = statsvalue
                break
            if row_list[0] == "GENOME_TERRITORY" or row_list[0] == "LIBRARY":
                headernames = row_list
                headercount = 0
                for headername in headernames:
                    statsdict[statstype][sampletype][headercount] = {} 
                    statsdict[statstype][sampletype][headercount]["colname"] = headername
                    headercount += 1
        return statsdict


def get_canvas_tumorinfo(canvasvcf):
    canvasdict = {}
    canvas_infofields = ["##OverallPloidy", "##DiploidCoverage", "##EstimatedTumorPurity", "##PurityModelFit", "##InterModelDistance", "##LocalSDmetric", "##EvennessScore", "##HeterogeneityProportion", "##EstimatedChromosomeCount"]
    with gzip.open(canvasvcf, 'rt') as vcf:
        for variant in vcf:
            variant = variant.rstrip('\n')
            variant_info = variant.split('\t')
            if variant_info[0].startswith('#'):
                if variant_info[0].split("=")[0] in canvas_infofields:
                    canvasfield = variant_info[0].split("=")[0]
                    canvasfield = canvasfield.replace("#", "")
                    canvasfield_value = variant_info[0].split("=")[1]
                    canvasfield_value = canvasfield_value.replace(".", ",")
                    canvasdict[canvasfield] = canvasfield_value
    return canvasdict


def read_tmb_file(filepath):
    tmb_dict = {}
    try:
        with open(filepath, 'r') as file:
            for line in file:
                key, value = line.strip().split('\t')
                try:
                    # Try to convert the value to a float or int if possible
                    if '.' in value:
                        tmb_dict[key] = float(value)
                    else:
                        tmb_dict[key] = int(value)
                except ValueError:
                    # If conversion fails, store the value as a string
                    tmb_dict[key] = value
    except FileNotFoundError:
        print(f"No file found at {filepath}")
    except Exception as e:
        print(f"An error occurred: {e}")
    return tmb_dict


def get_msi_info(msi, msi_red):
    msi_dict = {}
    try:
        with open(msi, 'r') as file:
            headers = file.readline().strip().split('\t')
            values = file.readline().strip().split('\t')
            for header, value in zip(headers, values):
                try:
                    if '.' in value:
                        msi_dict[f"msi_{header}"] = float(value)
                    else:
                        msi_dict[f"msi_{header}"] = int(value)
                except ValueError:
                    msi_dict[f"msi_{header}"] = value

        with open(msi_red, 'r') as file:
            headers = file.readline().strip().split('\t')
            values = file.readline().strip().split('\t')
            for header, value in zip(headers, values):
                try:
                    if '.' in value:
                        msi_dict[f"msi_filtered_{header}"] = float(value)
                    else:
                        msi_dict[f"msi_filtered_{header}"] = int(value)
                except ValueError:
                    msi_dict[f"msi_filtered_{header}"] = value

    except FileNotFoundError:
        print(f"No file found at {msi} or {msi_red}")
    except Exception as e:
        print(f"An error occurred: {e}")
    return msi_dict


def get_ascat_tumorinfo(ascatstats):
    """
    Reads ASCAT statistics from a file and returns a dictionary with relevant metrics.
    """
    ascatdict = {}
    ascat_filter_rows = ["ploidy", "purity", "goodness_of_fit"]
    with open(ascatstats, mode='r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue  # Skip empty lines and comments
            row = line.split("\t")
            if len(row) < 2:
                continue  # Skip rows that do not have at least two columns
            metric = row[0].strip()
            value = row[1].strip()
            if metric in ascat_filter_rows:
                try:
                    ascatdict[metric] = float(value)
                except ValueError:
                    ascatdict[metric] = value
    return ascatdict

def read_sample_purities(info_files):
    """
    Reads Sample_Purity values from multiple pileup_info.txt files.
    """
    purities = {}
    for info_file in info_files:
        try:
            with open(info_file, 'r') as file:
                ploidy = None
                purity = None
                for line in file:
                    line = line.strip()
                    if line.startswith("Output_Ploidy"):
                        ploidy = int(line.split("\t")[1])
                    elif line.startswith("Sample_Purity"):
                        purity = float(line.split("\t")[1])
                    # Add to dictionary when both ploidy and purity are found
                    if ploidy is not None and purity is not None:
                        if ploidy not in purities:
                            purities[ploidy] = []
                        purities[ploidy].append(purity)
                        ploidy = None  # Reset for the next entry
                        purity = None
        except FileNotFoundError:
            print(f"No file found at {info_file}")
        except Exception as e:
            print(f"An error occurred while reading {info_file}: {e}")
    return purities


def create_excel(statsdict, output, normalname='', tumorname='', match_dict={}, canvasdict={}, sex='', tmb_dict={}, msi_dict={}, ascatdict={}):
    current_date = time.strftime("%Y-%m-%d")
    excelfile = xlsxwriter.Workbook(output)
    worksheet = excelfile.add_worksheet("qc_stats")
    worksheet.set_column(0, 30, 20)
    cellformat = {}
    cellformat["header"] = excelfile.add_format({'bold': True, 'font_color': 'white', 'bg_color': 'black'})
    cellformat["section"] = excelfile.add_format({'bold': True, 'font_color': 'FFA764', 'bg_color': '878787'})
    cellformat["tumorname"] = excelfile.add_format({'bold': True, 'font_color': 'white', 'bg_color': 'BF000D'})
    cellformat["normalname"] = excelfile.add_format({'bold': True, 'font_color': 'white', 'bg_color': '00BF19'})
    cellformat["warning"] = excelfile.add_format({'bg_color': 'FF9A00'})
    cellformat["error"] = excelfile.add_format({'bg_color': 'FF0000'})
    cellformat["pass"] = excelfile.add_format({'bg_color': '95FF80'})
    cellformat["malesex"] = excelfile.add_format({'bold': True, 'bg_color': '89CFF0', 'font_size': 13})
    cellformat["femalesex"] = excelfile.add_format({'bold': True, 'bg_color': 'F4C2C2', 'font_size': 13})

    row = 1
    worksheet.write(row, 0, f"QC-report created: {current_date}")
    row += 1

    # Version numbers & tag
    worksheet.merge_range('A3:C3', f"{get_git_reponame()} tag: {get_git_tag()}, commit: {get_git_commit()}")

    print(f"STATSDICT: {statsdict}")

    # Input calculated sex
    if sex.lower() == "male":
        worksheet.merge_range('D3:E3', f"Computed sex of patient: {sex}", cellformat["malesex"])
    elif sex.lower() == "female":
        worksheet.merge_range('D3:E3', f"Computed sex of patient: {sex}", cellformat["femalesex"])
    else:
        worksheet.merge_range('D3:E3', f"Something went wrong with sex calculation: {sex}", cellformat["error"])

    # Specify sample used for sex calculation
    if "normal" in statsdict["coverage"].keys():
        worksheet.merge_range('D4:G4', f"Normal sample {normalname} used for sex calculation", cellformat["pass"])
    else:
        worksheet.merge_range('D4:G4', f"Warning: tumour sample {tumorname} used for sex calculation", cellformat["warning"])
    row += 2

    # Coverage stats
    for statstype in statsdict:
        header = False
        if statstype == "coverage":
            worksheet.write(row, 0, "COVERAGE STATS", cellformat["section"])
        else:
            worksheet.write(row, 0, "MAPPING STATS", cellformat["section"])
        row += 1
        for sampletype in statsdict[statstype]:
            if sampletype == "tumor":
                sname = tumorname
                nameformat = cellformat["tumorname"]
            else:
                sname = normalname
                nameformat = cellformat["normalname"]
            # loop through values in stats --> sample
            column_s = 1
            if not header:
                for column in statsdict[statstype][sampletype]:
                    write_col = column_s + column
                    worksheet.write(row, 0, "Samplename", cellformat["header"])
                    worksheet.write(row, write_col, statsdict[statstype][sampletype][column]["colname"], cellformat["header"])
                row += 1
                header = True
            for column in statsdict[statstype][sampletype]:
                write_col = column_s + column
                worksheet.write(row, 0, sname, nameformat)
                worksheet.write(row, write_col, statsdict[statstype][sampletype][column]["colvalue"])
            row += 1
        row += 1
    row += 1
    worksheet.write(row, 0, "TUMOR/NORMAL MATCH-CHECK", cellformat["section"])
    row += 1
    if match_dict:
        if "ERROR" in match_dict["match_status"]:
            worksheet.write(row, 0, "Error occurred", cellformat["error"])
        else:
            headerrow = row
            statrow = row + 1
            column_num = 0
            for stat in match_dict:
                worksheet.write(headerrow, column_num, stat, cellformat["header"])
                if stat == "match_status" or stat == "match_fraction":
                    if "warning" in match_dict["match_status"]:
                        style = "warning"
                    elif "error" in match_dict["match_status"]:
                        style = "error"
                    else:
                        style = "pass"
                    worksheet.write(statrow, column_num, match_dict[stat], cellformat[style])
                else:
                    worksheet.write(statrow, column_num, match_dict[stat])

                column_num += 1
    row += 4
    worksheet.write(row, 0, "CANVAS-STATS", cellformat["section"])
    worksheet.write(row, 1, tumorname, cellformat["tumorname"])
    row += 1
    if canvasdict:
        for key in canvasdict:
            worksheet.write(row, 0, key, cellformat["header"])
            worksheet.write(row, 1, canvasdict[key])
            row += 1

    row += 2
    worksheet.write(row, 0, "TMB-CALCULATION", cellformat["section"])
    worksheet.write(row, 1, tumorname, cellformat["tumorname"])
    row += 1
    if tmb_dict:
        for key in tmb_dict:
            worksheet.write(row, 0, key, cellformat["header"])
            worksheet.write(row, 1, tmb_dict[key])
            row += 1

    row += 2
    worksheet.write(row, 0, "MSI-STATS", cellformat["section"])
    worksheet.write(row, 1, tumorname, cellformat["tumorname"])
    row += 1
    if msi_dict:
        for key in msi_dict:
            worksheet.write(row, 0, key, cellformat["header"])
            worksheet.write(row, 1, msi_dict[key])
            row += 1

    row += 2
    worksheet.write(row, 0, "ASCAT-STATS", cellformat["section"])
    worksheet.write(row, 1, tumorname, cellformat["tumorname"])
    row += 1
    if ascatdict:
        for key in ascatdict:
            worksheet.write(row, 0, key, cellformat["header"])
            worksheet.write(row, 1, ascatdict[key])
            row += 1

    row += 2

    excelfile.close()


def create_excel_main(tumorcov='', ycov='', normalcov='', tumordedup='', normaldedup='', tumorvcf='', normalvcf='', canvasvcf='', tmb='', msi='', msi_red='', ascatstats='', tumor_info_files=[], output=''):
    statsdict = {}
    if tumorcov:
        tumorcovfile = os.path.basename(tumorcov)
        tumorname = tumorcovfile.replace("_WGScov.tsv", "")
        statsdict = extract_stats(tumorcov, "coverage",  "tumor", statsdict)
        statsdict = extract_stats(tumordedup, "dedup",  "tumor", statsdict)
        tmb_dict = read_tmb_file(tmb)
        ascatdict = get_ascat_tumorinfo(ascatstats)
    if normalcov:
        normalcovfile = os.path.basename(normalcov)
        normalname = normalcovfile.replace("_WGScov.tsv", "")
        statsdict = extract_stats(normalcov, "coverage", "normal", statsdict)
        statsdict = extract_stats(normaldedup, "dedup", "normal", statsdict)
    if tumorcov and normalcov:
        match_dict = determine_match(normalvcf, tumorvcf, 400000)
        canvas_dict = get_canvas_tumorinfo(canvasvcf)
        msi_dict = get_msi_info(msi, msi_red)

    if not output.endswith(".xlsx"):
        output = f"{output}.xlsx"

    # Determine which files to use to calculate sex
    if tumorcov:
        if normalcov:
            # Tumour + Normal
            calculated_sex = calc_sex(normalcov, ycov)
            create_excel(statsdict, output, normalname, tumorname, match_dict, canvas_dict, sex=calculated_sex, tmb_dict=tmb_dict, msi_dict=msi_dict, ascatdict=ascatdict)
        else:
            # Tumour only
            calculated_sex = calc_sex(tumorcov, ycov)
            create_excel(statsdict, output, tumorname=tumorname, sex=calculated_sex, tmb_dict=tmb_dict, ascatdict=ascatdict)
    else:
        # Normal only
        calculated_sex = calc_sex(normalcov, ycov)
        create_excel(statsdict, output, normalname, sex=calculated_sex)

def create_qc_toaggregate(tumorcov='', ycov='', normalcov='', tumordedup='', normaldedup='', tumorvcf='', normalvcf='', output='', tmb = '', tumorid = '', normalid = ''):

    statsdict = {}
    if tumorcov:
        statsdict = extract_stats(tumorcov, "coverage",  "tumor", statsdict)
        statsdict = extract_stats(tumordedup, "dedup",  "tumor", statsdict)
        tmb_dict = read_tmb_file(tmb)
        if normalcov:
            calculated_sex = calc_sex(normalcov, ycov)
            match_dict = determine_match(normalvcf, tumorvcf, 400000)
        else:
            calculated_sex = calc_sex(tumorcov, ycov)
            
    if normalcov:
        statsdict = extract_stats(normalcov, "coverage", "normal", statsdict)
        statsdict = extract_stats(normaldedup, "dedup", "normal", statsdict)
        if not tumorcov:
            calculated_sex = calc_sex(normalcov, ycov)

    # Prepare data for tumor and normal rows
    rows = []

    if tumorcov:
        tumor_row = {
            "Type": "Tumor",
            "SampleID": tumorid,
            "Sex": calculated_sex if normalcov else "N/A",
            "Coverage 10x": next((item["colvalue"] for item in statsdict["coverage"]["tumor"].values() if item["colname"] == "PCT_10X"), "N/A"),
            "Coverage 30x": next((item["colvalue"] for item in statsdict["coverage"]["tumor"].values() if item["colname"] == "PCT_30X"), "N/A"),
            "Coverage 60x": next((item["colvalue"] for item in statsdict["coverage"]["tumor"].values() if item["colname"] == "PCT_60X"), "N/A"),
            "Mean coverage": next((item["colvalue"] for item in statsdict["coverage"]["tumor"].values() if item["colname"] == "MEAN_COVERAGE"), "N/A"),
            "Median coverage": next((item["colvalue"] for item in statsdict["coverage"]["tumor"].values() if item["colname"] == "MEDIAN_COVERAGE"), "N/A"),
            "SD coverage": next((item["colvalue"] for item in statsdict["coverage"]["tumor"].values() if item["colname"] == "SD_COVERAGE"), "N/A"),
            "Match fraction": match_dict["match_fraction"] if normalcov else "N/A",
            "Match status": match_dict["match_status"] if normalcov else "N/A",
            "TMB": tmb_dict.get("TMB", "N/A"),
        }
        rows.append(tumor_row)

    if normalcov:
        normal_row = {
            "Type": "Normal",
            "SampleID": normalid,
            "Sex": calculated_sex,
            "Coverage 10x": next((item["colvalue"] for item in statsdict["coverage"]["normal"].values() if item["colname"] == "PCT_10X"), "N/A"),
            "Coverage 30x": next((item["colvalue"] for item in statsdict["coverage"]["normal"].values() if item["colname"] == "PCT_30X"), "N/A"),
            "Coverage 60x": next((item["colvalue"] for item in statsdict["coverage"]["normal"].values() if item["colname"] == "PCT_60X"), "N/A"),
            "Mean coverage": next((item["colvalue"] for item in statsdict["coverage"]["normal"].values() if item["colname"] == "MEAN_COVERAGE"), "N/A"),
            "Median coverage": next((item["colvalue"] for item in statsdict["coverage"]["normal"].values() if item["colname"] == "MEDIAN_COVERAGE"), "N/A"),
            "SD coverage": next((item["colvalue"] for item in statsdict["coverage"]["normal"].values() if item["colname"] == "SD_COVERAGE"), "N/A"),
            "Match fraction": match_dict["match_fraction"] if tumorcov else "N/A",
            "Match status": match_dict["match_status"] if tumorcov else "N/A",
            "TMB": "N/A",
        }
        rows.append(normal_row)

    if not output.endswith(".xlsx"):
        output = f"{output}.xlsx"
    # Create a DataFrame and save to Excel
    df = pd.DataFrame(rows, columns=["Type", "SampleID", "Sex", "Coverage 10x", "Coverage 30x", "Coverage 60x", "Mean coverage", "Median coverage", "SD coverage", "Match fraction", "Match status", "TMB"])
    df.to_excel(output, index=False)
    # Save to TSV
    tsv_output = output.replace(".xlsx", ".tsv")
    df.to_csv(tsv_output, sep="\t", index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-tc', '--tumorcov', nargs='?', help='Sentieon WGS cov file from tumorbam', required=False)
    parser.add_argument('-yc', '--ycov', nargs='?', help='Sentieon Y-cov file', required=False)
    parser.add_argument('-nc', '--normalcov', nargs='?', help='Sentieon WGS cov file from normalbam', required=True)
    parser.add_argument('-td', '--tumordedup', nargs='?', help='Sentieon dedup-stats for tumorbam', required=False)
    parser.add_argument('-nd', '--normaldedup', nargs='?', help='Sentieon dedup-stats for normalbam', required=True)
    parser.add_argument('-tv', '--tumorvcf', nargs='?', help='Tumor Germlinecalls', required=False)
    parser.add_argument('-nv', '--normalvcf', nargs='?', help='Normal Germlinecalls', required=True)
    parser.add_argument('-cv', '--canvasvcf', nargs='?', help='Somatic Canvas VCF', required=False)
    parser.add_argument('--tmb', nargs='?', help='TMB file', required=False)
    parser.add_argument('--msi', nargs='?', help='MSI result file', required=False)
    parser.add_argument('--msi_red', nargs='?', help='MSI filtered to bed file result file', required=False)
    parser.add_argument('-as', '--ascatstats', nargs='?', help='Somatic Ascat stats', required=False)
    parser.add_argument('--tumor_info_files', nargs='*', help='List of tumor pileup_info.txt files for different ploidies', required=False)
    parser.add_argument('-o', '--output', nargs='?', help='fullpath to file to be created (xlsx will be appended if not written)', required=True)
    args = parser.parse_args()
    create_excel_main(
        args.tumorcov, args.ycov, args.normalcov, args.tumordedup, args.normaldedup,
        args.tumorvcf, args.normalvcf, args.canvasvcf, args.tmb, args.msi, args.msi_red,
        args.ascatstats, args.tumor_info_files, args.output
    )
