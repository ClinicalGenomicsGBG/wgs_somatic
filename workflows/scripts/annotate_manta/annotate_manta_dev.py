#!/apps/bio/software/anaconda2/envs/mathias_general/bin/python3.6
import argparse
import xlsxwriter
import os
import re

def extract_variantlist(vcf):
    canvas_infofields = ["##OverallPloidy", "##DiploidCoverage", "##EstimatedTumorPurity", "##PurityModelFit", "##InterModelDistance", "##LocalSDmetric", "##EvennessScore", "##HeterogeneityProportion", "##EstimatedChromosomeCount"]
    canvas_info = []
    with open(vcf, 'r') as vcffile:
        variantlist = []
        for variant in vcffile:
            variant = variant.rstrip('\n')
            variant_info = variant.split('\t')
            if not variant_info[0].startswith('#'):
                variantlist.append(variant_info)
            else:
                if variant_info[0].startswith('#CHROM'):
                    vcf_header = variant_info
                if variant_info[0].split("=")[0] in canvas_infofields:
                    canvas_info.append(variant_info[0])

    return variantlist, vcf_header, canvas_info


def prepare_variantdict(variantlist, vcf_header):
    #columns = {
    #            "chrom": 0, 
    #            "pos": 1, 
    #            "varianttype": 2, 
    #            "reference": 3, 
    #            "variant": 4, 
    #            "filter": 6, 
    #            "info": 7,
    #            "format": 8,
    #            "normal": 9,
    #            "tumor": 10
    #            }
    variant_dict_list = []
    all_info_columns = []
    for variant in variantlist:
        variant_dict = {}
        for column_index, column_name in enumerate(vcf_header):
            variant_dict[column_name] = variant[column_index]
            # Collect all info-columnnames
            if column_name == "INFO":
                info_columns = [info_column.split("=")[0] for info_column in variant_dict["INFO"].split(";")]
                all_info_columns.extend(info_columns)
        variant_dict_list.append(variant_dict)

    # Remove all duplicate info-columnnames
    unique_info_columns = (list(set(all_info_columns)))

    # Replace variant info string with a variant info dict
    final_variant_dict_list = []
    for variant_dict in variant_dict_list:
        variant_info_dict = {}
        variant_info_list = [info_column.split("=") for info_column in variant_dict["INFO"].split(";")]
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
    with open(refseqgtf, 'r') as refseq:
        header = refseq.readline().split("\t")
        keep_cols = ["name", "chrom", "strand", "exonStarts", "exonEnds", "name2"]
        keep_cols_dict = {}
        for keepcol in keep_cols:
            for colnumber, column in enumerate(header):
                if column == keepcol:
                    keep_cols_dict.update({keepcol:colnumber})
    
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

            exon_start_list = transcript_info_list[keep_cols_dict["exonStarts"]].split(",")[:-1]
            exon_stop_list = transcript_info_list[keep_cols_dict["exonEnds"]].split(",")[:-1]

            refseqgene_dict[chrom][name2][name]["smallest"] = exon_start_list[0]
            refseqgene_dict[chrom][name2][name]["largest"] = exon_stop_list[-1]

            for exonnum, exon_start in enumerate(exon_start_list):
                real_exon_num = exonnum + 1
                if strand == "-":
                    real_exon_num = len(exon_start_list) - exonnum
                refseqgene_dict[chrom][name2][name]["exon"][real_exon_num] = {}
                exon_stop = exon_stop_list[exonnum]
                refseqgene_dict[chrom][name2][name]["exon"][real_exon_num]["start"] = exon_start
                refseqgene_dict[chrom][name2][name]["exon"][real_exon_num]["stop"] = exon_stop
        return refseqgene_dict

def find_overlapping_gene(refseq_dict, variant_chrom, variant_pos):
    overlapping_list = []
    for gene in refseq_dict[variant_chrom]:
        for transcript in refseq_dict[variant_chrom][gene]:        
            overlapping_dict = {}
            smallest = int(refseq_dict[variant_chrom][gene][transcript]["smallest"])
            largest = int(refseq_dict[variant_chrom][gene][transcript]["largest"])
            if variant_pos >= smallest and variant_pos <= largest:
                # variant is within transcript
                # within or between exons?
                for exon in refseq_dict[variant_chrom][gene][transcript]["exon"]:
                    start = int(refseq_dict[variant_chrom][gene][transcript]["exon"][exon]["start"])
                    stop = int(refseq_dict[variant_chrom][gene][transcript]["exon"][exon]["stop"])
                    if variant_pos >= start and variant_pos <= stop:
                        # variant is within exon
                        distance_start = abs(start - variant_pos)
                        distance_stop = abs(stop - variant_pos)
                        overlapping_dict["gene"] = gene
                        overlapping_dict["transcript"] = transcript
                        overlapping_dict["exon"] = f"exon{exon}"
                        overlapping_dict["minus_distance"] = distance_start
                        overlapping_dict["plus_distance"] = distance_stop
                        gene_info_string = f"{gene}:{transcript}:exon{exon} (-{distance_start}:+{distance_stop})"
                        overlapping_list.append(gene_info_string)    
                        break
                    else:
                        if variant_pos < start:
                            # variant is between this exon and the previous one
                            previous_exon_stop = int(refseq_dict[variant_chrom][gene][transcript]["exon"][previousexon]["stop"])
                            distance_upcoming_exon_start = abs(start - variant_pos)
                            distance_previous_exon_stop = abs(previous_exon_stop - variant_pos)
                            overlapping_dict["gene"] = gene
                            overlapping_dict["transcript"] = transcript
                            overlapping_dict["exon"] = f"exon{previousexon}-exon{exon}"
                            overlapping_dict["minus_distance"] = distance_previous_exon_stop
                            overlapping_dict["plus_distance"] = distance_upcoming_exon_start
                            gene_info_string = f"{gene}:{transcript}:exon{previousexon}(-{distance_previous_exon_stop})-exon{exon}(+{distance_upcoming_exon_start})"
                            overlapping_list.append(gene_info_string)
                            break
                        previousexon = exon            
    return overlapping_list

def find_nearby_gene(refseq_dict, variant_chrom, variant_pos):
    upstream_distances = {}
    downstream_distances = {}
    for gene in refseq_dict[variant_chrom]:
        # find smallest upstream and downstream distances
        for transcript in refseq_dict[variant_chrom][gene]:
            smallest = int(refseq_dict[variant_chrom][gene][transcript]["smallest"])
            largest = int(refseq_dict[variant_chrom][gene][transcript]["largest"])            
            if largest < variant_pos:
                # transcript is downstream of variant
                downstream_distances[f"{gene},{transcript}"] = abs(variant_pos - largest)
            else:
                # transcript is upstream of variant
                upstream_distances[f"{gene},{transcript}"] = abs(variant_pos - smallest)
    try:
        min_upstream = min(upstream_distances, key=upstream_distances.get)
        up_gene, up_transcript = min_upstream.split(",")
        up_distance = upstream_distances[min_upstream]
    except:
        up_gene, up_transcript, up_distance = "N/A", "N/A", "N/A"
        
    try:
        min_downstream = min(downstream_distances, key=downstream_distances.get)
        down_gene, down_transcript = min_downstream.split(",")
        down_distance = downstream_distances[min_downstream]
    except:
        down_gene, down_transcript, down_distance = "N/A", "N/A", "N/A"

    nearby_dict = {}
    nearby_dict["down_gene"] = down_gene
    nearby_dict["down_transcript"] = down_transcript
    nearby_dict["down_distance"] = down_distance
    nearby_dict["up_gene"] = up_gene
    nearby_dict["up_transcript"] = up_transcript
    nearby_dict["up_distance"] = up_distance
    gene_info_string = [f"{down_gene}:{down_transcript} (-{down_distance}) | {up_gene}:{up_transcript} (+{up_distance})"]

    return gene_info_string

def write_to_excel(output, vcfname, header, variantinfo, canvas_info):
    excelfile = xlsxwriter.Workbook(f"{output}/{vcfname}.xlsx")
    worksheet = excelfile.add_worksheet()
    worksheet.set_column(2, 2, 70)
    worksheet.set_column(4, 4, 70)
    row = 0
    column = 0
    header_format = excelfile.add_format({'bold': True, 'font_color': 'white', 'bg_color': 'black'})
    inside_gene_format = excelfile.add_format({'bg_color': 'yellow'})
    inside_exon_format = excelfile.add_format({'bg_color': 'lime'})
    for column_name in header:
        worksheet.write(row, column, column_name, header_format)
        column += 1
    
    column = 0
    for variant in variantinfo:
        # count overlapping transcripts in pos and end
        # column 2 and 4
        row += 1
        pos_overlap_count = len(variant[2])
        endpos_overlap_count = len(variant[4])
        # write everything except column 2 and 4
        for columnpos, info in enumerate(variant):
            if columnpos == 2 or columnpos == 4:
                continue
            worksheet.write(row, columnpos, info)
        for gene_info in variant[2]:
            if ":exon" in gene_info:
                if "-exon" in gene_info:
                    worksheet.write(row, 2, gene_info, inside_gene_format)
                else:
                    worksheet.write(row, 2, gene_info, inside_exon_format)
            else:
                worksheet.write(row, 2, gene_info)
            row += 1
        row -= pos_overlap_count
        for gene_info in variant[4]:
            if ":exon" in gene_info:
                if "-exon" in gene_info:
                    worksheet.write(row, 4, gene_info, inside_gene_format)
                else:
                    worksheet.write(row, 4, gene_info, inside_exon_format)
            else:
                worksheet.write(row, 4, gene_info)
            row += 1
        row -= endpos_overlap_count
        row += max(pos_overlap_count, endpos_overlap_count)
        
    row += 3
    worksheet.write(row, 2, "breakpoint between exons", inside_gene_format)
    row += 1
    worksheet.write(row, 2, "breakpoint inside exon", inside_exon_format)
    row += 1
    if canvas_info:
        for info in canvas_info:
            name, value = info.split("=")
            row += 1
            worksheet.write(row, 2, name, header_format)
            worksheet.write(row, 3, value, header_format)

    excelfile.close()

def annotate_vcf(vcf, refseq, output):
    
    # Prepare OutputNames
    vcfname = os.path.basename(vcf)
    if output.endswith("/"):
        output = output[:-1]
    # Prepare Dict of Variants
    variantlist, vcf_header, canvas_info = extract_variantlist(vcf)
    variant_dict_list, unique_info_columns = prepare_variantdict(variantlist, vcf_header)

    # Prepare Header for ExcelFile
    variant_write_table_header = ["Varianttype", "Breakpoint 1", "GeneInfo 1", "Breakpoint 2", "GeneInfo 2"]
    for columnname in variant_dict_list[0]:
        if columnname == "INFO":
            for info_columnname in unique_info_columns:
                variant_write_table_header.append(info_columnname)
        else:
                variant_write_table_header.append(columnname)

    # Prepare Dict of RefseqTranscripts
    refseq_dict = create_refseq_dict(refseq)
    # Chrom List
    chrom_list = []
    for chrom in refseq_dict:
        chrom_list.append(chrom)

    # Loop Through each Variant in dict and find overlapping / nearby genes
    variant_write_table = []
    for variant in variant_dict_list:
        variant_write = []
        variant_chrom = variant["#CHROM"]
        variant_pos = int(variant["POS"])
        variant_type = variant["ID"].split(":")[0]
        
        if variant_type == "MantaBND":
            # G]4:11470658] or [4:11470658[A
            regex = re.compile(r'[\[\]]([\w\.]*:[0-9]*)[\[\]]')
            end_location = regex.findall(variant["ALT"])[0]
            end_chrom = end_location.split(":")[0]
            end_pos = int(end_location.split(":")[1])
        else:
            end_chrom = variant_chrom
            end_pos = int(variant["INFO"]["END"])

        if variant_chrom not in chrom_list:
            pos_gene_info = ["Not available for chromosome"]
        else:
            pos_gene_info = find_overlapping_gene(refseq_dict, variant_chrom, variant_pos)
            if not pos_gene_info:
                pos_gene_info = find_nearby_gene(refseq_dict, variant_chrom, variant_pos)
                
        if not variant_type == "MantaINS":
            if end_chrom in chrom_list:
                endpos_gene_info = find_overlapping_gene(refseq_dict, end_chrom, end_pos)
                if not endpos_gene_info:
                    endpos_gene_info = find_nearby_gene(refseq_dict, end_chrom, end_pos)
            else:
                endpos_gene_info = ["Not available for chromosome"]
        else:
            endpos_gene_info = ["Not available for INS"]

        variant_write.extend([variant_type, f"{variant_chrom}:{variant_pos}", pos_gene_info, f"{end_chrom}:{end_pos}", endpos_gene_info])
        for stat in variant:
            if stat == "INFO":
                for unique_info_name in unique_info_columns:
                    variant_write.extend([variant[stat][unique_info_name]])
            else:
                variant_write.extend([variant[stat]])
        variant_write_table.append(variant_write)
    #######################################################
    # Write ExcelFile
    write_to_excel(output, vcfname, variant_write_table_header, variant_write_table, canvas_info) 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', nargs='?', help='Input MantaVCF to Annotate', required=True)
    parser.add_argument('-g', '--refseqgtf', nargs='?', help='Input Refseq GTF file', required=True)
    parser.add_argument('-o', '--output', nargs='?', help='location to output results', required=True)
    args = parser.parse_args()
    annotate_vcf(args.vcf, args.refseqgtf, args.output)
