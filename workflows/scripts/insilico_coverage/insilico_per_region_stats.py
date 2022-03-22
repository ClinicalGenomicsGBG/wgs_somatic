#!/usr/bin/env python
import os
import csv
import sys
import argparse
import statistics

def count_bases_in_threshold(coverage_list, threshold):
    total_positions = len(coverage_list)
    count_positions = 0
    for position in coverage_list:
        if position >= threshold:
            count_positions += 1
    percent_covered = float(count_positions*100) / total_positions
    percent_covered = float("{0:.3f}".format(percent_covered))
    return percent_covered

def calculate_mean(coverage_list):
    mean_calc = sum(coverage_list) / len(coverage_list) 
    mean = float("{0:.3f}".format(mean_calc))
    return mean

def calculate_median(coverage_list):
    median = statistics.median(coverage_list)
    return median

def calculate_stdev(coverage_list):
    stdev_calc = statistics.stdev(coverage_list)
    stdev = float("{0:.3f}".format(stdev_calc))
    return stdev

def group_region_coverage_positions(coveragefile, levels):
    region_dict = {}
    
    # Create Unique Annotation_list
    region_name_list = []
    for covpos in coveragefile:
        covpos = covpos.rstrip('\n')
        annotation_position = 2
        covpos = covpos.split(",")
        region_name = []
        for n in range(levels):
            annotation_position += 1
            try:
                region_name.append(covpos[annotation_position])
            except:
                region_name.append("Missing")
                print("Warning: desired annotation level does not exist in coveragefile at:")
                print(covpos)
        region_name_str = '_'.join(region_name)
        if region_dict.get(region_name_str):
            region_dict[region_name_str]["coveragelist"].append(int(covpos[2]))
        else:
            region_dict[region_name_str] = {}
            annotation_position = 2
            for n in range(levels):
                annotation_position += 1
                region_dict[region_name_str][n] = covpos[annotation_position]
            region_dict[region_name_str]["coveragelist"] = []
            region_dict[region_name_str]["coveragelist"].append(int(covpos[2]))
    return region_dict

def coverage_stats(level, coverage, thresholds, output):
    threshold_list = thresholds.split(",")
    level = int(level)
    if output.endswith('/'):
        output = output[:-1]
    with open(coverage, 'r') as coverage_file:
        # 4th column is annotation-region 1
        # 5th column is annotation-region 2 (optional)
        # 6th column is annotation-region 3 (optional)
        region_dict = group_region_coverage_positions(coverage_file, int(level))
        with open(f"{output}", 'w+') as csvout:
            fieldnames = []
            if level == 1:
                fieldnames.append("Gene/Transcript")
            if level == 2:
                fieldnames.append("Gene")
                fieldnames.append("Transcript")
            if level == 3:
                fieldnames.append("Gene")
                fieldnames.append("Transcript")
                fieldnames.append("Exon")

            for threshold in threshold_list:
                fieldnames.extend([f"%>={threshold}X"])
        
            fieldnames.extend(["median"])
            csv_writer = csv.writer(csvout)
            csv_write_list = []
            csv_write_list.append(fieldnames)

            for region in region_dict:
                    #region_values = region_dict[region]
                    region_write_values = []
                    for n in range(int(level)):
                        region_write_values.append(region_dict[region][n])
                    coverage_list = region_dict[region]["coveragelist"]
                    #region_dict[region]["mean"] = str(calculate_mean(coverage_list))
                    #region_dict[region]["stdev"] = str(calculate_stdev(coverage_list))
                    region_dict[region]["median"] = str(calculate_median(coverage_list))
                    for threshold in threshold_list:
                        region_dict[region][f"%>={threshold}X"] = str(count_bases_in_threshold(coverage_list, int(threshold)))
                        region_write_values.append(region_dict[region][f"%>={threshold}X"])                        
                    region_write_values.append(region_dict[region]["median"])
                    csv_write_list.append(region_write_values)
            csv_writer.writerows(csv_write_list)                    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--level', nargs='?', choices=['1','2','3'], help='on which level do you want stats? (1 tends to be gene, 2 transcript, 3 exon) ', required=True)
    parser.add_argument('-c', '--coverage', nargs='?', help='(from samtools mpileup -l $BED -Q 1 $BAM | cut -f1,2,4) annotated with gene,transcript,exon depending on level required', required=True)    
    parser.add_argument('-t', '--thresholds', nargs='?', help='list of coverage thresholds, 1,10,20,40,50 etc...', required=True)
    parser.add_argument('-o', '--output', nargs='?', help='(Output location of csv-file with per region coverage stats)', required=True)
    args = parser.parse_args()
    coverage_stats(args.level, args.coverage, args.thresholds, args.output)

