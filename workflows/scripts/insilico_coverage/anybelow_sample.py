#!/usr/bin/env python

import argparse
import csv
#import logging

# Takes input file and returns a csv dicreader object containing only the positions below coverage limit
def read_input(infile, limit):
    limit = int(limit)
    below_limit = []
    with open(infile, "r") as inp:
        fieldnames = ["chr", "pos", "cov", "gene", "transcript", "exon"]
        reader = csv.DictReader(inp, delimiter=",", fieldnames=fieldnames)
        for i in reader:
            if int(i["cov"]) >= limit:
                continue
            else:
                below_limit.append(i)
    return below_limit


# Takes the csv object containing only the regions below and merges adjacent positions into a single region
def merge_adjacent(indata):
    lastpos = 0
    mergedfile = []
    currentregion = []
    for h, i in enumerate(indata):
        try:        
            if int(i["pos"]) == lastpos + 1:
                currentregion.append(i)
            elif (int(i["pos"]) + 1 < int(indata[h + 1]["pos"]) or not (i["chr"] == indata[h + 1]["chr"])):
                i["start"] = i["pos"]
                i["stop"] = i["start"]
                mergedfile.append(i)
            else:
                if len(currentregion) > 0:
                    mergedfile.append(collapser(currentregion))                
                currentregion = []
        except IndexError as ie:
            i["start"] = i["pos"]
            i["stop"] = i["start"]
            mergedfile.append(i)
        lastpos = int(i["pos"])
    print("FINAL PRINT")
    for thing in mergedfile:
        print(thing)
    return mergedfile


# Collapses positions adjacent into single region
def collapser(regionlist):
    positions = []
    for i in regionlist:
        positions.append(int(i["pos"]))
    firstpos = min(positions)
    lastpos = max(positions)
    if firstpos == lastpos:
        firstpos -= 1
    result_dict = {"chr": regionlist[0]["chr"], "start": f"{firstpos}", "stop": f"{lastpos}", "gene": regionlist[0]["gene"], "transcript": regionlist[0]["transcript"],
                   "exon": regionlist[0]["exon"]}
    return result_dict

# Takes ready data and writes to file
def write_merged(merge_in, outdata):
    with open(outdata, "w+") as outfile:
            fieldnames = ["chr", "start", "stop", "gene", "transcript", "exon"]
            writer = csv.DictWriter(outfile, delimiter=",", fieldnames=fieldnames)
            for i in merge_in:
                writer.writerow({"chr": str(i["chr"]),
                                 "start": str(i["start"]),
                                 "stop": str(i["stop"]),
                                 "gene": str(i["gene"]),
                                 "transcript": str(i["transcript"]),
                                 "exon": str(i["exon"])})
    return outfile


#Function to be run from other scripts
def run(infile,  limits=[20, 10]):
    outs = []
    for i in limits:
        outs.append(merge_adjacent(read_input(infile, i)))
    return outs


def alt_main(csv_input, csv_output_10x, csv_output_20x):
    outfiles = [csv_output_10x, csv_output_20x]
    limits = [10, 20]
    for limit, outfile in zip(limits, outfiles):
        fileloc = write_merged(merge_adjacent(read_input(csv_input, limit)), outfile)
        print(f"The file for {limit}X is located at {fileloc.name}")

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("csv_input", nargs="?", type=str, help="Full path to input csv from previous script")
    parser.add_argument("csv_output_10x", nargs="?", type=str, help="Full path to csv output of below 10x regions")
    parser.add_argument("csv_output_20x", nargs="?", type=str, help="Full path to csv output of below 20x regions")

    args = parser.parse_args()

#    logging.getLogger().setLevel(logging.INFO)

    outfiles = [args.csv_output_10x, args.csv_output_20x]
    limits = [10, 20]

    for limit, outfile in zip(limits, outfiles):
        fileloc = write_merged(merge_adjacent(read_input(args.csv_input, limit)), outfile)
        print(f"The file for {limit}X is located at {fileloc.name}")


if __name__ == "__main__":
    main()
