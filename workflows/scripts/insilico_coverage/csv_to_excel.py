#!/usr/bin/env python
import argparse
import pandas as pd
import os

def csv_to_excel(raw_infile, my_outfile):
    my_dataframe = pd.read_csv(raw_infile, sep=',', dtype={"Chr": int, "Start": int, "Stop": int, "Gene": str, "Transcript": str, "Exon": str}, names=["Chromosome", "Start", "Stop", "Gene", "Transcript", "Exon"])
    my_dataframe.to_excel(my_outfile)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Full path to csv file")
    parser.add_argument("-o", "--outfile", help="Full path to result Excel-file")
    args = parser.parse_args()
    csv_to_excel(args.infile, args.outfile)

if __name__ == "__main__":
    main()
