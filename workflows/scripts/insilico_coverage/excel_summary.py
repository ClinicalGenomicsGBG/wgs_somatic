#!/usr/bin/env python
import pandas as pd
import argparse


def get_input(my_input):
    with open(my_input, "r") as my_open_file:
        my_data = pd.read_csv(my_open_file, sep=',', dtype={"Gene": str, "Transcript": str, "%>=10X": float, "%>=20X": float, "median": float})
        my_df = pd.DataFrame(my_data)
        return my_df


def write_excel(my_infile, outfile):
    indata = get_input(my_infile)
    below_10x = indata.loc[indata['%>=10X'] < 100.0]
    below_10x.to_excel(outfile)


def main():
    parser = argparse.ArgumentParser(prog='below10x convert to Excel')
    parser.add_argument("-i", "--infile", help="Full path to csv file")
    parser.add_argument("-o", "--outfile", help="Full path to result Excel-file")
    args = parser.parse_args()

    #raw_infile = "/medstore/results/wgs/KG/190920_A00687_0032_AHMVY7DSXX/DNA58618_190920_AHMVY7DSXX/new_insilico_panels/epilepsi.v1.0.csv"
    raw_infile = args.infile
    raw_outfile = args.outfile

    write_excel(raw_infile, raw_outfile)

if __name__ == "__main__":
    main()
