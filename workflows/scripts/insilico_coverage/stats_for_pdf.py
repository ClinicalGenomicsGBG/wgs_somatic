#!/usr/bin/env python

import argparse
import logging

def load_input():
    pass

def get_pct_cov():
    pass

def genes_and_transcripts():
    pass

def get_regions_covered_below():
    pass


main():
    parser = argparse.ArgumentParser()
    parser.add_argument("my_input", nargs="?", type=str, help="Full path to input")
    args = parser.parse_args()
    logging.getLogger().setLevel(logging.INFO)

    

if __name__ == "__main__":
    main()
