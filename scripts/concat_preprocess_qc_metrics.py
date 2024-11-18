#! /usr/bin/env python

#version = '1.0.0'

# import python modules
import pandas as pd

import sys
import argparse
import subprocess


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "--total_reads_R1_raw")
    parser.add_argument( "--total_reads_R2_raw")
    parser.add_argument( "--read_length_R1_raw")
    parser.add_argument( "--read_length_R2_raw")
    parser.add_argument( "--read_pairs_raw")
    parser.add_argument( "--total_reads_R1_cleaned")
    parser.add_argument( "--total_reads_R2_cleaned")
    parser.add_argument( "--read_length_R1_cleaned")
    parser.add_argument( "--read_length_R2_cleaned")
    parser.add_argument( "--read_pairs_cleaned")
    parser.add_argument( "--sample_name")
    options = parser.parse_args(args)
    return options


#### MAIN ####
if __name__ == '__main__':

    options = getOptions()
    sample_name = options.sample_name
    print(sample_name)


    total_reads_R1_raw = options.total_reads_R1_raw
    total_reads_R2_raw = options.total_reads_R2_raw
    read_length_R1_raw = options.read_length_R1_raw
    read_length_R2_raw = options.read_length_R2_raw
    read_pairs_raw = options.read_pairs_raw

    total_reads_R1_cleaned = options.total_reads_R1_cleaned
    total_reads_R2_cleaned = options.total_reads_R2_cleaned
    read_length_R1_cleaned = options.read_length_R1_cleaned
    read_length_R2_cleaned = options.read_length_R2_cleaned
    read_pairs_cleaned = options.read_pairs_cleaned


    df = pd.DataFrame()
    df.at[0, 'sample_name'] = sample_name

    df['read_diff'] = int(read_pairs_raw) - int(read_pairs_cleaned)

    df['read_length_R1_raw'] = read_length_R1_raw
    df['read_length_R2_raw'] = read_length_R2_raw
    df['reads_R1_raw'] = total_reads_R1_raw
    df['reads_R2_raw'] = total_reads_R1_raw
    df["read_pairs_raw"] = read_pairs_raw
    df["reads_raw"] = int(read_pairs_raw) * 2

    
    df['read_length_R1_cleaned'] = read_length_R1_cleaned
    df['read_length_R2_cleaned'] = read_length_R2_cleaned
    df['reads_R1_cleaned'] = total_reads_R1_cleaned
    df['reads_R2_cleaned'] = total_reads_R1_cleaned
    df["read_pairs_cleaned"] = read_pairs_cleaned
    df["reads_cleaned"] = int(read_pairs_cleaned) * 2
    


    outfile =  "%s_preprocess_qc_metrics.csv" % sample_name
    df.to_csv(outfile, index = False)