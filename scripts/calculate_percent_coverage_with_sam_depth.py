#! /usr/bin/env python

# import python modules
import pandas as pd
from datetime import date
# from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord

import sys
import argparse
import subprocess

# calcuate percent coverage using sam depth
# add up the number of base positons with >= 30x

def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "--sam_depth")
    parser.add_argument('--sample_name')
    parser.add_argument('--base_name') #{sample_name}_A_HA-H1 or {sample_name}_A_NP
    parser.add_argument('--segment') # HA, NP
    options = parser.parse_args(args)
    return options

def calc_percent_coverage(sam_depth):
    df = pd.read_csv(sam_depth, sep = '\t', header = None)
    df = df.rename(columns = {0: 'ref', 1: 'pos', 2: 'depth'})
    total_bases = df.shape[0]
    print(df.shape)

    df_30x = df[df.depth >= 30]
    print(df_30x.shape)
    total_30x_bases = df_30x.shape[0]
    df_30x

    percent_coverage = (total_30x_bases/total_bases) * 100
    print(percent_coverage)

    return percent_coverage, total_30x_bases, total_bases

def create_output(sample_name, segment, total_bases, total_30x_bases, percent_coverage):

    df = pd.DataFrame()
    description_list = ['total_bases', 'total_30x_bases', 'percent_coverage']
    value_list = [total_bases, total_30x_bases, percent_coverage]
    df['description'] = description_list
    df['value'] = value_list
    df['sample_name'] = sample_name
    df['segment'] = segment
    
    col_order = ['sample_name', 'segment', 'description', 'value']
    df = df[col_order]

    outfile=f'percent_coverage_results.csv' 
    df.to_csv(outfile, index = False)


### MAIN ####
if __name__ == '__main__':

    print('starting script')

    options = getOptions()
    sam_depth_file = options.sam_depth
    sample_name = options.sample_name
    base_name = options.base_name
    segment = options.segment


    percent_coverage, total_30x_bases, total_bases = calc_percent_coverage(sam_depth_file)

    create_output(sample_name = sample_name, 
                    segment = segment,
                    total_bases = total_bases, 
                    percent_coverage = percent_coverage, 
                    total_30x_bases = total_30x_bases)

    
    

