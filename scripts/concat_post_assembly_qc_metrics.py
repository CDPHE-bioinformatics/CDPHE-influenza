#! /usr/bin/env python

#version = '1.0.0'

# import python modules
import pandas as pd
from datetime import date
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

import argparse
import subprocess

# create a single file of the results for all gene segments
# reads in an array of bam results and percent coverage resutls
# outputs a csv

# what the inputs (mapped_reads_csv, percent_coverage_csv) looks like:
# headers: sample_name, basename, segment name, gene_name, description, value
# rows: 21000000, HA, mean_depth, 25.6
# rows: 21000000, HA, num_reads_mapped, 198
# rows: 21000000, HA, perc_cov, 98.876

### LISTS AND DICTIONARIES ###
segment_list = ['HA', 'NA', 'MP', 'NP', 'NS', 'PA', 'PB1', 'PB2']
metric_variables = ['per_cov', 'mean_depth', 'num_mapped_reads', 'seq_len', 'expected_len']

col_headers = ['sample_name', 'total_segments','total_flu_mapped_reads', 'average_per_cov', 'average_mean_depth']


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "--sample_name")
    parser.add_argument( "--mapped_reads_csv_array")
    parser.add_argument( "--percent_coverage_csv_array")
    options = parser.parse_args(args)
    return options

# def create_list_from_write_lines_input(write_lines_input):
#     list = []
#     with open(write_lines_input, 'r') as f:
#         for line in f:
#             list.append(line.strip())
#     return list

def create_col_headers(segment_list, metric_variables):
    header_list = ['sample_name', 'total_segments', 'total_flu_mapped_reads']
    for segment in segment_list:
        for metric in metric_variables:
            header_name = "%s_%s" % (segment, metric)
            header_list.append(header_name)
    return header_list


#### MAIN ####
if __name__ == '__main__':

    options = getOptions()
    sample_name = options.sample_name
    mapped_reads_csv_file_list = options.mapped_reads_csv_array
    percent_coverage_csv_file_list= options.percent_coverage_csv_array
    print(mapped_reads_csv_file_list)
    print()
    print(percent_coverage_csv_file_list)
    

    # set up the pandas dataframe
    header_list = create_col_headers(segment_list = segment_list, 
                                metric_variables = metric_variables)
    df = pd.DataFrame(columns = header_list)
    df.at[0, 'sample_name'] = sample_name

    # get a list of file paths
    # mapped_reads_csv_file_list = create_list_from_write_lines_input(write_lines_input = mapped_reads_csv_txt)
    # percent_coverage_csv_file_list = create_list_from_write_lines_input(write_lines_input = percent_coverage_csvs_txt)
    # mapped_reads_csv_file_list = mapped_reads_csv_txt
    
    # insert bam results into data frame
    # track number of gene segments assemblied and total mapped reads
    num_segs = 0
    total_mapped_reads = 0

    print('/n/n')
    for mapped_reads_csv_file in mapped_reads_csv_file_list:
        num_segs = num_segs + 1
        mapped_reads_df = pd.read_csv(mapped_reads_csv_file)
        # fill in "NAs" the NA gene is being read as NA
        mapped_reads_df = mapped_reads_df.fillna('NA')
        gene_name = mapped_reads_df.gene_name[0]
        for row in range(mapped_reads_df.shape[0]):
            
            description = mapped_reads_df.description[row]
            value = mapped_reads_df.value[row]

            if description == 'num_mapped_reads':
                total_mapped_reads = total_mapped_reads + value

            # get correct column header name
            col_name = "%s_%s" % (gene_name, description)
            df.at[0, col_name] = value

    # insert per cov results into data frame
    for percent_coverage_csv_file in percent_coverage_csv_file_list:
        percent_coverage_df = pd.read_csv(percent_coverage_csv_file)
        percent_coverage_df = percent_coverage_df.fillna("NA")
        gene_name = percent_coverage_df.gene_name[0]
        for row in range(percent_coverage_df.shape[0]):
            description = percent_coverage_df.description[row]
            value = percent_coverage_df.value[row]
            
            # get correct column header name
            col_name = "%s_%s" % (gene_name, description)
            df.at[0, col_name] = value


    # add in final columns
    df.at[0, 'total_segments'] = num_segs
    df.at[0, 'total_flu_mapped_reads'] = total_mapped_reads
    

    # save df
    outfile = f"{sample_name}_assembly_qc_metrics.csv"
    df.to_csv(outfile, index = False)
