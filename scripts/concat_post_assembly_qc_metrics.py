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

# what input looks like:
# headers: sample_name, segment, description, value
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
    parser.add_argument( "--bam_results")
    parser.add_argument( "--per_cov_results")
    parser.add_argument( "--ivar_parameters")
    options = parser.parse_args(args)
    return options

def create_list_from_write_lines_input(write_lines_input):
    list = []
    with open(write_lines_input, 'r') as f:
        for line in f:
            list.append(line.strip())
    return list

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
    bam_results_txt = options.bam_results
    per_cov_results_txt = options.per_cov_results
    ivar_parameters_path = options.ivar_parameters

    # set up the pandas dataframe
    header_list = create_col_headers(segment_list = segment_list, 
                                metric_variables = metric_variables)
    df = pd.DataFrame(columns = header_list)
    df.at[0, 'sample_name'] = sample_name

    # get a list of file paths
    bam_results_file_list = create_list_from_write_lines_input(write_lines_input = bam_results_txt)
    per_cov_results_file_list = create_list_from_write_lines_input(write_lines_input = per_cov_results_txt)
    
    # insert bam results into data frame
    # track number of gene segments assemblied and total mapped reads
    num_segs = 0
    total_mapped_reads = 0

    print('/n/n')
    for bam_result_file in bam_results_file_list:
        num_segs = num_segs + 1
        bam_df = pd.read_csv(bam_result_file)
        # fill in "NAs" the NA gene is being read as NA
        bam_df = bam_df.fillna('NA')
        gene_name = bam_df.gene_name[0]
        for row in range(bam_df.shape[0]):
            
            description = bam_df.description[row]
            value = bam_df.value[row]

            if description == 'num_mapped_reads':
                total_mapped_reads = total_mapped_reads + value

            # get correct column header name
            col_name = "%s_%s" % (gene_name, description)
            df.at[0, col_name] = value

    # insert per cov results into data frame
    for per_cov_result_file in per_cov_results_file_list:
        per_cov_df = pd.read_csv(per_cov_result_file)
        per_cov_df = per_cov_df.fillna("NA")
        gene_name = per_cov_df.gene_name[0]
        for row in range(per_cov_df.shape[0]):
            description = per_cov_df.description[row]
            value = per_cov_df.value[row]
            
            # get correct column header name
            col_name = "%s_%s" % (gene_name, description)
            df.at[0, col_name] = value


    # add in final columns
    df.at[0, 'total_segments'] = num_segs
    df.at[0, 'total_flu_mapped_reads'] = total_mapped_reads
    
    # add in ivar parameters
    ivar_df = pd.read_csv(ivar_parameters_path)
    df.at[0, "ivar_version"] = ivar_df.loc[0, 'ivar_version']
    df.at[0, 'ivar_docker'] = ivar_df.loc[0, 'ivar_docker']
    df.at[0, 'ivar_min_depth'] = ivar_df.loc[0, 'ivar_min_depth']
    df.at[0, 'ivar_min_freq'] = ivar_df.loc[0, 'ivar_min_freq']
    df.at[0, 'ivar_min_qual'] = ivar_df.loc[0, 'ivar_min_qual']

    # save df
    outfile = "%s_assembly_qc_metrics.csv" % sample_name
    df.to_csv(outfile, index = False)
