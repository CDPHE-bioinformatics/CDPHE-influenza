#! /usr/bin/env python

#version = '1.0.0'

# import python modules
import pandas as pd
from datetime import date
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import re

import argparse
import subprocess

# create a single file of the results for all gene segments
# reads in an array of bam results and percent coverage resutls
# outputs a csv

# what the inputs (bam_stats_csv, percent_coverage_csv) looks like:
# headers: sample_name, basename, segment name, gene_name, description, value
# rows: 21000000, HA, mean_depth, 25.6
# rows: 21000000, HA, num_reads_mapped, 198
# rows: 21000000, HA, perc_cov, 98.876

### LISTS AND DICTIONARIES ###
segment_list = ['HA', 'NA', 'MP', 'NP', 'NS', 'PA', 'PB1', 'PB2']
metric_variables = ['percent_coverage', 'mean_depth', 'mapped_reads', 'seq_len', 'expected_len']
# these will be determined for each segment

#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "--sample_name")
    parser.add_argument( "--bam_stats_csv_list")
    parser.add_argument( "--percent_coverage_csv_list")
    parser.add_argument( "--irma_read_counts")
    options = parser.parse_args(args)
    return options

def create_list_from_string_input(string_input):
    list = string_input.split(' ')
    return list

def create_col_headers(segment_list, metric_variables):
    header_list = ['sample_name']
    for segment in segment_list:
        for metric in metric_variables:
            header_name = "%s_%s" % (segment, metric)
            header_list.append(header_name)
    return header_list


#### MAIN ####
if __name__ == '__main__':

    options = getOptions()
    sample_name = options.sample_name
    bam_stats_files_string_input = options.bam_stats_csv_list
    percent_coverage_csv_file_string_input= options.percent_coverage_csv_list
    read_counts_txt = options.irma_read_counts

     # get a list of file paths
    bam_stats_csv_file_list = create_list_from_string_input(string_input = bam_stats_files_string_input)
    percent_coverage_csv_file_list = create_list_from_string_input(string_input = percent_coverage_csv_file_string_input)

    
    ##### set up the pandas dataframe
    print('Setting up final DF')
    header_list = create_col_headers(segment_list = segment_list, 
                                metric_variables = metric_variables)
    df = pd.DataFrame(columns = header_list)
    df.at[0, 'sample_name'] = sample_name


    #### READS MAPPED (FROM READ_COUNTS.txt)
    print('\n\n\n Pulling reads mapped from READ_COUNTS.txt')
    read_counts_df = pd.read_csv(read_counts_txt, sep = '\t')
    filtered_reads = read_counts_df[read_counts_df.Record == '1-initial'].Reads.iloc[0]
    mapped_reads = read_counts_df[read_counts_df.Record == '3-match'].Reads.iloc[0]
    if '3-altmatch' in read_counts_df.Record.to_list():
        alt_mapped_reads = read_counts_df[read_counts_df.Record == '3-altmatch'].Reads.iloc[0]
    else:
        alt_mapped_reads = 0

        print(f'filtered_reads: {filtered_reads}')
        print(f'alt_mapped_reads: {alt_mapped_reads}')
        print('')
    
    for row in range(read_counts_df.shape[0]):
        record = read_counts_df.Record[row]
        if re.search('4-', record):
            segment = record.split('-')[-1].split('_')[1]
            mapped_reads = read_counts_df[read_counts_df.Record == record ].Reads.iloc[0]

            # add to DF
            col_name = f'{segment}_mapped_reads'
            df.at[0, col_name] = mapped_reads

            print(f'{segment}')
            print(f'{record}')
            print(f'{col_name} : {mapped_reads}')
            print('')


    #### MEAN DEPTH FROM SAMTOOLS (FROM BAM STATS CSV):
    print('\n\nLooping through bam stats csv files')
    print('pulling out mean depth for each segment')
    for bam_stats_csv_file in bam_stats_csv_file_list:
        
        bam_stats_df = pd.read_csv(bam_stats_csv_file, 
                                   dtype = {'sample_name' : object},
                                    na_filter = False )
        segment_name = bam_stats_df.segment_name[0]
        print(f'\n{segment_name}')

        for row in range(bam_stats_df.shape[0]):
            description = bam_stats_df.description[row]
            if description == 'mean_depth':
                value = bam_stats_df.value[row]

                # get correct column header name
                col_name = f"{segment_name}_{description}"
                df.at[0, col_name] = value
                print(col_name, value)

    
    
    ##### PERCENT COVERAGE
    print('\n\nLooping through percent coverage files')
    print('adding up assembled, complete and total percent coverage')
    
    assembled_segments = 0
    complete_segments = 0
    percent_coverage_total = 0
    
    for percent_coverage_csv_file in percent_coverage_csv_file_list:
        assembled_segments += 1

        percent_coverage_df = pd.read_csv(percent_coverage_csv_file,
                                         dtype = {'sample_name' : object},
                                        na_filter = False )
        segment_name = percent_coverage_df.segment[0]
        print(f'\n{segment_name}')

        for row in range(percent_coverage_df.shape[0]):
            description = percent_coverage_df.description[row]
            value = percent_coverage_df.value[row]
            
            # get correct column header name
            col_name = f"{segment_name}_{description}"
            df.at[0, col_name] = value
            print(col_name, value)
            # iterate up if percent coverage = description
            # and percent coverage > 90%
            if description == 'percent_coverage':
                percent_coverage_total = percent_coverage_total + value
                if value >=90:
                    complete_segments +=1


    # add in final columns
    df.at[0, 'assembled_segments'] = assembled_segments
    df.at[0, 'complete_segments'] = complete_segments
    df.at[0, 'filtered_reads'] = filtered_reads # from READ_COUNTS.txt
    df.at[0, 'mapped_reads'] = mapped_reads # from READ_COUNTS.txt
    df.at[0, 'alt_mapped_reads'] = alt_mapped_reads # from READ_COUNTS.txt
    
    df.at[0, 'average_percent_coverage'] = percent_coverage_total/assembled_segments
    df.at[0, 'average_mean_depth'] = mapped_reads/assembled_segments

    # save df
    outfile = f"{sample_name}_assembly_qc_metrics.csv"
    df.to_csv(outfile, index = False)
