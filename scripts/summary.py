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
    parser.add_argument( "--sample_id")
    parser.add_argument( "--preprocess_qc_metrics")
    parser.add_argument( "--irma_typing")
    parser.add_argument( "--irma_assembly_qc_metrics")
    parser.add_argument( "--project_name")
    parser.add_argument( "--analysis_date")
    options = parser.parse_args(args)
    return options

def create_list_from_write_lines_input(write_lines_input):

    list = []
    with open(write_lines_input, 'r') as f:
        for line in f:
            list.append(line.strip())
    return list


#### MAIN ####
if __name__ == '__main__':

    options = getOptions()
    sample_id_txt = options.sample_id
    preprocess_qc_metrics_txt = options.preprocess_qc_metrics
    irma_typing_txt = options.irma_typing
    irma_qc_metrics_txt = options.irma_assembly_qc_metrics
    project_name = options.project_name
    analysis_date = options.analysis_date

    sample_id_list = create_list_from_write_lines_input(write_lines_input = sample_id_txt)
    preprocess_qc_metrics_list = create_list_from_write_lines_input(write_lines_input = preprocess_qc_metrics_txt)
    irma_typing_list = create_list_from_write_lines_input(write_lines_input = irma_typing_txt)
    irma_qc_metrics_list= create_list_from_write_lines_input(write_lines_input = irma_qc_metrics_txt)

    preprocess_qc_metrics_df_list = []
    for preprocess_qc_metrics in preprocess_qc_metrics_list:
        df = pd.read_csv(preprocess_qc_metrics, dtype = {'sample_id' : object})
        preprocess_qc_metrics_df_list.append(df)
    preprocess_qc_metrics_df = pd.concat(preprocess_qc_metrics_df_list).reset_index(drop = True)
    preprocess_qc_metrics_df = preprocess_qc_metrics_df.set_index('sample_id')

    irma_typing_df_list = []
    for irma_typing in irma_typing_list:
        df = pd.read_csv(irma_typing, dtype = {'sample_id' : object})
        irma_typing_df_list.append(df)
    irma_typing_df = pd.concat(irma_typing_df_list).reset_index(drop=True)
    irma_typing_df = irma_typing_df.set_index('sample_id')

    irma_qc_metrics_df_list = []
    first_item = irma_qc_metrics_list[0]
    if first_item != "" and len(irma_qc_metrics_list) != 1:
        for irma_qc_metrics in irma_qc_metrics_list:
            df = pd.read_csv(irma_qc_metrics, dtype = {'sample_id' : object})
            irma_qc_metrics_df_list.append(df)
        irma_qc_metrics_df = pd.concat(irma_qc_metrics_df_list).reset_index(drop = True)
        irma_qc_metrics_df = irma_qc_metrics_df.set_index('sample_id')
    else:
        # somehow need to create a fake sample so that the headers get included:
        adict = {'sample_id' : ['dummy'], 
        'total_segments' : [0], 
        'total_flu_mapped_reads':[0], 
        'HA_per_cov' :[0],
        'HA_mean_depth':[0], 
        'HA_num_mapped_reads':[0], 
        'HA_seq_len':[0], 
        'HA_expected_len':[0],
        'NA_per_cov':[0], 
        'NA_mean_depth':[0], 
        'NA_num_mapped_reads':[0], 
        'NA_seq_len':[0],
        'NA_expected_len':[0], 
        'MP_per_cov':[0], 
        'MP_mean_depth':[0], 
        'MP_num_mapped_reads':[0],
        'MP_seq_len':[0], 
        'MP_expected_len':[0], 
        'NP_per_cov':[0], 
        'NP_mean_depth':[0],
        'NP_num_mapped_reads':[0], 
        'NP_seq_len':[0], 
        'NP_expected_len':[0], 
        'NS_per_cov':[0],
        'NS_mean_depth':[0], 
        'NS_num_mapped_reads':[0], 
        'NS_seq_len':[0], 
        'NS_expected_len':[0],
        'PA_per_cov':[0], 
        'PA_mean_depth':[0], 
        'PA_num_mapped_reads':[0], 
        'PA_seq_len':[0],
        'PA_expected_len':[0], 
        'PB1_per_cov':[0], 
        'PB1_mean_depth':[0],
        'PB1_num_mapped_reads':[0], 
        'PB1_seq_len':[0], 
        'PB1_expected_len':[0],
        'PB2_per_cov':[0], 
        'PB2_mean_depth':[0], 
        'PB2_num_mapped_reads':[0], 
        'PB2_seq_len':[0],
        'PB2_expected_len':[0], 
        'ivar_version':[0], 
        'ivar_docker':[0], 
        'ivar_min_depth':[0],
        'ivar_min_freq':[0], 
        'ivar_min_qual':[0]}

        irma_qc_metrics_df = pd.DataFrame(adict)



    # join
    df = preprocess_qc_metrics_df.join(irma_typing_df, how = 'outer')
    df = df.join(irma_qc_metrics_df,how = 'outer')
    df = df.reset_index()
    df['percent_flu_mapped_reads'] = round((df.total_flu_mapped_reads / df.read_pairs_cleaned) * 100 , 2)

    # order columns
    col_order = ['sample_id', 'type', 'HA_subytpe', 'NA_subtype',
    'total_segments', 'total_flu_mapped_reads', 'percent_flu_mapped_reads',
    'read_pairs_cleaned',
    'HA_per_cov','HA_mean_depth', 'HA_num_mapped_reads', 'HA_seq_len', 'HA_expected_len',
    'NA_per_cov', 'NA_mean_depth', 'NA_num_mapped_reads', 'NA_seq_len',
    'NA_expected_len', 'MP_per_cov', 'MP_mean_depth', 'MP_num_mapped_reads',
    'MP_seq_len', 'MP_expected_len', 'NP_per_cov', 'NP_mean_depth',
    'NP_num_mapped_reads', 'NP_seq_len', 'NP_expected_len', 'NS_per_cov',
    'NS_mean_depth', 'NS_num_mapped_reads', 'NS_seq_len', 'NS_expected_len',
    'PA_per_cov', 'PA_mean_depth', 'PA_num_mapped_reads', 'PA_seq_len',
    'PA_expected_len', 'PB1_per_cov', 'PB1_mean_depth',
    'PB1_num_mapped_reads', 'PB1_seq_len', 'PB1_expected_len',
    'PB2_per_cov', 'PB2_mean_depth', 'PB2_num_mapped_reads', 'PB2_seq_len',
    'PB2_expected_len',
    'read_type', 'total_read_diff', 'read_length_R1_raw',
    'read_length_R2_raw', 'total_reads_R1_raw', 'total_reads_R2_raw',
    'read_pairs_raw', 'read_length_R1_cleaned', 'read_length_R2_cleaned',
    'total_reads_R1_cleaned', 'total_reads_R2_cleaned',
    'read_pairs_cleaned', 
    'fastqc_version', 'fastqc_docker',
    'seqyclean_version', 'seqyclean_docker',
    'irma_version', 'irma_docker', 'irma_module', 
    'ivar_version', 'ivar_docker', 'ivar_min_depth', 'ivar_min_freq', 'ivar_min_qual']

    df = df[col_order]
    df["analysis_date"] = analysis_date

    #drop dummy sample if exists:
    df = df[df.sample_id != "dummy"].reset_index(drop = True)

    #outfile
    outfile = '%s_sequencing_results.csv' % project_name
    df.to_csv(outfile, index = False)












