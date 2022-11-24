#! /usr/bin/env python

import pandas as pd
import re
from Bio import SeqIO
import argparse
import sys


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")

    parser.add_argument('--sample_name')
    parser.add_argument('--fasta_files_array')
    parser.add_argument('--num_fastas')
    parser.add_argument('--tables_directory')

    parser.add_argument('--total_reads_R1_raw')
    parser.add_argument('--total_reads_R2_raw')
    parser.add_argument('--read_length_R1_raw')
    parser.add_argument('--read_length_R2_raw')
    parser.add_argument('--read_pairs_raw')

    parser.add_argument('--total_reads_R1_cleaned')
    parser.add_argument('--total_reads_R2_cleaned')
    parser.add_argument('--read_length_R1_cleaned')
    parser.add_argument('--read_length_R2_cleaned')
    parser.add_argument('--read_pairs_cleaned')

    parser.add_argument('--fastqc_version')
    parser.add_argument('--fastqc_docker')
    parser.add_argument('--seqyclean_version')
    parser.add_argument('--seqyclean_docker')
    parser.add_argument('--irma_version')

    options = parser.parse_args(args)
    return options

def wdl_array_to_python_list (txt_file):
    python_list = []
    with open(txt_file, 'r') as f:
        for line in f:
            python_list.append(line.strip())
    return python_list

def get_subtype(fasta_files_list, flu_type):
    HA_subtype = ''
    NA_subtype = ''

    if flu_type == 'A':
        for fasta in fasta_files_list:
            fasta_name = fasta_files_list[0].split('/')[-1]
            if re.search('HA', fasta_name):
                HA_subtype = fasta_name.split('_')[-1]
            elif re.search('NA', fasta_name):
                NA_subtype = fasta_name.split('_')[-1]

    subtype = "%s%s" % (HA_subtype, NA_subtype)

    return subtype

def get_flu_type(fasta_files_list):
    flu_type = ''
    first_fasta_name = fasta_files_list[0].split('/')[-1]
    flu_type = first_fasta_name.split("_")[1]

    return flu_type


def calc_coverate(fasta_file, gene_segment):
    # gene_segment lenght dictionary
    ref_len_dict = {'A_MP': 982,'A_NP': 1497,'A_NS': 863,'A_PA': 2151,'A_PB1': 2274,'A_PB2': 2280,
    'A_HA_H1': 1704,'A_HA_H10': 1686,'A_HA_H11': 1698,'A_HA_H12': 1695,'A_HA_H13': 1701,'A_HA_H14': 1707,
    'A_HA_H15': 1713,'A_HA_H16': 1698,'A_HA_H2': 1689,'A_HA_H3': 1704,'A_HA_H4': 1695,'A_HA_H5': 1707,
    'A_HA_H6': 1704,'A_HA_H7': 1713,'A_HA_H8': 1701,'A_HA_H9': 1683,'A_NA_N1': 1413,'A_NA_N2': 1410,
    'A_NA_N3': 1410,'A_NA_N4': 1413,'A_NA_N5': 1422,'A_NA_N6': 1413,'A_NA_N7': 1416,
    'A_NA_N8': 1413,'A_NA_N9': 1413,'B_HA': 1758,'B_MP': 1139,'B_NA': 1408,'B_NP': 1683,
    'B_NS': 1034,'B_PA': 2181,'B_PB1': 2263,'B_PB2': 2313}

    record = SeqIO.read(fasta_file, 'fasta')
    consensus_seq_len = len(record.seq)
    expected_len = ref_len_dict[gene_segment]
    perc_cov = (seq_len/expected_len) * 100

    return {'consensus_seq_len':consesnus_seq_len,
            'expected_len':expected_len,
            'perc_cov':perc_cov}


def calc_av_depth(coverage_txt):
    cov_df = pd.read_csv(coverage_txt, sep = '\t')
    cov_df = cov_df.rename(columns = {'Coverage Depth': 'coverage_depth'})
    mean_depth = cov_df.coverage_depth.mean()

    return mean_depth

def get_gene_segment (fasta_file):
    fasta_name = fasta_files_list[0].split('/')[-1]
    gene_segment = fasta_name.split('_')[2]
    type = fasta_name.split('_')[1]
    combo = '%s_%s' % (type, gene_segment)

    return {'gene_segment' : gene_segment,'type' : type,'combo' : combo}

if __name__ == '__main__':

    options = getOptions()

    # note fasta_file is formated as /sample_name_A_HA_H1.fasta
    fasta_files_list = wdl_array_to_python_list(txt_file = options.fasta_files_array)
    tables_directory_list = wdl_array_to_python_list(txt_file = options.tables_directory)

    # set up dataframe
    col_headers = ['sample_id','flu_type','subtype',
    'HA_expected_len','HA_seq_len','HA_mean_depth','HA_percent_coverage',
    'NA_expected_len','NA_seq_len','NA_mean_depth','NA_percent_coverage',
    'MP_expected_len','MP_seq_len','MP_mean_depth','MP_percent_coverage',
    'NP_expected_len','NP_seq_len','NP_mean_depth','NP_percent_coverage',
    'NS_expected_len','NS_seq_len','NS_mean_depth','NS_percent_coverage',
    'PA_expected_len','PA_seq_len','PA_mean_depth','PA_percent_coverage',
    'PB1_expected_len','PB1_seq_len','PB1_mean_depth','PB1_percent_coverage',
    'PB2_expected_len','PB2_seq_len','PB2_mean_depth','PB2_percent_coverage']

    df = pd.DataFrame(columns = col_headers)
    df.at[0, 'sample_id'] = options.sample_name

    if options.num_fastas == 0:
        for col in df.columns:
            if col not in ['sample_id']:
                df.at[0, col] = ''

    else:
        # type
        flu_type = get_flu_type(fasta_files_list = fasta_files_list)
        df.at[0, 'flu_type'] = flu_type

        #subtype
        subtype = get_subtype(fasta_files_list = fasta_files_list, flu_type = flu_type)
        df.at[0, 'subtype'] = subtype

        # get coverage
        for fasta in fasta_files_list:
            seg_dict = get_gene_segment(fasta_file = fasta_file)
            gene_segment = seg_dict['gene_segment']
            gene_segment_combo= seg_dict['combo']

            coverage_dict = calc_coverate(fasta_file =fasta_file, gene_segment= gene_segment_combo)
            perc_cov = coverage_dict['perc_cov']
            expected_len = coverage_dict['expected_len']
            consensus_seq_len = coverage_dict['consensus_seq_len']

            col1 = '%s_seq_len' % gene_segment
            col2 = '%s_expected_len' % gene_segment

            df.at[0, col1] = consensus_seq_len
            df.at[0, col2] = expected_len

        # get depth
        for file in tables_directory_list:
            file_name = file.split('/')[-1]
            gene_segment = file_name.split('_')[2]
            if re.serach('coverage',file_name):
                depth= calc_av_depth(coverage_txt = file)

                col1 = '%s_mean_depth' % gene_segment
                df.at[0, col1] = depth

    # fill in the rest of the table
    df.at[0, 'num_gene_segments'] = options.num_fastas

    df.at[0, 'total_reads_R1_raw'] = options.total_reads_R1_raw
    df.at[0, 'total_reads_R2_raw'] = options.total_reads_R2_raw
    df.at[0, 'read_length_R1_raw'] = options.read_length_R1_raw
    df.at[0, 'read_length_R2_raw'] = options.read_length_R2_raw
    df.at[0, 'read_pairs_raw'] = options.read_pairs_raw

    df.at[0, 'total_reads_R1_cleaned'] = options.total_reads_R1_cleaned
    df.at[0, 'total_reads_R2_cleaned'] = options.total_reads_R2_cleaned
    df.at[0, 'read_length_R1_cleaned'] = options.read_length_R1_cleaned
    df.at[0, 'read_length_R2_cleaned'] = options.read_length_R2_cleaned
    df.at[0, 'read_pairs_cleaned'] = options.read_pairs_cleaned

    df.at[0, 'fastqc_version'] = options.fastqc_version
    df.at[0, 'fastqc_docker'] = options.fastqc_docker
    df.at[0, 'seqyclean_version'] = options.seqyclean_version
    df.at[0, 'seqyclean_docker'] = options.seqyclean_docker
    df.at[0, 'irma_version'] = options.irma_version

    # write outfile
    outfile_name = '%s_assembly_metrics_reordered.csv' % options.sample_name
    df.to_csv(outfile_name, index = False)

    # change the col_header order to have options for viewing
    col_order = ['sample_id','num_gene_segments','flu_type','subtype',
    'HA_percent_coverage', 'NA_percent_coverage', 'MP_percent_coverage',
    'NP_percent_coverage', 'NS_percent_coverage', 'PA_percent_coverage',
    'PB1_percent_coverage', 'PB2_percent_coverage',
    'HA_mean_depth', 'NA_mean_depth', 'MP_mean_depth',
    'NP_mean_depth', 'NS_mean_depth', 'PA_mean_depth',
    'PB1_mean_depth', 'PB2_mean_depth',
    'HA_seq_len', 'NA_seq_len', 'MP_seq_len',
    'NP_seq_len', 'NS_seq_len', 'PA_seq_len',
    'PB1_seq_len', 'PB2_seq_len',
    'num_gene_segments',
    'total_reads_R1_raw', 'total_reads_R2_raw', 'read_length_R1_raw', 'read_length_R2_raw', 'read_pairs_raw',
    'total_reads_R1_cleaned', 'total_reads_R2_cleaned', 'read_length_R1_cleaned', 'read_length_R2_cleaned','read_pairs_cleaned',
    'fastqc_version', 'fastqc_docker', 'seqyclean_version', 'seqyclean_docker', 'irma_version']

    df2 = df[col_order]
    outfile_name = '%s_assembly_metrics.csv' % options.sample_name
    df2.to_csv(outfile_name, index = False)
