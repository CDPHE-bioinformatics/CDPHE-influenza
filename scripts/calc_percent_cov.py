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


# calculate percent coverage for each gene segment
# read in single fasta file
# determine gene segment from file name
# calculates percent coverage
# output a single file with fasta file name and percent coverage

### Segement length dictionary
ref_len_dict = {'A_MP': 982,'A_NP': 1497,'A_NS': 863,'A_PA': 2151,'A_PB1': 2274,'A_PB2': 2280,
    'A_HA_H1': 1704,'A_HA_H10': 1686,'A_HA_H11': 1698,'A_HA_H12': 1695,'A_HA_H13': 1701,'A_HA_H14': 1707,
    'A_HA_H15': 1713,'A_HA_H16': 1698,'A_HA_H2': 1689,'A_HA_H3': 1704,'A_HA_H4': 1695,'A_HA_H5': 1707,
    'A_HA_H6': 1704,'A_HA_H7': 1713,'A_HA_H8': 1701,'A_HA_H9': 1683,'A_NA_N1': 1413,'A_NA_N2': 1410,
    'A_NA_N3': 1410,'A_NA_N4': 1413,'A_NA_N5': 1422,'A_NA_N6': 1413,'A_NA_N7': 1416,
    'A_NA_N8': 1413,'A_NA_N9': 1413,'B_HA': 1758,'B_MP': 1139,'B_NA': 1408,'B_NP': 1683,
    'B_NS': 1034,'B_PA': 2181,'B_PB1': 2263,'B_PB2': 2313}

#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "--fasta_file")
    parser.add_argument('--sample_name')
    options = parser.parse_args(args)
    return options


def get_fasta_file_basename(fasta_file_path):
    basename = fasta_file_path.split('/')[-1] # strip directories
    return basename

def get_segment_name(fasta_file_path, sample_name):
    basename = fasta_file_path.split('/')[-1] # strip directories
    segment_name = '_'.join(basename.split(sample_name)[1].split('_')[1:]).split('.')[0] # A_HA_H3, A_MP etc.

    return segment_name

def get_gene_name(fasta_file_path, sample_name):
    basename = fasta_file_path.split('/')[-1] # strip directories
    segment_name = '_'.join(basename.split(sample_name)[1].split('_')[1:]).split('.')[0]
    gene_name = segment_name.split('_')[1] #HA, MP etc
    
    return gene_name

# def get_sample_name(fasta_file_path):
#     basename = fasta_file_path.split('/')[-1] # strip directories
#     sample_name = basename.split('_')[0] # pull out sample id

#     return sample_name

def get_seq_len(fasta_file_path):
    # read in fasta file
    record = SeqIO.read(fasta_file_path, 'fasta')

    # get length of non ambigous bases
    seq = record.seq
    seq_len = (seq.count('A') + seq.count('C') + seq.count('G') + seq.count('T'))

    return seq_len

def calc_percent_cov(seq_len, ref_len_dict, segment_name):
    # calcuat per cov based on expected ref length
    expected_length = ref_len_dict[segment_name]
    per_cov = round(((seq_len/expected_length)*100), 2)

    return per_cov

def create_output(sample_name, basename, segment_name, gene_name, seq_len, per_cov, expected_len):

    df = pd.DataFrame()
    description_list = ['expected_len', 'seq_len', 'per_cov']
    value_list = [expected_len, seq_len, per_cov]
    df['description'] = description_list
    df['value'] = value_list
    df['sample_name'] = sample_name
    df['base_name'] = basename
    df['segment_name'] = segment_name
    df['gene_name'] = gene_name
    
    col_order = ['sample_name', 'base_name', 'segment_name', 'gene_name', 'description', 'value']
    df = df[col_order]

    outfile='perc_cov_results.csv' 
    df.to_csv(outfile, index = False)


#### MAIN ####
if __name__ == '__main__':

    options = getOptions()
    fasta_file_path = options.fasta_file
    sample_name = options.sample_name

    basename = get_fasta_file_basename(fasta_file_path = fasta_file_path)

    segment_name = get_segment_name(fasta_file_path=fasta_file_path, sample_name = sample_name)
    gene_name = get_gene_name(fasta_file_path=fasta_file_path, sample_name = sample_name)

    seq_len = get_seq_len(fasta_file_path=fasta_file_path)
    per_cov = calc_percent_cov(seq_len=seq_len, 
                                ref_len_dict=ref_len_dict, 
                                segment_name=segment_name)

    expected_len = ref_len_dict[segment_name]

    create_output(sample_name = sample_name, 
                    basename = basename, 
                    segment_name=segment_name, 
                    gene_name=gene_name,
                    seq_len=seq_len, 
                    per_cov=per_cov, 
                    expected_len=expected_len)

    
    
    











