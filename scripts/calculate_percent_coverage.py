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
    'A_HA-H1': 1704,'A_HA-H10': 1686,'A_HA-H11': 1698,'A_HA-H12': 1695,'A_HA-H13': 1701,'A_HA-H14': 1707,
    'A_HA-H15': 1713,'A_HA-H16': 1698,'A_HA-H2': 1689,'A_HA-H3': 1704,'A_HA-H4': 1695,'A_HA-H5': 1707,
    'A_HA-H6': 1704,'A_HA-H7': 1713,'A_HA-H8': 1701,'A_HA-H9': 1683,'A_NA-N1': 1413,'A_NA-N2': 1410,
    'A_NA-N3': 1410,'A_NA-N4': 1413,'A_NA-N5': 1422,'A_NA-N6': 1413,'A_NA-N7': 1416,
    'A_NA-N8': 1413,'A_NA-N9': 1413,'B_HA': 1758,'B_MP': 1139,'B_NA': 1408,'B_NP': 1683,
    'B_NS': 1034,'B_PA': 2181,'B_PB1': 2263,'B_PB2': 2313}

#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "--fasta_file")
    parser.add_argument('--sample_name')
    parser.add_argument('--base_name') #{sample_name}_A_HA-H1 or {sample_name}_A_NP
    parser.add_argument('--segment') # HA, NP
    options = parser.parse_args(args)
    return options


def get_ref_seg(base_name):
    print('base_name', base_name)
    #remove sample name from base_name
    ref_seg = base_name.replace(f'{sample_name}_', '')
    print('ref_seg', ref_seg)
    # ref seg will have type, segment and subytpe in it,
    # to match the dictionary of seq_lens
    return ref_seg

def get_seq_len(fasta_file_path):
    # read in fasta file
    record = SeqIO.read(fasta_file_path, 'fasta')

    # get length of non ambigous bases
    seq = record.seq
    seq_len = (seq.count('A') + seq.count('C') + seq.count('G') + seq.count('T'))
    return seq_len

def calc_percent_cov(seq_len, ref_len_dict, ref_seg):
    # calcuat per cov based on expected ref length
    expected_length = ref_len_dict[ref_seg]
    per_cov = round(((seq_len/expected_length)*100), 2)

    return per_cov

def create_output(sample_name, segment, seq_len, per_cov, expected_len):

    df = pd.DataFrame()
    description_list = ['expected_len', 'seq_len', 'percent_coverage']
    value_list = [expected_len, seq_len, per_cov]
    df['description'] = description_list
    df['value'] = value_list
    df['sample_name'] = sample_name
    df['segment'] = segment
    
    col_order = ['sample_name', 'segment', 'description', 'value']
    df = df[col_order]

    outfile=f'percent_coverage_results.csv' 
    df.to_csv(outfile, index = False)


#### MAIN ####
if __name__ == '__main__':

    print('starting script')

    options = getOptions()
    fasta_file_path = options.fasta_file
    sample_name = options.sample_name
    base_name = options.base_name
    segment = options.segment

    ref_seg = get_ref_seg(base_name = base_name)
    seq_len = get_seq_len(fasta_file_path=fasta_file_path)
    per_cov = calc_percent_cov(seq_len=seq_len, 
                                ref_len_dict=ref_len_dict, 
                                ref_seg = ref_seg)

    expected_len = ref_len_dict[ref_seg]

    create_output(sample_name = sample_name, 
                    segment = segment,
                    seq_len=seq_len, 
                    per_cov=per_cov, 
                    expected_len=expected_len)

    
    
    











