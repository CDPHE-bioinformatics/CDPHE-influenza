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

#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "--irma_assembled_gene_segments_csv")
    parser.add_argument( "--sample_name")
    parser.add_argument( "--irma_version")
    parser.add_argument( "--irma_docker")
    parser.add_argument( "--irma_module")
    options = parser.parse_args(args)
    return options

#### MAIN ####
if __name__ == '__main__':

    options = getOptions()
    irma_assembled_gene_segments_csv = options.irma_assembled_gene_segments_csv
    sample_name = options.sample_name
    irma_version = options.irma_version
    irma_docker = options.irma_docker
    irma_module = options.irma_module


    df = pd.read_csv(irma_assembled_gene_segments_csv, dtype = {'sample_name' : object})
    df = df.fillna('none')

    # check for mixed types
    TYPES = df.TYPE.unique().tolist()
    if len(TYPES) > 1 :
        TYPE = 'mixed'
    else:
        TYPE = TYPES[0]

    # pull out subtypes
    HA_subtype = 'none'
    NA_subtype = 'none'

    for row in range(df.shape[0]):
        gene_segment = df.gene_segement[row] 
        if gene_segment == 'HA':
            HA_subtype = df.subtype[row]
        elif gene_segment == 'NA':
            NA_subtype = df.subtype[row]

    # put all the info together

    summary_df = pd.DataFrame()
    summary_df['sample_name'] = sample_name
    summary_df['type'] = TYPE
    summary_df['HA_subtype'] = HA_subtype
    summary_df['NA_sbutype'] = NA_subtype
    summary_df['irma_module'] = irma_module
    summary_df['irma_docker'] = irma_docker
    summary_df['irma_version'] = irma_version

    outfile = f'{sample_name}_irma_typing.csv'
    summary_df.to_csv(outfile, index = False)

    #write some single line outfiles
    type_outfile = 'TYPE.txt'
    with open(type_outfile, 'w') as f:
        f.writelines(f'{TYPE}')

    HA_subtype_outfile = 'HA_SUBTYPE.txt'
    with open(HA_subtype_outfile, 'w') as f:
        f.writelines(f'HA_subtype')

    NA_subtype_outfile = 'NA_SUBTYPE.txt'
    with open(NA_subtype_outfile, 'w') as f:
        f.writelines(f'NA_subtype')



    


