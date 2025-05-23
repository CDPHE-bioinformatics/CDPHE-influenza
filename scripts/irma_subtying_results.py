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
    options = parser.parse_args(args)
    return options

#### MAIN ####
if __name__ == '__main__':

    options = getOptions()
    irma_assembled_gene_segments_csv = options.irma_assembled_gene_segments_csv
    sample_name = options.sample_name


    df = pd.read_csv(irma_assembled_gene_segments_csv, 
                     dtype = {'sample_name' : object},
                     na_filter = False)
    # use na_filter so that NA segment is not interpreted as Null value

    # check for mixed types
    TYPES = df.flu_type.unique().tolist()
    if len(TYPES) > 1 :
        TYPE = 'mixed'
    else:
        TYPE = TYPES[0]

    print(TYPE)

    # pull out subtypes
    subtyping_dict = dict(zip(df.gene_segment, df.subtype))
    if "HA" in subtyping_dict.keys():
        HA_subtype = subtyping_dict["HA"]
    else:
        HA_subtype = ""
    
    if "NA" in subtyping_dict.keys():
        NA_subtype = subtyping_dict["NA"]
    else:
        NA_subtype = ""

    
    # put all the info together

    summary_df = pd.DataFrame()
    summary_df['sample_name'] = [sample_name]
    summary_df['flu_type'] = [TYPE]
    summary_df['HA_subtype'] = [HA_subtype]
    summary_df['NA_subtype'] = [NA_subtype]

    outfile = f'{sample_name}_irma_typing.csv'
    summary_df.to_csv(outfile, index = False)

    print(TYPE, HA_subtype, NA_subtype)

    print(summary_df)
    
    #write some single line outfiles
    type_outfile = 'TYPE.txt'
    with open(type_outfile, 'w') as f:
        f.writelines(f'{TYPE}')

    HA_subtype_outfile = 'HA_SUBTYPE.txt'
    with open(HA_subtype_outfile, 'w') as f:
        f.writelines(f'{HA_subtype}')

    NA_subtype_outfile = 'NA_SUBTYPE.txt'
    with open(NA_subtype_outfile, 'w') as f:
        f.writelines(f'{NA_subtype}')



    


