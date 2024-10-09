#! /usr/bin/env python

# import python modules
import pandas as pd
import sys
import argparse
import subprocess
import re


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "--sample_name")
    parser.add_argument( "--preprocess_qc_metrics")
    parser.add_argument( "--irma_typing")
    parser.add_argument( "--assembly_qc_metrics")
    parser.add_argument( "--nextclade_tsv")
    parser.add_argument('--workflow_version')
    parser.add_argument( "--analysis_date")
    parser.add_argument( "--project_name")

    options = parser.parse_args(args)
    return options

def create_list_from_write_lines_input(write_lines_input):

    list = []
    with open(write_lines_input, 'r') as f:
        for line in f:
            line = line.strip()
            if line: 
                # ensure that only lines with values are added to the list
                list.append(line.strip())
    return list



#### MAIN ####
if __name__ == '__main__':

    options = getOptions()
    sample_name_txt = options.sample_name
    preprocess_qc_metrics_txt = options.preprocess_qc_metrics
    irma_typing_txt = options.irma_typing
    post_qc_metrics_txt = options.assembly_qc_metrics
    nextclade_tsv_txt = options.nextclade_tsv
    workflow_version = options.workflow_version
    project_name = options.project_name
    analysis_date = options.analysis_date

    sample_name_list = create_list_from_write_lines_input(write_lines_input = sample_name_txt)
    preprocess_qc_metrics_list = create_list_from_write_lines_input(write_lines_input = preprocess_qc_metrics_txt)
    irma_typing_list = create_list_from_write_lines_input(write_lines_input = irma_typing_txt)
    post_qc_metrics_list= create_list_from_write_lines_input(write_lines_input = post_qc_metrics_txt)
    nextclade_tsv_list = create_list_from_write_lines_input(write_lines_input = nextclade_tsv_txt)
    
   

    # create summary metrics file:
    # the issue I have is that if assemlby fails then the the post assembly
    # metrics won't be generated. Plus if the HA or NA segment fail then
    # we won't have the nextclade output. 

    # dropped seq_len and expected length from columns
    # nextclade columns to drop - totalMissing, totalNonACGTNs, totalUnknownAa,
    # aaSubstituttions, aaDeletions, aaInsertions, warning, errors
    # keep others for QC purposes.

    # set columns
    col_order = [
        'sample_name', 'project_name', 'analysis_date', 'flu_type', 
        'HA_subtype', 'NA_subtype',
        'HA_clade', 'HA_short-clade', 'HA_subclade',
        'NA_clade', 
        'complete_segments', 'assembled_segments', 
        'total_flu_mapped_reads', 'percent_flu_mapped_reads','total_reads_cleaned',
        'HA_percent_coverage','HA_mean_depth', 'HA_num_mapped_reads', 
        'NA_percent_coverage', 'NA_mean_depth', 'NA_num_mapped_reads', 
        'MP_percent_coverage', 'MP_mean_depth', 'MP_num_mapped_reads',
        'NP_percent_coverage', 'NP_mean_depth','NP_num_mapped_reads',  
        'NS_percent_coverage','NS_mean_depth', 'NS_num_mapped_reads', 
        'PA_percent_coverage', 'PA_mean_depth', 'PA_num_mapped_reads', 
        'PB1_percent_coverage', 'PB1_mean_depth','PB1_num_mapped_reads', 
        'PB2_percent_coverage', 'PB2_mean_depth', 'PB2_num_mapped_reads', 
        'total_read_diff',  'total_reads_R1_raw', 'total_reads_R2_raw', 'read_pairs_raw', 'total_reads_raw', 
        'read_length_R1_raw', 'read_length_R2_raw',
        'total_reads_R1_cleaned', 'total_reads_R2_cleaned', 'read_pairs_cleaned', 
        'read_length_R1_cleaned', 'read_length_R2_cleaned',
        'HA_totalSubstitutions', 'HA_totalDeletions', 'HA_totalInsertions', 'HA_totalFrameShifts',  
        'HA_totalAminoacidSubstitutions', 'HA_totalAminoacidDeletions', 'HA_totalAminoacidInsertions', 
        'NA_totalSubstitutions', 'NA_totalDeletions', 'NA_totalInsertions', 'NA_totalFrameShifts', 
        'NA_totalAminoacidSubstitutions', 'NA_totalAminoacidDeletions', 'NA_totalAminoacidInsertions'
        ]

    # preprocess
    preprocess_qc_metrics_df_list = []
    for preprocess_qc_metrics in preprocess_qc_metrics_list:
        df = pd.read_csv(preprocess_qc_metrics, dtype = {'sample_name' : object})
        # print()
        # print(df.columns)
        # print()
        preprocess_qc_metrics_df_list.append(df)
    preprocess_qc_metrics_df = pd.concat(preprocess_qc_metrics_df_list).set_index('sample_name')

    # irma subtyping
    irma_typing_df_list = []
    for irma_typing in irma_typing_list:
        df = pd.read_csv(irma_typing, dtype = {'sample_name' : object})
        irma_typing_df_list.append(df)
    irma_typing_df = pd.concat(irma_typing_df_list).set_index('sample_name')


    # nextclade
    na_nextclade_df_list = []
    ha_nextclade_df_list = []

    # check that files exist
    if len(nextclade_tsv_list) >= 1:
        for nextclade_tsv in nextclade_tsv_list:
            df = pd.read_csv(nextclade_tsv, sep ='\t')
            # determine if HA or NA segment using the file name
            # example file name "sample_A_HA-H1.tsv" or "sample_B_HA.tsv"
            segment = nextclade_tsv.split('.tsv')[0].split('_')[-1] # this gives you HA-H1 or HA depednig if A or B
            type = nextclade_tsv.split('.tsv')[0].split('_')[-2] # this give A or B
            
            if re.search("NA", segment):
                # pull sample_name from seqName
                seq_name = df.seqName[0]
                print('NA - seq_name: ', seq_name)
                sample_name = seq_name.split(f'{type}_{segment}')[0].rstrip('_')
                print('NA - sample_name: ', sample_name)

                # add and rename columns
                df['sample_name'] = sample_name
                df['nextclade_coverage'] = df['coverage']

                # reorder columns
                col_keep = ['sample_name', 'clade', 
                            'totalSubstitutions','totalDeletions', 'totalInsertions', 
                            'totalFrameShifts', 'totalMissing','totalNonACGTNs', 'totalAminoacidSubstitutions',
                            'totalAminoacidDeletions', 'totalAminoacidInsertions', 'totalUnknownAa', 
                            'nextclade_coverage','aaSubstitutions', 'aaDeletions', 'aaInsertions',
                            'warnings', 'errors']
                # add "na" prefix to all column headers
                rename_cols = {}
                for col in col_keep:
                    if col != 'sample_name':
                        new_column = f'NA_{col}'
                        rename_cols[col] = new_column

                df = df[col_keep]
                df = df.rename(columns = rename_cols)
                na_nextclade_df_list.append(df)

            elif re.search('HA', segment):
                # pull sample_name from seqName
                seq_name = df.seqName[0]
                print('HA - seq_name: ', seq_name)
                sample_name = seq_name.split(f'{type}_{segment}')[0].rstrip('_')
                print('HA - sample_name: ', sample_name)

                # add column and rename columns
                df['sample_name'] = sample_name
                df['nextclade_coverage'] = df['coverage']

                # add missing columns
                # HA: add short-clade (For bvic only)
                if "subclade" not in df.columns:
                    df['subclade'] = ""
                if "short-clade" not in df.columns:
                    df['short-clade'] = ""
                # reorder columns
                col_keep = ['sample_name', 'clade', 'short-clade', 'subclade', 
                            'totalSubstitutions','totalDeletions', 'totalInsertions', 
                            'totalFrameShifts', 'totalMissing','totalNonACGTNs', 'totalAminoacidSubstitutions',
                            'totalAminoacidDeletions', 'totalAminoacidInsertions', 'totalUnknownAa', 
                            'nextclade_coverage','aaSubstitutions', 'aaDeletions', 'aaInsertions',
                            'warnings', 'errors']
                # add "ha" prefix to all column headers
                rename_cols = {}
                for col in col_keep:
                    if col != 'sample_name':
                        new_column = f'HA_{col}'
                        rename_cols[col] = new_column

                df = df[col_keep]
                df = df.rename(columns = rename_cols)
                ha_nextclade_df_list.append(df)
    
    if len(ha_nextclade_df_list) > 0:
        ha_nextclade_df = pd.concat(ha_nextclade_df_list).set_index('sample_name')
    else:
        ha_nextclade_df = pd.DataFrame()
    if len(na_nextclade_df_list) > 0:
        na_nextclade_df = pd.concat(na_nextclade_df_list).set_index('sample_name')
    else:
        na_nextclade_df = pd.DataFrame()



    # post assembly qc metrics
    post_qc_metrics_df_list = []
    if len(post_qc_metrics_list) >= 1:
        for post_qc_metrics in post_qc_metrics_list:
            df = pd.read_csv(post_qc_metrics, dtype = {'sample_name' : object})
            # print()
            # print(df.columns)
            # print()
            post_qc_metrics_df_list.append(df)
        post_qc_metrics_df = pd.concat(post_qc_metrics_df_list).set_index('sample_name')
    else:
        # create empty df for joining downstream
        post_qc_metrics_df = pd.DataFrame()

    
    # join all the dfs together
    df = pd.DataFrame(sample_name_list, columns = ['sample_name']).set_index('sample_name')
    df = df.join(preprocess_qc_metrics_df, how = 'left')
    df = df.join(irma_typing_df, how = 'left')
    df = df.join(na_nextclade_df, how = 'left')
    df = df.join(ha_nextclade_df, how = 'left')
    df = df.join(post_qc_metrics_df, how = 'left')
    df = df.reset_index(drop = False)

    # add some columns and do a calcuation
    df["analysis_date"] = analysis_date
    df['percent_flu_mapped_reads'] = round((df.total_flu_mapped_reads / df.total_reads_cleaned) * 100 , 2)
    # TODO the percent flu mapped reads we are seeing is much lower than CDC's. 
    # I think they are somehow dividing by the total flu "classified" reads
    df['project_name'] = project_name
    df = df[col_order] 

    # check columns - if column doesn't exist then add column
    # for if assembly failed for all samples assemlby qc metrics or nextclade
    # df weren't every made; then this keeps the columns consisent regardless
    for column in col_order:
        if column not in df.columns:
            df[column] = pd.NA
   
    
    # outfile
    outfile = f'{project_name}_sequencing_results_{workflow_version}.csv' 
    df.to_csv(outfile, index = False)


    













