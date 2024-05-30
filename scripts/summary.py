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
    parser.add_argument( "--sample_name")
    parser.add_argument( "--preprocess_qc_metrics")
    parser.add_argument( "--irma_typing")
    parser.add_argument( "--assembly_qc_metrics")
    parser.add_argument( "--nextclade_na_tsv")
    parser.add_argument( "--nextclade_ha_tsv")
    parser.add_argument( "--version_catpure_file")
    parser.add_argument( "--workflow_version")
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

# create version capture file for assembly
def create_version_capture_file_for_assembly_workflow(
        version_capture_file_list,
        project_name,
        workflow_version) :
    
    '''
    this function concatenates and removes duplicates of the version capture
    of the illumina_pe_assembly workflow. Have to do this in the summary workflow
    because not all samples will make it through assembly and if the sample faisl, or
    if there not an ha/na segment then some version capture info will be missing
    for some samples. And the version capture would overwrite each sample in the
    assembly workflow. 
    '''
    
    version_capture_df_list = []
    for version_capture_file in version_capture_file_list:
        df = pd.read_csv(version_capture_file, dtype = {'sample_name' : object})
        version_capture_df_list.append(df)
    version_capture_df = pd.concat(version_capture_df).reset_index(drop = True)
    version_capture_df = version_capture_df.drop(columns = ['sample_name'])
    version_capture_df = version_capture_df.drop_duplicates(keep = 'first')
    version_capture_df = version_capture_df.sort(by = 'software')

    outfile = f'version_capture_influenza_illumina_pe_assembly_{project_name}_{workflow_version}.csv'
    version_capture_df.to_csv(outfile, index = False)

    return None

# create versioin capture file for summary workflow
def create_version_capture_file_for_summary_workflow(
        workflow_version,
        project_name,
        analysis_date):
    
    '''
    create version catpure file for summary workflow
    basically just captures workflow version
    '''

    df = pd.DataFrame()
    df['project_name'] = [project_name]
    df['analysis_date'] = [analysis_date]
    df['software'] = ['influenza_assembly_summary']
    df['associated_docker_container'] = ['']
    df['version'] = [workflow_version]

    outfile = f'version_capture_influenza_assembly_summary_{project_name}_{workflow_version}.csv'
    df.to_csv(outfile, index = False)

    return None

def rename_and_drop_columns(df):
    '''
    becasue not all dfs may be produced, we have to do 
    a clever (and hacky) method to join the dataframes
    together. To do this we loop and use the merge function.
    with the merge function, column conflicts (present in both
    dfs being merged) create _x and _y columns. Since the
    first df will be the empty df, _x columns will always be 
    associated with the column conflict of the empty df and 
    we can drop taht column because it will be empty. So this 
    functions will drop all _x columsn and rename all _y columsn 
    to the original name
    '''
    renamed_cols = {}
    for col in df.columns:
        if col.endswith('_y'):
            renamed_cols[col] = col[:-2]
        elif col.endswith('_x'):
            # because the first left df is the empty df
            # any column conflicts will be _x associated
            # with teh empty df
            del df[col]
    df.rename(columns=renamed_cols, inplace=True)
    return df

#### MAIN ####
if __name__ == '__main__':

    options = getOptions()
    sample_name_txt = options.sample_name
    preprocess_qc_metrics_txt = options.preprocess_qc_metrics
    irma_typing_txt = options.irma_typing
    post_qc_metrics_txt = options.assembly_qc_metrics
    nextclade_na_tsv_txt = options.nextclade_na_tsv
    nextclade_ha_tsv_txt = options.nextclade_ha_tsv
    project_name = options.project_name
    version_capture_file = options.version_capture_file
    workflow_version = options.workflow_version
    analysis_date = options.analysis_date

    sample_name_list = create_list_from_write_lines_input(write_lines_input = sample_name_txt)
    preprocess_qc_metrics_list = create_list_from_write_lines_input(write_lines_input = preprocess_qc_metrics_txt)
    irma_typing_list = create_list_from_write_lines_input(write_lines_input = irma_typing_txt)
    post_qc_metrics_list= create_list_from_write_lines_input(write_lines_input = post_qc_metrics_txt)
    nextclade_na_tsv_list = create_list_from_write_lines_input(write_lines_input = nextclade_na_tsv_txt)
    nextclade_ha_tsv_list = create_list_from_write_lines_input(write_lines_input = nextclade_ha_tsv_txt)
    version_capture_file_list = create_list_from_write_lines_input(write_lines_input = version_capture_file)

    # create version capture file for assembly workflow
    create_version_capture_file_for_assembly_workflow(
            version_capture_file_list = version_capture_file_list,
            project_name = project_name,
            workflow_version = workflow_version)
    
    # crete version capture file for summary workflow
    create_version_capture_file_for_summary_workflow(
        workflow_versioin = workflow_version,
        project_name = project_name,
        analysis_date = analysis_date)


    # create summary metrics file:
    # the issue I have is that if assemlby fails then the the post assembly
    # metrics won't be generated. Plus if the HA or NA segment fail then
    # we won't ahve the nextclade output. So I need to create an empty dataframe
    # that has all the columns so that all the columns will always be output
    # when I'm joining together the dataframes

    # set columns
    col_order = [
        'hsn', 'sample_name', 'project_name', 'analysis_date', 'flu_type', 
        'HA_subtype', 'NA_subtype',
        'HA_clade', 'HA_short-clade', 'HA_subclade',
        'NA_clade', 'NA_short-clade', 'NA_subclade',
        'total_segments', 'total_flu_mapped_reads', 'percent_flu_mapped_reads','total_reads_cleaned',
        'HA_per_cov','HA_mean_depth', 'HA_num_mapped_reads', 'HA_seq_len', 'HA_expected_len',
        'NA_per_cov', 'NA_mean_depth', 'NA_num_mapped_reads', 'NA_seq_len','NA_expected_len', 
        'MP_per_cov', 'MP_mean_depth', 'MP_num_mapped_reads','MP_seq_len', 'MP_expected_len', 
        'NP_per_cov', 'NP_mean_depth','NP_num_mapped_reads', 'NP_seq_len', 'NP_expected_len', 
        'NS_per_cov','NS_mean_depth', 'NS_num_mapped_reads', 'NS_seq_len', 'NS_expected_len',
        'PA_per_cov', 'PA_mean_depth', 'PA_num_mapped_reads', 'PA_seq_len','PA_expected_len', 
        'PB1_per_cov', 'PB1_mean_depth','PB1_num_mapped_reads', 'PB1_seq_len', 'PB1_expected_len',
        'PB2_per_cov', 'PB2_mean_depth', 'PB2_num_mapped_reads', 'PB2_seq_len','PB2_expected_len',
        'total_read_diff',  'total_reads_R1_raw', 'total_reads_R2_raw', 'read_pairs_raw', 'total_reads_raw', 
        'read_length_R1_raw', 'read_length_R2_raw',
        'total_reads_R1_cleaned', 'total_reads_R2_cleaned', 'read_pairs_cleaned', 
        'read_length_R1_cleaned', 'read_length_R2_cleaned',
        'HA_totalSubstitutions', 'HA_totalDeletions', 'HA_totalInsertions', 'HA_totalFrameShifts', 'HA_totalMissing', 'HA_totalNonACGTNs', 'HA_totalAminoacidSubstitutions', 'HA_totalAminoacidDeletions', 'HA_totalAminoacidInsertions', 'HA_totalUnknownAa', 'HA_nextclade_coverage', 'HA_aaSubstitutions', 'HA_aaDeletions', 'HA_aaInsertions', 'HA_warnings', 'HA_errors',
        'NA_totalSubstitutions', 'NA_totalDeletions', 'NA_totalInsertions', 'NA_totalFrameShifts', 'NA_totalMissing', 'NA_totalNonACGTNs', 'NA_totalAminoacidSubstitutions', 'NA_totalAminoacidDeletions', 'NA_totalAminoacidInsertions', 'NA_totalUnknownAa', 'NA_nextclade_coverage', 'NA_aaSubstitutions', 'NA_aaDeletions', 'NA_aaInsertions', 'NA_warnings', 'NA_errors'
        ]

    empty_df = pd.DataFrame(columns = col_order)
    
    to_merge_dfs = [empty_df]

    # preprocess
    preprocess_qc_metrics_df_list = []
    for preprocess_qc_metrics in preprocess_qc_metrics_list:
        df = pd.read_csv(preprocess_qc_metrics, dtype = {'sample_name' : object})
        preprocess_qc_metrics_df_list.append(df)
    preprocess_qc_metrics_df = pd.concat(preprocess_qc_metrics_df_list).reset_index(drop = True)
    to_merge_dfs.append(preprocess_qc_metrics_df)

    # irma subtyping
    irma_typing_df_list = []
    for irma_typing in irma_typing_list:
        df = pd.read_csv(irma_typing, dtype = {'sample_name' : object})
        irma_typing_df_list.append(df)
    irma_typing_df = pd.concat(irma_typing_df_list).reset_index(drop=True)
    to_merge_dfs.append(irma_typing_df)

    # na nextclade
    nextclade_na_df_list = []
    if len(nextclade_na_df_list) >= 1:
        for nextclade_tsv in nextclade_na_tsv_list:
            sample_name = nextclade_tsv.split('_na')[0]
            df = pd.read_csv(nextclade_tsv, sep ='\t')
            df['sample_name'] = sample_name
            df['nextclade_coverage'] = df['coverage']
            # add missing columns
            # NA: add subclade, short-clade (for Bvic, H1N1, and H3N2)
            # HA: add short-clade (For bvic only)
            if "subclade" not in df.columns:
                print('DNE')
                df['subclade'] = ""
            if "short-clade" not in df.columns:
                print('DNE')
                df['short-clade'] = ""
            # reorder columns
            col_keep = ['sample_name', 'clade', 'short-clade', 'subclade', 
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
            nextclade_na_df_list.append(df)
        nextclade_na_df = pd.concat(nextclade_na_df_list).rest_index(drop = True)
        to_merge_dfs.append(nextclade_na_df)

    # ha nextclade
    nextclade_ha_df_list = []
    if len(nextclade_ha_df_list) >= 1:
        for nextclade_tsv in nextclade_ha_tsv_list:
            sample_name = nextclade_tsv.split('_ha')[0]
            df = pd.read_csv(nextclade_tsv, sep ='\t')
            df['sample_name'] = sample_name
            df['nextclade_coverage'] = df['coverage']
            # add missing columns
            # NA: add subclade, short-clade (for Bvic, H1N1, and H3N2)
            # HA: add short-clade (For bvic only)
            if "subclade" not in df.columns:
                print('DNE')
                df['subclade'] = ""
            if "short-clade" not in df.columns:
                print('DNE')
                df['short-clade'] = ""
            # reorder columns
            col_keep = ['sample_name', 'clade', 'short-clade', 'subclade', 
                        'totalSubstitutions','totalDeletions', 'totalInsertions', 
                        'totalFrameShifts', 'totalMissing','totalNonACGTNs', 'totalAminoacidSubstitutions',
                        'totalAminoacidDeletions', 'totalAminoacidInsertions', 'totalUnknownAa', 
                        'nextclade_coverage','aaSubstitutions', 'aaDeletions', 'aaInsertions',
                        'warnings', 'errors']
            # add "nha" prefix to all column headers
            rename_cols = {}
            for col in col_keep:
                if col != 'sample_name':
                    new_column = f'HA_{col}'
                    rename_cols[col] = new_column

            df = df[col_keep]
            df = df.rename(columns = rename_cols)
            nextclade_ha_df_list.append(df)
        nextclade_ha_df = pd.concat(nextclade_ha_df_list).rest_index(drop = True)
        to_merge_dfs.append(nextclade_ha_df)

    # post assembly qc metrics
    post_qc_metrics_df_list = []
    if len(post_qc_metrics_df_list) >= 1:
        for post_qc_metrics in post_qc_metrics_list:
            df = pd.read_csv(post_qc_metrics, dtype = {'sample_name' : object})
            post_qc_metrics_df_list.append(df)
        post_qc_metrics_df = pd.concat(post_qc_metrics_df_list).reset_index(drop = True)
        to_merge_dfs.append(post_qc_metrics_df)
    
    # combinig everything together
    # merge by looping through all existing dfs
    result = to_merge_dfs[0] # this is the empty df
    for df in to_merge_dfs[1:]:
        result = pd.merge(result, df, on = 'col2', how = 'outer')
    result = rename_and_drop_columns(df = result)
    result = result.reset_index()
    result["analysis_date"] = analysis_date
    result['percent_flu_mapped_reads'] = round((df.total_flu_mapped_reads / df.total_reads_cleaned) * 100 , 2)
    result['project_name'] = project_name
    
    # ordering colums while making sure I didn't forget any
    # all columns in not in col_order will be sorted alphabetically at the end of the df
    columns = result.columns.tolist()
    columns.sort() # sort alphabetically

    for n, column in enumerate(col_order):
        print(column)
        columns.remove(column)
        columns.insert(n, column)

    result = result[columns]
    
    #outfile
    outfile = f'{project_name}_sequencing_results.csv' 
    result.to_csv(outfile, index = False)


    













