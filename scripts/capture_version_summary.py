#! /usr/bin/env python

# import python modules
import pandas as pd
import sys
import argparse
import subprocess


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "--version_capture_file")
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
    version_capture_df = pd.concat(version_capture_df_list).reset_index(drop = True)
    version_capture_df = version_capture_df.drop(columns = ['sample_name'])
    version_capture_df = version_capture_df.drop_duplicates(keep = 'first')
    version_capture_df = version_capture_df.sort_values(by = 'software')

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

    outfile = f'version_capture_influenza_illumina_pe_assembly_summary_{project_name}_{workflow_version}.csv'
    df.to_csv(outfile, index = False)

    return None


#### MAIN ####
if __name__ == '__main__':

    options = getOptions()
    project_name = options.project_name
    version_capture_file = options.version_capture_file
    workflow_version = options.workflow_version
    analysis_date = options.analysis_date

    version_capture_file_list = create_list_from_write_lines_input(write_lines_input = version_capture_file)

    # create version capture file for assembly workflow
    create_version_capture_file_for_assembly_workflow(
            version_capture_file_list = version_capture_file_list,
            project_name = project_name,
            workflow_version = workflow_version)
    
    # crete version capture file for summary workflow
    create_version_capture_file_for_summary_workflow(
        workflow_version = workflow_version,
        project_name = project_name,
        analysis_date = analysis_date)


    


    













