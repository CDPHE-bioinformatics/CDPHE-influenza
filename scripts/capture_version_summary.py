#! /usr/bin/env python

# import python modules
import pandas as pd
import sys
import argparse
import subprocess


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "--workflow_version")
    parser.add_argument('--workflow_name')
    parser.add_argument( "--analysis_date")
    parser.add_argument( "--project_name")

    options = parser.parse_args(args)
    return options

# create versioin capture file for summary workflow
def create_version_capture_file_for_summary_workflow(
        workflow_version,
        workflow_name,
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

    outfile = f'version_capture_{workflow_name}_{project_name}.csv'
    df.to_csv(outfile, index = False)

    return None


#### MAIN ####
if __name__ == '__main__':

    options = getOptions()
    project_name = options.project_name
    workflow_name = options.workflow_name
    workflow_version = options.workflow_version
    analysis_date = options.analysis_date

    
    # crete version capture file for summary workflow
    create_version_capture_file_for_summary_workflow(
        workflow_version = workflow_version,
        workflow_name = workflow_name,
        project_name = project_name,
        analysis_date = analysis_date)