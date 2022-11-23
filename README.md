# CDPHE-influenza

## In development
This repository is in active development. This document describes the Colorado Department of Public Health and Environment's workflow for the assembly and anlaysis of whole genome sequenicng data of influenza A and B ulitizin ght Terra.bio platform. 

## Workflow overview
Currently assembly is preformed using IRMA and contained in a WDL workflow that can be run on Terra.bio. THe results from IRMA are then summarized using a juypter notebook. Hopefully soon the code from the jupyter notebook will be integrated in the WDL workflow.

## IRMA overview
Iterative Refinement Meta-Assembler (IRMA) was developed by the CDC. More information about IRMA can be found here: https://wonder.cdc.gov/amd/flu/irma/. IRMA performs genome assembly and variant calling of flu. Illumina paired-end, Illumina single-end, and ONT data can be used with IRMA; however our workflow uses illumina pair-end data.

## WDL workflow: cdphe-influnenza
The workflow performs the following steps:
1. **FastQC.** FastQC is performed on the raw reads and on the cleaned reads generated after step 2, Seqyclean. Output files for both the forward and reverse reads are generated. The total number of reads and the read length for the forward and reverse reads are recorded from the "fastq_data.txt" file. And finally the total number of read pairs is recorded.

2. **Seqyclean.** Seqyclean is run on the raw reads with the parameters min length set to 70bp and quality set to 30 30. The adapters and contaminats fasta file can be found in the workspace data directory and should be linked to your workspace data in terra.

3. **IRMA.** The 8 gene sgements are assembled using IRMA (see IRMA overview above). We use the "FLU" module. The fasta, bam, bam.bai, and data table files are renamed to include the sample name in the file name and are collected from the output. We are currently using the fasta files from the main directory and not the fasta files from the fasta files from the "amended_consensus" directory. The difference is that the fasta files in the "amended_consensus" directory use IUPAC ambiguity codes. Note: the fastq files must follow the illumina style fastq header. Fastq files downloaded from SRA will need to be modified prior to running the workflow.   

4. **Collect run parameters.** This task pulls the versioning information from the FastaQC, Seqyclean, and IRMA tasks and the read quality data from the FastQC task and concatates it into a csv file.

5. **Transfer intermediate files to GCP bucket.** This task transfers output from FastQC, Seqyclean, IRMA and the csv file with run parameters, to a GCP bucket of your choosing.


### Running cdphe-influenza on Terra

#### Import Workflow from DockStore
#### Setting up your data table
The datatable should look like the following and be saved as a tsv or txt file:

| entity:sample_id   |  fastq_R1   | fastq_R2 | gcp_out_bucket_path |
| :------------- | :------------- |
| sample_name     | gs://path_to_fastq_R1      | gs://path_to_fastq_R2 | gs://path_to_transfer_output


#### Setting up your workspace data

## Python Jupyter notebook
This notebook uses the consensus sequences, and -coverage.txt files generated from IRMA to determine the type, the subtype (if flu A), calculate the percent coverage of each gene segment and the average sequencing depth for each gene segement (i.e the assembly metrics). Then the notebook combines the assembly metrics with the run parameters data to generate a single summary csv file.  

Specially, this notebook, reads in the terra datatable to generate a list of sample names. Then the script pulls down the the consensus sequences (fasta files), -coverage.txt files, and the run_parameters.csv file from the gcp bucket. Next it uses the names of the fasta files to determine the type and subytpe (if A). Next it calculates the percnet coverage of each gene segement using the conensus sequence compared to the length of the refernce gene segmenet (the reference gene segement sequences are those listed in reference file IRMA). Next the coverage.txt file is used to calculate the average sequencing depth for each gene segmenet. Finally the information for the run_parameters.csv file is combmined withe assembly metrics and a summary file as a csv is generated. Examples of output can be found in example_data directory.  
