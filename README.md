# CDPHE-influenza

## In development
This repository is in active development. This document describes the Colorado Department of Public Health and Environment's workflow for the assembly and anlaysis of whole genome sequenicng data of influenza A and B ulitizing the Terra.bio platform. The influenza_assembly.wdl workflow was adapted from Thieagen Genomics's theacov_pe workflow for influenza. 

## Workflow overview
Currently assembly is preformed using IRMA and contained in a WDL workflow that can be run on Terra.bio. 

## IRMA overview
Iterative Refinement Meta-Assembler (IRMA) was developed by the CDC. More information about IRMA can be found here: https://wonder.cdc.gov/amd/flu/irma/. IRMA performs genome assembly and variant calling of flu. Illumina paired-end, Illumina single-end, and ONT data can be used with IRMA; however our workflow uses illumina pair-end data. We use the default configuration file for the IRMA FLU module whch mean we use ALIGN_PROG=SAM and DEL_TYPE=''. 

## WDL workflow: influenza_assembly.wdl
The workflow performs the following steps:

1. **FastQC.** FastQC is performed on the raw reads and on the cleaned reads generated after step 2, Seqyclean. Output files for both the forward and reverse reads are generated. The total number of reads and the read length for the forward and reverse reads are recorded from the "fastq_data.txt" file. And finally the total number of read pairs is recorded.

2. **Seqyclean.** Seqyclean is run on the raw reads with the parameters min length set to 70bp and quality set to 30 30. The adapters and contaminats fasta file can be found in the workspace data directory and should be linked to your workspace data in terra.

3. **IRMA.** The 8 gene sgements are assembled using IRMA (see IRMA overview above). We use the "FLU" module. The fasta header is renamed to include the sample id (e.g. ">{sample_id}_A_HA_H1"). The fasta, bam, and vcf files are renamed to include the sample name in the file name (e.g. "{sample_id}_A_HA_H1.fasta").  We are currently using the fasta files from the main directory and not the fasta files from the "amended_consensus" directory. The difference is that the fasta files in the "amended_consensus" directory use IUPAC ambiguity codes. Note: the fastq files must follow the illumina style fastq header. Fastq files downloaded from SRA will need to be modified prior to running the workflow.   

4. **Post Assembly QC Metrics.** For each gene segment this step calculates the average depth and number of mapped reads using samtools and calcuates the percent coverage and determined the lenght using a custom python script. The results are then writed to a summary table ({sample_id}_qc_metrics.csv).

5. **Transfer.** This task transfers output from FastQC, Seqyclean, IRMA and Post Assmebly QC Metrics to the GCP bucket specified in the "out_dir" column in the terra data table. The structure is as follows:

    * gs://{out_dir}/fastqc_raw/
    * gs://{out_dir}/fastqc_clean/
    * gs://{out_dir}/seqyclean/
    * gs://{out_dir}/irma/{sample_id}/assemblies/
    * gs://{out_dir}/irma/{sample_id}/bam_files/
    * gs://{out_dir}/irma/{sample_id}/vcfs/
    * gs://{out_dir}/irma/{sample_id}/

### Running influenza_assembly.wdl on Terra

#### Import Workflow from DockStore
#### Setting up your data table
The datatable should look like the following and be saved as a tsv or txt file:

| entity:sample_id   |  fastq_R1   | fastq_R2 | out_dir |
|-------------------|-------------|-----------|---------------------|
| sample_name     | gs://path_to_fastq_R1 | gs://path_to_fastq_R2 | gs://path_to_transfer_output


#### Setting up your workspace data


