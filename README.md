# CDPHE-influenza

## In development
This repository is in active development. This document describes the Colorado Department of Public Health and Environment's workflow for the assembly and anlaysis of whole genome sequenicng data of influenza A and B ulitizing the Terra.bio platform. The influenza_assembly.wdl workflow was adapted from and influenced by Thieagen Genomics's wf_theiacov_illumina_pe workflow for influenza. 


## Workflow overview
This workflow can handle paired end (PE) and single end (SE) illumina short read data. Currently assembly is preformed using IRMA and contained in a WDL workflow that can be run on Terra.bio. The influenza_assembly.wdl is run on the sample and the influenza_assembly_summary.wdl is run on a sample_set. 

## IRMA overview
Iterative Refinement Meta-Assembler (IRMA) was developed by the CDC. More information about IRMA can be found here: https://wonder.cdc.gov/amd/flu/irma/. IRMA performs genome assembly and variant calling of flu. Illumina paired-end, Illumina single-end, and ONT data can be used with IRMA; however our workflow uses illumina pair-end data. We use the default configuration file for the IRMA FLU module whch mean we use ALIGN_PROG=SAM and DEL_TYPE=''. 

## WDL workflow: influenza_assembly.wdl
(insert drop down)
(insert workflow diagram)

The workflow performs the following steps:

1. **FastQC.** FastQC is performed on the raw reads and on the cleaned reads generated after step 2, Seqyclean. Output files for both the forward and reverse reads are generated. The total number of reads and the read length for the forward and reverse reads are recorded from the "fastq_data.txt" file. And finally the total number of read pairs is recorded.

2. **Seqyclean.** Seqyclean is run on the raw reads with the parameters min length set to 70bp and quality set to 30 30. The adapters and contaminats fasta file can be found in the workspace data directory and should be linked to your workspace data in terra.

3. **IRMA.** The 8 gene sgements are assembled using IRMA (see IRMA overview above). We use the "FLU" module. The fasta header is renamed to include the sample id (e.g. ">{sample_id}_A_HA_H1"). The fasta, bam, and vcf files are renamed to include the sample name in the file name (e.g. "{sample_id}_A_HA_H1.fasta").  We are currently using the fasta files from the main directory and not the fasta files from the "amended_consensus" directory. The difference is that the fasta files in the "amended_consensus" directory use IUPAC ambiguity codes. Note: the fastq files must follow the illumina style fastq header. Fastq files downloaded from SRA will need to be modified prior to running the workflow.   

4. **Post Assembly QC Metrics.** For each gene segment this step calculates the average depth and number of mapped reads using samtools and calcuates the percent coverage and determined the lenght using a custom python script. The results are then writed to a summary table ({sample_id}_qc_metrics.csv).

5. **Transfer.** This task transfers output from FastQC, Seqyclean, IRMA and Post Assmebly QC Metrics to the GCP bucket specified in the "out_dir" column in the terra data table. The structure is as follows:

    * gs://{out_dir}/fastqc_raw/
    * gs://{out_dir}/fastqc_clean/
    * gs://{out_dir}/seqyclean/
    * gs://{out_dir}/preprocess_qc_metrics/
    * gs://{out_dir}/irma/{sample_id}/assemblies/
    * gs://{out_dir}/irma/{sample_id}/bam_files/
    * gs://{out_dir}/irma/{sample_id}/vcfs/
    * gs://{out_dir}/irma/{sample_id}/

### Running influenza_assembly.wdl on Terra

#### Setting up your data table
The datatable should look like the following and be saved as a tsv or txt file:

| entity:sample_id   | read_type |  fastq_R1   | fastq_R2 | out_dir |
|-------------------|-----|-------------|-----------|---------------------|
| sample_name     | paired | gs://path_to_fastq_R1 | gs://path_to_fastq_R2 | gs://path_to_transfer_output

* note if you are using SE data, the fastq_R2 column should be blank and read_type should be "single". 

#### Workflow inputs
- to see an example of inputs as we have it set up, see the influenza_assembly_input.json
- optional inputs

| WDL variable | task | default| options|
|----|----|----|---|
|flu_module| irma | "FLU" | "FLU",  |



#### Setting up your workspace data
We have our workflow setup so that the following data files are stored in our workspace data. You can find these files in the scripts and references directory within this repo.

| Worksapce variable name | WDL variable name | File name |
| ------------------- | ----------- | ----------- |
| adapters_and_contaminants | adapters_and_contaminants | Adapters_plus_PhiX_174.fasta |
| flu_concat_preprocess_qc_metrics_py  |concat_preprocess_qc_metrics_py | concat_preprocess_qc_metrics.py|
| flu_calc_percent_cov_py | calc_percent_cov_py | calculate_percnet_cov.py
| flu_concat_post_assemby_qc_metrics_py | concat_post_assembly_qc_py | concat_post_assembly_qc.py| 





### Outputs

**Preprocessing**

|WDL Output variable name | File Name | Description |
|-------|------|------------|
| fastqc_version | N/A | version of fastqc |
| fastqc_docker | N/A | docker used for fastqc |
| fastqc1_html_raw | {sample_id}_R1_001_fastqc.html | |
| fastqc1_zip_raw | {sample_id}_R1_001_fastqc.zip| |
| fastqc2_html_raw | {sample_id}_R2_001_fastqc.html | empty if read_type == "single"|
| fastqc2_zip_raw | {sample_id}_R1_001_fastqc.zip| empty if read_type == "single"|
| fastqc1_html_cleaned | {sample_id}_R1_001_fastqc.html| |
| fastqc1_zip_clenaed | {sample_id}_R1_001_fastqc.zip| |
| fastqc2_html_clenaed | {sample_id}_R2_001_fastqc.html| empty if read_type == "single"|
| fastqc2_zip_cleaned | {sample_id}_R2_001_fastqc.zip| empty if read_type == "single"|
| seqyclean_version | N/A | version of seqyclean |
| seqyclean_docker | N/A | docker used for seqyclean | 
| seqyclean_summary | {smaple_id}_clean_SummaryStatistics.tsv | |



**IRMA Assembly**

|WDL Output variable name | File Name | Description |
|-------|------|------------|
| irma_verison | N/A | version of IRMA |
| irma_module | N/A | IRMA module used; default "FLU"|
| irma_docker| N/A | docker used for IRMA |
| irma_type | N/A | influenza type called by IRMA; options A, B, N/A|
| irma_ha_subtype | N/A | if influenza type == "A" then it is the influenza subtype for the HA gene called by IRMA; commonly "H1" or "H3"|
| irma_na_subtype | N/A | if influenza type == "A" then it is the influenza subtype for the NA gene called by IRMA; commonly "N1" or "N2" |
|irma_typing| {sample_id}_irma_typing.csv | csv file with the sample id, irma type, irma ha subytpe and irma na subtype listed in a tabluar format|
|irma_assemblies| {sample_id}_{flu_type}\_{gene_segment}.fasta | array of consensus assembly fasta files. Each assembled gene segment has a fasta file. The fasta header is formatted as : ">{sample_id}_{flu_type}\_{gene_segment}"|
|irma_bam_files| {sample_id}_{flu_type}\_{gene_semgnet}.bam | Array of bam files. Each assembled gene segment has a bam file. The reference sequence is the final iterative plurality consensus |
|irma_vcfs | {sample_id}_{flu_type}\_{gene_semgnet}.vcf | Array of vcf files. Each assembled gene segment has a vcf file. The reference sequence is the final iterative plurality consensus. |

 

**Post Assembly QC**
|WDL Output variable name | File Name | Description |
|-------|------|------------|
|irma_bam_results| bam_results.csv| Array of bam_results files. Each assembied gene segment has a bam_resutl.csv file. Contains the segment name, number of reads mapped and the mean depth for that segment in a tabular format. Produced only if assembly is successful.|
|irma_per_cov_results| | Array of per_cov_results files. Each assembied gene segment has a per_cov_results.csv file. Contains the segment name, percent of genome covered (percent_coverage), the expected gene segmenet length (based on the seed reference segement size IRMA usees), and the assemblied gene segement length for that segment in a tabular format.Produced only if assembly is successful.|
|irma_qc_metrics | {sample_id}_qc_metrics.csv | A tablular formated file that contains the totally of the perc_cov_resutls and the bam_resutls files for all gene segments successfully assemblied. Produced only if assembly is successful.|


