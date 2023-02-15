version 1.0

# import tasks
import "../tasks/preprocess_tasks.wdl" as fastq_preprocess
import "../tasks/irma_task.wdl" as irma
import "../tasks/ivar_task.wdl" as ivar
import "../tasks/post_assembly_tasks.wdl" as post_assembly_qc
import "../tasks/transfer_tasks.wdl" as transfer

# in order to check right import
# begin workflow
workflow influenza_assembly {

    input {
        String sample_id
        String read_type
        File fastq_R1
        File? fastq_R2
        File adapters_and_contaminants
        String bucket_path

        # python scripts
        File concat_preprocess_qc_metrics_py
        File calc_percent_cov_py
        File concat_post_assembly_qc_metrics_py
    }

    # 1 - Preprocess QC raw fastq files
    call fastq_preprocess.fastqc as fastqc_raw {
        input:
            sample_id = sample_id,
            fastq_R1 = fastq_R1,
            fastq_R2 = fastq_R2,
            read_type = read_type
    }
    call fastq_preprocess.seqyclean as seqyclean {
        input:
            sample_id = sample_id,
            read_type = read_type,
            fastq_R1 = fastq_R1,
            fastq_R2 = fastq_R2,
            adapters_and_contaminants = adapters_and_contaminants  
    }
    call fastq_preprocess.fastqc as fastqc_cleaned {
        input:
            sample_id = sample_id,
            fastq_R1 = seqyclean.fastq_R1_cleaned,
            fastq_R2 = seqyclean.fastq_R2_cleaned,
            read_type = read_type
    }
    call fastq_preprocess.concat_preprocess_qc_metrics as concat_preprocess_qc_metrics{
        input:
            python_script = concat_preprocess_qc_metrics_py,
            sample_id = sample_id,
            read_type = read_type,

            fastqc_version = fastqc_raw.fastqc_version,
            fastqc_docker = fastqc_raw.fastqc_docker,

            total_reads_R1_raw = fastqc_raw.total_reads_R1,
            total_reads_R2_raw = fastqc_raw.total_reads_R2,
            read_length_R1_raw = fastqc_raw.read_length_R1,
            read_length_R2_raw = fastqc_raw.read_length_R2,
            read_pairs_raw = fastqc_raw.read_pairs,

            total_reads_R1_cleaned = fastqc_cleaned.total_reads_R1,
            total_reads_R2_cleaned = fastqc_cleaned.total_reads_R2,
            read_length_R1_cleaned = fastqc_cleaned.read_length_R1,
            read_length_R2_cleaned = fastqc_cleaned.read_length_R2,
            read_pairs_cleaned = fastqc_cleaned.read_pairs,

            seqyclean_version = seqyclean.seqyclean_version,
            seqyclean_docker = seqyclean.seqyclean_docker


    }

    # 2- run irma
    call irma.irma as irma {
        input:
            sample_id = sample_id,
            read_type = read_type,
            fastq_R1 = seqyclean.fastq_R1_cleaned,
            fastq_R2 = seqyclean.fastq_R2_cleaned
    }

    # proceed with post assembly QC metrics if irma assembly was successful
    if (irma.irma_type != "no IRMA assembly generated") {
        # 3- post assembly QC metrics
        ## scatter over bam file array

        scatter (bam_file in select_all(irma.irma_bam_files)) {
            call post_assembly_qc.samtools_mapped_reads as irma_samtools_mapped_reads { 
                input:
                    bam_file = bam_file
            }

            call ivar.ivar_consensus as ivar_consensus {
                input:
                bam_file = bam_file
            }
        }
        
        ## scatter over fasta file array from ivar
        scatter (fasta_file in select_all(ivar_consensus.ivar_consensus_fasta)) {
            call post_assembly_qc.calc_percent_coverage as irma_percent_coverage { 
                input:
                    fasta_file  = fasta_file,
                    python_script = calc_percent_cov_py
                    
            }
        }

        # input depth and cov arrays to create a single outfile for all gene segments
        call post_assembly_qc.concat_post_qc_metrics as irma_concat_post_qc_metrics{
            input:
                python_script = concat_post_assembly_qc_metrics_py,
                sample_id = sample_id,
                bam_results_array = irma_samtools_mapped_reads.bam_results,
                per_cov_results_array = irma_percent_coverage.perc_cov_results,
                ivar_parameters = ivar_consensus.ivar_parameters

        }
    
    }
    # 5 - Transfer
    
    call transfer.transfer_assembly_wdl as transfer_assembly_wdl {
        input:
            sample_id = sample_id, 
            bucket_path = bucket_path,
            fastqc1_html_raw = fastqc_raw.fastqc1_html,
            fastqc1_zip_raw = fastqc_raw.fastqc1_zip,
            fastqc2_html_raw = fastqc_raw.fastqc2_html,
            fastqc2_zip_raw = fastqc_raw.fastqc2_zip,

            seqyclean_summary = seqyclean.seqyclean_summary,

            preprocess_qc_metrics = concat_preprocess_qc_metrics.preprocess_qc_metrics,

            fastqc1_html_cleaned = fastqc_cleaned.fastqc1_html,
            fastqc1_zip_cleaned = fastqc_cleaned.fastqc1_zip,
            fastqc2_html_cleaned = fastqc_cleaned.fastqc2_html,
            fastqc2_zip_cleaned = fastqc_cleaned.fastqc2_zip,

            irma_assemblies = irma.irma_assemblies,
            irma_bam_files = irma.irma_bam_files,
            irma_vcfs = irma.irma_vcfs,

            irma_sorted_bams = irma_samtools_mapped_reads.sorted_bam,

            ivar_assemblies = ivar_consensus.ivar_consensus_fasta,
            ivar_outputs = ivar_consensus.ivar_output,

            irma_qc_metrics = irma_concat_post_qc_metrics.qc_metrics_summary


    }


    output {
        # output from preprocess
        String fastqc_version = fastqc_raw.fastqc_version
        String fastqc_docker = fastqc_raw.fastqc_docker

        File fastqc1_html_raw = fastqc_raw.fastqc1_html
        File fastqc1_zip_raw = fastqc_raw.fastqc1_zip
        File? fastqc2_html_raw = fastqc_raw.fastqc2_html
        File? fastqc2_zip_raw = fastqc_raw.fastqc2_zip

        String seqyclean_version = seqyclean.seqyclean_version
        String seqyclean_docker = seqyclean.seqyclean_docker
        File seqyclean_summary = seqyclean.seqyclean_summary

        File fastqc1_html_cleaned = fastqc_cleaned.fastqc1_html
        File fastqc1_zip_cleaned = fastqc_cleaned.fastqc1_zip
        File? fastqc2_html_cleaned = fastqc_cleaned.fastqc2_html
        File? fastqc2_zip_cleaned = fastqc_raw.fastqc2_zip

        File preprocess_qc_metrics = concat_preprocess_qc_metrics.preprocess_qc_metrics

        # output from irma
        String irma_type = irma.irma_type
        String irma_ha_subtype = irma.irma_ha_subtype
        String irma_na_subtype = irma.irma_na_subtype
        File irma_typing = irma.irma_typing
        # Array[String] segment_array = irma.segment_array
        Array[File] irma_assemblies = irma.irma_assemblies
        Array[File] irma_bam_files = irma.irma_bam_files
        Array[File] irma_vcfs = irma.irma_vcfs
        String irma_version = irma.irma_version
        String irma_docker = irma.irma_docker

        # output from ivar_consensus and samtools
        Array[File]? irma_sorted_bams = irma_samtools_mapped_reads.sorted_bam
        Array[File]? ivar_assemblies = ivar_consensus.ivar_consensus_fasta
        Array[File]? ivar_outputs = ivar_consensus.ivar_output

        # output from post assembly QC metrics
        Array[File]? irma_bam_results = irma_samtools_mapped_reads.bam_results
        Array[File]? irma_per_cov_results = irma_percent_coverage.perc_cov_results
        File? irma_assembly_qc_metrics = irma_concat_post_qc_metrics.qc_metrics_summary
        
        # output from transfer
        String transfer_date=transfer_assembly_wdl.transfer_date
    }
}
