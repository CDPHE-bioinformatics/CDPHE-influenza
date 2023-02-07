version 1.0

# import tasks
import "/home/molly_hetheringtonrauth/sandbox_influenza/pipeline_development/tasks/preprocess_tasks.wdl" as fastq_preprocess
import "/home/molly_hetheringtonrauth/sandbox_influenza/pipeline_development/tasks/irma_task.wdl" as irma
import "/home/molly_hetheringtonrauth/sandbox_influenza/pipeline_development/tasks/post_assembly_tasks.wdl" as post_assembly_qc
import "/home/molly_hetheringtonrauth/sandbox_influenza/pipeline_development/tasks/transfer_task.wdl" as transfer

# begin workflow
workflow influenza_assembly {

    input {
        String sample_id
        String read_type
        File fastq_R1
        File fastq_R2
        File adapters_and_contaminants
        String bucket_path

        # python scripts
        File calc_percent_cov_py
        File concat_qc_metrics_py
    }

    # 1 - Preprocess QC raw fastq files
    call fastq_preprocess.fastqc as fastqc_raw {
        input:
            fastq_R1 = fastq_R1,
            fastq_R2 = fastq_R2
    }
    call fastq_preprocess.seqyclean as seqyclean {
        input:
            sample_id = sample_id,
            fastq_R1 = fastq_R1,
            fastq_R2 = fastq_R2,
            adapters_and_contaminants = adapters_and_contaminants  
    }
    call fastq_preprocess.fastqc as fastqc_cleaned {
        input:
            fastq_R1 = seqyclean.fastq_R1_cleaned,
            fastq_R2 = seqyclean.fastq_R2_cleaned
    }

    # 2- run irma
    call irma.irma as irma {
        input:
            sample_id = sample_id,
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
        }
        
        ## scatter over fasta file array
        scatter (fasta_file in select_all(irma.irma_assemblies)) {
            call post_assembly_qc.calc_percent_coverage as irma_percent_coverage { 
                input:
                    fasta_file  = fasta_file,
                    python_script = calc_percent_cov_py
            }
        }

        # input depth and cov arrays to create a single outfile for all gene segments
        call post_assembly_qc.concat_post_qc_metrics as irma_concat_post_qc_metrics{
            input:
                python_script = concat_qc_metrics_py,
                sample_id = sample_id,
                bam_results_array = irma_samtools_mapped_reads.bam_results,
                per_cov_results_array = irma_percent_coverage.perc_cov_results

        }
    
    }

    #4 - Write summary outputs for each sample so just have to concatenate
    ## files for results wdl+
    5 - Transfer
    call transfer.transfer_assembly_wdl as transfer {
        input:
            sample_id = sample_id, 
            bucket_path = bucket_path,
            fastqc1_html_raw = fastqc_raw.fastqc1_html,
            fastqc1_zip_raw = fastqc_raw.fastqc1_zip,
            fastqc2_html_raw = fastqc_raw.fastqc2_html,
            fastqc2_zip_raw = fastqc_raw.fastqc2_zip,

            seqyclean_summary = seqyclean.seqyclean_summary,

            fastqc1_html_cleaned = fastqc_cleaned.fastqc1_html,
            fastqc1_zip_cleaned = fastqc_cleaned.fastqc1_zip,
            fastqc2_html_cleaned = fastqc_cleaned.fastqc2_html,
            fastqc2_zip_cleaned = fastqc_raw.fastqc2_zip,

            irma_assemblies = irma.irma_assemblies,
            irma_bam_files = irma.irma_bam_files,
            irma_vcfs = irma.irma_vcfs,

            irma_qc_metrics = irma_concat_post_qc_metrics.qc_metrics_summary


    # }


    output {
        # output from preprocess
        String fastqc_version = fastqc_raw.fastqc_version
        String fastqc_docker = fastqc_raw.fastqc_docker

        File fastqc1_html_raw = fastqc_raw.fastqc1_html
        File fastqc1_zip_raw = fastqc_raw.fastqc1_zip
        File fastqc2_html_raw = fastqc_raw.fastqc2_html
        File fastqc2_zip_raw = fastqc_raw.fastqc2_zip

        String seqyclean_version = seqyclean.seqyclean_version
        String seqyclean_docker = seqyclean.seqyclean_docker
        File seqyclean_summary = seqyclean.seqyclean_summary

        File fastqc1_html_cleaned = fastqc_cleaned.fastqc1_html
        File fastqc1_zip_cleaned = fastqc_cleaned.fastqc1_zip
        File fastqc2_html_cleaned = fastqc_cleaned.fastqc2_html
        File fastqc2_zip_cleaned = fastqc_raw.fastqc2_zip

        # output from irma
        String irma_type = irma.irma_type
        String irma_ha_subtype = irma.irma_ha_subtype
        String irma_na_subtype = irma.irma_na_subtype
        # Array[String] segment_array = irma.segment_array
        Array[File] irma_assemblies = irma.irma_assemblies
        Array[File] irma_bam_files = irma.irma_bam_files
        Array[File] irma_vcfs = irma.irma_vcfs
        String irma_version = irma.irma_version
        String irma_docker = irma.irma_docker

        # output from post assembly QC metrics
        Array[File]? irma_bam_results = irma_samtools_mapped_reads.bam_results
        Array[File]? irma_per_cov_results = irma_percent_coverage.perc_cov_results
        File? irma_qc_metrics = irma_concat_post_qc_metrics.qc_metrics_summary
        
        # output from transfer
        String transfer=transfer.transfer_date
    }
}
