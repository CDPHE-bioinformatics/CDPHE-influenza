version 1.0

import "../tasks/summary_task.wdl" as summary 
import "../tasks/transfer_tasks.wdl" as transfer

workflow influenza_assembly_summary{
    
    input {
        Array[String] sample_name
        Array[String] project_name_array
        Array[File] preprocess_qc_metrics
        Array[File] irma_typing
        Array[File] irma_assembly_qc_metrics
        Array[File] nextclade_na_tsv
        Array[File] nextclade_ha_tsv
        String out_bucket_path
        Array[String] analysis_date_array

        # python scripts
        File summary_py
    }

    # String bucket_path = select_first(bucket_path_array)
    String project_name = select_first(project_name_array)
    String analysis_date = select_first(analysis_date_array)

    call summary.summary as summary {
        input:
            sample_name = sample_name,
            preprocess_qc_metrics = preprocess_qc_metrics,
            irma_typing = irma_typing,
            irma_assembly_qc_metrics = irma_assembly_qc_metrics,
            nextclade_na_tsv = nextclade_na_tsv,
            nextclade_ha_tsv = nextclade_ha_tsv,

            
            python_script = summary_py,
            project_name = project_name,
            analysis_date = analysis_date
        
    }


     call version_capture.workflow_version_capture  as workflow_version_capture{
        input:
    }

    Array[VersionInfo] version_array = [
        seqyclean.seqyclean_version_info,
        fastqc_cleaned.fastqc_version_info,
        align_reads.bwa_version_info,
        align_reads.samtools_version_info,
        ivar_consensus.ivar_version_info,
        ivar_consensus.samtools_version_info,
        bam_stats.samtools_version_info
    ]
    if (scrub_reads) {
        Array[VersionInfo] version_array_with_hostile = flatten([version_array, select_all([hostile.hostile_version_info])])
    }
    call version_capture.task_version_capture as task_version_capture {
        input:
            version_array = select_first([version_array_with_hostile, version_array]),
            workflow_name = "SC2_illumina_pe_assembly",
            workflow_version = workflow_version_capture.workflow_version,
            project_name = project_name,
            analysis_date = workflow_version_capture.analysis_date,
            version_capture_py = version_capture_py
    }

    call transfer.transfer_assembly_summary_wdl as summary_transfer {
        input:
            sequencing_results_csv = summary.sequencing_results_csv,
            bucket_path = out_bucket_path
    }

    output {

        File sequencing_results_csv = summary.sequencing_results_csv
        String transfer_date = summary_transfer.transfer_date

    }

}