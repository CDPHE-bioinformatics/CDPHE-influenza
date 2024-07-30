version 1.0

import "../tasks/summary_task.wdl" as summary 
import "../tasks/transfer_tasks.wdl" as transfer
import "../tasks/capture_version_tasks.wdl" as capture_version

workflow influenza_assembly_summary{
    
    input {
        Array[String] sample_name
        Array[String] project_name_array
        Array[File] preprocess_qc_metrics
        Array[File] irma_typing
        Array[File] assembly_qc_metrics
        Array[File] na_nextclade_tsv
        Array[File] ha_nextclade_tsv
        String out_bucket_path
        Array[File] version_capture_file

        # python scripts
        File summary_py
    }

    String project_name = select_first(project_name_array)



    call capture_version.capture_workflow_version  as capture_workflow_version{
        input:
    }


    call summary.summary as summary {
        input:
            sample_name = sample_name,
            preprocess_qc_metrics = preprocess_qc_metrics,
            irma_typing = irma_typing,
            assembly_qc_metrics = assembly_qc_metrics,
            na_nextclade_tsv = na_nextclade_tsv,
            ha_nextclade_tsv = ha_nextclade_tsv,
            version_capture_file = version_capture_file,
            workflow_version = capture_workflow_version.workflow_version,
            analysis_date = capture_workflow_version.analysis_date,
            python_script = summary_py,
            project_name = project_name
        
    }


    call transfer.transfer_assembly_summary_wdl as summary_transfer {
        input:
            sequencing_results_csv = summary.sequencing_results_csv,
            version_capture_influenza_assembly_csv = summary.version_capture_influenza_assembly_csv,
            version_capture_influenza_assembly_summary_csv = summary.version_capture_influenza_assembly_summary_csv,
            bucket_path = out_bucket_path
    }

    output {

        File sequencing_results_csv = summary.sequencing_results_csv
        File version_capture_influenza_assembly_csv = summary.version_capture_influenza_assembly_csv
        File version_capture_influenza_assembly_summary_csv = summary.version_capture_influenza_assembly_summary_csv
        String transfer_date = summary_transfer.transfer_date

    }

}