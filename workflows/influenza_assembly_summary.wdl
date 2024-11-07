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
        Array[Array[File]] nextclade_tsv
        String out_bucket_path

        # python scripts
        File results_summary_py
        File capture_version_summary_py
    }

    String project_name = select_first(project_name_array)
    Array[File] nextclade_tsv_flatten = flatten(select_all(nextclade_tsv))



    call capture_version.capture_workflow_version  as capture_workflow_version{
        input:
    }


    call summary.results_summary as results_summary {
        input:
            sample_name = sample_name,
            preprocess_qc_metrics = preprocess_qc_metrics,
            irma_typing = irma_typing,
            assembly_qc_metrics = assembly_qc_metrics,
            nextclade_tsv_flatten = nextclade_tsv_flatten,
            analysis_date = capture_workflow_version.analysis_date,
            # workflow_version = capture_workflow_version.workflow_version,
            python_script = results_summary_py,
            project_name = project_name
        
    }

    call summary.capture_version_summary as capture_version_summary {
    input:
        workflow_version = capture_workflow_version.workflow_version,
        workflow_name = "influenza_assembly_summary",
        analysis_date = capture_workflow_version.analysis_date,
        python_script = capture_version_summary_py,
        project_name = project_name
        
    }

    call transfer.transfer_assembly_summary_wdl as summary_transfer {
        input:
            workflow_version = capture_workflow_version.workflow_version,
            sequencing_results_csv = results_summary.sequencing_results_csv,
            version_capture_influenza_assembly_summary_csv = capture_version_summary.version_capture_influenza_assembly_summary_csv,
            bucket_path = out_bucket_path
    }

    output {

        File sequencing_results_csv = results_summary.sequencing_results_csv
        File version_capture_influenza_assembly_summary_csv = capture_version_summary.version_capture_influenza_assembly_summary_csv
        String transfer_date = summary_transfer.transfer_date

    }

}