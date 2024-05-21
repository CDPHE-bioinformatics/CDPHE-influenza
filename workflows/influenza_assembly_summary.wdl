version 1.0

import "../tasks/summary_task.wdl" as summary 
import "../tasks/transfer_tasks.wdl" as transfer

workflow influenza_assembly_summary{
    
    input {
        Array[String] sample_name
        Array[String] project_name_array
        Array[File] preprocess_qc_metrics
        Array[File] irma_typing
        Array[File] post_assembly_qc_metrics
        Array[File] nextclade_na_tsv
        Array[File] nextclade_ha_tsv
        String out_bucket_path
        Array[File] version_capture_file

        # python scripts
        File summary_py
    }

    String project_name = select_first(project_name_array)



    call version_capture.workflow_version_capture  as workflow_version_capture{
        input:
    }


    call summary.summary as summary {
        input:
            sample_name = sample_name,
            preprocess_qc_metrics = preprocess_qc_metrics,
            irma_typing = irma_typing,
            post_assembly_qc_metrics = post_assembly_qc_metrics,
            nextclade_na_tsv = nextclade_na_tsv,
            nextclade_ha_tsv = nextclade_ha_tsv,
            version_capture_file = version_capture_file,
            workflow_version = workflow_version_capture.workflow_version,
            analysis_date = workflow_version_capture.analysis_date,
            python_script = summary_py,
            project_name = project_name
        
    }



    call transfer.transfer_assembly_summary_wdl as summary_transfer {
        input:
            sequencing_results_csv = summary.sequencing_results_csv,
            version_capture_influenza_assembly_csv = summary.version_capture_influenza_assembly_csv
            version_capture_influena_assembly_summary_csv = summary.version_capture_influenza_summary_csv
            bucket_path = out_bucket_path
    }

    output {

        File sequencing_results_csv = summary.sequencing_results_csv
        File version_capture_influenza_assembly_csv = summary.version_capture_influenza_assembly_csv
        File version_capture_influena_assembly_summary_csv = summary.version_capture_influena_assembly_summary_csv
        String transfer_date = summary_transfer.transfer_date

    }

}