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
        Array[String] bucket_path_array
        Array[String] analysis_date_array
        Array[File] workbook_path_array

        File summary_py
    }

    String bucket_path = select_first(bucket_path_array)
    String project_name = select_first(project_name_array)
    String analysis_date = select_first(analysis_date_array)
    String workbook_path = select_first(workbook_path_array)

    call summary.summary as summary {
        input:
            sample_name = sample_name,
            preprocess_qc_metrics = preprocess_qc_metrics,
            irma_typing = irma_typing,
            irma_assembly_qc_metrics = irma_assembly_qc_metrics,
            python_script = summary_py,
            project_name = project_name,
            analysis_date = analysis_date,
            workbook_path = workbook_path
    }

    call transfer.transfer_assembly_summary_wdl as summary_transfer {
        input:
            sequencing_results_csv = summary.sequencing_results_csv,
            bucket_path = bucket_path
    }

    output {

        File sequencing_results_csv = summary.sequencing_results_csv
        String transfer_date = summary_transfer.transfer_date

    }

}