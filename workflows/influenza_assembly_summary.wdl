version 1.0

import "../tasks/summary_task.wdl" as summary 
import "../tasks/transfer_task.wdl" as transfer

workflow influenza_assembly_summary{
    
    input {
        Array[String] sample_id
        Array[String] project_name_array
        Array[File] preprocess_qc_metrics
        Array[File] irma_typing
        Array[File] irma_qc_metrics
        Array[String] bucket_path_array

        File summary_py
    }

    String bucket_path = select_first(bucket_path_array)
    String project_name = select_first(project_name_array)

    call summary.summary as summary {
        input:
            sample_id = sample_id,
            preprocess_qc_metrics = preprocess_qc_metrics,
            irma_typing = irma_typing,
            irma_qc_metrics = irma_qc_metrics,
            python_script = summary_py,
            project_name = project_name
    }

    call transfer.summary_transfer as summary_transfer {
        input:
            summary_file = summary_file,
            bucket_path = bucket_path
    }

    output {

        File assembly_results_csv = summary.assembly_results_csv
        String transfer_date = summary_transfer.transfer_date

    }

}