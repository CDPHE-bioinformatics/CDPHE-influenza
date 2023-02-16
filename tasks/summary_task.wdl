version 1.0

task summary {
    meta{
        description: "concat all important metrics across samples into a single file"
    }

    input{
        Array[String] sample_id
        Array[File] preprocess_qc_metrics
        Array[File] irma_typing
        Array[File] irma_assembly_qc_metrics
        File python_script
        String project_name
        String run_date
    }

    command <<<

    python ~{python_script} \
        --sample_id ~{write_lines(sample_id)} \
        --preprocess_qc_metrics ~{write_lines(preprocess_qc_metrics)} \
        --irma_typing ~{write_lines(irma_typing)} \
        --irma_qc_metrics ~{write_lines(irma_assembly_qc_metrics)} \
        --project_name ~{project_name} \
        --transfer_date ~{run_date}

    >>>

    output {
        File sequencing_results_csv = "~{project_name}_sequencing_results.csv" 
    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}