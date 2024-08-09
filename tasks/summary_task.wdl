version 1.0

task results_summary {
    meta{
        description: "concat all important metrics across samples into a single file"
    }

    input{
        Array[String] sample_name
        Array[File] preprocess_qc_metrics
        Array[File] irma_typing
        Array[File] assembly_qc_metrics
        Array[File] na_nextclade_tsv
        Array[File] ha_nextclade_tsv
        String workflow_version
        String analysis_date
        File python_script
        String project_name

    }

    command <<<

    python ~{python_script} \
        --sample_name ~{write_lines(sample_name)} \
        --preprocess_qc_metrics ~{write_lines(preprocess_qc_metrics)} \
        --irma_typing ~{write_lines(irma_typing)} \
        --assembly_qc_metrics ~{write_lines(assembly_qc_metrics)} \
        --na_nextclade_tsv ~{write_lines(na_nextclade_tsv)} \
        --ha_nextclade_tsv ~{write_lines(ha_nextclade_tsv)} \
        --workflow_version ~{workflow_version}\
        --project_name "~{project_name}" \
        --analysis_date "~{analysis_date}" 

    >>>

    output {
        File sequencing_results_csv = "~{project_name}_sequencing_results_~{workflow_version}.csv" 

    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

task capture_version_summary {
    meta{
        description: "generate version capture files"
    }

    input{
        String workflow_name
        String workflow_version
        String analysis_date
        File python_script
        String project_name

    }

    command <<<

    python ~{python_script} \
        --workflow_version "~{workflow_version}" \
        --workflow_name "~{workflow_name}" \
        --project_name "~{project_name}" \
        --analysis_date "~{analysis_date}" 

    >>>

    output {
        
        File version_capture_influenza_assembly_summary_csv = "version_capture_~{workflow_name}_~{project_name}.csv"
    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}