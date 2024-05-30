version 1.0

task summary {
    meta{
        description: "concat all important metrics across samples into a single file; and generate version capture files"
    }

    input{
        Array[String] sample_name
        Array[File] preprocess_qc_metrics
        Array[File] irma_typing
        Array[File] assembly_qc_metrics
        Array[File] nextclade_na_tsv
        Array[File] nextclade_ha_tsv
        Array[File] version_capture_file
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
        --nextclade_na_tsv ~{write_lines(nextclade_na_tsv)} \
        --nextclade_ha_tsv ~{write_lines(nextclade_ha_tsv)} \
        --version_capture_file ~{write_lines(version_capture_file)} \
        --workflow_verion "~{workflow_version}" \
        --project_name "~{project_name}" \
        --analysis_date "~{analysis_date}" 

    >>>

    output {
        File sequencing_results_csv = "~{project_name}_sequencing_results.csv" 
        File version_capture_influenza_assembly_csv = "version_capture_influenza_assembly_~{project_name}_~{workflow_version}.csv"
        File version_capture_influenza_assembly_summary_csv = "version_capture_influenza_assembly_summary_~{project_name}_~{workflow_version}.csv"
    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}