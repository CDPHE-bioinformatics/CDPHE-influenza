version 1.0

task transfer_assembly_wdl{
    meta{
        description: "transfer fastqc, fasta, bam, qc_metrics files to gcp"
    }

    input{
        String workflow_version
        String sample_name
        String bucket_path

        # pre-process outputs
        File fastqc1_html_raw
        File fastqc1_zip_raw
        File fastqc2_html_raw
        File fastqc2_zip_raw

        File seqyclean_summary

        File fastqc1_html_cleaned
        File fastqc1_zip_cleaned
        File fastqc2_html_cleaned
        File fastqc2_zip_cleaned

        # irma assembly outputs
        File? irma_read_counts
        File? irma_run_info
        File irma_assembled_gene_segments_csv
        File? irma_multifasta
        Array[File]? irma_fasta_array
        Array[File]? irma_bam_array
        Array[File]? irma_vcf_array

        # ivar & sorted bams
        Array[File]? sam_coverage_array
        Array[File]? sam_depth_array
        Array[File]? sorted_bam_array
        Array[File]? sorted_bai_array

        # nextclade
        # Array[File]? nextclade_json_array 
        # Array[File]? nextclade_tsv_array 
        File? nextclade_HA_json
        File? nextclade_NA_json
        File? nextclade_HA_tsv
        File? nextclade_NA_tsv
        File? nextclade_SigPep_translation_fasta
        File? nextclade_HA1_translation_fasta
        File? nextclade_HA2_translation_fasta
        File? nextclade_HA_translation_fasta
        File? nextclade_NA_translation_fasta

        # version
        File? version_capture_file


    }
    
    String out_path1 = sub(bucket_path, "/$", "") # fix if have a / at end
    String out_path = "~{out_path1}/~{workflow_version}"

    command <<<
        # transfer version capture file
        gsutil -m cp ~{version_capture_file} ~{out_path}/version_capture/

        # transfer fastqc raw
        gsutil -m cp ~{fastqc1_html_raw} ~{out_path}/fastqc_raw/
        gsutil -m cp ~{fastqc1_zip_raw} ~{out_path}/fastqc_raw/
        gsutil -m cp ~{fastqc2_html_raw} ~{out_path}/fastqc_raw/
        gsutil -m cp ~{fastqc2_zip_raw} ~{out_path}/fastqc_raw/

        # transfter seqyclean
        gsutil -m cp ~{seqyclean_summary} ~{out_path}/seqyclean/

        # transfer fastqc clean
        gsutil -m cp ~{fastqc1_html_cleaned} ~{out_path}/fastqc_cleaned/
        gsutil -m cp ~{fastqc1_zip_cleaned} ~{out_path}/fastqc_cleaned/
        gsutil -m cp ~{fastqc2_html_cleaned} ~{out_path}/fastqc_cleaned/
        gsutil -m cp ~{fastqc2_zip_cleaned} ~{out_path}/fastqc_cleaned/


        # transfer irma
        gsutil -m cp ~{irma_read_counts} ~{out_path}/irma_logs/
        gsutil -m cp ~{irma_run_info} ~{out_path}/irma_logs/
        gsutil -m cp ~{irma_assembled_gene_segments_csv} ~{out_path}/irma_assembly_results/
        gsutil -m cp ~{irma_multifasta} ~{out_path}/irma_assembly_multifastas/
        gsutil -m cp ~{sep = " " irma_fasta_array} ~{out_path}/irma_assemblies/~{sample_name}/
        gsutil -m cp ~{sep = " " irma_bam_array} ~{out_path}/irma_alignments/~{sample_name}/
        gsutil -m cp ~{sep = " " irma_vcf_array} ~{out_path}/irma_vcfs/~{sample_name}/

         # transfer sorted bams 
        gsutil -m cp ~{sep = " " sorted_bam_array} ~{out_path}/sorted_bams/~{sample_name}/
        gsutil -m cp ~{sep = " " sorted_bai_array} ~{out_path}/sorted_bams/~{sample_name}/
        gsutil -m cp ~{sep = " " sam_coverage_array} ~{out_path}/bam_stats/~{sample_name}/
        gsutil -m cp ~{sep = " " sam_depth_array} ~{out_path}/bam_stats/~{sample_name}/

        # transfer nextclade
        gsutil -m cp ~{nextclade_HA_json} ~{out_path}/nextclade_out/~{sample_name}/
        gsutil -m cp ~{nextclade_NA_json} ~{out_path}/nextclade_out/~{sample_name}/
        gsutil -m cp ~{nextclade_HA_tsv} ~{out_path}/nextclade_out/~{sample_name}/
        gsutil -m cp ~{nextclade_NA_tsv} ~{out_path}/nextclade_out/~{sample_name}/
        gsutil -m cp ~{nextclade_HA_translation_fasta} ~{out_path}/nextclade_out/~{sample_name}/
        gsutil -m cp ~{nextclade_HA1_translation_fasta} ~{out_path}/nextclade_out/~{sample_name}/
        gsutil -m cp ~{nextclade_HA2_translation_fasta} ~{out_path}/nextclade_out/~{sample_name}/
        gsutil -m cp ~{nextclade_SigPep_translation_fasta} ~{out_path}/nextclade_out/~{sample_name}/
        gsutil -m cp ~{nextclade_NA_translation_fasta} ~{out_path}/nextclade_out/~{sample_name}/

        # transfer date
        transferdate=`date`
        echo $transferdate | tee TRANSFERDATE

    >>>
    output {
        String transfer_date = read_string("TRANSFERDATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "16 GiB"
        cpu: 4
        disks: "local-disk 50 SSD"
        preemptible: 0
    }

}

task transfer_assembly_summary_wdl {
    meta {
        description: ""
    }

    input {

        String workflow_version
        String bucket_path
        File sequencing_results_csv
        File version_capture_influenza_assembly_summary_csv
    }

    String out_path1 = sub(bucket_path, "/$", "") # fix if have a / at end
    String out_path = "~{out_path1}/~{workflow_version}"


    command <<< 
        gsutil -m cp ~{sequencing_results_csv} ~{out_path}/summary_results/
        gsutil -m cp ~{version_capture_influenza_assembly_summary_csv} ~{out_path}/version_capture/
        
        # transfer date
        transferdate=`date`
        echo $transferdate | tee TRANSFERDATE
    >>>

    output {
        String transfer_date = read_string("TRANSFERDATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "16 GiB"
        cpu: 4
        disks: "local-disk 50 SSD"
        preemptible: 0
    }

    
}


