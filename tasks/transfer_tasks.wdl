version 1.0

task transfer_assembly_wdl{
    meta{
        description: "transfer fastqc, fasta, bam, qc_metrics files to gcp"
    }

    input{
        String sample_name
        String bucket_path

        # pre-process outputs
        File fastqc1_html_raw
        File fastqc1_zip_raw
        File fastqc2_html_raw
        File fastqc2_zip_raw

        File seqyclean_summary

        File preprocess_qc_metrics

        File fastqc1_html_cleaned
        File fastqc1_zip_cleaned
        File fastqc2_html_cleaned
        File fastqc2_zip_cleaned

        # irma assembly outputs
        File irma_assembled_gene_segments_csv
        File? irma_all_assembled_segments_fasta
        Array[File]? irma_fasta_array
        Array[File]? irma_bam_array
        Array[File]? irma_vcf_array

        # ivar & sorted bams
        Array[File]? ivar_fasta_array
        Array[File]? sorted_bam_array

        # post assembly qc outputs
        File? assembly_qc_metrics

        # nextclade
        File? na_nextclade_json
        File? na_nextclade_tsv
        File? na_translation_fasta

        File? ha_nextclade_json
        File? ha_nextclade_tsv
        File? ha_HA1_translation_fasta
        File? ha_HA2_translation_fasta
        File? ha_SigPep_translation_fasta


    }
    
    String out_path = sub(bucket_path, "/$", "") # fix if have a / at end

    command <<<
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

        # transfer preprocess qc metrics summary
        gsutil -m cp ~{preprocess_qc_metrics} ~{out_path}/preprocess_qc_metrics/

        # transfer irma
        gsutil -m cp ~{irma_assembled_gene_segments_csv} ~{out_path}/irma_assembly_results/
        gsutil -m cp ~{irma_all_assembled_segments_fasta} ~{out_path}/irma_assembly_multifastas/
        gsutil -m cp ~{sep = " " irma_fasta_array} ~{out_path}/irma_assemblies/~{sample_name}/
        gsutil -m cp ~{sep = " " irma_bam_array} ~{out_path}/irma_alignments/~{sample_name}/
        gsutil -m cp ~{sep = " " irma_vcf_array} ~{out_path}/irma_vcfs/~{sample_name}/

        # transfer ivar assemblies and sorted bams 
        gsutil -m cp ~{sep = " " sorted_bam_array} ~{out_path}/sorted_bams/~{sample_name}/
        gsutil -m cp ~{sep = " " ivar_fasta_array} ~{out_path}/ivar_assemblies/~{sample_name}
         

        # transfer post assembly qc
        gsutil -m cp ~assembly_qc_metrics ~{out_path}/assembly_qc_metrics/

        # transfer nextclade
        gustil -m cp ~{na_nextclade_json} ~{out_path}/nextclade_out/~{sample_name}/
        gustil -m cp ~{na_nextclade_tsv} ~{out_path}/nextclade_out/~{sample_name}/
        gustil -m cp ~{na_translation_fasta} ~{out_path}/nextclade_out/~{sample_name}/
        gustil -m cp ~{ha_nextclade_json} ~{out_path}/nextclade_out/~{sample_name}/
        gustil -m cp ~{ha_nextclade_tsv} ~{out_path}/nextclade_out/~{sample_name}/
        gustil -m cp ~{ha_HA1_translation_fasta} ~{out_path}/nextclade_out/~{sample_name}/
        gustil -m cp ~{ha_HA2_translation_fasta} ~{out_path}/nextclade_out/~{sample_name}/
        gustil -m cp ~{ha_SigPep_translation_fasta} ~{out_path}/nextclade_out/~{sample_name}/
        


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

        String bucket_path
        File sequencing_results_csv
        File version_capture_influenza_assembly_csv
        File version_capture_influenza_assembly_summary_csv
    }

    String out_path = sub(bucket_path, "/$", "") # fix if have a / at end


    command <<< 
        gsutil -m cp ~{sequencing_results_csv} ~{out_path}/summary_files/
        gsutil -m cp ~{version_capture_influenza_assembly_csv} ~{out_path}/summary_files/
        gsutil -m cp ~{version_capture_influenza_assembly_summary_csv} ~{out_path}/summary_files/
        
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


