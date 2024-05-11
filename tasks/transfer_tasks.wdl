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
        File? fastqc2_html_raw
        File? fastqc2_zip_raw

        File seqyclean_summary

        File preprocess_qc_metrics

        File fastqc1_html_cleaned
        File fastqc1_zip_cleaned
        File? fastqc2_html_cleaned
        File? fastqc2_zip_cleaned

        # irma assembly outputs
        File irma_assembled_gene_segments_csv
        Array[File] irma_assemblies
        Array[File] irma_bam_files
        Array[File] irma_vcfs

        # ivar_assemblies (from irma assembler) and samtools tools
        Array[File]? irma_sorted_bams
        Array[File]? irma_ivar_assemblies
        Array[File]? irma_ivar_outputs

        # post assembly qc outputs
        File? irma_qc_metrics

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
        gsutil -m cp ~{irma_assembled_gene_segments_csv} ~{out_path}/irma/~{sample_name}/
        gsutil -m cp ~{sep = " " irma_assemblies} ~{out_path}/irma/~{sample_name}/assemblies/
        gsutil -m cp ~{sep = " " irma_bam_files} ~{out_path}/irma/~{sample_name}/bam_files/
        gsutil -m cp ~{sep = " " irma_vcfs} ~{out_path}/irma/~{sample_name}/vcf_files/

        # transfer ivar assemblies and sorted bams 
        gsutil -m cp ~{sep = " " irma_sorted_bams} ~{out_path}/irma/~{sample_name}/sorted_bam_files/
         gsutil -m cp ~{sep = " " irma_ivar_assemblies} ~{out_path}/irma/~{sample_name}/irma_ivar_consensus/
         gsutil -m cp ~{sep = " " irma_ivar_outputs} ~{out_path}/irma/~{sample_name}/irma_ivar_outputs/

        # transfer post assembly qc
        gsutil -m cp ~{irma_qc_metrics} ~{out_path}/irma/~{sample_name}/

        # transfer nextclade
        gustil -m cp ~{na_nextclade_json} ~{out_path}/irma/~{sample_name}/nextclade/
        gustil -m cp ~{na_nextclade_tsv} ~{out_path}/irma/~{sample_name}/nextclade/
        gustil -m cp ~{na_translation_fasta} ~{out_path}/irma/~{sample_name}/nextclade/
        gustil -m cp ~{ha_nextclade_json} ~{out_path}/irma/~{sample_name}/nextclade/
        gustil -m cp ~{ha_nextclade_tsv} ~{out_path}/irma/~{sample_name}/nextclade/
        gustil -m cp ~{ha_HA1_translation_fasta} ~{out_path}/irma/~{sample_name}/nextclade/
        gustil -m cp ~{ha_HA2_translation_fasta} ~{out_path}/irma/~{sample_name}/nextclade/
        gustil -m cp ~{ha_SigPep_translation_fasta} ~{out_path}/irma/~{sample_name}/nextclade/
        
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
    }

    String out_path = sub(bucket_path, "/$", "") # fix if have a / at end


    command <<< 
         gsutil -m cp ~{sequencing_results_csv} ~{out_path}/summary_files/

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


