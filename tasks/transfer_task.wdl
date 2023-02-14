version 1.0

task transfer_assembly_wdl{
    meta{
        description: "transfer fastqc, fasta, bam, qc_metrics files to gcp"
    }

    input{
        String sample_id
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
        Array[File] irma_assemblies
        Array[File] irma_bam_files
        Array[File] irma_vcfs

        # ivar_assemblies (from irma assembler)
        Array[File]? ivar_assemblies
        Array[File]? ivar_outputs

        # post assembly qc outputs
        File? irma_qc_metrics

    }
    
    String out_path = sub(~{bucket_path}, "/$", "") # fix if have a / at end

    command <<<
        # transfer fastqc raw
        gsutil -m cp ~{fastqc1_html_raw} ~{out_path}/prefastqc_raw/
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
        gsutil -m cp ~{sep = " " irma_assemblies} ~{out_path}/irma/~{sample_id}/assemblies/
        gsutil -m cp ~{sep = " " irma_bam_files} ~{out_path}/irma/~{sample_id}/bam_files/
        gsutil -m cp ~{sep = " " irma_vcfs} ~{out_path}/irma/~{sample_id}/vcf_files/

        # transfer ivar assemblies
         gsutil -m cp ~{sep = " " irma_assemblies} ~{out_path}/irma/~{sample_id}/ivar_consensus/
         gsutil -m cp ~{sep = " " irma_assemblies} ~{out_path}/irma/~{sample_id}/ivar_outputs/

        # transfer post assemlby qc
        gsutil -m cp ~{irma_qc_metrics} ~{out_path}/irma/~{sample_id}/

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

