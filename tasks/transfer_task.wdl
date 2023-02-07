version 1.0

task transfer_assembly_wdl{
    meta{
        description: "transfer fastqc, fasta, bam, qc_metrics files to gcp"
    }

    input{
        String sample_id
        String bucket_path
        String out_path = sub(bucket_path, "/$", "") # fix if have a / at end

        # pre-process outputs
        File fastqc1_raw_html
        File fastqc1_raw_zip
        File fastqc2_raw_html
        File fastqc2_raw_zip

        File seqyclean_summary

        File fastqc1_cleaned_html
        File fastqc1_cleaned_zip
        File fastqc2_cleaned_html
        File fastqc2_cleaned_zip

        # irma assembly outputs
        File[Array] irma_assemblies
        File[Array] irma_bam_files
        File[Array] irma_vcfs

        # post assembly qc outputs
        File irma_qc_metrics

    }

    command <<<
        # transfer fastqc raw
        gsutil -m cp ~{fastqc1_raw_html} ~{out_path}/fastqc_raw/
        gsutil -m cp ~{fastqc1_raw_zip} ~{out_path}/fastqc_raw/
        gsutil -m cp ~{fastqc2_raw_html} ~{out_path}/fastqc_raw/
        gsutil -m cp ~{fastqc2_raw_zip} ~{out_path}/fastqc_raw/

        # transfter seqyclean
        gsutil -m cp ~{seqyclean_summary} ~{out_path}/seqyclean/

        # transfer fastqc clean
        gsutil -m cp ~{fastqc1_cleaned_html} ~{out_path}/fastqc_cleaned/
        gsutil -m cp ~{fastqc1_cleaned_zip} ~{out_path}/fastqc_cleaned/
        gsutil -m cp ~{fastqc2_cleaned_html} ~{out_path}/fastqc_cleaned/
        gsutil -m cp ~{fastqc2_cleaned_zip} ~{out_path}/fastqc_cleaned/

        # transfer irma
        gsutil -m cp ~{sep = " " irma_assemblies} ~{out_path}/irma/~{sample_name}/assemblies/
        gsutil -m cp ~{sep = " " irma_bam_files} ~{out_path}/irma/~{sample_name}/bam_files/
        gsutil -m cp ~{sep = " " irma_vcfs} ~{out_path}/irma/~{sample_name}/vcf_files/

        # transfer post assemlby qc
        gsutil -m cp ~{sep = " " irma_qc_metrics} ~{out_path}/irma/~{sample_name}/

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