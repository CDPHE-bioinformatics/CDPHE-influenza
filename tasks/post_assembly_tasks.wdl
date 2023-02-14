version 1.0

task samtools_mapped_reads {
    meta {
        description: "use samtools to calc depth metrics"
    }

    input {
        File bam_file
    }

    command <<<

    # get segment name from bam file name
    base_name=$(basename ~{bam_file})
    sample_id=$(basename ~{bam_file} | cut -d "_" -f 1)
    segment_name=$(basename ~{bam_file} | cut -d "." -f 1 | cut -d "_" -f 2-)
    gene_name=$(basename ~{bam_file} | cut -d "." -f 1 | cut -d "_" -f 3)

    # count mapped flu reads
    samtools view -c -F 260 ~{bam_file} > num_mapped_reads.txt
    samtools coverage ~{bam_file} | tail -1 | cut -f 7 > mean_depth.txt

    echo "sample_id,base_name,segment_name,gene_name,description,value" > bam_results.csv
    echo "${sample_id},${base_name},${segment_name},${gene_name},num_mapped_reads,$(cat num_mapped_reads.txt)" >> bam_results.csv
    echo "${sample_id},${base_name},${segment_name},${gene_name},mean_depth,$(cat mean_depth.txt)" >> bam_results.csv


    >>>

    output {
        File bam_results = "bam_results.csv"
        # String mapped_reads = read_string("num_mapped_reads")
        # Striing mean_depth = read_string('mean_depth')
    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/samtools:1.10"

    }
}

task calc_percent_coverage{

    meta{
        description: "use python script to calculate the percent coverage of each gene segment"
    }
    
    input{
        File python_script
        File fasta_file

    }

    command <<<

    # get segment name from fasta file name
    # segment_name=$(basename ~{fasta_file} | cut -d "." -f 1 | cut -d "_" -f 2-)


    python ~{python_script} --fasta_file ~{fasta_file}
    
    >>>

    output{
        File perc_cov_results = "perc_cov_results.csv"
    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }

}

task concat_post_qc_metrics{
    meta{
        description:" pull the alignment (depth results from \samtools results) and the assembly (percent coverage \resutls from the biopython script) into a single file \for the sample"
    }

    input{
        File python_script
        String sample_id
        Array[File] bam_results_array
        Array[File] per_cov_results_array
        Array[File] ivar_parameters

    }
    File ivar_parameters_file = select_first(ivar_parameters)

    command <<<

    python ~{python_script} \
        --sample_id ~{sample_id} \
        --bam_results ~{write_lines(bam_results_array)} \
        --per_cov_results ~{write_lines(per_cov_results_array)} \
        --ivar_parameters ~{ivar_parameters_file}

    >>>

    output{
        File? qc_metrics_summary = "~{sample_id}_qc_metrics.csv"

    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }

}