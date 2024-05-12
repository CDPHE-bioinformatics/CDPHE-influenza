version 1.0

task samtools_mapped_reads {
    meta {
        description: "output sorted bams and use samtools to calc depth metrics"
    }

    input {
        File? bam_file
        String sample_name
    }

    command <<<

        # create name for sorted bam file
        prefix=$(basename ~{bam_file} | cut -d "." -f 1)
        sorted_bam=$(echo ${prefix}.sorted.bam)
        
        # # pull sample id, segment name, and gene name from original ba file
        # these are the same because I don't want to input the type and subtpye into
        # and downstream that info is not pulled from this file
        segment_name=$(echo "${prefix/~{sample_name}/*}" | cut -d "_" -f 2)
        gene_name=$(echo "${prefix/~{sample_name}/*}" | cut -d "_" -f 2)

        # create sorted bam file
        samtools sort ~{bam_file} -o ${sorted_bam}

        # use sorted bam file to get number mapped reads and mean depth
        samtools view -c -F 260 ${sorted_bam} > num_mapped_reads.txt
        samtools coverage ${sorted_bam} | tail -1 | cut -f 7 > mean_depth.txt

        # create output file
        # echo "sample_name,file_name,segment_name,gene_name,description,value" > mapped_reads.csv
        # echo "~{sample_name},${sorted_bam},${segment_name},${gene_name},num_mapped_reads,$(cat num_mapped_reads.txt)" >> mapped_reads.csv
        # echo "~{sample_name},${sorted_bam},${segment_name},${gene_name},mean_depth,$(cat mean_depth.txt)" >> mapped_reads.csv

        echo "sample_name,file_name,segment_name,gene_name,description,value" > mapped_reads.csv
        echo "~{sample_name},${sorted_bam},${segment_name},${gene_name},num_mapped_reads,$(cat num_mapped_reads.txt)" >> mapped_reads.csv
        echo "~{sample_name},${sorted_bam},${segment_name},${gene_name},mean_depth,$(cat mean_depth.txt)" >> mapped_reads.csv

    >>>

    output {
        File mapped_reads_csv = "mapped_reads.csv"
        File sorted_bam = select_first(glob("*.sorted.bam"))
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
        String sample_name

    }

    command <<<

    python ~{python_script} \
    --fasta_file "~{fasta_file}" \
    --sample_name "~{sample_name}"
    
    >>>

    output{
        File percent_coverage_csv = "percent_coverage_results.csv"
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
        description: "concatenate all post assembly qc metrics (depth and coverage) into a single file"
    }

    input{
        File python_script
        String sample_name
        Array[File?] mapped_reads_csv_array
        Array[File?] percent_coverage_csv_array
        
    }
    


    command <<<

    python ~{python_script} \
        --sample_name "~{sample_name}" \
        --mapped_reads_csv_array "[~{sep = ", " mapped_reads_csv_array}]" \
        --percent_coverage_csv_array "[~sep = "," percent_coverage_csv_array}]"

    >>>

    output{
        File? qc_metrics_summary = "~{sample_name}_assembly_qc_metrics.csv"
    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }

}