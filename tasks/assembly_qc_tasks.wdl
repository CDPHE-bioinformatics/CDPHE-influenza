version 1.0

# define structure
struct VersionInfo {
  String software
  String docker
  String version
}

task calc_bam_stats_samtools {
    meta {
        description: "output sorted bams and use samtools to calc depth metrics"
    }

    input {
        File? bam_file
        String sample_name
    }

    String docker = "staphb/samtools:1.10"

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
        echo "sample_name,file_name,segment_name,gene_name,description,value" > bam_stats.csv
        echo "~{sample_name},${sorted_bam},${segment_name},${gene_name},num_mapped_reads,$(cat num_mapped_reads.txt)" >> bam_stats.csv
        echo "~{sample_name},${sorted_bam},${segment_name},${gene_name},mean_depth,$(cat mean_depth.txt)" >> bam_stats.csv

        samtools --version | awk '/samtools/ {print $2}' | tee VERSION

    >>>

    output {
        File bam_stats_csv = "bam_stats.csv"
        File sorted_bam = select_first(glob("*.sorted.bam"))

        VersionInfo samtools_version_info = object{
            software: "samtools",
            docker: docker,
            version: read_string("VERSION")
        }
    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker

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

task concat_assembly_qc_metrics{
    meta{
        description: "concatenate all post assembly qc metrics (depth and coverage) into a single file"
    }

    input{
        File python_script
        String sample_name
        Array[File] bam_stats_csv_array
        Array[File] percent_coverage_csv_array
        
    }
    


    command <<<

    python ~{python_script} \
        --sample_name "~{sample_name}" \
        --bam_stats_csv_list "~{sep= " " bam_stats_csv_array}" \
        --percent_coverage_csv_list "~{sep = " " percent_coverage_csv_array}"

    >>>

    output{
        File? assembly_qc_metrics_summary = "~{sample_name}_assembly_qc_metrics.csv"
    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }

}


task make_multifasta {
    meta {
        description: "create a mulitfasta of the consensus assemblies"
    }

    input {
        Array[File] fasta_array
        String sample_name
    }

    command <<<
        
        # Concatenate all the FASTA files into a single file
        cat ~{sep=' ' fasta_array} > ~{sample_name}_ivar.fasta


    >>>

    output {
        File multifasta = "~{sample_name}_ivar.fasta"

    }
        runtime {
        docker: "theiagen/utility:1.0"
        memory: "16 GiB"
        cpu: 4
        disks: "local-disk 50 SSD"
        preemptible: 0
    }
}