version 1.0

# define structure
# struct VersionInfo {
#   String software
#   String docker
#   String version
# }
import "../tasks/capture_version_tasks.wdl" as capture_version

task calc_bam_stats_samtools {
    meta {
        description: "output sorted bams and use samtools to calc depth metrics"
    }

    input {
        File? bam_file
        String sample_name
        String segment_name
        String base_name
    }

    String sorted_bam_fn = "~{base_name}.sorted.bam"
    String sorted_bai_fn = "~{base_name}.sorted.bam.bai"
    String sam_coverage_fn = '~{base_name}_coverage.txt'
    String sam_depth_fn = '~{base_name}_depth.txt'
    String docker = "staphb/samtools:1.10"

    command <<<

        echo "base_name"
        echo ~{base_name}
        echo "segment_name"
        echo ~{segment_name}
        echo "bam_file"
        echo ~{bam_file}
        echo "sorted_bam_fn"
        echo ~{sorted_bam_fn}
        echo "sorted_bai_fn"
        echo ~{sorted_bai_fn}
        echo "sam_coverage"
        echo ~{sam_coverage_fn}
        echo "sam_depth"
        echo ~{sam_depth_fn}

        # create sorted bam file
        samtools sort ~{bam_file} -o ~{sorted_bam_fn}
        samtools index ~{sorted_bam_fn}

        # get coverage and depth files as outputs 
        # will use the depth file to calcuate the percent coverage at 30x.

        samtools coverage -o ~{sam_coverage_fn} ~{bam_file}
        samtools depth -a -o ~{sam_depth_fn} ~{bam_file}


        # use sorted bam file to get number mapped reads and mean depth
        samtools view -c -F 260 ~{sorted_bam_fn} > num_mapped_reads.txt
        samtools coverage ~{sorted_bam_fn} | tail -1 | cut -f 7 > mean_depth.txt

        # create output file
        # why is the file set up like this? why not have num_mapped_reads and mean_depth be their own columns?
        # it has to do with the way the summary script formats headers
        # it uses {segment_name}_{description} -- "HA- num_mapped_reads"
        # so it was just easier to loop through the data with a description column and segment name column
        
        echo "sample_name,segment_name,description,value" > bam_stats.csv
        echo "~{sample_name},~{segment_name},num_mapped_reads,$(cat num_mapped_reads.txt)" >> bam_stats.csv
        echo "~{sample_name},~{segment_name},mean_depth,$(cat mean_depth.txt)" >> bam_stats.csv

        samtools --version | awk '/samtools/ {print $2}' | tee VERSION
    >>>

    output {
        File bam_stats_csv = "bam_stats.csv"
        File sam_coverage = sam_coverage_fn
        File sam_depth = sam_depth_fn
        File sorted_bam = sorted_bam_fn
        File sorted_bai = sorted_bai_fn

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
        String base_name
        String segment
    }

    command <<<
        python ~{python_script}  --fasta_file "~{fasta_file}" \
        --sample_name "~{sample_name}" \
        --base_name "~{base_name}" \
        --segment "~{segment}"
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
        File irma_read_counts_txt
        Array[File] bam_stats_csv_array
        Array[File] percent_coverage_csv_array
        
    }
    
    command <<<
        python ~{python_script} \
            --sample_name "~{sample_name}" \
            --irma_read_counts_txt "~{irma_read_counts_txt}"
            --bam_stats_csv_list "~{sep= " " bam_stats_csv_array}" \
            --percent_coverage_csv_list "~{sep = " " percent_coverage_csv_array}"
    >>>

    output{
        File assembly_qc_metrics_summary = "~{sample_name}_assembly_qc_metrics.csv"
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
        # Is this different than what IRMA outputs?
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