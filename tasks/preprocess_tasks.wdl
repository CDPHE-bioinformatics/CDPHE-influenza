version 1.0

# define structure
import "../tasks/capture_version_tasks.wdl" as capture_version
# struct VersionInfo {
#   String software
#   String docker
#   String version
# }

# begin tasks
task assess_quality_fastqc {
    meta {
      description: "this task uses fastqc to evaluate the quality of reads in the fastq files. The output is an html file with various quality statistics, from which the task pulls the number of reads. Modified from Theiagen Genomics- public helath viral genomics."
    }

    input {
        String sample_name
        File fastq_R1
        File fastq_R2
        
    }
    String docker = 'staphb/fastqc:0.11.9'

    command <<<
        # grab version
        fastqc --version | tee VERSION

        # get basename of fastq file
        fastq_R1_name=$(basename ~{fastq_R1} | cut -d "." -f 1 | cut -d "." -f 1)
        fastq_R2_name=$(basename ~{fastq_R2} | cut -d "." -f 1 | cut -d "." -f 1)

        # run fastqc
        fastqc --outdir $PWD ~{fastq_R1} ~{fastq_R2}

        # pull some info from the zip file regarding number of reads and read length
        unzip -p ${fastq_R1_name}_fastqc.zip */fastqc_data.txt | grep "Sequence length" | cut -f 2 | tee READ1_LEN
        unzip -p ${fastq_R2_name}_fastqc.zip */fastqc_data.txt | grep "Sequence length" | cut -f 2 | tee READ2_LEN

        unzip -p ${fastq_R1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ1_SEQS
        unzip -p ${fastq_R2_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ2_SEQS

        READ1_SEQS=$(unzip -p ${fastq_R1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )
        READ2_SEQS=$(unzip -p ${fastq_R2_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )

        if [ $READ1_SEQS == $READ2_SEQS ]; then
            read_pairs=$READ1_SEQS
        else
            read_pairs="Uneven pairs: R1=$READ1_SEQS, R2=$READ2_SEQS"
        fi

        echo $read_pairs | tee READ_PAIRS

        # # calculate the total number of reads(i.e. READ_PAIRS * 2 for paired data; else *1 for SE)
        # total_reads=$((2*READ_PAIRS))
        # echo $total_reads | tee TOTAL_READS

        # rename files 
        if [ ${fastq_R1_name} != ~{sample_name}_R1 ]; then 
            mv ${fastq_R1_name}_fastqc.html ~{sample_name}_R1_fastqc.html
            mv ${fastq_R1_name}_fastqc.zip ~{sample_name}_R1_fastqc.zip
            mv ${fastq_R2_name}_fastqc.html ~{sample_name}_R2_fastqc.html
            mv ${fastq_R2_name}_fastqc.zip ~{sample_name}_R2_fastqc.zip
        fi
        
        

    >>>
    output {
        File fastqc1_html = "~{sample_name}_R1_fastqc.html"
        File fastqc1_zip = "~{sample_name}_R1_fastqc.zip"
        File fastqc2_html = "~{sample_name}_R2_fastqc.html"
        File fastqc2_zip = "~{sample_name}_R2_fastqc.zip"
 
        # these outputs will go into the concatenated preprocess qc metrics file
        Int total_reads_R1 = read_string("READ1_SEQS")
        Int total_reads_R2 = read_string("READ2_SEQS")
        
        String read_length_R1 = read_string('READ1_LEN')
        String read_length_R2 = read_string('READ2_LEN')
        
        String read_pairs = read_string("READ_PAIRS")

        VersionInfo fastqc_version_info = object{
            software: "fastqc",
            docker: docker,
            version: read_string("VERSION")
        }
    } 

    runtime {
      docker: docker
      memory: "1 GiB"
      cpu: 2
      disks: "local-disk 100 SSD"
      preemptible: 0
      maxRetries: 3
    }
}

task filter_reads_seqyclean {
    meta{
        description : "This task uses seqyclean to remove containments and adapaters and then filters on min len and quality."
    }
    input {
        File contam_fasta
        String sample_name
        File fastq_R1
        File fastq_R2
        
    }

    String docker = "staphb/seqyclean:1.10.09"

    command <<<

        # run seqyclean
        seqyclean -minlen 70 -qual 30 30 -gz -1 ~{fastq_R1} -2 ~{fastq_R2} -c ~{contam_fasta} -o ~{sample_name}_clean
       
        # pull version out of the summary file
        awk 'NR==2 {print $1}' ~{sample_name}_clean_SummaryStatistics.tsv | tee VERSION
    >>>

    output {
        # String seqyclean_version = read_string("VERSION")
        # String seqyclean_docker = docker
        File fastq_R1_cleaned = "${sample_name}_clean_PE1.fastq.gz"
        File fastq_R2_cleaned = "${sample_name}_clean_PE2.fastq.gz"
        File seqyclean_summary = "${sample_name}_clean_SummaryStatistics.tsv"

        VersionInfo seqyclean_version_info = object{
            software: "seqyclean",
            docker: docker,
            version: read_string("VERSION")
        }
    }
    runtime {
        docker: docker
        memory: "2 GiB"
        cpu: 2
        disks: "local-disk 100 SSD"
        preemptible: 0
        maxRetries: 3
    }
}

task concat_preprocess_qc_metrics {
    meta {
        description: "creates a single csv file with all preprocess qc metrics"
    }

    input {
        File python_script
        String sample_name

        Int total_reads_R1_raw
        Int total_reads_R2_raw
        String read_length_R1_raw 
        String read_length_R2_raw 
        String read_pairs_raw


        Int total_reads_R1_cleaned
        Int total_reads_R2_cleaned
        String read_length_R1_cleaned
        String read_length_R2_cleaned
        String read_pairs_cleaned

        
    }

    command <<<

        python ~{python_script} \
            --total_reads_R1_raw "~{total_reads_R1_raw}" \
            --total_reads_R2_raw "~{total_reads_R2_raw}" \
            --read_length_R1_raw "~{read_length_R1_raw}" \
            --read_length_R2_raw "~{read_length_R2_raw}" \
            --read_pairs_raw "~{read_pairs_raw}" \
            --total_reads_R1_cleaned "~{total_reads_R1_cleaned}" \
            --total_reads_R2_cleaned "~{total_reads_R2_cleaned}" \
            --read_length_R1_cleaned "~{read_length_R1_cleaned}" \
            --read_length_R2_cleaned "~{read_length_R2_cleaned}" \
            --read_pairs_cleaned "~{read_pairs_cleaned}" \
            --sample_name "~{sample_name}"
        

    >>>

    output {
        File preprocess_qc_metrics = "~{sample_name}_preprocess_qc_metrics.csv"
    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }  
}

