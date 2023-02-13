version 1.0

# begin tasks
task fastqc {
    meta {
      description: "this task uses fastqc to evaluate the quality of reads in the fastq files. The output is an html file with various quality statistics, from which the task pulls the number of reads. Modified from Theiagen Genomics- public helath viral genomics."
    }

    input {
        String sample_id
        File fastq_R1
        File fastq_R2 = "default.fastq.gz" # need a dummy value so can use the basename on fastq_R2
        String read_type
        String docker = 'staphb/fastqc:0.11.9'
    }

    String fastq_R1_name = basename(basename(basename(fastq_R1, ".gz"), ".fastq"), ".fq")
    String fastq_R2_name = basename(basename(basename(fastq_R2, ".gz"), ".fastq"), ".fq")

    command <<<
        # capture date and version
        fastqc --version | tee VERSION

        if [ ~{read_type} == "paired" ]; then
            # get the base name of the fastq files
            # fastq_R1_name=$(basename ~{fastq_R1} | cut -d "." -f 1 | cut -d "." -f 1)
            # fastq_R2_name=$(basename ~{fastq_R2} | cut -d "." -f 1 | cut -d "." -f 1)

            # run fastqc
            fastqc --outdir $PWD ~{fastq_R1} ~{fastq_R2}

            # # rename outputs
            # mv ${fastq_R1_name}_fastqc.html ~{sample_id}_R1_fastqc.html
            # mv ${fastq_R1_name}_fastqc.zip ~{sample_id}_R1_fastqc.zip
            # mv ${fastq_R2_name}_fastqc.html ~{sample_id}_R2_fastqc.html
            # mv ${fastq_R2_name}_fastqc.zip ~{sample_id}_R2_fastqc.zip

            # pull some info from the zip file regarding number of reads and read length
            unzip -p ~{fastq_R1_name}_fastqc.zip */fastqc_data.txt | grep "Sequence length" | cut -f 2 | tee READ1_LEN
            unzip -p ~{fastq_R2_name}_fastqc.zip */fastqc_data.txt | grep "Sequence length" | cut -f 2 | tee READ2_LEN

            unzip -p ~{fastq_R1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ1_SEQS
            unzip -p ~{fastq_R2_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ2_SEQS

            READ1_SEQS=$(unzip -p ~{fastq_R1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )
            READ2_SEQS=$(unzip -p ~{fastq_R2_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )

            if [ $READ1_SEQS == $READ2_SEQS ]; then
                read_pairs=$READ1_SEQS
            else
                read_pairs="Uneven pairs: R1=$READ1_SEQS, R2=$READ2_SEQS"
            fi

            echo $read_pairs | tee READ_PAIRS
        
        elif [ ~{read_type} == "single" ]; then
            # get base name
            # fastq_R1_name=$(basename ~{fastq_R1} | cut -d "." -f 1 | cut -d "." -f 1)


            # run fastqc
            fastqc --outdir $PWD ~{fastq_R1}

            # # rename outputs
            # mv ${fastq_R1_name}_fastqc.html ~{sample_id}_R1_fastqc.html
            # mv ${fastq_R1_name}_fastqc.zip ~{sample_id}_R1_fastqc.zip

            # pull some info from the zip file regarding number of reads and read length
            unzip -p ~{fastq_R1_name}_fastqc.zip */fastqc_data.txt | grep "Sequence length" | cut -f 2 | tee READ1_LEN
            # unzip -p ~{fastq_R2_name}_fastqc.zip */fastqc_data.txt | grep "Sequence length" | cut -f 2 | tee READ2_LEN
            
            unzip -p ~{fastq_R1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ1_SEQS
            # unzip -p ~{fastq_R2_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ2_SEQS

            READ1_SEQS=$(unzip -p ~{fastq_R1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )
            # READ2_SEQS=$(unzip -p ~{fastq_R2_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )
    
            echo $READ1_SEQS| tee READ_PAIRS
            
            # create dumby variables
            echo 0 | tee READ2_SEQS
            echo 0 | tee READ2_LEN
        fi

    >>>
    output {
        File fastqc1_html = "~{fastq_R1_name}_fastqc.html"
        File fastqc1_zip = "~{fastq_R1_name}_fastqc.zip"
        File? fastqc2_html = "~{fastq_R2_name}_fastqc.html"
        File? fastqc2_zip = "~{fastq_R2_name}_fastqc.zip"

        Int total_reads_R1 = read_string("READ1_SEQS")
        Int total_reads_R2 = read_string("READ2_SEQS")
        
        String read_length_R1 = read_string('READ1_LEN')
        String read_length_R2 = read_string('READ2_LEN')
        
        String read_pairs = read_string("READ_PAIRS")
        String fastqc_version = read_string("VERSION")
        String fastqc_docker = "~{docker}"
    } 
    runtime {
      docker: "~{docker}"
      memory: "1 GiB"
      cpu: 2
      disks: "local-disk 100 SSD"
      preemptible: 0
      maxRetries: 3
    }
}

task seqyclean {
    meta{
        description : "This task reads in a pair of fastq files (R1 and R2), uses seqyclean to remove containments and adapaters and then filters on min len and quality. The output is two clean fastq files (R1 and R2)."
    }
    input {
        File adapters_and_contaminants
        String sample_id
        File fastq_R1
        File? fastq_R2
        String read_type
        String docker = "staphb/seqyclean:1.10.09"
    }
    command <<<

        # get version number of seqyclean
        if [ ~{read_type} == "paired" ]; then
            seqyclean -minlen 70 -qual 30 30 -gz -1 ~{fastq_R1} -2 ~{fastq_R2} -c ~{adapters_and_contaminants} -o ~{sample_id}_clean
        elif [ ~{read_type} == "single" ]; then
            seqyclean -minlen 70 -qual 30 30 -gz -U ~{fastq_R1} -c ~{adapters_and_contaminants} -o ~{sample_id}_clean
            mv ~{sample_id}_clean_SE.fastq.gz ~{sample_id}_clean_PE1.fastq.gz # change name so matches output
        fi
        
        # pull version out of the summary file
        awk 'NR==2 {print $1}' ~{sample_id}_clean_SummaryStatistics.tsv | tee VERSION
    >>>

    output {
        String seqyclean_version = read_string("VERSION")
        String seqyclean_docker = "~{docker}"
        File fastq_R1_cleaned = "${sample_id}_clean_PE1.fastq.gz"
        File? fastq_R2_cleaned = "${sample_id}_clean_PE2.fastq.gz"
        File seqyclean_summary = "${sample_id}_clean_SummaryStatistics.tsv"
    }
    runtime {
        docker: "~{docker}"
        memory: "2 GiB"
        cpu: 2
        disks: "local-disk 100 SSD"
        preemptible: 0
        maxRetries: 3
    }
}

task concat_preprocess_qc_metrics {
    meta {
        description: "creating a single csv file with all pre=processing qc"
    }

    input {
        File python_script
        String sample_id
        String read_type
        
        String fastqc_version
        String fastqc_docker

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

        String seqyclean_version
        String seqyclean_docker
        
    }

    command <<<

        if [ ~{read_type} == "paired" ]; then 
            python ~{python_script} \
                --fastqc_version "~{fastqc_version}" \
                --fastqc_docker "~{fastqc_docker}" \
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
                --seqyclean_version "~{seqyclean_version}" \
                --seqyclean_docker "~{seqyclean_docker}" \
                --read_type "~{read_type}" \
                --sample_id "~{sample_id}"
        
        elif [ ~{read_type} == "single" ]; then
            python ~{python_script} \
                --fastqc_version "~{fastqc_version}" \
                --fastqc_docker "~{fastqc_docker}" \
                --total_reads_R1_raw "~{total_reads_R1_raw}" \
                --read_length_R1_raw "~{read_length_R1_raw}" \
                --read_pairs_raw "~{read_pairs_raw}" \
                --total_reads_R1_cleaned "~{total_reads_R1_cleaned}" \
                --read_length_R1_cleaned "~{read_length_R1_cleaned}" \
                --read_pairs_cleaned "~{read_pairs_cleaned}" \
                --seqyclean_version "~{seqyclean_version}" \
                --seqyclean_docker "~{seqyclean_docker}" \
                --read_type "~{read_type}" \
                --sample_id "~{sample_id}"
        
        fi

    >>>

    output {
        File preprocess_qc_metrics = "~{sample_id}_preprocess_qc_metrics.csv"
    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }  
}

