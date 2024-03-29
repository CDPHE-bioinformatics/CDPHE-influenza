version 1.0

# begin tasks
task fastqc {
    meta {
      description: "this task uses fastqc to evaluate the quality of reads in the fastq files. The output is an html file with various quality statistics, from which the task pulls the number of reads. Modified from Theiagen Genomics- public helath viral genomics."
    }

    input {
        String sample_name
        File fastq_R1
        File? fastq_R2
        String read_type
        String docker = 'staphb/fastqc:0.11.9'
    }

    command <<<
        # grab version
        fastqc --version | tee VERSION

        if [ ~{read_type} == "paired" ]; then

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
        
        elif [ ~{read_type} == "single" ]; then
           # get basename of fastq file
            fastq_R1_name=$(basename ~{fastq_R1} | cut -d "." -f 1 | cut -d "." -f 1)

            # run fastqc
            fastqc --outdir $PWD ~{fastq_R1}

            # pull some info from the zip file regarding number of reads and read length
            unzip -p ${fastq_R1_name}_fastqc.zip */fastqc_data.txt | grep "Sequence length" | cut -f 2 | tee READ1_LEN

            unzip -p ${fastq_R1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ1_SEQS

            READ1_SEQS=$(unzip -p ${fastq_R1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )
    
            echo $READ1_SEQS| tee READ_PAIRS

            # # calculate the total number of reads(i.e. for SE data it's the same as READ1_SEQS)
            # echo $READ1_SEQS | tee TOTAL_READS
            
            # create dummy variables for second read so WDL is happy
            echo 0 | tee READ2_SEQS
            echo 0 | tee READ2_LEN

            # rename files 
             if [ ${fastq_R1_name} != ~{sample_name}_R1 ]; then
                mv ${fastq_R1_name}_fastqc.html ~{sample_name}_R1_fastqc.html
                mv ${fastq_R1_name}_fastqc.zip ~{sample_name}_R1_fastqc.zip
            fi
        fi



    >>>
    output {
        File fastqc1_html = "~{sample_name}_R1_fastqc.html"
        File fastqc1_zip = "~{sample_name}_R1_fastqc.zip"
        File? fastqc2_html = "~{sample_name}_R2_fastqc.html"
        File? fastqc2_zip = "~{sample_name}_R2_fastqc.zip"
 
        # these outputs will go into the concatenated preprocess qc metrics file
        Int total_reads_R1 = read_string("READ1_SEQS")
        Int total_reads_R2 = read_string("READ2_SEQS")
        
        String read_length_R1 = read_string('READ1_LEN')
        String read_length_R2 = read_string('READ2_LEN')
        
        String read_pairs = read_string("READ_PAIRS")
        # String total_reads = read_string("TOTAL_READS")
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
        description : "This task uses seqyclean to remove containments and adapaters and then filters on min len and quality."
    }
    input {
        File adapters_and_contaminants
        String sample_name
        File fastq_R1
        File? fastq_R2
        String read_type
        String docker = "staphb/seqyclean:1.10.09"
    }

    command <<<

        # run seqyclean
        if [ ~{read_type} == "paired" ]; then
            seqyclean -minlen 70 -qual 30 30 -gz -1 ~{fastq_R1} -2 ~{fastq_R2} -c ~{adapters_and_contaminants} -o ~{sample_name}_clean
        elif [ ~{read_type} == "single" ]; then
            seqyclean -minlen 70 -qual 30 30 -gz -U ~{fastq_R1} -c ~{adapters_and_contaminants} -o ~{sample_name}_clean
            mv ~{sample_name}_clean_SE.fastq.gz ~{sample_name}_clean_PE1.fastq.gz # change name so matches output
        fi
        
        # pull version out of the summary file
        awk 'NR==2 {print $1}' ~{sample_name}_clean_SummaryStatistics.tsv | tee VERSION
    >>>

    output {
        String seqyclean_version = read_string("VERSION")
        String seqyclean_docker = "~{docker}"
        File fastq_R1_cleaned = "${sample_name}_clean_PE1.fastq.gz"
        File? fastq_R2_cleaned = "${sample_name}_clean_PE2.fastq.gz"
        File seqyclean_summary = "${sample_name}_clean_SummaryStatistics.tsv"
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
        description: "creates a single csv file with all preprocess qc metrics"
    }

    input {
        File python_script
        String sample_name
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
                --sample_name "~{sample_name}"
        
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
                --sample_name "~{sample_name}"
        
        fi

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

