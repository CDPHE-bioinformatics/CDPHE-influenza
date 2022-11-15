version 1.0

task fastq_scan {
  meta {
    description: "this task reads in a pair of fastq files (R1 and R2),uses fastq-scan to count hte reads in each file, and outputs the total number of reads"
  }
  input {
    String sample_name
    File fastq_R1
    File fastq_R2
    String docker = 'staphb/fastq-scan:0.4.4'
    Int cpu = 2
    Int memory = 2

  }
  command <<<

  # get fastq scan VERSION
  fastq-scan -v | tee VERSION
  zcat ~{fastq_R1} | fastq-scan | jq .qc_stats.read_total > total_reads_R1
  zcat ~{fastq_R1} | fastq-scan | jq .qc_stats.read_total > total_reads_R2
  >>>
  output{
    String version = read_string("VERSION")
    String fastq_scan_docker = "~{docker}"
    String total_reads_R1 = read_string('total_reads_R1')
    String total_reads_R2 = read_string('total_reads_R2')
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: "~{cpu}"
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}

task fastqc {
  meta {
    description : "this task uses fastqc to evaluate the quality of readsin the fastq files. the output is an html file with various quality statistics. Modified from Theiagen Genomics-public health viral genomics"
  }
  input {
    File fastq_R1
    File fastq_R2
    String fastq_R1_file_name = basename(basename(basename(fastq_R1, ".gz"), ".fastq"), ".fq")
    String fastq_R2_file_name = basename(basename(basename(fastq_R2, ".gz"), ".fastq"), ".fq")
    String docker = 'staphb/fastqc:0.11.9'
    Int cpu = 1
    Int memory = 2
  }
  command <<<
    # capture date and version
    fastqc --version | tee VERSION

    fastqc --outdir $PWD ~{fastq_R1} ~{fastq_R2}

    # pull some info from the zip file regarding number of reads and read length
    unzip -p ~{fastq_R1_file_name}_fastqc.zip */fastqc_data.txt | grep "Seqeunce Length" | cut -f 2 | tee READ1_SEQS
    unzip -p ~{fastq_R2_file_name}_fastqc.zip */fastqc_data.txt | grep "Sequence Length" | cut -f 2 | tee READ2_SEQS

    unzip -p ~{fastq_R1_file_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ1_SEQS
    unzip -p ~{fastq_R2_file_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ2_SEQS

    READ1_SEQS=$(unzip -p ~{fastq_R1_file_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )
    READ2_SEQS=$(unzip -p ~{fastq_R2_file_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )

    if [ $READ1_SEQS == $READ2_SEQS ]; then
      read_pairs=$READ1_SEQS
    else
      read_pairs="Uneven pairs: R1=$READ1_SEQS, R2=$READ2_SEQS"
    fi
    echo $read_pairs | tee READ_PAIRS
  >>>
  output {
    File fastqc1_html = "~{fastq_R1_file_name}_fastqc.html"
    File fastqc1_zip = "~{fastq_R1_file_name}_fastqc.zip"
    File fastqc2_html = "~{fastq_R2_file_name}_fastqc.html"
    File fastqc2_zip = "~{fastq_R2_file_name}_fastqc.zip"
    Int total_reads_R1 = read_string("READ1_SEQS")
    Int total_reads_R2 = read_string("READ2_SEQS")
    Int read_length_R1 = read_string('READ1_LEN')
    Int read_length_R2 = read_string('READ2_LEN')
    String read_pairs = read_string("READ_PAIRS")
    String fastqc_version = read_string("VERSION")
    String fastqc_docker = "~{docker}"
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GiB"
    cpu: "~{cpu}"
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}

task seqyclean {
  meta {
    description: " this tasks reads in a pair of fastq files (R1 and R2), uses seqyclean to remove containments, adpators and filters on min length and quality. The output is two clean fastq files (R1 and R2)"
  }
  input {
    File adapters_and_contaminants
    String sample_name
    File fastq_R1
    File fastq_R2
    String docker = 'staphb/seqyclean:1.10.09'
    Int cpu = 2
    Int memory = 2
  }
  command <<<
  # get version number of seqyclean

  seqyclean -minlen 70 -qual 30 30 -gz -1 ~{fastq_R1} -2 ~{fastq_R1} -c ~{adapters_and_contaminants} -o ~{sample_name}_clean

  # pull version out of the summary file
  awk 'NR==2 {print $1}' ~{sample_name}_clean_SummaryStatistics.tsv | tee VERSION
  >>>

  output {
    String seqyclean_version = read_string("VERSION")
    String seqyclean_docker = "~{docker}"
    File fastq_R1_cleaned = "${sample_name}_clean_PE1.fastq.gz"
    File fastq_R2_cleaned = "${sample_name}_clean_PE2.fastq.gz"
    File seqyclean_summary = "${sample_name}_clean_SummaryStatistics.tsv"
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GiB"
    cpu: "~{cpu}"
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }

}
