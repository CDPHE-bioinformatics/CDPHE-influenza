version 1.0

# begin workflow
workflow irma_beta_wf{

  input{
    String sample_name
    File fastq_R1
    File fastq_R2
    File adapters_and_contaminants
    String gcp_out_bucket_path
  }

  call fastqc as fastqc_raw {
    input:
      fastq_R1 = fastq_R1,
      fastq_R2 = fastq_R2
  }

  call seqyclean as seqyclean{
    input:
      adapters_and_contaminants = adapters_and_contaminants,
      sample_name = sample_name,
      fastq_R1 = fastq_R1,
      fastq_R2 = fastq_R2

  }
  call fastqc as fastqc_cleaned {
    input:
      fastq_R1 = seqyclean.fastq_R1_cleaned,
      fastq_R2 = seqyclean.fastq_R2_cleaned
  }

  call cdc_IRMA as irma {
    input:
      sample_name = sample_name,
      fastq_R1 = seqyclean.fastq_R1_cleaned,
      fastq_R2 = seqyclean.fastq_R2_cleaned
  }

  call get_run_parameters as run_parameters {
    input:
      sample_name = sample_name,
      fastqc_version = fastqc_raw.fastqc_version,
      fastqc_docker = fastqc_raw.fastqc_docker,

      total_reads_R1_raw = fastqc_raw.total_reads_R1,
      total_reads_R2_raw = fastqc_raw.total_reads_R2,
      read_length_R1_raw = fastqc_raw.read_length_R1,
      read_length_R2_raw = fastqc_raw.read_length_R2,
      read_pairs_raw = fastqc_raw.read_pairs,

      seqyclean_version = seqyclean.seqyclean_version,
      seqyclean_docker = seqyclean.seqyclean_docker,

      total_reads_R1_cleaned = fastqc_cleaned.total_reads_R1,
      total_reads_R2_cleaned  = fastqc_cleaned.total_reads_R2,
      read_length_R1_cleaned  = fastqc_cleaned.read_length_R1,
      read_length_R2_cleaned  = fastqc_cleaned.read_length_R2,
      read_pairs_cleaned  = fastqc_cleaned.read_pairs,

      irma_version = irma.irma_version

  }
  call transfer_IRMA as transfer_irma {
    input:
      sample_name = sample_name,
      gcp_out_bucket_path = gcp_out_bucket_path,
      fastqc1_raw_html = fastqc_raw.fastqc1_html,
      fastqc1_raw_zip = fastqc_raw.fastqc1_zip,
      fastqc2_raw_html = fastqc_raw.fastqc2_html,
      fastqc2_raw_zip = fastqc_raw.fastqc2_zip,

      seqyclean_summary = seqyclean.seqyclean_summary,

      fastqc1_cleaned_html = fastqc_cleaned.fastqc1_html,
      fastqc1_cleaned_zip = fastqc_cleaned.fastqc1_zip,
      fastqc2_cleaned_html = fastqc_cleaned.fastqc2_html,
      fastqc2_cleaned_zip = fastqc_cleaned.fastqc2_zip,

      consensus_fastas_array = irma.consensus_fastas_array,
      bam_files_array = irma.bam_files_array,
      bai_files_array = irma.bai_files_array,
      vcf_files_array = irma.vcf_files_array,
      tables_directory = irma.tables_directory,
      parameters_and_outputs_tsv = run_parameters.parameters_and_outputs_tsv
}
  output {
    # output from preprcoess
    String fastqc_version = fastqc_raw.fastqc_version
    String fastqc_docker = fastqc_raw.fastqc_docker

    File fastqc1_html_raw = fastqc_raw.fastqc1_html
    File fastqc1_zip_raw = fastqc_raw.fastqc1_zip
    File fastqc2_html_raw = fastqc_raw.fastqc2_html
    File fastqc2_zip_raw = fastqc_raw.fastqc2_zip

    String seqyclean_version = seqyclean.seqyclean_version
    String seqyclean_docker = seqyclean.seqyclean_docker

    File fastqc1_html_cleaned = fastqc_cleaned.fastqc1_html
    File fastqc1_zip_cleaned = fastqc_cleaned.fastqc1_zip
    File fastqc2_html_cleaned = fastqc_cleaned.fastqc2_html
    File fastqc2_zip_cleaned = fastqc_raw.fastqc2_zip


    # output from irma
    String irma_version = irma.irma_version
    Array[File] consensus_fastas_array = irma.consensus_fastas_array
    Array[File] bam_files_array = irma.bam_files_array
    Array[File] vcf_files_array = irma.vcf_files_array
    Array[File] tables_directory = irma.tables_directory

    # output from run_parameters
    File parameters_and_outputs_tsv = run_parameters.parameters_and_outputs_tsv


  }

}


# read in tasks
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
   unzip -p ~{fastq_R1_file_name}_fastqc.zip */fastqc_data.txt | grep "Sequence length" | cut -f 2 | tee READ1_LEN
   unzip -p ~{fastq_R2_file_name}_fastqc.zip */fastqc_data.txt | grep "Sequence length" | cut -f 2 | tee READ2_LEN

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
   String read_length_R1 = read_string('READ1_LEN')
   String read_length_R2 = read_string('READ2_LEN')
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

task cdc_IRMA{
  meta {
    description : 'run irma'
  }

  input {
    String sample_name
    File fastq_R1
    File fastq_R2
    String docker = "staphb/irma:latest"
    Int cpu = 2
    Int memory = 8
  }

  command <<<

  IRMA --version | grep -o v[[:digit:]].[[:digit:]].[[:digit:]] | tee VERSION

  IRMA FLU ~{fastq_R1} ~{fastq_R2} ~{sample_name}

  # rename fasta files and fasta header
  cd ~{sample_name}/

  if ls | grep -q ".fasta"
  then
  for fasta in *.fasta
  do
    sed "s/>\(.*\)/>~{sample_name}_\1/g" $fasta > "~{sample_name}_${fasta%}"
    #mv -- "$fasta" "~{sample_name}_${fasta%}"
  done
  fi

  # rename bam files
  for bam in *.bam*
  do
    mv -- "$bam" "~{sample_name}_${bam%}"
  done

  # rename vcf files
  for vcf in *.vcf
  do
    mv -- "$vcf" "~{sample_name}_${vcf%}"
  done

  # rename coverage files
  cd tables/
  for file in *.txt
  do
    mv -- "$file" "~{sample_name}_${file%}"
  done

  >>>

  output{
    #Array[File] irma_output_directory = glob('~{sample_name}/*')
    Array[File] consensus_fastas_array = glob('~{sample_name}/~{sample_name}*.fasta')
    Array[File] bam_files_array = glob('~{sample_name}/*.bam*')
    Array[File] bai_files_array = glob('~{sample_name}/*.bam.bai')
    Array[File] vcf_files_array = glob('~{sample_name}/*.vcf')
    Array[File] tables_directory = glob('~{sample_name}/tables/*')
    String irma_version = read_string("VERSION")
  }

  runtime {
    docker: "~{docker}"
    memory: "~{memory} GiB"
    cpu: "~{cpu}"
    disks: "local-disk 50 SSD"
    preemptible: 0
  }

}



task get_run_parameters {
  meta {
    description : "create a csv file of all the run parameter info that was output into the terra data table"
  }

  input {
    String sample_name

    # pull string and int outputs from preprocess tasks
    String fastqc_version
    String fastqc_docker
    Int total_reads_R1_raw
    Int total_reads_R2_raw
    String read_length_R1_raw
    String read_length_R2_raw
    String read_pairs_raw

    String seqyclean_version
    String seqyclean_docker

    Int total_reads_R1_cleaned
    Int total_reads_R2_cleaned
    String read_length_R1_cleaned
    String read_length_R2_cleaned
    String read_pairs_cleaned

    # pull irma string outputs
    String irma_version

    String docker = "mchether/py3-bio:v1"
    Int cpu = 2
    Int memory = 8

  }

  command <<<
    python3 <<CODE

    import pandas as pd

    # create dictionary and then dataframe from dictionary
    parameter_dict = {'sample_name' :['~{sample_name}'],
      'fastqc_version' : ['~{fastqc_version}'],
      'fastqc_docker' : ['~{fastqc_docker}'],
      'seqyclean_version' : ['~{seqyclean_version}'],
      'seqyclean_docker' : ['~{seqyclean_docker}'],
      'irma_version': ['~{irma_version}'],
      'total_reads_R1_raw' : ['~{total_reads_R1_raw}'],
      'total_reads_R2_raw' : ['~{total_reads_R2_raw}'],
      'read_length_R1_raw' : ['~{read_length_R1_raw}'],
      'read_length_R2_raw' : ['~{read_length_R2_raw}'],
      'read_pairs_raw' : ['~{read_pairs_raw}'],
      'total_reads_R1_cleaned': ['~{total_reads_R1_cleaned}'],
      'total_reads_R2_cleaned' : ['~{total_reads_R2_cleaned}'],
      'read_length_R1_cleaned': ['~{read_length_R1_cleaned}'],
      'read_length_R2_cleaned': ['~{read_length_R2_cleaned}'],
      'read_pairs_cleaned': ['~{read_pairs_cleaned}']
    }

    df = pd.DataFrame.from_dict(parameter_dict)
    outfile_name = '~{sample_name}_parameters_and_outputs.tsv'
    df.to_csv(outfile_name, sep = '\t', index = False)

    CODE
  >>>

  output {
    File parameters_and_outputs_tsv = '${sample_name}_parameters_and_outputs.tsv'
  }

  runtime {
    docker: "~{docker}"
    memory: "~{memory} GiB"
    cpu: "~{cpu}"
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}

task transfer_IRMA {
  meta {
    description : "transfer IRMA output to gcp bucket"
  }

  input {
    String sample_name
    String gcp_out_bucket_path

    # fastq pre-processing outputs
    File fastqc1_raw_html
    File fastqc1_raw_zip
    File fastqc2_raw_html
    File fastqc2_raw_zip

    File seqyclean_summary

    File fastqc1_cleaned_html
    File fastqc1_cleaned_zip
    File fastqc2_cleaned_html
    File fastqc2_cleaned_zip

    # irma results
    Array[File] consensus_fastas_array
    Array[File] bam_files_array
    Array[File] bai_files_array
    Array[File] vcf_files_array
    Array[File] tables_directory

    # parameters and outputs
    File parameters_and_outputs_tsv

    # run time info
    String docker = "theiagen/utility:1.0"
    Int cpu = 4
    Int memory = 16
  }

  # remove last backslash incase it exists
  String out_path = sub(gcp_out_bucket_path, "/$", "")

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
  gsutil -m cp ~{sep = " " consensus_fastas_array} ~{out_path}/irma/~{sample_name}/consesnus_fastas/
  gsutil -m cp ~{sep = " " bam_files_array} ~{out_path}/irma/~{sample_name}/bam_files/
  gsutil -m cp ~{sep = " " bai_files_array} ~{out_path}/irma/~{sample_name}/bam_files/
  gsutil -m cp ~{sep = " " vcf_files_array} ~{out_path}/irma/~{sample_name}/vcf_files/
  gsutil -m cp ~{sep = " " tables_directory} ~{out_path}/irma/~{sample_name}/tables/

  # transfer parameters and outputs
  gsutil -m cp ~{parameters_and_outputs_tsv} ~{out_path}/summary_stats/

  # transfer date
  transferdate=`date`
  echo $transferdate | tee TRANSFERDATE

  >>>

  output {
    String transfer_date = read_string("TRANSFERDATE")
  }

  runtime {
    docker: "~{docker}"
    memory: "~{memory} GiB"
    cpu: "~{cpu}"
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}
