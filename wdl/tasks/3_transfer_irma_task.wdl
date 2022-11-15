version 1.0

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
    Int read_length_R1_raw
    Int read_length_R2_raw
    String read_pairs_raw

    String seqyclean_version
    String seqyclean_docker

    Int total_reads_R1_cleaned
    Int total_reads_R2_cleaned
    Int read_length_R1_cleaned
    Int read_length_R2_cleaned
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
  dict = {'sample_name' : ~{sample_name},

          'fastqc_version' : ~{fastqc_version},
          'fastqc_docker' : ~{fastqc_docker},

          'seqyclean_version' : ~{seqyclean_version},
          'seqyclean_docker' : ~{seqyclean_docker},

          'irma_version': ~{irma_version},

          'total_reads_R1_raw' : ~{total_reads_R1_raw},
          'total_reads_R2_raw' : ~{total_reads_R2_raw},
          'read_length_R1_raw' : ~{read_length_R1_raw},
          'read_length_R2_raw' : ~{read_length_R2_raw},
          'read_pairs_raw' : ~{read_pairs_raw},

          'total_reads_R1_cleaned': ~{total_reads_R1_cleaned},
          'total_reads_R2_cleaned' : ~{total_reads_R2_cleaned},
          'read_length_R1_cleaned': ~{read_length_R1_cleaned},
          'read_length_R2_cleaned': ~{read_length_R2_cleaned},
          'read_pairs_cleaned': ~{read_paris_cleaned}
          }

  df = pd.DataFrame.from_dict(dict)
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
    File fastqc2_raw_html_
    File fastqc2_raw_zip

    File seqyclean_summary

    File fastqc1_cleaned_html
    File fastqc1_cleaned_zip
    File fastqc2_cleaned_html
    File fastqc2_cleaned_zip

    # irma results
    Array[File] irma_output_directory

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
  gsutil -m cp ~{seqyclean_summary}/seqyclean/

  # transfer fastqc clean
  gsutil -m cp ~{fastqc1_cleaned_html} ~{out_path}/fastqc_cleaned/
  gsutil -m cp ~{fastqc1_cleaned_zip} ~{out_path}/fastqc_cleaned/
  gsutil -m cp ~{fastqc2_cleaned_html} ~{out_path}/fastqc_cleaned/
  gsutil -m cp ~{fastqc2_cleaned_zip} ~{out_path}/fastqc_cleaned/

  # transfer irma
  gsutil -m cp ~{sep = " " irma_output_directory} ~{out_path}/irma/

  # transfer parameters and outputs
  gsutil -m cp ~{parameters_and_outputs_tsv} ~{output}/summary_stats/

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
