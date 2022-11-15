version 1.0

# import tasks
import "../tasks/fastq_preprocess_tasks.wdl" as preprocess
import "../tasks/3_irma_task.wdl" as irma_task
import "../tasks/3_transfer_irma_task" as transfer

# begin workflow
workflow irma_beta_wf{

  input{
    String sample_name
    File fastq_R1
    File fastq_R2
    File adapters_and_contaminants
    String gcp_out_bucket_path
  }

  call preprocess.fastqc as fastqc_raw {
    input:
      fastq_R1 = fastq_R1,
      fastq_R2 = fastq_R2
  }

  call preprocess.seqyclean as seqyclean{
    input:
      adapters_and_contaminants = adapters_and_contaminants,
      sample_name = sample_name,
      fastq_R1 = fastq_R1,
      fastq_R2 = fastq_R2

  }
  call preprocess.fastqc as fastqc_cleaned {
    input:
      fastq_R1 = seqyclean.fastq_R1_cleaned,
      fastq_R2 = seqyclean.fastq_R2_cleaned
  }

  call irma_task.cdc_IRMA as irma {
    input:
      sample_name = sample_name
      fastq_R1 = seqyclean.fastq_R1_cleaned,
      fastq_R2 = seqyclean.fastq_R2_cleaned
  }

  call transfer.get_run_parameters as run_parameters {
    input:
      sample_name = sample_name,
      fastqc_version = fastqc_raw.fastqc_version,
      fastqc_docker = fastqc_raw.fastqc_docker,

      total_reads_R1_raw = fastqc_raw.total_reads_R1,
      total_reads_R2_raw = fastqc_raw.total_reads_R2,
      read_length_R1_raw = fastqc_raw.read_length_R1,
      read_length_R2_raw = fastqc_raw.read_length_R2,
      read_pairs_raw = fastqc_raw.read_pairs_raw,

      seqyclean_version = seqyclean.seqyclean_version,
      seqyclean_docker = seqyclean.seqyclean_docker,

      total_reads_R1_cleaned = fastqc_cleaned .total_reads_R1,
      total_reads_R2_cleaned  = fastqc_cleaned .total_reads_R2,
      read_length_R1_cleaned  = fastqc_cleaned .read_length_R1,
      read_length_R2_cleaned  = fastqc_cleaned .read_length_R2,
      read_pairs_cleaned  = fastqc_cleaned.read_pairs_raw,

      irma_version = irma.version

  }
  call transfer.transfer_IRMA as transfer_irma
    input:
      sample_name = sample_name,
      gcp_out_bucket_path = gcp_out_bucket_path,
      fastqc1_raw_html = fastq_raw.fastqc1_html,
      fastqc1_raw_zip = fastqc_raw.fastqc_zip,
      fastqc2_raw_html = fastqc_raw.fastqc2_html,
      fastqc2_raw_zip = fastqc_raw.fastqc2_zip,

      seqyclean_summary = seqyclean.seqyclean_summary,

      fastqc1_cleaned_html = fastqc_cleaned.fastqc1_html,
      fastqc1_cleaned_zip = fastqc_cleaned.fastqc1_zip,
      fastqc2_cleaned_html = fastqc_cleaned.fastqc2_html,
      fastqc2_cleaned_zip = fastqc_cleaned.fastqc2_zip,

      irma_output_directory = irma.irma_version,
      parameters_and_outpus_tsv = run_parameters.parameters_and_outputs_tsv

  output {
    # output from preprcoess
    String fastqc_version = fastqc_raw.fastqc_version,
    String fastqc_docker = fastqc_raw.fastqc_docker,

    File fastqc1_html_raw = fastqc_raw.fastqc1_html,
    File fastqc1_zip_raw = fastqc_raw.fastqc1_zip,
    File fastqc2_html_raw = fastqc_raw.fastqc2_html,
    File fastqc2_zip_raw = fastqc_raw.fastqc2_zip,

    String seqyclean_version =
    String seqyclean_docker =

    File fastqc1_html_cleaned = fastqc_cleaned.fastqc1_html,
    File fastqc1_zip_cleaned = fastqc_cleaned.fastqc1_zip,
    File fastqc2_html_cleaned = fastqc_cleaned.fastqc2_html,
    File fastqc2_zip_cleaned = fastqc_raw.fastqc2_zip,


    # output from irma
    String irma_version = irma.irma_version,
    Array[File] irma_out_dir = irma.irma_output_directory,

    # output from run_parameters
    File parameters_and_outputs_tsv = run_parameters.parameters_and_outputs_tsv


  }

}
