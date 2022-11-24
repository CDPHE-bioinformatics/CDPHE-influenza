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

  call get_assembly_metrics as assembly_metrics {
    input:
      sample_name = sample_name,
      consensus_fastas_array = irma.consensus_fastas_array,
      tables_directory = irma.tables_directory,

      total_reads_R1_raw = fastqc_raw.total_reads_R1,
      total_reads_R2_raw = fastqc_raw.total_reads_R2,
      read_length_R1_raw = fastqc_raw.read_length_R1,
      read_length_R2_raw = fastqc_raw.read_length_R2,
      read_pairs_raw = fastqc_raw.read_pairs,

      total_reads_R1_cleaned = fastqc_cleaned.total_reads_R1,
      total_reads_R2_cleaned  = fastqc_cleaned.total_reads_R2,
      read_length_R1_cleaned  = fastqc_cleaned.read_length_R1,
      read_length_R2_cleaned  = fastqc_cleaned.read_length_R2,
      read_pairs_cleaned  = fastqc_cleaned.read_pairs

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

      parameters_and_outputs_tsv = run_parameters.parameters_and_outputs_tsv,
      assembly_metrics_tsv = assembly_metrics.assembly_metrics_tsv
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

    # output from assembly_metrics
    File assembly_metrics_tsv = assembly_metrics.assembly_metrics_tsv
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
   unzip -p ~{fastq_R1_file_name}_fastqc.zip */fastqc_data.txt | grep "Seqeunce length" | cut -f 2 | tee READ1_LEN
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
  for fasta in *.fasta
  do
    sed "s/>\(.*\)/>~{sample_name}_\1/g" $fasta > "~{sample_name}_${fasta%}"
    #mv -- "$fasta" "~{sample_name}_${fasta%}"
  done

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

    # asembly metrics
    File assembly_metrics_tsv
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

  # transfer assembly metrics outputs
  gsutil -m cp ~{assembly_metrics_tsv} ~{out_path}/summary_stats/
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


task get_assembly_metrics {

  meta {
    descripiton : "pull together assembly metrics and perform 'subtyping'. The output will be strings which I will read into the task parameters"
  }

  input {
    String sample_name
    Array[File] consensus_fastas_array
    Int num_fastas = length(consensus_fastas_array)
    Array[File] tables_directory

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

    String docker = "mchether/py3-bio:v1"
    Int cpu = 2
    Int memory = 1

  }

  command <<<
    python3 <<CODE

    import pandas as pd
    import re
    from Bio import SeqIO

    ######################################################
    #### load some reference data about gene segements ###
    # gene segments
    flu_gene_segs = ['HA',  'NA', 'MP','NP', 'NS', 'PA', 'PB1', 'PB2']

    # expected length of each gene segment
    ref_len_dict = {'A_MP': 982,'A_NP': 1497,'A_NS': 863,'A_PA': 2151,'A_PB1': 2274,'A_PB2': 2280,
   'A_HA_H1': 1704,'A_HA_H10': 1686,'A_HA_H11': 1698,'A_HA_H12': 1695,'A_HA_H13': 1701,'A_HA_H14': 1707,
   'A_HA_H15': 1713,'A_HA_H16': 1698,'A_HA_H2': 1689,'A_HA_H3': 1704,'A_HA_H4': 1695,'A_HA_H5': 1707,
   'A_HA_H6': 1704,'A_HA_H7': 1713,'A_HA_H8': 1701,'A_HA_H9': 1683,'A_NA_N1': 1413,'A_NA_N2': 1410,
   'A_NA_N3': 1410,'A_NA_N4': 1413,'A_NA_N5': 1422,'A_NA_N6': 1413,'A_NA_N7': 1416,
   'A_NA_N8': 1413,'A_NA_N9': 1413,'B_HA': 1758,'B_MP': 1139,'B_NA': 1408,'B_NP': 1683,
   'B_NS': 1034,'B_PA': 2181,'B_PB1': 2263,'B_PB2': 2313}

   #######################
   #### prep dataframe ###
   # set empty dataframe (will add some other column for fastqc stuff later)
   col_headers = ['sample_id','flu_type','subtype',
   'HA_expected_len','HA_seq_len','HA_mean_depth','HA_percent_coverage',
   'NA_expected_len','NA_seq_len','NA_mean_depth','NA_percent_coverage',
   'MP_expected_len','MP_seq_len','MP_mean_depth','MP_percent_coverage',
   'NP_expected_len','NP_seq_len','NP_mean_depth','NP_percent_coverage',
   'NS_expected_len','NS_seq_len','NS_mean_depth','NS_percent_coverage',
   'PA_expected_len','PA_seq_len','PA_mean_depth','PA_percent_coverage',
   'PB1_expected_len','PB1_seq_len','PB1_mean_depth','PB1_percent_coverage',
   'PB2_expected_len','PB2_seq_len','PB2_mean_depth','PB2_percent_coverage']

   df = pd.DataFrame(columns = col_headers)
   df['sample_name'] = ["~{sample_name}"]

   ######################################################################
   #### check to see if fasta array is empty, if empty then get zeros ###
   if num_fastas == 0:
      df.at[0, 'flu_type'] = ''
      df.at[0, 'subtype'] = ''
      for column in df.columns:
        if column not in ['sample_name', 'flu_type', 'subtype']:
          df.at[0, column] = 0

   #####################################################################
   ### now determine subytpe, coverage and depth to fill in
   else:
      ################################################################
      ### get list of fasta files so can see which genes got sequenced
      sequenced_gene_segs = []
      for fasta in "~{sep = "" consensus_fastas_array}":
         # fasta file format = "{sample_name}_A_HA_H1.fasta" or "{sample_name}_B_PB1.fasta" for example
         #gene_seg format = "A_HA_H1" or "B_PB1" for example
         gene_seg = file.split('/')[-1].split('.fasta')[0].split('~{sample_name}_')[-1]
         sequenced_gene_segs.append(gene_seg)

       ######################################
       ### get type and subtype
      flu_type = ''
      subtype = ''
      flu_type = sequenced_gene_segs[0].split("_")[0]
      # if flu A determine subtype
      if flu_type == "A":
         HA_sub = ''
         NA_sub = ''
         for gene_seg in sequenced_gene_segs:
            if re.search('HA', gene_seg):
               HA_sub = gene_seg.split('A_HA_')[1]
            elif res.earch('NA', gene_seg):
               NA_sub = gene_seg.split('A_NA_')[1]
         subtype = '%s%s' % (HA_sub, NA_sub)
         df.at[0, 'flu_type'] = flu_type
         df.at[0, 'subtype'] = subtype

      #####################################
      ### get coverage of each gene seg
      for gene_seg in flu_gene_segs:
         # gene_seg_name format = 'HA' or "PB1" for example
         gene_seg_name= gene_seg.split('_')[1] # pulls out gene seg from A_HA_H1 format
         col_name = '%s_expected_len' % gene_seg_name
         expected_len = ref_len_dict[gene_seg]
         df.at[0, col_name] = expected_len

      # calc percent coverage
      for fasta in "~{sep = "" consensus_fastas_array}":
         # get gene seg name
         gene_seg = file.split('/')[-1].split('.fasta')[0].split('~{sample_name}_')[-1]
         gene_seg_name = gene_seg.split('_')[1]

         # get expected gene segment length
         col_name = '%s_expected_len' % gene_seg_name
         expected_len = ref_len_dict[gene_seg]
         df.at[0, col_name] = expected_len

         # read fasta file and pull out seq lenght
         record = SeqIO.read(fasta, 'fasta')
         seq_len = len(record.seq)
         percent_coverage = (seq_len/expected_len) * 100

         col_name = '%s_seq_len' % gene_seg
         df.at[0, col_name] = [seq_len]

         col_name= '%s_percent_coverage' % col_gene_seg
         df.at[0, col_name] = [percent_coverage]


       ######################################
       ### get average depth at each gene segment
       for file in "~{sep = "" tables_directory}":
          if re.search('coverage') in file:
            gene_seg_name = coverage_file.split('/')[-1].split('_')[2]

            cov_df = pd.read_csv(file, sep = '\t')
            cov_df = cov_df.rename(columns = {'Coverage Depth' : 'coverage_depth'})
            mean_depth = cov_df.coverage_depth.mean()

            # write mean deapth for each gene seg to df
            col_name = '%s_mean_depth' % gene_seg_name
            df.at[0, col_name] = mean_depth

    ##########################################################
    ### add in other columns for read qc stuff
    df.at[0, 'total_reads_R1_raw'] = "~{total_reads_R1_raw}"
    df.at[0, 'total_reads_R2_raw'] = '~{total_reads_R2_raw}'
    df.at[0, 'read_length_R1_raw'] = '~{read_length_R1_raw}'
    df.at[0, 'read_length_R2_raw'] = '~{read_length_R2_raw}'
    df.at[0, 'read_pairs_raw'] = '~{read_pairs_raw}'

    df.at[0, 'total_reads_R1_cleaned'] = "~{total_reads_R1_cleaned}"
    df.at[0, 'total_reads_R2_cleaned'] = '~{total_reads_R2_cleaned}'
    df.at[0, 'read_length_R1_cleaned'] = '~{read_length_R1_cleaned}'
    df.at[0, 'read_length_R2_cleaned'] = '~{read_length_R2_cleaned}'
    df.at[0, 'read_pairs_cleaned'] = '~{read_pairs_cleaned}'

    outfile = '~{sample_name}_assembly_metrics.tsv'
    df.to_csv(outfile, sep = '\t', index = False)
    print('written outfile')
    CODE
  >>>

  output {
    File assembly_metrics_tsv = "~{sample_name}_assembly_metrics.tsv"

  }

    runtime {
      docker: "~{docker}"
      memory: "~{memory} GiB"
      cpu: "~{cpu}"
      disks: "local-disk 50 SSD"
      preemptible: 0
    }
}
