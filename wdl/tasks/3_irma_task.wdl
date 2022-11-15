verison 1.0

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
    sed "s/>\(.*\)/>~{sample_name}_\1/g" $fasta
    mv -- "$fasta" "~{sample_name}_${fasta%}"

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
  for file in *-coverage.txt
  do
    mv -- "$file" "~{sample_name}${file%}"
  done

  >>>

  output{
    Array[File] irma_output_directory = glob('~{sample_name}/*')
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
