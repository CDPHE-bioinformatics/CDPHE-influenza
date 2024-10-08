version 1.0

# define structure
struct VersionInfo {
  String software
  String docker
  String version
}

# begin tasks
task call_consensus_ivar {
    meta {
        description : "generate consensus sequence from sorted bam files input"
    }

    input {
        File? bam_file
        String sample_name
        String base_name
    }

    String docker = "andersenlabapps/ivar:1.3.1"
    # String prefix = sub(base_name, "_irma", "")

    #record parameters
    Int ivar_min_depth = 25
    Float ivar_min_freq = 0.5
    Int ivar_min_qual = 20

    command <<<

        
        # generate consensus; first sort bam file
        # sorted bam and bai files are made and transfered from
        # the calc samtools stats task
        echo "base_name"
        echo ~{base_name}
        
        samtools sort ~{bam_file} -o sorted.bam
        samtools mpileup -A --a -B -Q ~{ivar_min_qual} sorted.bam | \
        ivar consensus -p ~{base_name} -q ~{ivar_min_qual} -t ~{ivar_min_freq} -m ~{ivar_min_depth} | tee ~{base_name}_ivar_screen_capture.txt
        
        # fasta will be named base_name.fa out of ivar; rename with .fasta ending
        mv ~{base_name}.fa ~{base_name}.fasta
        cat ~{base_name}.fasta # for troubleshooting purposes print fasta contents to screen
        
        # rename fasta header
        sed -i "s/>.*/>~{base_name}/" ~{base_name}.fasta
        # for sed, -i means edit file in place, s means substitution
        # s/regular expression/replacement/

        # for troubleshooting-print fasta contents to screen
        cat ~{base_name}.fasta 

        # output ivar parameters and version
        ivar version | awk '/version/ {print $3}' | tee VERSION_ivar
        echo "ivar_min_depth,ivar_min_freq,ivar_min_qual" > ivar_parameters.csv
        echo "~{ivar_min_depth},~{ivar_min_freq},~{ivar_min_qual}" >> ivar_parameters.csv

        # samtools version
        samtools --version | awk '/samtools/ {print $2}' | tee VERSION_samtools
    >>>

    output {
        File ivar_consensus_fasta = "~{base_name}.fasta"
        File ivar_output = "~{base_name}_ivar_screen_capture.txt"# currently not an output or transfered
        File ivar_parameters = "ivar_parameters.csv"

        VersionInfo ivar_version_info = object{
            software: "ivar",
            docker: docker,
            version: read_string("VERSION_ivar")
        }

        VersionInfo samtools_version_info = object{
            software: "samtools",
            docker: docker,
            version: read_string("VERSION_samtools")
        }

    }

        runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 10 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker: docker
    }

}
