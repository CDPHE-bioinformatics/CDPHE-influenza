version 1.0

# begin tasks
task ivar_consensus {
    meta {
        description : "generate consesnus sequnce from bam files input"
    }

    input {
        File bam_file
        String docker = "andersenlabapps/ivar:1.3.1"
    }

    Int ivar_min_depth = 10
    Float ivar_min_freq = 0.5
    Int ivar_min_qual = 20
    

    command <<<

    # pull sample id and segment info from bam file
    sample_id=$(basename ~{bam_file} | cut -d "_" -f 1)
    segment_name=$(basename ~{bam_file} | cut -d "." -f 1 | cut -d "_" -f 2-)
    ivar_prefix=$(basename ~{bam_file} | cut -d "." -f 1)
    
    # generate consensus; first sort bam file
    samtools sort ~{bam_file} -o sorted.bam
    samtools mpileup -A --a -B -Q 20 sorted.bam | \
    ivar consensus -p ${ivar_prefix} -q 20 -t 0.6 -m 10
    # fasta will be named prefix.fa
    # create txt file with the name of the fasta file
    echo ${ivar_prefix}.fa | tee fasta_file_name.txt
    cat ${ivar_prefix}.fa

    # rename consesnus header
    header_name=$(echo ${sample_id}_${segment_name})
    sed -i "s/>.*/>$header_name/" ${ivar_prefix}.fa


    # output ivar parameters
    ivar version | head -n1 | cut -d " " -f 3 | tee ivar_version.txt
    ivar_version=$(cat ivar_version.txt)
    echo "ivar_version,ivar_docker,ivar_min_depth,ivar_min_freq,ivar_min_qual" > ivar_parameters.csv
    echo "${ivar_version},~{docker},~{ivar_min_depth},~{ivar_min_freq},~{ivar_min_qual}" >> ivar_parameters.csv


    >>>

    output {
        File ivar_consensus_fasta = read_string("fasta_file_name.txt")
        # Int ivar_min_depth = ~{ivar_min_depth}
        # Float ivar_min_freq = ~{ivar_min_freq}
        # Int ivar_min_qual = ~{ivar_min_qual}
        String ivar_docker = "~{docker}"
        String ivar_version = read_string("ivar_version.txt")
        File ivar_parameters = "ivar_parameters.csv"

    }

        runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 10 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker: "~{docker}"
    }

}