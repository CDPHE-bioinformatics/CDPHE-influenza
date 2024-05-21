version 1.0

# begin tasks
task ivar_consensus {
    meta {
        description : "generate consesnus sequnce from sorted bam files input"
    }

    input {
        File? bam_file
        String sample_name
        String irma_type
        String irma_na_subtype
        String irma_ha_subtype
        String docker = "andersenlabapps/ivar:1.3.1"
    }

    #record parameters
    Int ivar_min_depth = 25
    Float ivar_min_freq = 0.5
    Int ivar_min_qual = 20
    

    command <<<

    # create name for sorted bam file (34345_HA.bam)
    prefix=$(basename ~{bam_file} | cut -d "." -f 1)
    # bam_file = some_directory/2400000_HA.bam
    # this line of code returns 240000_HA (so strips the upper directories and the file suffix)
    # some_directory/2400000_HA.bam ---> 24000000_HA.bam ---> 2400000_HA
    
    # pull sample id, segment name from original ba file
    segment_name=$(echo "${prefix/~{sample_name}/}" | cut -d "_" -f 2-) 
    # shell replacement syntax: varname/pattern/replacement - to replace teh first occurance of teh pattern
    # 2400000_HA --> _HA ---> HA
    
    # generate consensus; first sort bam file
    samtools sort ~{bam_file} -o sorted.bam
    samtools mpileup -A --a -B -Q ~{ivar_min_qual} sorted.bam | \
    ivar consensus -p ${prefix} -q ~{ivar_min_qual} -t ~{ivar_min_freq} -m ~{ivar_min_depth} | tee ${prefix}_ivar_screen_capture.txt
    
    # fasta will be named prefix.fa
    cat ${prefix}.fa # for troubleshooting purposes print fasta contents to screen

    # rename consensus header
    # header name = 24000000_A_HA_H1, 24000000_A_NP, etc.
    if [ "~{irma_type}" == "A" ]; then
        if [ "$segment_name" == "HA" ]; then
            subtype="~{irma_ha_subtype}"
        elif [ "$segment_name" == "NA" ]; then
            subtype="~{irma_na_subtype}"
        else
            subtype=""  # Default if segment is neither "HA" nor "NA"
        fi
    elif [ "$irma_type" == "B" ]; then
        subtype=""  # Set subtype to empty string if irma_type is "B"

    fi
    
    if [${subtype} == ""]; then
        header_name=$(echo ~{sample_name}_~{irma_type}_${segment_name})
    else
        header_name=$(echo ~{sample_name}_~{irma_type}_${segment_name}_${subtype})
    fi

    echo ${header_name} # print for troubleshooting purposes
    
    sed -i "s/>.*/>${header_name}/" ${prefix}.fa
    # for sed, -i means edit file in place, s means substitution
    # s/regular expression/replacement/

    cat ${prefix}.fa # for troubleshooting purposes print fasta contents to screen

    # output ivar parameters
    ivar version | head -n1 | cut -d " " -f 3 | tee VERSION
    ivar_version=$(cat ivar_version.txt)
    echo "ivar_min_depth,ivar_min_freq,ivar_min_qual" > ivar_parameters.csv
    echo "~{ivar_min_depth},~{ivar_min_freq},~{ivar_min_qual}" >> ivar_parameters.csv

    >>>

    output {
        File ivar_consensus_fasta = select_first(glob("*.fa"))
        File? ivar_seg_ha_fasta = "~{sample_name}_HA*fa"
        File? ivar_seg_na_fasta = "~{sample_name}_NA*fa"
        File ivar_output = select_first(glob("*_ivar_screen_capture.txt")) # currently not an output or transfered
        File ivar_parameters = "ivar_parameters.csv"

        VersionInfo ivar_version_info = object{
            sofware: 'ivar',
            docker: "~{docker}",
            version: read_string("VERSION")
        }

    }

        runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 10 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker: "andersenlabapps/ivar:1.3.1"
    }

}