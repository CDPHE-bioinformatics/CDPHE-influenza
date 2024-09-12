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
        String irma_type
        String irma_subtype
        String base_name
    }

    String docker = "andersenlabapps/ivar:1.3.1"

    #record parameters
    Int ivar_min_depth = 25
    Float ivar_min_freq = 0.5
    Int ivar_min_qual = 20

    command <<<

        # create name for sorted bam file (34345_A_HA-H1.bam or 3455_A_NP.bam)
        # bam_file = some_directory/2400000_A_HA-H1.bam
        # this line of code returns 240000_A_HA-H1 (so strips the upper directories and the file suffix)
        # some_directory/2400000_A_HA-H1.bam ---> 24000000_A_HA-H1.bam ---> 2400000_A_HA-H1
        
        # pull sample id, segment base_name/~{sample_name}/}" | cut -d "_" -f 2-) 
        # shell replacement syntax: varname/pattern/replacement - to replace teh first occurance of teh pattern
        # 2400000_HA --> _HA ---> HA
        
        # generate consensus; first sort bam file
        samtools sort ~{bam_file} -o sorted.bam
        samtools mpileup -A --a -B -Q ~{ivar_min_qual} sorted.bam | \
        ivar consensus -p ${base_name} -q ~{ivar_min_qual} -t ~{ivar_min_freq} -m ~{ivar_min_depth} | tee ${base_name}_ivar_screen_capture.txt
        
        # fasta will be named base_name.fa
        cat ${base_name}.fa # for troubleshooting purposes print fasta contents to screen

        # rename consensus header
        # header_name = $(echo ${base_name})

        # rename consensus header
        # header name = 24000000_A_HA-H1, 24000000_A_NP, etc.
        # if [ "~{irma_type}" == "A" ]; then
        #     if [ "$segment_name" == "HA" ]; then
        #         subtype="~{irma_ha_subtype}"
        #     elif [ "$segment_name" == "NA" ]; then
        #         subtype="~{irma_na_subtype}"
        #     else
        #         subtype=""  # Default if segment is neither "HA" nor "NA"
        #     fi
        # elif [ "$irma_type" == "B" ]; then
        #     subtype=""  # Set subtype to empty string if irma_type is "B"

        # fi
        
        # if [${subtype} == ""]; then
        #     header_name=$(echo ~{sample_name}_~{irma_type}_${segment_name})
        # else
        #     header_name=$(echo ~{sample_name}_~{irma_type}_${segment_name}_${subtype})
        # fi

        # echo ${header_name} # print for troubleshooting purposes
        
        sed -i "s/>.*/>${base_name}/" ${base_name}.fa
        # for sed, -i means edit file in place, s means substitution
        # s/regular expression/replacement/

        cat ${base_name}.fa # for troubleshooting purposes print fasta contents to screen

        # output ivar parameters and version
        ivar version | awk '/version/ {print $3}' | tee VERSION_ivar
        echo "ivar_min_depth,ivar_min_freq,ivar_min_qual" > ivar_parameters.csv
        echo "~{ivar_min_depth},~{ivar_min_freq},~{ivar_min_qual}" >> ivar_parameters.csv

        # samtools version
        samtools --version | awk '/samtools/ {print $2}' | tee VERSION_samtools
    >>>

    output {
        File ivar_consensus_fasta = select_first(glob("*.fa"))
        File? ivar_seg_ha_fasta = "~{sample_name}_HA*fa"
        File? ivar_seg_na_fasta = "~{sample_name}_NA*fa"
        File ivar_output = select_first(glob("*_ivar_screen_capture.txt")) # currently not an output or transfered
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
