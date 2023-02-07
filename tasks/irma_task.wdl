version 1.0

task irma {
    meta {
        description: "runs CDC's IRMA for FLU. Form more info on IRMA visit: https://wonder.cdc.gov/amd/flu/irma/configuration.html. Modified from theiagen genomics - publich helath viral genomics. This task runs the default config of IRMA meaning it uses: \ALIGN_PROG=SAM \DEL_TYPE=''"
    }

    input {
        String sample_id
        File fastq_R1
        File fastq_R2
        String irma_module = "FLU"
        String docker = "staphb/irma:latest"

    }
    command <<<
        # potentially adjust config file to match theiagen
        IRMA | head -n1 | awk -F' ' '{ print "IRMA " $5 }' | tee VERSION

        # run IRMA
        IRMA ~{irma_module} ~{fastq_R1} ~{fastq_R2} ~{sample_id}

        # determine if assemly was successful
        if compgen -G "~{sample_id}/*.fasta"; then
            proceed=yes
        else 
            TYPE="no IRMA assembly generated"
            HA_SUBTYPE=""
            NA_SUBTYPE=""
        fi

        # if succesful then determine type; if type == A then grab subtype
        if [ $proceed == 'yes' ]; then
            TYPE=$(basename $(find ~{sample_id}/*.fasta | head -n 1 ) | cut -d "_" -f 1)

            ## Type == A then grab subtype if exists
            if [ $TYPE == "A"]; then
                if compgen -G "~{sample_id}/*_HA*.fasta"; then
                    HA_SUBTYPE=$(basename $(find ~{sample_id}/*HA*.fasta | head -n 1 ) | cut -d "_" -f 3 | cut -d "." -f1)
                else
                    HA_SUBTYPE=""
                fi

                if compgen -G "~{sample_id}/*_NA*.fasta"; then
                    NA_SUBTYPE=$(basename $(find ~{sample_id}/*NA*.fasta | head -n 1 ) | cut -d "_" -f 3 | cut -d "." -f1)
                else
                    NA_SUBTYPE=""
                fi
            else
                HA_SUBTYPE=""
                NA_SUBTYPE=""
            fi

            # rename header and file name for fasta
            ## also create an array of the segment names
            for file in ~{sample_id}/*.fasta; do
                # grab base name and drop .fasta
                segment=$(basename ${file%.*})
                # echo $segement >> segment_list.txt
                header_name=$(echo ~{sample_id}_${segment})
                sed -i "s/>.*/>${header_name}/" ${file}

                # rename file
                new_name=$(echo ~{sample_id}_$(basename $file))
                mv "${file}" "${new_name}"
            done

            # rename bam and vcf files
            for file in ~{sample_id}/*{.vcf,.bam,.bai}; do
                new_name=$(echo ~{sample_id}_$(basename $file))
                mv "${file}" "${new_name}"
            done

        fi

        # print variables to files so can read output
        echo ${TYPE} > TYPE.txt
        echo ${HA_SUBTYPE} > HA_SUBTYPE.txt
        echo ${NA_SUBTYPE} > NA_SUBTYPE.txt
    >>>

    output {
        String irma_type = read_string("TYPE.txt")
        String irma_ha_subtype = read_string("HA_SUBTYPE.txt")
        String irma_na_subtype = read_string("NA_SUBTYPE.txt")
        # Array[String] segment_array = read_lines("segment_list.txt")
        Array[File] irma_assemblies = glob("~{sample_id}*.fasta")
        Array[File] irma_bam_files = glob("~{sample_id}*.bam")
        Array[File] irma_vcfs = glob("~{sample_id}*.vcf")
        String irma_version = read_string("VERSION")
        String irma_docker = "~{docker}"
    }

    runtime {
        docker: "staphb/irma:latest"
        memory: "8 GiB"
        cpu: 2
        disks: "local-disk 50 SSD"
        preemptible: 0
  }
}