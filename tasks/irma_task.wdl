version 1.0

# define structure
import "../tasks/capture_version_tasks.wdl" as capture_version
# struct VersionInfo {
#   String software
#   String docker
#   String version
# }

task perform_assembly_irma {
    meta {
        description: "runs CDC's IRMA for FLU. For more info on IRMA visit: https://wonder.cdc.gov/amd/flu/irma/configuration.html. Modified from theiagen genomics - public helath viral genomics (theiaCov PE workflow). This task runs the default config of IRMA meaning it uses: \ALIGN_PROG=SAM \DEL_TYPE=''"
    }

    input {
        String sample_name
        File fastq_R1
        File? fastq_R2
    }

    String docker = "cdcgov/irma:v1.2.1"

    command <<<
        # grab version
        IRMA | head -n1 | awk -F' ' '{ print "IRMA " $5 }' | tee VERSION
        version=$(cat VERSION)

        # set config file
        touch irma_config.sh 
        echo 'MIN_CONS_SUPPORT="50"' >> irma_config.sh
        # 50 appears to be what MIRA uses
        # any base with less than 50x depth will be called an N
        # the fasta files in the amended_consensus directory will have the MIN_CONS_SUPPORT added
        # The fasta files in the amended_consensus directory will also have IUPAC for mixed based calls
        # I will change IUPAC letters to Ns
        
        echo 'DEL_TYPE="DEL"' >> irma_config.sh
       
        # echo 'MIN_CONS_QUALITY="20"' >> irma_config.sh
        # MIRA appears to leave this as the default value at 0, so comment this out for v1.0.0
       
        echo 'MIN_LEN="70"' >> irma_config.sh


        # run IRMA
        IRMA FLU ~{fastq_R1} ~{fastq_R2} ~{sample_name} --external-config irma_config.sh

        # declare associative arrays for segment numbers
        # declare formatted name assoicate array which will be [seg_num] = [A_HA-H1] or [seg_num] = [B_MP]
        # and will be filled in during the loop
        # formatted_name_dict: [segment number] = header name
        declare -A FluA=(["PB2"]="1" ["PB1"]="2" ["PA"]="3" ["HA"]="4" ["NP"]="5" ["NA"]="6" ["MP"]="7" ["NS"]="8" )
        declare -A FluB=(["PB1"]="1" ["PB2"]="2" ["PA"]="3" ["HA"]="4" ["NP"]="5" ["NA"]="6" ["MP"]="7" ["NS"]="8" )      
        declare -A formatted_name_dict

        # IUPAC bases to replace in amended fasta files
        IUPAC=( "B" "D" "H" "K" "M" "N" "R" "S" "V" "W" "Y" )

        echo ""
        echo "IRMA DONE"
        # determine if assembly was successful
        if compgen -G "~{sample_name}/*.fasta"; then
            echo "irma assembly pass" | tee irma_qc.txt
            echo "sample_name,flu_type,gene_segment,subtype" > ~{sample_name}_irma_assembled_gene_segments.csv
            
            # first want ot make my assembled dataframe
            # this is old, with exception that I'm now tracking subtype and TYPE
            echo -e '\n\n\n'
            echo "LOOPING THROUGH FASTA FILES"
            for file in ~{sample_name}/*.fasta; do
                echo -e '\n'
                echo $file
                full=$(basename ${file%.*} | cut -d "." -f 1) # A_HA_H1 or A_NP
                TYPE=$(echo ${full} | cut -d "_" -f 1) # A
                segment=$(echo ${full} | cut -d "_" -f 2) # HA or NP
                segment_subtype=$(echo ${full} | cut -d "_" -f 2-) # HA_H1 or NP
                segment_subtype=${segment_subtype//_/-} # HA-H1 or NP
                subtype=$(echo ${full} | cut -d "_" -f 3) # H1 or none
                
                header_name=$(echo ~{sample_name}_${TYPE}_${segment_subtype})
                echo "header name: $header_name"
                
                # add to assembled_gene_segments.csv
                echo "~{sample_name},${TYPE},${segment},${subtype}" >> ~{sample_name}_irma_assembled_gene_segments.csv

                # this if statement won't work if the PB1 and PB2 are mixed types
                # one will be overwritten in the associative array
                # but also the files would be overwritten in the ammended fasta
                # this should be so rare; but wanted to make a note in the failed logic
                if [ $TYPE == "A" ]; then
                    segment_num=${FluA[$segment]}
                    formatted_name_dict+=( [$segment_num]=$header_name )
                    echo "segment number: $segment_num"
                    echo "fasta header - in formatted_name_dict: ${formatted_name_dict[$segment_num]}"
                elif [ $TYPE == "B" ]; then
                    segment_num=${FluB[$segment]}
                    formatted_name_dict+=( [$segment_num]=$header_name )
                    echo "segment number: $segment_num"
                    echo "fasta header - in formatted_name_dict: ${formatted_name_dict[$segment_num]}"
            
                fi

            done

            # use amended fastas because they ahve the 50x cut off
            # this is new
            echo -e '\n\n\n'
            echo "LOOPING THROUGH AMENDED CONSENSUS FASTAS"
            for file in ~{sample_name}/amended_consensus/*.fa; do
                echo -e '\n'
                echo ${file}

                # grab segment number
                BFN=$(basename ${file%.*})
                segment_number=$(echo $BFN | grep -o '[^_]*$')

                # use associative array to get the formatted name
                header_name=${formatted_name_dict[$segment_number]}
                echo "DEBUG: checking header name pulled form the formatted_name_dict"
                echo "formatted_name_dict: ${formatted_name_dict[$segment_number]}"
                echo "segment number: $segment_number"
                echo "header: $header_name"

                # replace header
                sed -i "s/>.*/>${header_name}/" ${file}

                # replace IUPAC bases with Ns
                for base in ${IUPAC[@]}; do
                    sed -i "/^>/! s/${base}/N/g" $file
                done

                # remove "-" since these represent gaps relative to refernece
                # replace periods with Ns
                sed -i "/^>/! s/-//g" $file
                sed -i "/^>/! s/\./N/g" $file

                # rename file
                new_name=$(echo ${header_name}_irma.fasta)
                mv "${file}" "${new_name}"

                # if HA or NA rename with generic name to use for nextclade
                if [ $segment_number == 4 ]; then 
                    echo "generating HA fasta and HA basename file for nextclade inputs"
                    echo "value in HA_basename.txt: ${header_name}"
                    new_name_HA="HA.fasta"
                    cp "${new_name}" "${new_name_HA}"
                    echo $header_name > "HA_basename.txt"
                fi

                if [ $segment_number == 6 ]; then
                    echo "generating NA fasta and NA basename file for nextclade inputs"
                    echo "value in NA_basename.txt: ${header_name}"
                    new_name_NA="NA.fasta"
                    cp "${new_name}" "${new_name_NA}"
                    echo $header_name > "NA_basename.txt"
                fi


                echo "DEBUG: print contents of final fasta file"
                echo "fasta file name: $new_name"
                cat $new_name

                # add file contents to concatenated fasta file
                cat ${new_name} >> ~{sample_name}_irma_multi.fasta

            done

            echo "\n\n\n"
            echo "creating dummy files if needed"
            # if NA_basename and HA_basename txt files are not created...
            # create dummy file
            if [ ! -f HA_basename.txt ]; then
                echo "creating dummy HA_basename.txt file"
                echo "no HA fasta generated" > "HA_basename.txt"
            fi

            if [ ! -f NA_basename.txt ]; then
                echo "creating dummy NA_basename.txt file"
                echo "no NA fasta generated" > "NA_basename.txt"
            fi

            echo -e '\n\n\n'
            # rename bam and vcf files
            echo "RENAMING BAM AND VCF FILES"
            for file in ~{sample_name}/*{.vcf,.bam}; do
                echo $file
                base_name=$(basename ${file%.*}) # to grab the extenstion

                full=$(basename ${file} | cut -d "." -f 1) # A_HA_H1 or A_NP
                TYPE=$(echo ${full} | cut -d "_" -f 1) # A
                segment_subtype=$(echo ${full} | cut -d "_" -f 2-) # HA_H1 or NP
                segment_subtype=${segment_subtype//_/-} # HA-H1 or NP
                
                extension="${file##*.}"
                
                new_name=$(echo ~{sample_name}_${TYPE}_${segment_subtype}.${extension})
                mv "${file}" "${new_name}"
            done

        else 
            echo "irma assembly fail" | tee irma_qc.txt
            echo "sample_name,flu_type,gene_segment,subtype" > ~{sample_name}_irma_assembled_gene_segments.csv
            echo "~{sample_name},no IRMA assembly generated,none,none" >> ~{sample_name}_irma_assembled_gene_segments.csv

        fi 
        echo -e '\n\n\n'

        echo "RENAMING TABLES AND LOGS"
        # copy read_counts file: path = sample_name/tables/READ_COUNTS.txt
        # copy run_info.tx file: path = sample_name/logs/run_info.txt
        # copy NR counts log: pat = sample_name/logs/NR_COUNTS_log.txt
        # rename with sample name in the file name
        read_counts_fn='~{sample_name}/tables/READ_COUNTS.txt'
        echo "read_counts.txt:"
        cat $read_counts_fn
        echo ""
        new_fn="~{sample_name}_READ_COUNTS.txt"
        mv ${read_counts_fn} ${new_fn}

        echo "read_counts.txt moved:"
        cat $new_fn
        echo ""

        run_info_fn='~{sample_name}/logs/run_info.txt'
        echo "run_info.txt: "
        cat $run_info_fn
        echo ""
        new_fn="~{sample_name}_run_info.txt"
        mv ${run_info_fn} ${new_fn}

        echo "run_info.txt moved: "
        cat $new_fn
        echo ""

        echo -e '\n\n\n'

    >>>

    output {

        # want some of the irma output files
        File? irma_read_counts = "~{sample_name}_READ_COUNTS.txt"
        File? irma_run_info = "~{sample_name}_run_info.txt"
    
        File irma_assembled_gene_segments_csv = "~{sample_name}_irma_assembled_gene_segments.csv"
        # Added '_multi' to file name to differentiate from segment fastas
        File? irma_multifasta = "~{sample_name}_irma_multi.fasta"
        File? HA_fasta = "HA.fasta" # for nextclade
        String? HA_basename_txt = read_string("HA_basename.txt") # for nextclade
        File? NA_fasta = "NA.fasta" # for nextclade
        String? NA_basename_txt = read_string("NA_basename.txt") # for nextclade
        # globs are ordered, so if the diffierent file types all have the same names, these should all be in the same order
        # However this is dependent on all three files being created for every segment and subtype- does that
        # ever not happen? If not, the logic would need to be changed but I don't think it would be difficult
        Array[File] assemblies = glob("*_irma.fasta")
        Array[File] alignments = glob("*.bam")
        Array[File] vcfs = glob("*.vcf")
        String irma_assembly_qc = read_string("irma_qc.txt")

        VersionInfo IRMA_version_info = object{
            software: "IRMA",
            docker: docker,
            version: read_string("VERSION")
        }
        
    }

    runtime {
        docker: docker
        memory: "8 GiB"
        cpu: 2
        disks: "local-disk 50 SSD"
        preemptible: 0
  }
}


task grab_segment_info {
    meta {
        description: "create assembled segment structs"
    }

    input {
        String sample_name
        File fasta
    }

    String base_name = sub(basename(fasta, ".fasta"), "~{sample_name}_", "") # A_HA-H1_irma or A_NP_irma

    command <<<
        echo "base_name"
        echo ~{base_name}

        echo "TYPE"
        echo ~{base_name} | cut -d "_" -f 1 | tee TYPE # A
        echo ""
        echo "Segment"
        echo ~{base_name} | cut -d "_" -f 2 | cut -d "-" -f 1 | tee SEGMENT # HA or 
        echo ""
        echo "segment variable"
        # capture segment as variable to use in if statement
        segment=$(echo ~{base_name} | cut -d "_" -f 2 | cut -d "-" -f 1) 
        echo $segment
        echo ""

        # Does IRMA ever output things with hyphens?
        echo "subtype"

        if [[ $segment == "HA" ]] || [[ $segment == "NA" ]]; then
            echo "yes this HA or NA!"
            echo ~{base_name} | cut -d "_" -f 2 | cut -d "-" -f 2 | tee SUBTYPE # H1 or N1
        else
            echo "" | tee SUBTYPE
        fi
    >>>

    output {
        String type = read_string('TYPE')
        String segment = read_string('SEGMENT')
        String subtype = read_string('SUBTYPE')
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "4 GiB"
        cpu: 4
        disks: "local-disk 50 SSD"
        preemptible: 0
  }
}



task get_irma_subtyping_results {
    meta {
        description: "taking the assembled gene segments info to pull out the type and subtype; added this task to account for potentially mixed types"
    }

    input {
        File irma_assembled_gene_segments_csv
        String sample_name

        File python_script
    }

    command <<<
        python ~{python_script} \
            --irma_assembled_gene_segments_csv "~{irma_assembled_gene_segments_csv}" \
            --sample_name "~{sample_name}" 
    >>>

    output {
        File irma_typing = "~{sample_name}_irma_typing.csv"
        String irma_type = read_string("TYPE.txt")
        String irma_ha_subtype = read_string("HA_SUBTYPE.txt")
        String irma_na_subtype = read_string("NA_SUBTYPE.txt")

    }
        
    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}
