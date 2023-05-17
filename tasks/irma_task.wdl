version 1.0

task irma {
    meta {
        description: "runs CDC's IRMA for FLU. For more info on IRMA visit: https://wonder.cdc.gov/amd/flu/irma/configuration.html. Modified from theiagen genomics - public helath viral genomics (theiaCov PE workflow). This task runs the default config of IRMA meaning it uses: \ALIGN_PROG=SAM \DEL_TYPE=''"
    }

    input {
        String sample_name
        File fastq_R1
        File? fastq_R2
        String read_type
        String module = "FLU"
        String docker = "staphb/irma:1.0.3"

    }
    command <<<
        # grab version
        IRMA | head -n1 | awk -F' ' '{ print "IRMA " $5 }' | tee VERSION
        version=$(cat VERSION)

        # run IRMA
        if [ ~{read_type} == 'paired' ]; then
            IRMA ~{module} ~{fastq_R1} ~{fastq_R2} ~{sample_name}
        elif [ ~{read_type} == 'single' ]; then
            IRMA ~{module} ~{fastq_R1} ~{sample_name}
        fi
        
        # determine if assemly was successful
        if compgen -G "~{sample_name}/*.fasta"; then
            echo "sample_name,flu_type,gene_segment,subtype" > ~{sample_name}_irma_assembled_gene_segments.csv
            for file in ~{sample_name}/*.fasta; do
                # grab type
                TYPE=$(head -n 1 $file | cut -d "_" -f 1 | cut -d ">" -f 2)

                #grab gene segment
                gene_segment=$(head -n 1 $file | cut -d "_" -f 2 )

                # grab subtype (ok if doesn't exist)
                subtype=$(head -n 1 $file | cut -d "_" -f 3 )

                echo "~{sample_name},${TYPE},${gene_segment},${subtype}" >> ~{sample_name}_irma_assembled_gene_segments.csv
            done

            # rename header and file name for fasta
            ## also create an array of the segment names
            for file in ~{sample_name}/*.fasta; do
                # grab base name and drop .fasta
                segment=$(basename ${file%.*})
                # echo $segement >> segment_list.txt
                header_name=$(echo ~{sample_name}_${segment})
                sed -i "s/>.*/>${header_name}/" ${file}

                # rename file
                new_name=$(echo ~{sample_name}_$(basename $file)_irma)
                mv "${file}" "${new_name}"
            done

            # rename bam and vcf files
            for file in ~{sample_name}/*{.vcf,.bam,.bai}; do
                new_name=$(echo ~{sample_name}_$(basename $file))
                mv "${file}" "${new_name}"
            done

        
        else 
            echo "sample_name,flu_type,gene_segment,subtype" > ~{sample_name}_irma_assembled_gene_segments.csv
            echo "~{sample_name},no IRMA assembly generated,none,none" >> ~{sample_name}_irma_assembled_gene_segments.csv

        fi 

        # create an output file with all the irma info
        echo "irma_version,irma_module,irma_docker" > 'irma_runtime.csv'
        echo "${version},~{module},~{docker}" >> 'irma_runtime.csv'

    >>>

    output {

        File irma_assembled_gene_segments_csv = "~{sample_name}_irma_assembled_gene_segments.csv"
        Array[File] irma_assemblies = glob("~{sample_name}*.fasta")
        Array[File] irma_bam_files = glob("~{sample_name}*.bam")
        Array[File] irma_vcfs = glob("~{sample_name}*.vcf")
        File irma_runtime_csv = "irma_runtime.csv"
        String irma_version = read_string("VERSION")
        String irma_docker = "~{docker}"
        String irma_module = "~{module}"
        
    }

    runtime {
        docker: "~{docker}"
        memory: "8 GiB"
        cpu: 2
        disks: "local-disk 50 SSD"
        preemptible: 0
  }
}

task irma_subtyping_results {
    meta {
        description: "taking the assembled gene segments info to pull out the type and subtype; added this task to account for potentially mixed types"
    }

    input {
        File irma_assembled_gene_segments_csv
        String sample_name
        File irma_runtime_csv

        File python_script
    }

    command <<<
        python ~{python_script} \
            --irma_assembled_gene_segments_csv "~{irma_assembled_gene_segments_csv}" \
            --sample_name "~{sample_name}" \
            --irma_runtime_csv "~{irma_runtime_csv}"
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
