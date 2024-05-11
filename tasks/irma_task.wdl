version 1.0

task irma {
    meta {
        description: "runs CDC's IRMA for FLU. For more info on IRMA visit: https://wonder.cdc.gov/amd/flu/irma/configuration.html. Modified from theiagen genomics - public helath viral genomics (theiaCov PE workflow). This task runs the default config of IRMA meaning it uses: \ALIGN_PROG=SAM \DEL_TYPE=''"
    }

    input {
        String sample_name
        File fastq_R1
        File? fastq_R2
        # String module = "FLU"
        # String docker = "staphb/irma:1.0.3"

    }
    command <<<
        # grab version
        IRMA | head -n1 | awk -F' ' '{ print "IRMA " $5 }' | tee VERSION
        version=$(cat VERSION)

        # run IRMA
        IRMA FLU ~{fastq_R1} ~{fastq_R2} ~{sample_name}
        
        # determine if assembly was successful
        if compgen -G "~{sample_name}/*.fasta"; then
            echo "sample_name,flu_type,gene_segment,subtype" > ~{sample_name}_irma_assembled_gene_segments.csv
            for file in ~{sample_name}/*.fasta; do
            echo ${file}
                # grab type
                # read in header of fasta file and grab type (e.g. A,B), 
                # gene_segment (e.g. HA, NA, PB1) and subtype (e.g. N2, H1) (ok if subtype doesn't exist)

                segment=$(basename ${file%.*} | cut -d "." -f 1)
                TYPE=$(echo ${segment} | cut -d "_" -f 1)
                gene_segment=$(echo ${segment} | cut -d "_" -f 2)
                subtype=$(echo ${segment} | cut -d "_" -f 3)
                   
                echo "~{sample_name},${TYPE},${gene_segment},${subtype}" >> ~{sample_name}_irma_assembled_gene_segments.csv
            done

            # rename header and file name for fasta
            ## also create an array of the segment names
            for file in ~{sample_name}/*.fasta; do
                # grab base name and drop .fasta
                segment=$(basename ${file} | cut -d "." -f 1)
                # echo $segement >> segment_list.txt
                header_name=$(echo ~{sample_name}_${segment})
                sed -i "s/>.*/>${header_name}/" ${file}

                # add file contents to concatenated fasta file
                cat ${fiile} >> ~{sample_name}_all_assembled_segments.fasta

                # rename file
                new_name=$(echo ~{sample_name}_${segment}_irma.fasta)
                mv "${file}" "${new_name}"

            
            done

            # rename bam and vcf files
            for file in ~{sample_name}/*{.vcf,.bam,.bai}; do
                base_name=$(basename ${file%.*})
                extension="${file##*.}"
                new_name=$(echo ~{sample_name}_${base_name}.${extension})
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
        File? irma_all_assembled_segments_fasta = "~{sample_name}_all_assembled_segments.fasta"
        
        # assemblies
        File? irma_seg_ha_fasta = select_first(glob("~{sample_name}_HA*.fasta"))
        File? irma_seg_na_fasta = select_first(glob("~{sample_name}_HA*.fasta"))
        File? irma_seg_pb1_fasta = select_first(glob("~{sample_name}_HA*.fasta"))
        File? irma_seg_pb2_fasta = select_first(glob("~{sample_name}_HA*.fasta"))
        File? irma_seg_np_fasta = select_first(glob("~{sample_name}_HA*.fasta"))
        File? irma_seg_pa_fasta = select_first(glob("~{sample_name}_HA*.fasta"))
        File? irma_seg_ns_fasta = select_first(glob("~{sample_name}_HA*.fasta"))
        File? irma_seg_mp_fasta = select_first(glob("~{sample_name}_HA*.fasta"))

        # alignments
        File? irma_seg_ha_bam = select_first(glob("~{sample_name}_HA*.bam"))
        File? irma_seg_na_bam = select_first(glob("~{sample_name}_HA*.bam"))
        File? irma_seg_pb1_bam = select_first(glob("~{sample_name}_HA*.bam"))
        File? irma_seg_pb2_bam = select_first(glob("~{sample_name}_HA*.bam"))
        File? irma_seg_np_bam = select_first(glob("~{sample_name}_HA*.bam"))
        File? irma_seg_pa_bam = select_first(glob("~{sample_name}_HA*.bam"))
        File? irma_seg_ns_bam = select_first(glob("~{sample_name}_HA*.bam"))
        File? irma_seg_mp_bam = select_first(glob("~{sample_name}_HA*.bam"))

        # vcfs
        File? irma_seg_ha_vcf = select_first(glob("~{sample_name}_HA*.vcf"))
        File? irma_seg_na_vcf = select_first(glob("~{sample_name}_HA*.vcf"))
        File? irma_seg_pb1_vcf = select_first(glob("~{sample_name}_HA*.vcf"))
        File? irma_seg_pb2_vcf = select_first(glob("~{sample_name}_HA*.vcf"))
        File? irma_seg_np_vcf = select_first(glob("~{sample_name}_HA*.vcf"))
        File? irma_seg_pa_vcf = select_first(glob("~{sample_name}_HA*.vcf"))
        File? irma_seg_ns_vcf = select_first(glob("~{sample_name}_HA*.vcf"))
        File? irma_seg_mp_vcf = select_first(glob("~{sample_name}_HA*.vcf"))

        # runtime
        File irma_runtime_csv = "irma_runtime.csv"
        String irma_version = read_string("VERSION")
        String irma_docker = "~{docker}"
        String irma_module = "~{module}"
        
    }

    runtime {
        docker: "staphb/irma:1.0.3"
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
