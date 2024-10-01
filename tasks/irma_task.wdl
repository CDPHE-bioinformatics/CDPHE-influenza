version 1.0

# define structure
struct VersionInfo {
  String software
  String docker
  String version
}

task perform_assembly_irma {
    meta {
        description: "runs CDC's IRMA for FLU. For more info on IRMA visit: https://wonder.cdc.gov/amd/flu/irma/configuration.html. Modified from theiagen genomics - public helath viral genomics (theiaCov PE workflow). This task runs the default config of IRMA meaning it uses: \ALIGN_PROG=SAM \DEL_TYPE=''"
    }

    input {
        String sample_name
        File fastq_R1
        File? fastq_R2
    }

    String docker = "cdcgov/irma:v1.1.5"

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

                full=$(basename ${file%.*} | cut -d "." -f 1) # A_HA_H1 or A_NP
                TYPE=$(echo ${full} | cut -d "_" -f 1) # A
                segment=$(echo ${full} | cut -d "_" -f 2) # HA or NP
                subtype=$(echo ${full} | cut -d "_" -f 3) # H1 or none
                echo "assembled_gene_segments.csv"
                echo $full $TYPE $segment $subtype 
                   
                echo "~{sample_name},${TYPE},${segment},${subtype}" >> ~{sample_name}_irma_assembled_gene_segments.csv
            done

            # rename header and file name for fasta
            ## also create an array of the segment names
            for file in ~{sample_name}/*.fasta; do
                # grab base name and drop .fasta
                full=$(basename ${file} | cut -d "." -f 1) # A_HA_H1 or A_NP
                TYPE=$(echo ${full} | cut -d "_" -f 1) # A
                segment_subtype=$(echo ${full} | cut -d "_" -f 2-) # HA_H1 or NP
                segment_subtype=${segment_subtype//_/-} # HA-H1 or NP
                # echo $segement >> segment_list.txt
                header_name=$(echo ~{sample_name}_${TYPE}_${segment_subtype})

                echo "fasta file variables"
                echo $full $TYPE $segment_subtype 
                echo $header_name

                sed -i "s/>.*/>${header_name}/" ${file}
                # add file contents to concatenated fasta file
                cat ${file} >> ~{sample_name}_irma_multi.fasta

                # rename file
                new_name=$(echo ~{sample_name}_${TYPE}_${segment_subtype}_irma.fasta)
                mv "${file}" "${new_name}"
                echo $new_name

            
            done

            # rename bam and vcf files
            for file in ~{sample_name}/*{.vcf,.bam}; do
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
            echo "sample_name,flu_type,gene_segment,subtype" > ~{sample_name}_irma_assembled_gene_segments.csv
            echo "~{sample_name},no IRMA assembly generated,none,none" >> ~{sample_name}_irma_assembled_gene_segments.csv

        fi 
    >>>

    output {

        File irma_assembled_gene_segments_csv = "~{sample_name}_irma_assembled_gene_segments.csv"
        # Added '_multi' to file name to differentiate from segment fastas
        File? irma_multifasta = "~{sample_name}_irma_multi.fasta"
        # globs are ordered, so if the different file types all have the same names, these should all be in the same order
        # However this is dependent on all three files being created for every segment and subtype- does that
        # ever not happen? If not, the logic would need to be changed but I don't think it would be difficult
        Array[File] assemblies = glob("*_irma.fasta")
        Array[File] alignments = glob("*.bam")
        Array[File] vcfs = glob("*.vcf")

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

    String base_name = sub(basename(fasta, ".fasta"), "~{sample_name}_", "")

    command <<<
        # base file name looks like
        # {sample_name}_A_HA-H1 --> A_HA-H1
        # {sample_name}_A_NP --> A_NP

        echo ~{base_name} | cut -d "_" -f 1 | tee TYPE # A
        echo ~{base_name} | cut -d "_" -f 2 | tee SEGMENT # HA or NP
        # Does IRMA ever output things with hyphens?
        echo ~{base_name} | cut -d "-" -f 3 | tee SUBTYPE # H1 or N1
    >>>

    output {
        String type = read_lines('TYPE')
        String segment = read_lines('SEGMENT')
        String subtype = read_lines('SUBTYPE')
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
