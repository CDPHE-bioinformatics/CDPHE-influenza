version 1.0

# define structure
struct VersionInfo {
  String software
  String docker
  String version
}

task nextclade {
    meta {
        description: "run nextclade"
    }

    input {
        File ivar_seg_fasta
        String segment_name
        String subtype_name
        String sample_name 
    }
    
    String docker = "nextstrain/nextclade:3.8.2"
    Map[String, String] a_dict = {
        "H1": "flu_h1n1pdm_ha",
        "H3": "flu_h3n2_ha",
        "H5": "community/moncla-lab/iav-h5/ha/all-clades"
    }
    String base_name = "~{sample_name}_nextclade_~{irma_subtype}"
    
    command <<<
        # Check the value of irma_type and assign dataset accordingly
        if [ "~{irma_type}" = "A" ]; then
            dataset = a_dict["~{irma_type}"]
        elif [ "~{irma_type}" = "B" ]; then
            dataset="flu_vic_ha"
        else
            echo "Invalid irma_type: $irma_type"
            exit 1  # Exit the script with an error status
        fi

        # run nextclade:
        # 0- capture nextclade version
        nextclade --version | tee VERSION

        # 1- download the dataset
        nextclade dataset get --name ${dataset} --output-dir "data/flu"

        #2- run nextclade
        nextclade run --input-dataset data/flu --output-all ~{ivar_seg_fasta} --output-basename ~{base_name}
    >>>

    output {
        File nextclade_json = "~{base_name}.json"
        File nextclade_tsv = "~{base_name}.tsv"
        File? nextclade_translation_fasta = "~{base_name}.cds_translation.~{irma_subtype}.fasta"
        File? nextclade_HA1_translation_fasta = "~{base_name}.cds_translation.HA1.fasta"
        File? nextclade_HA2_translation_fasta = "~{base_name}.cds_translation.HA2.fasta"
        File? nextclade_SigPep_translation_fasta = "~{base_name}.cds_translation.SigPep.fasta"

        VersionInfo ha_nextclade_version_info = object{
            software: "nextclade",
            docker: docker,
            version: read_string("VERSION")
        }
    }

    runtime {
        cpu:    4
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker

    }
}

