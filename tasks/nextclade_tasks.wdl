version 1.0

# define structure
# struct VersionInfo {
#   String software
#   String docker
#   String version
# }

import "../tasks/capture_version_tasks.wdl" as capture_version

task nextclade_HA {
    meta {
        description: "run nextclade"
    }

    input {
        File fasta
        String type
        String segment
        String subtype
        String sample_name
        String base_name
    }
    
    String docker = "nextstrain/nextclade:3.8.2"

    # use subtype_name to determine flu b dataset
    Map[String, String] a_dict = {
        "H1": "flu_h1n1pdm_ha",
        "H3": "flu_h3n2_ha",
        "H5": "community/moncla-lab/iav-h5/ha/all-clades",
        "N1" : "flu_h1n1pdm_na",
        "N2" : "flu_h3n2_na"
    }

    # use segment_name to determine flu b dataset
    Map[String, String] b_dict = {
        "NA" : "flu_vic_na",
        "HA" : "flu_vic_ha"
    }

    # # select dataset
    # if ("~{type}" == "A" ) {
    #        String dataset = a_dict["~{subtype}"]
    # }

    # if ("~{type}" == 'B') {
    #     String dataset = b_dict["~{segment}"]
    # }
    String dataset = if "~{type}" == "A" then a_dict["~{subtype}"] else b_dict["~{subtype}"]

    # String base_name = "~{sample_name}_~{type}_~{segment}-~{subtype}"
    
    command <<<
        # # Check the value of irma_type and assign dataset accordingly
        # if [[ "~{type}" == "A" ]]; then
        #     dataset = a_dict["~{subtype}"]
        # elif [[ "~{type}" == "B" ]]; then
        #     dataset = b_dict["~{segment}"]
        # else
        #     echo "Invalid irma_type: $type"
        #     exit 1  # Exit the script with an error status
        # fi

        echo "datasets selected"
        echo ~{dataset}


        # run nextclade:
        # 0- capture nextclade version
        nextclade --version | tee VERSION

        # 1- download the dataset
        nextclade dataset get --name "~{dataset}" --output-dir "data/flu"

        echo "got dataset"

        #2- run nextclade
        nextclade run --input-dataset data/flu --output-all=output/ --output-basename "~{base_name}" "~{fasta}"
        echo "ran nextclade"
    >>>

    output {
        File nextclade_HA_json = "output/~{base_name}.json"
        File nextclade_HA_tsv = "output/~{base_name}.tsv"
        File? nextclade_HA_translation_fasta = "output/~{base_name}.cds_translation.HA.fasta" # H5 only
        File? nextclade_HA1_translation_fasta = "output/~{base_name}.cds_translation.HA1.fasta" # H1, H3
        File? nextclade_HA2_translation_fasta = "output/~{base_name}.cds_translation.HA2.fasta" # H1, H3
        File? nextclade_SigPep_translation_fasta = "output/~{base_name}.cds_translation.SigPep.fasta" # H1, H3
        # File? nextclade_NA_translation_fasta = "output/~{base_name}.cds_translation.NA.fasta" # N1, N3

        VersionInfo nextclade_version_info = object{
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


task nextclade_NA {
    meta {
        description: "run nextclade"
    }

    input {
        File fasta
        String type
        String segment
        String subtype
        String sample_name
        String base_name
    }
    
    String docker = "nextstrain/nextclade:3.8.2"

    # use subtype_name to determine flu b dataset
    Map[String, String] a_dict = {
        "H1": "flu_h1n1pdm_ha",
        "H3": "flu_h3n2_ha",
        "H5": "community/moncla-lab/iav-h5/ha/all-clades",
        "N1" : "flu_h1n1pdm_na",
        "N2" : "flu_h3n2_na"
    }

    # use segment_name to determine flu b dataset
    Map[String, String] b_dict = {
        "NA" : "flu_vic_na",
        "HA" : "flu_vic_ha"
    }

    # # select dataset
    # if ("~{type}" == "A" ) {
    #        String dataset = a_dict["~{subtype}"]
    # }

    # if ("~{type}" == 'B') {
    #     String dataset = b_dict["~{segment}"]
    # }
    String dataset = if "~{type}" == "A" then a_dict["~{subtype}"] else b_dict["~{subtype}"]

    # String base_name = "~{sample_name}_~{type}_~{segment}-~{subtype}"
    
    command <<<
        # # Check the value of irma_type and assign dataset accordingly
        # if [[ "~{type}" == "A" ]]; then
        #     dataset = a_dict["~{subtype}"]
        # elif [[ "~{type}" == "B" ]]; then
        #     dataset = b_dict["~{segment}"]
        # else
        #     echo "Invalid irma_type: $type"
        #     exit 1  # Exit the script with an error status
        # fi

        echo "datasets selected"
        echo ~{dataset}


        # run nextclade:
        # 0- capture nextclade version
        nextclade --version | tee VERSION

        # 1- download the dataset
        nextclade dataset get --name "~{dataset}" --output-dir "data/flu"

        echo "got dataset"

        #2- run nextclade
        nextclade run --input-dataset data/flu --output-all=output/ --output-basename "~{base_name}" "~{fasta}"
        echo "ran nextclade"
    >>>

    output {
        File nextclade_NA_json = "output/~{base_name}.json"
        File nextclade_NA_tsv = "output/~{base_name}.tsv"
        # File? nextclade_HA_translation_fasta = "output/~{base_name}.cds_translation.HA.fasta" # H5 only
        # File? nextclade_HA1_translation_fasta = "output/~{base_name}.cds_translation.HA1.fasta" # H1, H3
        # File? nextclade_HA2_translation_fasta = "output/~{base_name}.cds_translation.HA2.fasta" # H1, H3
        # File? nextclade_SigPep_translation_fasta = "output/~{base_name}.cds_translation.SigPep.fasta" # H1, H3
        File? nextclade_NA_translation_fasta = "output/~{base_name}.cds_translation.NA.fasta" # N1, N3

        VersionInfo nextclade_version_info = object{
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

