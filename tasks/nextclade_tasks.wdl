version 1.0

import "../tasks/capture_version_tasks.wdl" as capture_version

task nextclade_HA {
    meta {
        description: "run nextclade"
    }

    input {
        File? fasta
        String sample_name
        String? base_name
    }
    
    String docker = "nextstrain/nextclade:3.8.2"
    
    command <<<
        
        # figure out the correct dataset to use
        # grab base_name, type, segment and subtype from fasta header
        declare -A HA_datasets=(['A_HA-H1']="flu_h1n1pdm_ha" ['A_HA-H3']="flu_h3n2_ha" ['A_HA-H5']="community/moncla-lab/iav-h5/ha/all-clades" ['B_HA']="flu_vic_ha")

        # base_name=$(basename ${file} | cut -d "." -f 1) # sample_A_HA-H1 or sample_A_NP
        echo "base_name string: ~{base_name}"
        echo ~{base_name} > temp.txt
        sed -i "s/~{sample_name}_//g" temp.txt
        base=$(sed -n "1p" temp.txt) # A_HA-H3 B_HA etc.
        echo "key for dataset selection: ${base}"
        
        dataset=${HA_datasets[${base}]}

        echo "datasets selected: ${dataset}"


        # run nextclade:
        # 0- capture nextclade version
        nextclade --version | tee VERSION

        # 1- download the dataset
        nextclade dataset get --name "${dataset}" --output-dir "data/flu"

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
        File? fasta
        String sample_name
        String? base_name
    }
    
    String docker = "nextstrain/nextclade:3.8.2"

    
    command <<<

        # figure out the correct dataset to use
        # grab base_name, type, segment and subtype from fasta header
        declare -A NA_datasets=(['A_NA-N1']="flu_h1n1pdm_na" ['A_NA-N2']="flu_h3n2_na" ['B_NA']="flu_vic_na")

        # base_name=$(basename ${file} | cut -d "." -f 1) # sample_A_HA-H1 or sample_A_NP
        echo "base_name string: ~{base_name}"
        echo ~{base_name} > temp.txt
        sed -i "s/~{sample_name}_//g" temp.txt
        base=$(sed -n "1p" temp.txt) # A_HA-H3 B_HA etc.
        echo "key for dataset selection: ${base}"
        
        dataset=${NA_datasets[${base}]}

        echo "datasets selected: ${dataset}"


        # run nextclade:
        # 0- capture nextclade version
        nextclade --version | tee VERSION

        # 1- download the dataset
        nextclade dataset get --name "${dataset}" --output-dir "data/flu"

        echo "got dataset"

        #2- run nextclade
        nextclade run --input-dataset data/flu --output-all=output/ --output-basename "~{base_name}" "~{fasta}"
        echo "ran nextclade"
    >>>

    output {
        File nextclade_NA_json = "output/~{base_name}.json"
        File nextclade_NA_tsv = "output/~{base_name}.tsv"
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

