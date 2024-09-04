version 1.0

# define structure
struct VersionInfo {
  String software
  String docker
  String version
}

task ha_nextclade {
    meta {
        description: "run nextclade"
    }

    input {
        File ivar_seg_ha_fasta
        String irma_type
        String irma_ha_subtype
        String sample_name 

    }
    
    String docker = "nextstrain/nextclade:3.8.2"


    command <<<

    # Check the value of irma_type and assign dataset accordingly
    if [ "~{irma_type}" = "A" ]; then
        # Check the value of irma_ha_subtype if irma_type is "A"
        if [ "~{irma_ha_subtype}" = "H1" ]; then
            dataset="flu_h1n1pdm_ha"
        elif [ "~{irma_ha_subtype}" = "H3" ]; then
            dataset="flu_h3n2_ha"
        elif [ "~{irma_ha_subtype}" = "H5" ]; then
            dataset="community/moncla-lab/iav-h5/ha/all-clades"
        else
            echo 'Invalid irma_ha_subtype: "~{irma_ha_subtype}"'
            exit 1  # Exit the script with an error status
        fi
    elif [ "~{irma_type}" = "B" ]; then
        dataset="flu_vic_ha"
    else
        echo "Invalid irma_type: $irma_type"
        exit 1  # Exit the script with an error status
    fi
        

    # run nextclade:
    # 0- caputre nextclade version
    nextclade --version | tee VERSION

    # 1- download the dataset
    nextclade dataset get --name ${dataset} --output-dir "data/flu_ha"

    #2- run nextclade
    nextclade run \
        --input-dataset data/flu_ha --output-all=. "~{ivar_seg_ha_fasta}"

    # 3- rename files
    if [ "~{irma_ha_subtype}" = "H5" ]; then
        mv nextclade.json "~{sample_name}_nextclade_HA.json"
        mv nextclade.csv "~{sample_name}_nextclade_HA.csv"
        mv nextclade.tsv "~{sample_name}_nextclade_HA.tsv"
        mv nextclade.cds_translation.HA.fasta "~{sample_name}_nextcalde.cds_translation.HA.fasta"
    else
        mv nextclade.json "~{sample_name}_nextclade_HA.json"
        mv nextclade.csv "~{sample_name}_nextclade_HA.csv"
        mv nextclade.tsv "~{sample_name}_nextclade_HA.tsv"
        mv nextclade.cds_translation.HA1.fasta "~{sample_name}_nextclade_HA.cds_translation.HA1.fasta"
        mv nextclade.cds_translation.HA2.fasta "~{sample_name}_nextclade_HA.cds_translation.HA2.fasta"
        mv nextclade.cds_translation.SigPep.fasta "~{sample_name}_nextclade_HA.cds_translation.SigPep.fasta"
    fi


    >>>

    output {
        File ha_nextclade_json = "~{sample_name}_nextclade_HA.json"
        File ha_nextclade_tsv = "~{sample_name}_nextclade_HA.tsv"
        File? ha_nextclade_HA1_translation_fasta = "~{sample_name}_nextclade_HA.cds_translation.HA1.fasta"
        File? ha_nextclade_HA2_translation_fasta = "~{sample_name}_nextclade_HA.cds_translation.HA2.fasta"
        File? ha_nextclade_SigPep_translation_fasta = "~{sample_name}_nextclade_HA.cds_translation.SigPep.fasta"
        File? ha_nextclade_translation_fasta = "~{sample_name}_nextcalde.cds_translation.HA.fasta"

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


task na_nextclade {
    meta {
        description: "run nextclade"
    }

    input {
        File ivar_seg_na_fasta
        String irma_type
        String irma_na_subtype
        String sample_name  

    }
    String docker = "nextstrain/nextclade:3.8.2"
    

    command <<<

    # Check the value of irma_type and assign dataset accordingly
    if [ "~{irma_type}" = "A" ]; then
        # Check the value of irma_ha_subtype if irma_type is "A"
        if [ "~{irma_na_subtype}" = "N1" ]; then
            dataset="flu_h1n1pdm_na"
        elif [ "~{irma_na_subtype}" = "N2" ]; then
            dataset="flu_h3n2_na"
        else
            echo 'Invalid irma_na_subtype: "~{irma_na_subtype}"'
            exit 1  # Exit the script with an error status

        fi
    elif [ "~{irma_type}" = "B" ]; then
        dataset="flu_vic_na"
    else
        echo "Invalid irma_type: ~{irma_type}"
        exit 1  # Exit the script with an error status
    fi
        

    # run nextclade:
    # 0- caputre nextclade version
    nextclade --version | tee VERSION

    # 1- download the dataset
    nextclade dataset get --name ${dataset} --output-dir "data/flu_na"

    #2- run nextclade
    nextclade run --input-dataset data/flu_na --output-all=. "~{ivar_seg_na_fasta}"

    # 3- rename files
    mv nextclade.json "~{sample_name}_nextclade_NA.json"
    mv nextclade.csv "~{sample_name}_nextclade_NA.csv"
    mv nextclade.tsv "~{sample_name}_nextclade_NA.tsv"
    mv nextclade.cds_translation.NA.fasta "~{sample_name}_nextclade_NA.cds_translation.NA.fasta"
    


    >>>

    output {
        File na_nextclade_json = "~{sample_name}_nextclade_NA.json"
        File na_nextclade_tsv = "~{sample_name}_nextclade_NA.tsv"
        File na_nextclade_translation_fasta = "~{sample_name}_nextclade_NA.cds_translation.NA.fasta"

        VersionInfo na_nextclade_version_info = object{
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