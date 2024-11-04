# define struct
struct VersionInfo {
  String software
  String docker
  String version
}

# begin workflow
workflow influenza_assembly{

    input {
        File fastq1 
    }


    # perform assembly
    call assembly {
        input: 
            fastq1 = fastq1
    }

    # if assembly passes move on 
    if (assembly.assembly_qc == 'pass'){

        # scatter over genes
        scatter (gene in genes) { 
            
            # perform bam_stats on all genes
            call bam_stats{}

            # if gene HA assembled then perform nextclade
            if (gene == 'HA' || gene == 'NA') {

                call nextclade{}
            }

        }
    }

    # here is where the problems start
    # There will always be 1 assembly_version_info,
    # Assumming there were some assembled genes then there should be
    # at least one bam_stats_version_info struct array
    # Assiming HA and NA was assembled then there would be a nextclade_version_info struct array 
    # however both bam_stats_version_info and nextclade_version_info may not exist
    
    # select all must only have one data type, in this case a VersionInfo Struct
    # so arrays created through the scatter must be "de-arrayed" using something like 
    # and select_frist() doesn't accept Array[Struct]?
    Array[VersionInfo] version_array = select_all(
        assembly.assembly_version_info,
        select_first(bam_stats.bam_stats_version_info),
        select_first(nextclade.nextclade_verion_info)
    )

    task write_versions{
        input:
            version_array = select_all(version_array)
    }

    output {

    }
}

###########

# example task
# every task outputs a VersionInfo Struct
task assembly {
    ....

    output{
        Array[String?] assembled_genes = read_lines("assembled_genes.txt")
        # this array would like like ['HA', 'NA', 'MP'] 
        # or could look like ['MP', 'PB1'] 
        # or could be empty
        # just depends on what is successuflly assembled

        VersionInfo assembly_version_info = object{
            software: "assembly_program",
            docker: "assembly_docker",
            version: "assembly_version"
        }
    }

}