version 1.0

# import tasks
import "../tasks/preprocess_tasks.wdl" as fastq_preprocess
import "../tasks/irma_task.wdl" as irma_task
import "../tasks/ivar_task.wdl" as ivar
import "../tasks/assembly_qc_tasks.wdl" as assembly_qc
import "../tasks/transfer_tasks.wdl" as transfer
import "../tasks/nextclade_tasks.wdl" as nextclade_tasks
import "../tasks/capture_version_tasks.wdl" as capture_version

# define struct
struct VersionInfo {
  String software
  String docker
  String version
}

# begin workflow
workflow influenza_assembly {

    input {
        String sample_name
        String project_name
        File fastq_R1
        File fastq_R2
        File contam_fasta
        String out_bucket_path

        # python scripts
        File concat_preprocess_qc_metrics_py
        File irma_subtyping_results_py
        File calc_percent_coverage_py
        File concat_assembly_qc_metrics_py
        File capture_version_py
    }

    # 1 - Preprocess QC raw fastq files
    call fastq_preprocess.assess_quality_fastqc as fastqc_raw {
        input:
            sample_name = sample_name,
            fastq_R1 = fastq_R1,
            fastq_R2 = fastq_R2
    }

    call fastq_preprocess.filter_reads_seqyclean as seqyclean {
        input:
            sample_name = sample_name,
            fastq_R1 = fastq_R1,
            fastq_R2 = fastq_R2,
            contam_fasta = contam_fasta  
    }

    call fastq_preprocess.assess_quality_fastqc as fastqc_cleaned {
        input:
            sample_name = sample_name,
            fastq_R1 = seqyclean.fastq_R1_cleaned,
            fastq_R2 = seqyclean.fastq_R2_cleaned
    }

    # concatenate all preprocess qc metrics into single file
    call fastq_preprocess.concat_preprocess_qc_metrics as concat_preprocess_qc_metrics{
        input:
            python_script = concat_preprocess_qc_metrics_py,
            sample_name = sample_name,

            total_reads_R1_raw = fastqc_raw.total_reads_R1,
            total_reads_R2_raw = fastqc_raw.total_reads_R2,
            read_length_R1_raw = fastqc_raw.read_length_R1,
            read_length_R2_raw = fastqc_raw.read_length_R2,
            read_pairs_raw = fastqc_raw.read_pairs,

            total_reads_R1_cleaned = fastqc_cleaned.total_reads_R1,
            total_reads_R2_cleaned = fastqc_cleaned.total_reads_R2,
            read_length_R1_cleaned = fastqc_cleaned.read_length_R1,
            read_length_R2_cleaned = fastqc_cleaned.read_length_R2,
            read_pairs_cleaned = fastqc_cleaned.read_pairs
    }

    # 2- run irma
    call irma_task.perform_assembly_irma as irma {
        input:
            sample_name = sample_name,
            fastq_R1 = seqyclean.fastq_R1_cleaned,
            fastq_R2 = seqyclean.fastq_R2_cleaned
    }

    call irma_task.get_irma_subtyping_results as irma_subtyping_results {
        input:
            irma_assembled_gene_segments_csv = irma.irma_assembled_gene_segments_csv,
            sample_name = sample_name,
            python_script = irma_subtyping_results_py
    }

    # add in if irma assembly was successful

    if (irma.irma_assembly_qc == "irma assembly fail") {
        String exit_reason = "no IRMA assemlby generated"
        
        call exit_wdl {
            input:
            exit_reason = exit_reason
        }
    }

    # only allow workflow to proceed if irma assembly is pass
    if (irma.irma_assembly_qc == "irma assembly pass") {
    # Scatter each segment assembled, run samtools, ivar consensus, calculate percent_coverage 
    # if HA or NA run nextclade

        Array[Int] indexes = range(length(irma.assemblies))
        scatter (idx in indexes) {
            File fasta = irma.assemblies[idx]
            File bam = irma.alignments[idx]
            File vcf = irma.vcfs[idx]
            String base_name = sub(basename(fasta, ".fasta"), "_irma", "")
            
            call irma_task.grab_segment_info as grab_segment_info {
                input:
                    sample_name = sample_name,
                    fasta = fasta
            }

            
            # # grab subtyping info
            # if (grab_segment_info.segment == 'HA') {
            #     String? HA_subtype = grab_segment_info.subtype
            # }
            # if (grab_segment_info.segment == 'NA') {
            #     String? NA_subtype = grab_segment_info.subytpe
            # }

            call assembly_qc.calc_bam_stats_samtools as bam_stats {
                input:
                    sample_name = sample_name,
                    bam_file = bam, 
                    segment_name = grab_segment_info.segment,
                    base_name = base_name
            }

            call ivar.call_consensus_ivar as ivar_consensus {
                input:
                    bam_file = bam,
                    sample_name = sample_name,
                    base_name = base_name
            }

            call assembly_qc.calc_percent_coverage as calc_percent_coverage{
                input:
                    fasta_file = fasta,
                    python_script = calc_percent_coverage_py,
                    sample_name = sample_name,
                    segment = grab_segment_info.segment,
                    base_name = base_name
            }

            # create boolean for ha and na to see if that works
            Boolean is_ha = grab_segment_info.segment == "HA"
            Boolean is_na = grab_segment_info.segment == "NA"

            if (is_ha || is_na) {
                call nextclade_tasks.nextclade as nextclade {
                    input:
                        ivar_fasta = ivar_consensus.ivar_consensus_fasta,
                        type = grab_segment_info.type,
                        segment = grab_segment_info.segment,
                        subtype = grab_segment_info.subtype,
                        sample_name = sample_name, 
                        base_name = base_name
                }
            }
        }
        # create arrays to better handle groups of files

        # IRMA - fasta, bam, vcf
        Array[File] irma_fasta_array = irma.assemblies
        Array[File] irma_bam_array = irma.alignments
        Array[File] irma_vcf_array = irma.vcfs

        # IVAR - fasta
        Array[File] ivar_fasta_array = ivar_consensus.ivar_consensus_fasta
        File ivar_parameters_file = select_first(ivar_consensus.ivar_parameters)

        # Samtools - mapped reads csv, sorted bam
        Array[File] bam_stats_csv_array = bam_stats.bam_stats_csv
        Array[File] sorted_bam_array = bam_stats.sorted_bam
        Array[File] sorted_bai_array = bam_stats.sorted_bai

        # percent coverage - percent coverage csv
        Array[File] percent_coverage_csv_array = calc_percent_coverage.percent_coverage_csv

        # nextclade
        Array[File?] nextclade_json_array = select_all(nextclade.nextclade_json)
        Array[File?] nextclade_tsv_array = select_all(nextclade.nextclade_tsv)
        Array[File?] nextclade_HA1_translation_fasta = select_all(nextclade.nextclade_HA1_translation_fasta)
        Array[File?] nextclade_HA2_translation_fasta = select_all(nextclade.nextclade_HA2_translation_fasta)
        Array[File?] nextclade_SigPep_translation_fasta = select_all(nextclade.nextclade_SigPep_translation_fasta)
        Array[File?] nextclade_HA_translation_fasta = select_all(nextclade.nextclade_HA_translation_fasta)
        Array[File?] nextclade_NA_translation_fasta = select_all(nextclade.nextclade_NA_translation_fasta)

        # concantenate post assembly qc metrics (coverage, depth) into a single file
        call assembly_qc.concat_assembly_qc_metrics as concat_assembly_qc_metrics{
            input:
                python_script = concat_assembly_qc_metrics_py,
                sample_name = sample_name,
                percent_coverage_csv_array = percent_coverage_csv_array,
                bam_stats_csv_array = bam_stats_csv_array
        }

        call assembly_qc.make_multifasta as make_ivar_multifasta{
            input:
                fasta_array = ivar_fasta_array,
                sample_name = sample_name
        }
        

        # create array of structs
        Array[VersionInfo] version_array = select_all([
            fastqc_raw.fastqc_version_info,
            seqyclean.seqyclean_version_info,
            irma.IRMA_version_info,
            select_first(bam_stats.samtools_version_info),
            select_first(ivar_consensus.ivar_version_info),
            select_first(ivar_consensus.samtools_version_info),
            select_all(nextclade.nextclade_version_info)[0]
        ])
        

        # capture version
        call capture_version.capture_workflow_version as capture_workflow_version{
            input:
        }

        call capture_version.capture_task_version as capture_task_version {
            input:
                version_array = version_array,
                workflow_name = "influenza_assembly",
                workflow_version = capture_workflow_version.workflow_version,
                project_name = project_name,
                sample_name = sample_name,
                analysis_date = capture_workflow_version.analysis_date,
                capture_version_py = capture_version_py
        }
        
        
        # 6 - Transfer some intermediate files and all final files to gcp bucket
        call transfer.transfer_assembly_wdl as transfer_assembly_wdl {
            input:
                workflow_version = capture_workflow_version.workflow_version,
                version_capture_file = capture_task_version.version_capture_file,
                sample_name = sample_name, 
                bucket_path = out_bucket_path,

                # preprocess and preprocesss qc metrics files
                fastqc1_html_raw = fastqc_raw.fastqc1_html,
                fastqc1_zip_raw = fastqc_raw.fastqc1_zip,
                fastqc2_html_raw = fastqc_raw.fastqc2_html,
                fastqc2_zip_raw = fastqc_raw.fastqc2_zip,

                seqyclean_summary = seqyclean.seqyclean_summary,

                fastqc1_html_cleaned = fastqc_cleaned.fastqc1_html,
                fastqc1_zip_cleaned = fastqc_cleaned.fastqc1_zip,
                fastqc2_html_cleaned = fastqc_cleaned.fastqc2_html,
                fastqc2_zip_cleaned = fastqc_cleaned.fastqc2_zip,

                preprocess_qc_metrics = concat_preprocess_qc_metrics.preprocess_qc_metrics,

                # irma
                irma_assembled_gene_segments_csv = irma.irma_assembled_gene_segments_csv,
                irma_multifasta = irma.irma_multifasta,
                irma_fasta_array = irma_fasta_array,
                irma_bam_array = irma_bam_array,
                irma_vcf_array = irma_vcf_array,

                # ivar
                ivar_fasta_array = ivar_fasta_array,
                ivar_multifasta = make_ivar_multifasta.multifasta,
                ivar_parameters = ivar_parameters_file,


                # from samtools - sorted bams
                sorted_bam_array = sorted_bam_array,
                sorted_bai_array = sorted_bai_array,
                assembly_qc_metrics = concat_assembly_qc_metrics.assembly_qc_metrics_summary,

                # nextclade
                nextclade_json_array  = nextclade_json_array,
                nextclade_tsv_array  = nextclade_tsv_array,
                nextclade_HA_translation_fasta  = nextclade_HA_translation_fasta,
                nextclade_HA1_translation_fasta  = nextclade_HA1_translation_fasta ,
                nextclade_HA2_translation_fasta  = nextclade_HA2_translation_fasta ,
                nextclade_SigPep_translation_fasta  = nextclade_SigPep_translation_fasta,
                nextclade_NA_translation_fasta = nextclade_NA_translation_fasta
        }
    }

    output {
        # output from preprocess
        File fastqc1_html_raw = fastqc_raw.fastqc1_html
        File fastqc1_zip_raw = fastqc_raw.fastqc1_zip
        File fastqc2_html_raw = fastqc_raw.fastqc2_html
        File fastqc2_zip_raw = fastqc_raw.fastqc2_zip

        File seqyclean_summary = seqyclean.seqyclean_summary

        File fastqc1_html_cleaned = fastqc_cleaned.fastqc1_html
        File fastqc1_zip_cleaned = fastqc_cleaned.fastqc1_zip
        File fastqc2_html_cleaned = fastqc_cleaned.fastqc2_html
        File fastqc2_zip_cleaned = fastqc_cleaned.fastqc2_zip

        File preprocess_qc_metrics = concat_preprocess_qc_metrics.preprocess_qc_metrics

        # output from irma
        File? irma_assembled_gene_segments_csv = irma.irma_assembled_gene_segments_csv
        File? irma_multifasta = irma.irma_multifasta
        Array[File]? irma_fasta_array_out = irma_fasta_array
        Array[File]? irma_bam_array_out = irma_bam_array
        Array[File]? irma_vcf_array_out = irma_vcf_array


        # output from irma_subtyping_results
        File irma_typing = irma_subtyping_results.irma_typing
        String irma_type = irma_subtyping_results.irma_type
        String irma_ha_subtype = irma_subtyping_results.irma_ha_subtype
        String irma_na_subtype = irma_subtyping_results.irma_na_subtype

        # output from post assembly
        Array[File?]? ivar_fasta_array_out = ivar_fasta_array
        File? ivar_parameters = ivar_parameters_file
        File? ivar_multifasta = make_ivar_multifasta.multifasta
        Array[File?]? percent_coverage_csv_array_out = percent_coverage_csv_array
        Array[File?]? sorted_bam_array_out = sorted_bam_array
        Array[File?]? sorted_bai_array_out = sorted_bai_array
        File? assembly_qc_metrics = concat_assembly_qc_metrics.assembly_qc_metrics_summary


        # output from nextclade
        Array[File?]? nextclade_json = nextclade_json_array
        Array[File?]? nextclade_tsv = nextclade_tsv_array
        Array[File?]? nextclade_HA1_translation_fasta1 = nextclade_HA1_translation_fasta 
        Array[File?]? nextclade_HA2_translation_fasta1 = nextclade_HA2_translation_fasta
        Array[File?]? nextclade_SigPep_translation_fasta1 = nextclade_SigPep_translation_fasta
        Array[File?]? nextclade_HA_translation_fasta1 = nextclade_HA_translation_fasta
        Array[File?]? nextclade_NA_translation_fasta1 = nextclade_NA_translation_fasta
        
        # version capture
        File? version_capture_file = capture_task_version.version_capture_file

        # output from transfer
        String? transfer_date=transfer_assembly_wdl.transfer_date
    }
}

task exit_wdl {
    input {
        String exit_reason
    }
    
    command <<<
        echo "~{exit_reason}" 1>&2
        
    >>>

    runtime {
        return_codes: 0
        docker: "ubuntu:latest"
    }
}
