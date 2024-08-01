version 1.0

# import tasks
import "../tasks/preprocess_tasks.wdl" as fastq_preprocess
import "../tasks/irma_task.wdl" as irma_task
import "../tasks/ivar_task.wdl" as ivar
import "../tasks/assembly_qc_tasks.wdl" as assembly_qc
import "../tasks/transfer_tasks.wdl" as transfer
import "../tasks/nextclade_tasks.wdl" as nextclade
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
            read_pairs_cleaned = fastqc_cleaned.read_pairs,



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


    if (irma_subtyping_results.irma_type != 'no IRMA assembly generated') {
        # for each successfully assembled gene segment run-
        # 1- samtools, 2 - ivar concensus, 3 - calcualte percent_coverage, 4 - if HA or NA run nextclade

        ####### 1 - HA ########
        if (defined(irma.irma_seg_ha_bam)) {
            call assembly_qc.calc_bam_stats_samtools as ha_bam_stats {
                input:
                    bam_file = irma.irma_seg_ha_bam,
                    sample_name = sample_name
            }

            call ivar.call_consensus_ivar as ha_ivar_consensus {
                input:
                    bam_file = irma.irma_seg_ha_bam,
                    sample_name = sample_name,
                    irma_type = irma_subtyping_results.irma_type,
                    irma_na_subtype = irma_subtyping_results.irma_na_subtype,
                    irma_ha_subtype = irma_subtyping_results.irma_ha_subtype
            }

            call assembly_qc.calc_percent_coverage as ha_calc_percent_coverage{
                input:
                    fasta_file = ha_ivar_consensus.ivar_consensus_fasta,
                    python_script = calc_percent_coverage_py,
                    sample_name = sample_name

            }

            call nextclade.ha_nextclade as ha_nextclade{
                input:
                    ivar_seg_ha_fasta = ha_ivar_consensus.ivar_consensus_fasta,
                    irma_type = irma_subtyping_results.irma_type,
                    irma_ha_subtype = irma_subtyping_results.irma_ha_subtype,
                    sample_name = sample_name
            }
        }

        ####### 2 - NA ########
        if (defined(irma.irma_seg_na_bam)) {
            call assembly_qc.calc_bam_stats_samtools as na_bam_stats {
                input:
                    bam_file = irma.irma_seg_na_bam,
                    sample_name = sample_name
            }

            call ivar.call_consensus_ivar as na_ivar_consensus {
                input:
                    bam_file = irma.irma_seg_na_bam,
                    sample_name = sample_name,
                    irma_type = irma_subtyping_results.irma_type,
                    irma_na_subtype = irma_subtyping_results.irma_na_subtype,
                    irma_ha_subtype = irma_subtyping_results.irma_ha_subtype
            }

            call assembly_qc.calc_percent_coverage as na_calc_percent_coverage{
                input:
                    fasta_file = na_ivar_consensus.ivar_consensus_fasta,
                    python_script = calc_percent_coverage_py,
                    sample_name = sample_name

            }

            call nextclade.na_nextclade as na_nextclade{
                input:
                    ivar_seg_na_fasta = na_ivar_consensus.ivar_consensus_fasta,
                    irma_type = irma_subtyping_results.irma_type,
                    irma_na_subtype = irma_subtyping_results.irma_na_subtype,
                    sample_name = sample_name
            }
        }    

        ####### 3 - PB1 ########
        if (defined(irma.irma_seg_pb1_bam)) {
            call assembly_qc.calc_bam_stats_samtools as pb1_bam_stats {
                input:
                    bam_file = irma.irma_seg_pb1_bam,
                    sample_name = sample_name
            }

            call ivar.call_consensus_ivar as pb1_ivar_consensus {
                input:
                    bam_file = irma.irma_seg_pb1_bam,
                    sample_name = sample_name,
                    irma_type = irma_subtyping_results.irma_type,
                    irma_na_subtype = irma_subtyping_results.irma_na_subtype,
                    irma_ha_subtype = irma_subtyping_results.irma_ha_subtype
            }

            call assembly_qc.calc_percent_coverage as pb1_calc_percent_coverage{
                input:
                    fasta_file = pb1_ivar_consensus.ivar_consensus_fasta,
                    python_script = calc_percent_coverage_py,
                    sample_name = sample_name

            }

        }  
        ####### 4 - PB2 ########
        if (defined(irma.irma_seg_pb2_bam)) {
            call assembly_qc.calc_bam_stats_samtools as pb2_bam_stats {
                input:
                    bam_file = irma.irma_seg_pb2_bam,
                    sample_name = sample_name
            }

            call ivar.call_consensus_ivar as pb2_ivar_consensus {
                input:
                    bam_file = irma.irma_seg_pb2_bam,
                    sample_name = sample_name,
                    irma_type = irma_subtyping_results.irma_type,
                    irma_na_subtype = irma_subtyping_results.irma_na_subtype,
                    irma_ha_subtype = irma_subtyping_results.irma_ha_subtype
            }

            call assembly_qc.calc_percent_coverage as pb2_calc_percent_coverage{
                input:
                    fasta_file = pb2_ivar_consensus.ivar_consensus_fasta,
                    python_script = calc_percent_coverage_py,
                    sample_name = sample_name

            }

        }  

        ####### 5 - NP ########
        if (defined(irma.irma_seg_np_bam)) {
            call assembly_qc.calc_bam_stats_samtools as np_bam_stats {
                input:
                    bam_file = irma.irma_seg_np_bam,
                    sample_name = sample_name
            }

            call ivar.call_consensus_ivar as np_ivar_consensus {
                input:
                    bam_file = irma.irma_seg_np_bam,
                    sample_name = sample_name,
                    irma_type = irma_subtyping_results.irma_type,
                    irma_na_subtype = irma_subtyping_results.irma_na_subtype,
                    irma_ha_subtype = irma_subtyping_results.irma_ha_subtype
            }

            call assembly_qc.calc_percent_coverage as np_calc_percent_coverage{
                input:
                    fasta_file = np_ivar_consensus.ivar_consensus_fasta,
                    python_script = calc_percent_coverage_py,
                    sample_name = sample_name

            }

        }  

        ####### 6 - PA ########
        if (defined(irma.irma_seg_pa_bam)) {
            call assembly_qc.calc_bam_stats_samtools as pa_bam_stats {
                input:
                    bam_file = irma.irma_seg_pa_bam,
                    sample_name = sample_name
            }

            call ivar.call_consensus_ivar as pa_ivar_consensus {
                input:
                    bam_file = irma.irma_seg_pa_bam,
                    sample_name = sample_name,
                    irma_type = irma_subtyping_results.irma_type,
                    irma_na_subtype = irma_subtyping_results.irma_na_subtype,
                    irma_ha_subtype = irma_subtyping_results.irma_ha_subtype
            }

            call assembly_qc.calc_percent_coverage as pa_calc_percent_coverage{
                input:
                    fasta_file = pa_ivar_consensus.ivar_consensus_fasta,
                    python_script = calc_percent_coverage_py,
                    sample_name = sample_name

            }

        } 
        ####### 7 - NS ########
        if (defined(irma.irma_seg_ns_bam)) {
            call assembly_qc.calc_bam_stats_samtools as ns_bam_stats {
                input:
                    bam_file = irma.irma_seg_ns_bam,
                    sample_name = sample_name
            }

            call ivar.call_consensus_ivar as ns_ivar_consensus {
                input:
                    bam_file = irma.irma_seg_ns_bam,
                    sample_name = sample_name,
                    irma_type = irma_subtyping_results.irma_type,
                    irma_na_subtype = irma_subtyping_results.irma_na_subtype,
                    irma_ha_subtype = irma_subtyping_results.irma_ha_subtype
            }

            call assembly_qc.calc_percent_coverage as ns_calc_percent_coverage{
                input:
                    fasta_file = ns_ivar_consensus.ivar_consensus_fasta,
                    python_script = calc_percent_coverage_py,
                    sample_name = sample_name

            }

        } 

        ####### 8 - MP ########
        if (defined(irma.irma_seg_mp_bam)) {
            call assembly_qc.calc_bam_stats_samtools as mp_bam_stats {
                input:
                    bam_file = irma.irma_seg_mp_bam,
                    sample_name = sample_name
            }

            call ivar.call_consensus_ivar as mp_ivar_consensus {
                input:
                    bam_file = irma.irma_seg_mp_bam,
                    sample_name = sample_name,
                    irma_type = irma_subtyping_results.irma_type,
                    irma_na_subtype = irma_subtyping_results.irma_na_subtype,
                    irma_ha_subtype = irma_subtyping_results.irma_ha_subtype
            }

            call assembly_qc.calc_percent_coverage as mp_calc_percent_coverage{
                input:
                    fasta_file = mp_ivar_consensus.ivar_consensus_fasta,
                    python_script = calc_percent_coverage_py,
                    sample_name = sample_name

            }

        } 
        # create arrays to better handle groups of files
        # IRMA - fasta, bam, vcf
        Array[File] irma_fasta_array = select_all([ irma.irma_seg_ha_fasta,
                                                    irma.irma_seg_na_fasta,
                                                    irma.irma_seg_pb1_fasta,
                                                    irma.irma_seg_pb2_fasta,
                                                    irma.irma_seg_np_fasta,
                                                    irma.irma_seg_pa_fasta,
                                                    irma.irma_seg_ns_fasta,
                                                    irma.irma_seg_mp_fasta])

        Array[File] irma_bam_array = select_all([irma.irma_seg_ha_bam,
                                                    irma.irma_seg_na_bam,
                                                    irma.irma_seg_pb1_bam,
                                                    irma.irma_seg_pb2_bam,
                                                    irma.irma_seg_np_bam,
                                                    irma.irma_seg_pa_bam,
                                                    irma.irma_seg_ns_bam,
                                                    irma.irma_seg_mp_bam])

        Array[File] irma_vcf_array = select_all([irma.irma_seg_ha_vcf,
                                                    irma.irma_seg_na_vcf,
                                                    irma.irma_seg_pb1_vcf,
                                                    irma.irma_seg_pb2_vcf,
                                                    irma.irma_seg_np_vcf,
                                                    irma.irma_seg_pa_vcf,
                                                    irma.irma_seg_ns_vcf,
                                                    irma.irma_seg_mp_vcf])

        # IVAR - fasta
        Array[File] ivar_fasta_array = select_all([ha_ivar_consensus.ivar_consensus_fasta,
                                                    na_ivar_consensus.ivar_consensus_fasta,
                                                    pb1_ivar_consensus.ivar_consensus_fasta,
                                                    pb2_ivar_consensus.ivar_consensus_fasta,
                                                    np_ivar_consensus.ivar_consensus_fasta,
                                                    pa_ivar_consensus.ivar_consensus_fasta,
                                                    ns_ivar_consensus.ivar_consensus_fasta,
                                                    mp_ivar_consensus.ivar_consensus_fasta])

        # Samtools - mapped reads csv, sorted bam
        Array[File] bam_stats_csv_array = select_all([ha_bam_stats.bam_stats_csv,
                                                            na_bam_stats.bam_stats_csv,
                                                            pb1_bam_stats.bam_stats_csv,
                                                            pb2_bam_stats.bam_stats_csv,
                                                            np_bam_stats.bam_stats_csv,
                                                            pa_bam_stats.bam_stats_csv,
                                                            ns_bam_stats.bam_stats_csv,
                                                            mp_bam_stats.bam_stats_csv])

        Array[File] sorted_bam_array = select_all([ha_bam_stats.sorted_bam,
                                                        na_bam_stats.sorted_bam,
                                                        pb1_bam_stats.sorted_bam,
                                                        pb2_bam_stats.sorted_bam,
                                                        np_bam_stats.sorted_bam,
                                                        pa_bam_stats.sorted_bam,
                                                        ns_bam_stats.sorted_bam,
                                                        mp_bam_stats.sorted_bam])

        # percent coverage - percent coverage csv
        Array[File] percent_coverage_csv_array = select_all([ha_calc_percent_coverage.percent_coverage_csv,
                                                                na_calc_percent_coverage.percent_coverage_csv,
                                                                pb1_calc_percent_coverage.percent_coverage_csv,
                                                                pb2_calc_percent_coverage.percent_coverage_csv,
                                                                np_calc_percent_coverage.percent_coverage_csv,
                                                                pa_calc_percent_coverage.percent_coverage_csv,
                                                                ns_calc_percent_coverage.percent_coverage_csv,
                                                                mp_calc_percent_coverage.percent_coverage_csv])

        
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
    }

    # 5 - Version capture
    call capture_version.capture_workflow_version  as capture_workflow_version{
        input:
    }
    
    # create array of structs
    Array[VersionInfo] version_array = select_all([
        fastqc_raw.fastqc_version_info,
        seqyclean.seqyclean_version_info,
        irma.IRMA_version_info,
        ha_bam_stats.samtools_version_info,
        ha_ivar_consensus.ivar_version_info,
        ha_ivar_consensus.samtools_version_info,
        ha_nextclade.ha_nextclade_version_info,
        na_bam_stats.samtools_version_info,
        na_ivar_consensus.ivar_version_info,
        na_ivar_consensus.samtools_version_info,
        na_nextclade.na_nextclade_version_info,
        pb1_bam_stats.samtools_version_info,
        pb1_ivar_consensus.ivar_version_info,
        pb1_ivar_consensus.samtools_version_info,
        pb2_bam_stats.samtools_version_info,
        pb2_ivar_consensus.ivar_version_info,
        pb2_ivar_consensus.samtools_version_info,
        np_bam_stats.samtools_version_info,
        np_ivar_consensus.ivar_version_info,
        np_ivar_consensus.samtools_version_info,
        pa_bam_stats.samtools_version_info,
        pa_ivar_consensus.ivar_version_info,
        pa_ivar_consensus.samtools_version_info,
        ns_bam_stats.samtools_version_info,
        ns_ivar_consensus.ivar_version_info,
        ns_ivar_consensus.samtools_version_info,
        mp_bam_stats.samtools_version_info,
        mp_ivar_consensus.ivar_version_info,
        mp_ivar_consensus.samtools_version_info

    ])

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

            # from samtoosls - sorted bams
            sorted_bam_array = sorted_bam_array,

            assembly_qc_metrics = concat_assembly_qc_metrics.assembly_qc_metrics_summary,

            # nextclade
            na_nextclade_json = na_nextclade.na_nextclade_json,
            na_nextclade_tsv = na_nextclade.na_nextclade_tsv,
            na_nextclade_translation_fasta = na_nextclade.na_nextclade_translation_fasta,

            ha_nextclade_json = ha_nextclade.ha_nextclade_json,
            ha_nextclade_tsv = ha_nextclade.ha_nextclade_tsv,
            ha_nextclade_HA1_translation_fasta = ha_nextclade.ha_nextclade_HA1_translation_fasta,
            ha_nextclade_HA2_translation_fasta = ha_nextclade.ha_nextclade_HA2_translation_fasta,
            ha_nextclade_SigPep_translation_fasta = ha_nextclade.ha_nextclade_SigPep_translation_fasta,


    
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
        File irma_assembled_gene_segments_csv = irma.irma_assembled_gene_segments_csv
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
        Array[File]? ivar_fasta_array_out = ivar_fasta_array
        File? ivar_multifasta = make_ivar_multifasta.multifasta
        Array[File]? percent_coverage_csv_array_out = percent_coverage_csv_array
        Array[File]? sorted_bam_array_out = sorted_bam_array
        File? assembly_qc_metrics = concat_assembly_qc_metrics.assembly_qc_metrics_summary


        # output from nextclade
        File? na_nextclade_json = na_nextclade.na_nextclade_json
        File? na_nextclade_tsv = na_nextclade.na_nextclade_tsv
        File? na_translation_fasta = na_nextclade.na_nextclade_translation_fasta

        File? ha_nextclade_json = ha_nextclade.ha_nextclade_json
        File? ha_nextclade_tsv = ha_nextclade.ha_nextclade_tsv
        File? ha_nextclade_HA1_translation_fasta = ha_nextclade.ha_nextclade_HA1_translation_fasta
        File? ha_nextclade_HA2_translation_fasta = ha_nextclade.ha_nextclade_HA2_translation_fasta
        File? ha_nextclade_SigPep_translation_fasta = ha_nextclade.ha_nextclade_SigPep_translation_fasta
        
        # version capture
        File version_capture_file = capture_task_version.version_capture_file

        # output from transfer
        String transfer_date=transfer_assembly_wdl.transfer_date
    }
}
