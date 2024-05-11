version 1.0

# import tasks
import "../tasks/preprocess_tasks.wdl" as fastq_preprocess
import "../tasks/irma_task.wdl" as irma_task
import "../tasks/ivar_task.wdl" as ivar
import "../tasks/post_assembly_tasks.wdl" as post_assembly_qc
import "../tasks/transfer_tasks.wdl" as transfer
import "../tasks/flu_nextclade_tasks.wdl" as nextclade

# begin workflow
workflow influenza_assembly {

    input {
        String sample_name
        File fastq_R1
        File fastq_R2
        File adapters_and_contaminants
        String bucket_path

        # python scripts
        File concat_preprocess_qc_metrics_py
        File irma_subtyping_results_py
        File calc_percent_coverage_py
        File concat_post_assembly_qc_metrics_py
    }

    # 1 - Preprocess QC raw fastq files
    call fastq_preprocess.fastqc as fastqc_raw {
        input:
            sample_name = sample_name,
            fastq_R1 = fastq_R1,
            fastq_R2 = fastq_R2
    }
    call fastq_preprocess.seqyclean as seqyclean {
        input:
            sample_name = sample_name,
            fastq_R1 = fastq_R1,
            fastq_R2 = fastq_R2,
            adapters_and_contaminants = adapters_and_contaminants  
    }
    call fastq_preprocess.fastqc as fastqc_cleaned {
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

            fastqc_version = fastqc_raw.fastqc_version,
            fastqc_docker = fastqc_raw.fastqc_docker,

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

            seqyclean_version = seqyclean.seqyclean_version,
            seqyclean_docker = seqyclean.seqyclean_docker


    }

    # 2- run irma
    call irma_task.irma as irma {
        input:
            sample_name = sample_name,
            read_type = read_type,
            fastq_R1 = seqyclean.fastq_R1_cleaned,
            fastq_R2 = seqyclean.fastq_R2_cleaned
    }

    call irma_task.irma_subtyping_results as irma_subtyping_results {
        input:
            irma_assembled_gene_segments_csv = irma.irma_assembled_gene_segments_csv,
            sample_name = sample_name,
            irma_runtime_csv = irma.irma_runtime_csv,
            python_script = irma_subtyping_results_py

    }

    # for each successfully assembled gene segment -
    # 1- samtools
    # 2 - ivar concensus
    # 3 - calcualte percent_coverage
    # 4 - if HA or NA run nextclade


    ####### 1 - HA ########
    if (defined(irma.irma_seg_ha_bam)) {
        call post_assembly_qc.samtools_mapped_reads as ha_mapped_reads {
            input:
                bam_file = irma.irma_seg_ha_bam,
                sample_name = sample_name
        }

        call ivar.ivar_consensus as ha_ivar_consensus {
            input:
                bam_file = irma_seg_ha_bam,
                sample_name = sample_name
        }

        call post_assembly_qc.calc_percent_coverage as ha_calc_percent_coverage{
            input:
                fasta_file = ha_ivar_consensus.ivar_consensus_fasta,
                python_script = calc_percent_coverage_py,
                sample_name = sample_name

        }

        call nextclade.nextclade_ha as nextclade_ha{
            input:
                ivar_seg_ha_fasta = ha_ivar_consensus.ivar_consensus_fasta,
                irma_type = irma_subtyping_results.irma_type,
                irma_ha_subtype = irma_subtyping_results.irma_ha_subtype,
                sample_name = sample_name
        }
    }

    ####### 2 - NA ########
    if (defined(irma.irma_seg_na_bam)) {
        call post_assembly_qc.samtools_mapped_reads as na_mapped_reads {
            input:
                bam_file = irma.irma_seg_na_bam,
                sample_name = sample_name
        }

        call ivar.ivar_consensus as na_ivar_consensus {
            input:
                bam_file = irma_seg_na_bam,
                sample_name = sample_name
        }

        call post_assembly_qc.calc_percent_coverage as na_calc_percent_coverage{
            input:
                fasta_file = na_ivar_consensus.ivar_consensus_fasta,
                python_script = calc_percent_coverage_py,
                sample_name = sample_name

        }

        call nextclade.nextclade_na as nextclade_na{
            input:
                ivar_seg_na_fasta = na_ivar_consensus.ivar_consensus_fasta,
                irma_type = irma_subtyping_results.irma_type,
                irma_ha_subtype = irma_subtyping_results.irma_ha_subtype,
                sample_name = sample_name
        }
    }    

    ####### 3 - PB1 ########
    if (defined(irma.irma_seg_pb1_bam)) {
        call post_assembly_qc.samtools_mapped_reads as pb1_mapped_reads {
            input:
                bam_file = irma.irma_seg_pb1_bam,
                sample_name = sample_name
        }

        call ivar.ivar_consensus as pb1_ivar_consensus {
            input:
                bam_file = irma_seg_pb1_bam,
                sample_name = sample_name
        }

        call post_assembly_qc.calc_percent_coverage as pb1_calc_percent_coverage{
            input:
                fasta_file = pb1_ivar_consensus.ivar_consensus_fasta,
                python_script = calc_percent_coverage_py,
                sample_name = sample_name

        }

    }  
    ####### 4 - PB2 ########
    if (defined(irma.irma_seg_pb2_bam)) {
        call post_assembly_qc.samtools_mapped_reads as pb2_mapped_reads {
            input:
                bam_file = irma.irma_seg_pb2_bam,
                sample_name = sample_name
        }

        call ivar.ivar_consensus as pb12_ivar_consensus {
            input:
                bam_file = irma_seg_pb2_bam,
                sample_name = sample_name
        }

        call post_assembly_qc.calc_percent_coverage as pb2_calc_percent_coverage{
            input:
                fasta_file = pb2_ivar_consensus.ivar_consensus_fasta,
                python_script = calc_percent_coverage_py,
                sample_name = sample_name

        }

    }  

    ####### 5 - NP ########
    if (defined(irma.irma_seg_np_bam)) {
        call post_assembly_qc.samtools_mapped_reads as np_mapped_reads {
            input:
                bam_file = irma.irma_seg_np_bam,
                sample_name = sample_name
        }

        call ivar.ivar_consensus as np_ivar_consensus {
            input:
                bam_file = irma_seg_np_bam,
                sample_name = sample_name
        }

        call post_assembly_qc.calc_percent_coverage as np_calc_percent_coverage{
            input:
                fasta_file = np_ivar_consensus.ivar_consensus_fasta,
                python_script = calc_percent_coverage_py,
                sample_name = sample_name

        }

    }  

    ####### 6 - PA ########
    if (defined(irma.irma_seg_pa_bam)) {
        call post_assembly_qc.samtools_mapped_reads as pa_mapped_reads {
            input:
                bam_file = irma.irma_seg_pa_bam,
                sample_name = sample_name
        }

        call ivar.ivar_consensus as pa_ivar_consensus {
            input:
                bam_file = irma_seg_pa_bam,
                sample_name = sample_name
        }

        call post_assembly_qc.calc_percent_coverage as pa_calc_percent_coverage{
            input:
                fasta_file = pa_ivar_consensus.ivar_consensus_fasta,
                python_script = calc_percent_coverage_py,
                sample_name = sample_name

        }

    } 
    ####### 7 - NS ########
    if (defined(irma.irma_seg_ns_bam)) {
        call post_assembly_qc.samtools_mapped_reads as nsmapped_reads {
            input:
                bam_file = irma.irma_seg_ns_bam,
                sample_name = sample_name
        }

        call ivar.ivar_consensus as ns_ivar_consensus {
            input:
                bam_file = irma_seg_ns_bam,
                sample_name = sample_name
        }

        call post_assembly_qc.calc_percent_coverage as ns_calc_percent_coverage{
            input:
                fasta_file = ns_ivar_consensus.ivar_consensus_fasta,
                python_script = calc_percent_coverage_py,
                sample_name = sample_name

        }

    } 

    ####### 8 - MP ########
    if (defined(irma.irma_seg_mp_bam)) {
        call post_assembly_qc.samtools_mapped_reads as mp_mapped_reads {
            input:
                bam_file = irma.irma_seg_mp_bam,
                sample_name = sample_name
        }

        call ivar.ivar_consensus as mp_ivar_consensus {
            input:
                bam_file = irma_seg_mp_bam,
                sample_name = sample_name
        }

        call post_assembly_qc.calc_percent_coverage as mp_calc_percent_coverage{
            input:
                fasta_file = mp_ivar_consensus.ivar_consensus_fasta,
                python_script = calc_percent_coverage_py,
                sample_name = sample_name

        }

    } 


        # concantenate post assembly qc metrics (coverage, depth) into a single file
        call post_assembly_qc.concat_post_qc_metrics as irma_concat_post_qc_metrics{
            input:
                python_script = concat_post_assembly_qc_metrics_py,
                sample_name = sample_name,

                # percent coverage results
                seg_ha_percent_coverage_csv = ha_calc_percent_coverage.percent_coverage_csv,
                seg_na_percent_coverage_csv = na_calc_percent_coverage.percent_coverage_csv,
                seg_pb1_percent_coverage_csv = pb1_calc_percent_coverage.percent_coverage_csv,
                seg_pb2_percent_coverage_csv = pb2_calc_percent_coverage.percent_coverage_csv,
                seg_np_percent_coverage_csv = np_calc_percent_coverage.percent_coverage_csv,
                seg_pa_percent_coverage_csv = pa_calc_percent_coverage.percent_coverage_csv,
                seg_ns_percent_coverage_csv = ns_calc_percent_coverage.percent_coverage_csv,
                seg_mp_percent_coverage_csv = mp_calc_percent_coverage.percent_coverage_csv,

                # mapped_reads_csv
                File? seg_ha_mapped_reads_csv = ha_mapped_reads.mapped_reads_csv,
                File? seg_na_mapped_reads_csv = na_mapped_reads.mapped_reads_csv,
                File? seg_pb1_mapped_reads_csv = pb1_mapped_reads.mapped_reads_csv,
                File? seg_pb2_mapped_reads_csv = pb2_mapped_reads.mapped_reads_csv,
                File? seg_np_mapped_reads_csv = np_mapped_reads.mapped_reads_csv,
                File? seg_pa_mapped_reads_csv = pa_mapped_reads.mapped_reads_csv,
                File? seg_ns_mapped_reads_csv = ns_mapped_reads.mapped_reads_csv,
                File? seg_mp_mapped_reads_csv = mp_mapped_reads.mapped_reads_csv


        }

       
    
    }
    # 5 - Transfer some intermediate files and all final files to gcp bucket
    call transfer.transfer_assembly_wdl as transfer_assembly_wdl {
        input:
            sample_name = sample_name, 
            bucket_path = bucket_path,

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

            # irma - assemblies
            irma_seg_ha_fasta = irma.irma_seg_ha_fasta,
            irma_seg_na_fasta = irma.irma_seg_na_fasta,
            irma_seg_pb1_fasta = irma.irma_seg_pb1_fasta,
            irma_seg_pb2_fasta = irma.irma_seg_pb2_fasta,
            irma_seg_np_fasta =  irma.irma_seg_np_fasta,
            irma_seg_pa_fasta = irma.irma_seg_pa_fasta,
            irma_seg_ns_fasta = irma.irma_seg_ns_fasta,
            irma_seg_mp_fasta = irma.irma_seg_mp_fasta,
            irma_all_assembled_segments_fasta = irma.irma_all_assembled_segments_fasta

            # irma - alignments
            irma_seg_ha_bam = irma.irma_seg_ha_bam,
            irma_seg_na_bam = irma.irma_seg_na_bam,
            irma_seg_pb1_bam = irma.irma_seg_pb1_bam,
            irma_seg_pb2_bam = irma.irma_seg_pb2_bam,
            irma_seg_np_bam = irma.irma_seg_np_bam,
            irma_seg_pa_bam = irma.irma_seg_pa_bam,
            irma_seg_ns_bam = irma.irma_seg_ns_bam,
            irma_seg_mp_bam  = irma.irma_seg_mp_bam,

            # irma - vcf
            irma_seg_ha_vcf = irma.irma_seg_ha_vcf,
            irma_seg_na_vcf = irma.irma_seg_na_vcf,
            irma_seg_pb1_vcf = irma.irma_seg_pb1_vcf,
            irma_seg_pb2_vcf = irma.irma_seg_b2_vcf,
            irma_seg_np_vcf = irma.irma_seg_np_vcf,
            irma_seg_pa_vcf = irma.irma_seg_pa_vcf,
            irma_seg_ns_vcf = irma.irma_seg_ns_vcf,
            irma_seg_mp_vcf = irma.irma_seg_mp_vcf,

            # ivar - assemblies
            ivar_seg_ha_fasta = ha_ivar_consensus.ivar_consensus_fasta,
            ivar_seg_na_fasta = na_ivar_consensus.ivar_consensus_fasta,
            ivar_seg_pb1_fasta = pb1_ivar_consensus.ivar_consensus_fasta,
            ivar_seg_pb2_fasta = pb2_ivar_consensus.ivar_consensus_fasta,
            ivar_seg_np_fasta = np_ivar_consensus.ivar_consensus_fasta,
            ivar_seg_pa_fasta = pa_ivar_consensus.ivar_consensus_fasta,
            ivar_seg_ns_fasta = ns_ivar_consensus.ivar_consensus_fasta,
            ivar_seg_mp_fasta = mp_ivar_consensus.ivar_consensus_fasta,

            # sorted bams from samtools
            irma_seg_ha_bam_sorted = ha_mapped_reads.sorted_bam,
            irma_seg_na_bam_sorted = na_mapped_reads.sorted_bam,
            irma_seg_pb1_bam_sorted = pb1_mapped_reads.sorted_bam,
            irma_seg_pb2_bam_sorted = pb2_mapped_reads.sorted_bam,
            irma_seg_np_bam_sorted = np_mapped_reads.sorted_bam,
            irma_seg_pa_bam_sorted = pa_mapped_reads.sorted_bam,
            irma_seg_ns_bam_sorted = ns_mapped_reads.sorted_bam,
            irma_seg_mp_bam_sorted = mp_mapped_reads.sorted_bam,

        


            irma_qc_metrics = irma_concat_post_qc_metrics.qc_metrics_summary,

            # nextclade
            na_nextclade_json = nextclade_na.na_nextclade_json,
            na_nextclade_tsv = nextclade_na.na_nextclade_tsv,
            na_translation_fasta = nextclade_na.na_translation_fasta,

            ha_nextclade_json = nextclade_ha.ha_nextclade_json,
            ha_nextclade_tsv = nextclade_ha.ha_nextclade_tsv,
            ha_HA1_translation_fasta = nextclade_ha.ha_HA1_translation_fasta,
            ha_HA2_translation_fasta = nextclade_ha.ha_HA2_translation_fasta,
            ha_SigPep_translation_fasta = nextclade_ha.ha_SigPep_translation_fasta
    }


    output {
        # output from preprocess
        String fastqc_version = fastqc_raw.fastqc_version
        String fastqc_docker = fastqc_raw.fastqc_docker

        File fastqc1_html_raw = fastqc_raw.fastqc1_html
        File fastqc1_zip_raw = fastqc_raw.fastqc1_zip
        File? fastqc2_html_raw = fastqc_raw.fastqc2_html
        File? fastqc2_zip_raw = fastqc_raw.fastqc2_zip

        String seqyclean_version = seqyclean.seqyclean_version
        String seqyclean_docker = seqyclean.seqyclean_docker
        File seqyclean_summary = seqyclean.seqyclean_summary

        File fastqc1_html_cleaned = fastqc_cleaned.fastqc1_html
        File fastqc1_zip_cleaned = fastqc_cleaned.fastqc1_zip
        File? fastqc2_html_cleaned = fastqc_cleaned.fastqc2_html
        File? fastqc2_zip_cleaned = fastqc_cleaned.fastqc2_zip

        File preprocess_qc_metrics = concat_preprocess_qc_metrics.preprocess_qc_metrics

        # output from irma
        File irma_assembled_gene_segments_csv = irma.irma_assembled_gene_segments_csv
        File? irma_all_assembled_segments_fasta = irma.irma_all_assembled_segments_fasta
        
        # assemblies
        File? irma_seg_ha_fasta = irma.irma_seg_ha_fasta
        File? irma_seg_na_fasta = irma.irma_seg_na_fasta
        File? irma_seg_pb1_fasta = irma.irma_seg_pb1_fasta
        File? irma_seg_pb2_fasta = irma.irma_seg_pb2_fasta
        File? irma_seg_np_fasta = irma.irma_seg_np_fasta
        File? irma_seg_pa_fasta = irma.irma_seg_pa_fasta
        File? irma_seg_ns_fasta = irma.irma_seg_ns_fasta
        File? irma_seg_mp_fasta = irma.irma_seg_mp_fasta

        # alignments
        File? irma_seg_ha_bam = irma.irma_seg_ha_bam
        File? irma_seg_na_bam = irma.irma_seg_na_bam
        File? irma_seg_pb1_bam = irma.irma_seg_pb1_bam
        File? irma_seg_pb2_bam = irma_.irma_seg_pb2_bam
        File? irma_seg_np_bam = irma_.irma_seg_np_bam
        File? irma_seg_pa_bam = irma.irma_seg_pa_bam
        File? irma_seg_ns_bam = irma.irma_seg_ns_bam
        File? irma_seg_mp_bam = irma.irma_seg_mp_bam

        # vcfs
        File? irma_seg_ha_vcf = irma.irma_seg_ha_vcf
        File? irma_seg_na_vcf = irma.irma_seg_na_vcf
        File? irma_seg_pb1_vcf = irma.irma_seg_pb1_vcf
        File? irma_seg_pb2_vcf = irma.irma_seg_pb2_vcf
        File? irma_seg_np_vcf = irma.irma_seg_np_vcf
        File? irma_seg_pa_vcf = irma.irma_seg_pa_vcf
        File? irma_seg_ns_vcf = irma.irma_seg_ns_vcf
        File? irma_seg_mp_vcf = irma.irma_seg_mp_vcf
        
        String irma_version = irma.irma_version
        String irma_docker = irma.irma_docker
        String irma_module = irma.irma_module

        # output from irma_subtyping_results
        File irma_typing = irma_subtyping_results.irma_typing
        String irma_type = irma_subtyping_results.irma_type
        String irma_ha_subtype = irma_subtyping_results.irma_ha_subtype
        String irma_na_subtype = irma_subtyping_results.irma_na_subtype

        # output from ivar_consensus
        File? ivar_seg_ha_fasta = ha_ivar_consensus.ivar_consensus_fasta
        File? ivar_seg_na_fasta = na_ivar_consensus.ivar_consensus_fasta
        File? ivar_seg_pb1_fasta = pb1_ivar_consensus.ivar_consensus_fasta
        File? ivar_seg_pb2_fasta = pb2_ivar_consensus.ivar_consensus_fasta
        File? ivar_seg_np_fasta = np_ivar_consensus.ivar_consensus_fasta
        File? ivar_seg_pa_fasta = pa_ivar_consensus.ivar_consensus_fasta
        File? ivar_seg_ns_fasta = ns_ivar_consensus.ivar_consensus_fasta
        File? ivar_seg_mp_fasta = mp_ivar_consensus.ivar_consensus_fasta


        # output from calculate percent coverage
        File? seg_ha_percent_coverage_csv = ha_calc_percent_coverage.percent_coverage_csv
        File? seg_na_percent_coverage_csv = na_calc_percent_coverage.percent_coverage_csv
        File? seg_pb1_percent_coverage_csv = pb1_calc_percent_coverage.percent_coverage_csv
        File? seg_pb2_percent_coverage_csv = pb2_calc_percent_coverage.percent_coverage_csv
        File? seg_np_percent_coverage_csv = np_calc_percent_coverage.percent_coverage_csv
        File? seg_pa_percent_coverage_csv = pa_calc_percent_coverage.percent_coverage_csv
        File? seg_ns_percent_coverage_csv = ns_calc_percent_coverage.percent_coverage_csv
        File? seg_mp_percent_coverage_csv = mp_calc_percent_coverage.percent_coverage_csv


        # sorted bams from samtools 
        File? irma_seg_ha_bam_sorted = ha_mapped_reads.sorted_bam
        File? irma_seg_na_bam_sorted = na_mapped_reads.sorted_bam
        File? irma_seg_pb1_bam_sorted = pb1_mapped_reads.sorted_bam
        File? irma_seg_pb2_bam_sorted = pb2_mapped_reads.sorted_bam
        File? irma_seg_np_bam_sorted = np_mapped_reads.sorted_bam
        File? irma_seg_pa_bam_sorted = pa_mapped_reads.sorted_bam
        File? irma_seg_ns_bam_sorted = ns_mapped_reads.sorted_bam
        File? irma_seg_mp_bam_sorted = mp_mapped_reads.sorted_bam
        
        # mapped reads from samtools

       

        # output from post assembly QC metrics
        Array[File]? irma_bam_results = irma_samtools_mapped_reads.bam_results
        Array[File]? irma_per_cov_results = irma_percent_coverage.perc_cov_results
        File? irma_assembly_qc_metrics = irma_concat_post_qc_metrics.qc_metrics_summary

        # output from nextclade
        File? na_nextclade_json = nextclade_na.na_nextclade_json
        File? na_nextclade_tsv = nextclade_na.na_nextclade_tsv
        File? na_translation_fasta = nextclade_na.na_translation_fasta

        File? ha_nextclade_json = nextclade_ha.ha_nextclade_json
        File? ha_nextclade_tsv = nextclade_ha.ha_nextclade_tsv
        File? ha_HA1_translation_fasta = nextclade_ha.ha_HA1_translation_fasta
        File? ha_HA2_translation_fasta = nextclade_ha.ha_HA2_translation_fasta
        File? ha_SigPep_translation_fasta = nextclade_ha.ha_SigPep_translation_fasta
        
        # output from transfer
        String transfer_date=transfer_assembly_wdl.transfer_date
    }
}
