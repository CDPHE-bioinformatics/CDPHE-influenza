 #1/bin/bash

bam_file="/home/molly_hetheringtonrauth/sandbox_influenza/pipeline_development/test_python/input/2103754030_A_HA_H1.bam"

 # get segment name from bam file name
base_name=$(basename ${bam_file})
sample_id=$(basename ${bam_file} | cut -d "_" -f 1)
segment_name=$(basename ${bam_file} | cut -d "." -f 1 | cut -d "_" -f 2-)
gene_name=$(basename ${bam_file} | cut -d "_" -f 3)

# count mapped flu reads
samtools view -c -F 260 ${bam_file} > num_mapped_reads.txt
samtools coverage ${bam_file} | tail -1 | cut -f 7 > mean_depth.txt

echo "sample_id,base_name,segment_name,gene_name,description,value" > ${segment_name}_bam_results.csv
echo "${sample_id},${base_name},${segment_name},${gene_name},num_mapped_reads,$(cat num_mapped_reads.txt)" >> ${segment_name}_bam_results.csv
echo "${sample_id},${base_name},${segment_name},${gene_name},mean_depth,$(cat mean_depth.txt)" >> ${segment_name}_bam_results.csv