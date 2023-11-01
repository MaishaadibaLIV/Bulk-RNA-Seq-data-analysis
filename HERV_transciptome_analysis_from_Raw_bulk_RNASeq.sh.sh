#!/bin/bash

#create a file named sample_ids.txt containing the headings of all sample in directory
# Read the sample IDs from sample_ids.txt and process each pair of R1 and R2 files
while IFS= read -r sample_id; do
    # Define R1 and R2 file names
    file1="${sample_id}_R1_UMI_removed.fastq"				# if trimming performed beforehand, replace the extensions as _R1_val_1.fq and _R2_val_2.fq
    file2="${sample_id}_R2.fastq"

    # Extract the sample name and read number
    sample_name="${sample_id%_R*}"
    read_number="${sample_id##*_}"

    echo "Processing $file1 and $file2"

#Source data - trimmed fastq file both, ref folder with bowtie2 index, HERV_rmsk.hg38.v2.2.gtf

#Alignment

    conda activate singleloc_align
    bowtie2 -p 12 -k 100 --very-sensitive-local --score-min "L,0,1.6" --rg-id sample_01 -x refs/hg38 -q -1 $file1 -2 $file2  -S ${sample_name}_multi.sam  2>&1 | tee ${sample_name}_multi.log
    conda deactivate

    samtools flagstat ${sample_name}_multi.sam > ${sample_name}_alignment_flagstat_report.txt

#Telescope run


  conda activate telescope_env
  telescope assign --theta_prior 200000 --max_iter 200 --updated_sam ${sample_name}_multi.sam '/media/ag72/Local Data 2/Telescope/refs/HERV_rmsk.hg38.v2.2.gtf' 2>&1 | tee ${sample_name}_telescope.log
  conda deactivate

#Rename telescope report files

  mv telescope-telescope_report.tsv ${sample_name}_telescope-telescope_report.tsv


  echo "Finished processing $file1 and $file2"
done < sample_ids.txt
