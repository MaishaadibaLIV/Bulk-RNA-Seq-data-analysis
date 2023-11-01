##Script for gene experssion count from RNAseq data
#!/bin/bash





####################################################################################
####################################################################################
#                                                                                  #
#                                                                                  #
#                     initial code for ref index file generation                   #
#           Run these lines manually once if analysing for the first time          #
#            This will create a ref index from ref genome and annotation           #
#                                                                                  #
#                                                                                  #
####################################################################################
####################################################################################

#Download ref files and build genome index
#here unmasked ref file was used of GRCh38.108
#input file name at 34,35,82,83


#Download gtf file
#wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz
#Download ref fasta
#wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#unzip both files
#gunzip *.gz

#Build genome index

#STAR --runMode genomeGenerate --genomeDir RNAseq_gene_ref/ --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.108.gtf --runThreadN 16

#the index files will be created, here, genomeDir is RNAseq_gene_ref2

####################################################################################
####################################################################################
#                                                                                  #
#                                                                                  #
#                                                                                  #
#                                 Main loop                                        #
#                                                                                  #
#                                                                                  #
#                                                                                  #
####################################################################################
####################################################################################


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


## mapping
## from https://www.reneshbedre.com/blog/star-aligner.html#mapping-reads-to-genome

#ref index = RNAseq_gene_ref2

STAR --runMode alignReads --runThreadN 12 --readFilesIn $file1 $file2 --genomeDir RNAseq_gene_ref2 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix file_Coorsorted_aligned  --outSAMunmapped Within

#cell statistics after conversion

samtools flagstat file_Coorsorted_alignedAligned.sortedByCoord.out.bam > ${sample_name}_Coorsorted_aligned_flagstat.txt

#cell sorting bam fie for names

samtools sort -n file_Coorsorted_alignedAligned.sortedByCoord.out.bam -o file_Coorsorted_aligned_name_sorted.bam

#make count tables

python -m HTSeq.scripts.count -f bam -r name file_Coorsorted_aligned_name_sorted.bam -s no -t exon "/media/ag72/Local Data 2/ERVmap workfile/Ref/RNAseq_gene_ref/Homo_sapiens.GRCh38.108.gtf" > ${sample_name}_htseq.cnt 2> ${sample_name}_htseq_log


#index coor sorted bam file

samtools index file_Coorsorted_alignedAligned.sortedByCoord.out.bam
#cell generate bigwig from coordinate sorted bam file

bamCoverage -b file_Coorsorted_alignedAligned.sortedByCoord.out.bam -o ${sample_name}_Coorsorted_alignedAligned.sortedByCoord.out.bw


  echo "Finished processing $file1 and $file2"
done < sample_ids.txt
