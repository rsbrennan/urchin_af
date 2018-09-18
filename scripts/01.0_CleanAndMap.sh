### ead cleaning (Trimmomatic) & mapping (bwamem)

#################################################################################################
### BASH SCRIPT FOR CLEANING READS VIA TRIMMOMATIC-0.36

#!/bin/bash

for f1 in *_R1.fastq.gz

do

 f2=${f1%%_R1.fastq.gz}"_R2.fastq.gz"



java -jar /data/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
        -threads 10 \
        -phred33 \
         "$f1" \
         "$f2" \
         /data/users/a/m/amakukho/RESEARCH/OASV2/data/clippedtrimmed/"$f1"_left_clean_paired.fq \
         /data/users/a/m/amakukho/RESEARCH/OASV2/data/clippedtrimmed/"$f1"_left_clean_unpaired.fq \
         /data/users/a/m/amakukho/RESEARCH/OASV2/data/clippedtrimmed/"$f2"_right_clean_paired.fq \
         /data/users/a/m/amakukho/RESEARCH/OASV2/data/clippedtrimmed/"$f2"_right_clean_unpaired.fq \
        ILLUMINACLIP:/users/a/m/amakukho/RESEARCH/OASV2/scripts/OASV2CaptureSeqAdapters.fa:2:30:10 \
        SLIDINGWINDOW:20:2 \
        LEADING:2 \
        TRAILING:2 \
        MINLEN:35

done

#################################################################################################
### MAPPING CLEANED READS TO S. PURPURATUS GENOME v 3.0 (build 7) VIA BWAMEM

# Example for sample OASV2_DNA_D7_8_0_S_15
/data/programs/bwa/bwa mem -t 10 -k 5 -R '@RG\tID:OASV2_DNA_D7_8_0_S_15\tSM:OASV2_DNA_D7_8_0_S_15\tPL:Illumina' ~/RESEARCH/OASV2/data/bwa_mapped/ref/ref ~/RESEARCH/OASV2/data/clippedtrimmed-2/OASV2_DNA_D7_8_0_S_15_R1_clean_paired.fa ~/RESEARCH/OASV2/data/clippedtrimmed-2/OASV2_DNA_D7_8_0_S_15_R2_clean_paired.fa > OASV2_DNA_D7_8_0_S_15_bwamem.sam

#################################################################################################
### MERGE, SORT, INDEX

#!/bin/bash

cd ~/RESEARCH/OASV2/data/bwa_mem

# Merge all the bam files into one
/data/programs/samtools merge bwamem_final.merged.bam *.bam

# NOTE: Need to specify a different tmp directory or it will fill the rhel-root...
/data/programs/sambamba_v0.6.0 sort -m 8G -t 8 -p --tmpdir=/users/a/m/amakukho/RESEARCH/OASV2/data/bwa_mem/tmp bwamem_final.merged.bam bwamem_final.merged.sorted.bam

# Run samtools index
/data/programs/samtools index bwamem_final.merged.sorted.bam






