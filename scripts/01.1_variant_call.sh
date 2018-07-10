#!/bin/bash

# call variants w/ samtools. run in parallel, 10 jobs

# the probe bed file is just being used to break the variant calling into scaffolds that contain capture probes. This should make things run more quickly. Will merge variants after.

tail -n +2 /data/OASV2/RG_7601_Probe.extended.bed | cut -f 1 | sort | uniq| grep -v "#scaffold" | xargs -I {} -n 1 -P 10 sh -c "/data/programs/samtools-1.4.1/samtools mpileup -u --skip-indels -d 10000  -t DP,AD,INFO/AD -f /data/OASV2/Spur_3.1.LinearScaffold.fa -r {} /data/OASV2/bwamem/fullymapped/bwamem_final.merged.sorted.bam | /data/programs/bcftools/bin/bcftools call -mv - > ~/urchin_af/variants/temp/tmp.{}.vcf"
