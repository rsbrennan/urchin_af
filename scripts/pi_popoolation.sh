#sort the bam file
# D1 bam is D1 ph 8.0

~/bin/samtools-1.6/samtools view -q 20 -b ~/urchin_af/variants/D1.bam | ~/bin/samtools-1.6/samtools sort - -o ~/urchin_af/variants/D1.sort.bam

~/bin/samtools-1.6/samtools view -q 20 -b ~/urchin_af/variants/D7_8.bam | ~/bin/samtools-1.6/samtools sort -  -o ~/urchin_af/variants/D7_8.sort.bam

~/bin/samtools-1.6/samtools view -q 20 -b ~/urchin_af/variants/D7_7.bam | ~/bin/samtools-1.6/samtools sort - -o ~/urchin_af/variants/D7_7.sort.bam

#convert to pileup:
~/bin/samtools-1.6/samtools mpileup ~/urchin_af/variants/D1.sort.bam > ~/urchin_af/variants/D1_8.sort.pileup
~/bin/samtools-1.6/samtools mpileup ~/urchin_af/variants/D7_8.sort.bam > ~/urchin_af/variants/D7_8.sort.pileup
~/bin/samtools-1.6/samtools mpileup ~/urchin_af/variants/D7_7.sort.bam > ~/urchin_af/variants/D7_7.sort.pileup

# calc pi
perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/urchin_af/variants/D1_8.sort.pileup --output ~/urchin_af/variants/D1_8.pi --measure pi --window-size 400 --step-size 200 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 20  --min-covered-fraction 0.5 --pool-size 1000

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/urchin_af/variants/D7_7.sort.pileup --output ~/urchin_af/variants/D7_7.pi --measure pi --window-size 400 --step-size 200 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 20 --min-covered-fraction 0.5 --pool-size 1000

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/urchin_af/variants/D7_8.sort.pileup --output ~/urchin_af/variants/D7_8.pi --measure pi --window-size 400 --step-size 200 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 20 --min-covered-fraction 0.5 --pool-size 1000


# tajima's D

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/urchin_af/variants/D1_8.sort.pileup --output ~/urchin_af/variants/D1_8.TD --measure D --window-size 400 --step-size 200 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 20 --min-covered-fraction 0.5 --pool-size 1000

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/urchin_af/variants/D7_7.sort.pileup --output ~/urchin_af/variants/D7_7.TD --measure D --window-size 400 --step-size 200 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 20 --min-covered-fraction 0.5 --pool-size 1000
