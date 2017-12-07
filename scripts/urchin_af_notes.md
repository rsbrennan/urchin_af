# urchin selection experiment

- Experiment run by April. Not sure when.
- all data is capture seq, exons + promoter regions. 
- 25 indivs were spawned. 10 females, 15 males? maybe
- low and ambient pH conditions
- genotyping at day 1 and day 7 in both conditions. Looking for shifts in AF.


Issues/approaches

- working on method to identify shifts in AF
	+ probably CMH
- validate this approach
	+ is there extended LD around selected regions
		* use LDx
	+ simulations to assess power?
- can we detect balancing selection?	

### variant calling

`nohup bash variant_call.sh 2> ~/urchin_af/log_out/variants.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/variants.stdout_$(date +"%F_%R").txt &`

### merge and filter variants:

`nohup bash vcf_filter.sh 2> ~/urchin_af/log_out/filter_stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/filter_stdout_$(date +"%F_%R").txt &`

parameters based on: https://www.nature.com/articles/srep33735


filter: 
- maf 0.025
- qual 75
- bi-allelic only
- min avg depth 10x
- not doing max, can filter that later

After filtering, kept 529,296 out of a possible 9835932 Sites

100x: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0140462

65x: http://onlinelibrary.wiley.com/doi/10.1111/mec.12360/full

### pull out allele freqs/pop

```bash
zcat ~/urchin_af/variants/urchin_dupsincl.vcf.gz | grep -v '^##' | cut -f 1-8 | sed 's/#CHROM/CHROM/g' > ~/urchin_af/analysis/af.info.txt

zcat ~/urchin_af/variants/urchin_dupsincl.vcf.gz | grep -v '^##' | cut -f 10- > ~/urchin_af/analysis/af.out.txt

```

`filter.R`

### filter varaints by depth

`nohup bash filter_final.sh 2> ~/urchin_af/log_out/final_filter.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/final_filter.stdout_$(date +"%F_%R").txt &`

274,180,614 before dups removed
130,427,229 after dups removed

48% of reads remain

### calculate changes in AF

`AF_change.R`

- initially used cmh, D1 to D7, in control and treatment. 
- but, this might not be best approach
	+ see: http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12810/full
- Argue that GLM approach probably more powerful
	+ previous simulation studies didn't look at these methods
- also add:
	+ quasibinomial GLM
	+ CMH comparing D7 control and low pH
	+ look at overlap between all approaches



### estimate LD

Using LDx

need to subsample sam file to "control" and "treatment"


want to compare ld in day 1 with ld in day 7, ph8 and ph7

isolate each group for comparisons. get bam out. these should be around 15g total.

```bash
samtools view -H /data/OASV2/merged.fixmate.sorted.bam > ~/urchin_af/variants/sam.header 

samtools view /data/OASV2/merged.fixmate.sorted.bam   | grep 'OASV2_DNA_D1.*' | cat ~/urchin_af/variants/sam.header - | samtools view -Sb > ~/urchin_af/variants/D1.bam

samtools view /data/OASV2/merged.fixmate.sorted.bam   | grep 'OASV2_DNA_D7_8_0.*' | cat ~/urchin_af/variants/sam.header - | samtools view -Sb > ~/urchin_af/variants/D7_8.bam

samtools view /data/OASV2/merged.fixmate.sorted.bam   | grep 'OASV2_DNA_D7_7_5.*' | cat ~/urchin_af/variants/sam.header - | samtools view -Sb > ~/urchin_af/variants/D7_7.bam

```


calculate ld for each of the 3 groups.

currently just changing 'rep' to `D1`, `D7_8`, `D7_7`, and running each sep. 

`nohup bash ld.sh 2> ~/urchin_af/log_out/ld_stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/ld_stdout_$(date +"%F_%R").txt &`


to plot, etc:
`ld.R`


### overlap of sig loci with genes

make bed file of sig loci:

```r

cmh <-read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE) 

selected <- cmh[which(cmh$pH_response == TRUE & cmh$pH_selection_pval < cut_off),]

#write bed file:
# remember, 1-2 specifies position 2. 0 based
# chr, start, stop

out <- cbind(as.character(selected$CHROM), selected$POS-1, selected$POS)

write.table(out,"~/urchin_af/analysis/cmh.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

```

bed file with annotations: /data/OASV2/RG_7601_Probe.extended.bed


```bash

sort -k1,1 -k2,2n ~/urchin_af/analysis/cmh.bed > ~/urchin_af/analysis/cmh.sorted.bed
sort -k1,1 -k2,2n /data/OASV2/RG_7601_Probe.extended.bed > ~/reference/urchin_probes.sorted.bed
sort -k1,1 -k2,2n ~/reference/SpBase3.1_build8.gff3/Transcriptome.gff3 > ~/reference/SpBase3.1_build8.gff3/Transcriptome.sorted.gff3

bedtools closest -a ~/urchin_af/analysis/cmh.sorted.bed \
-b ~/reference/urchin_probes.sorted.bed \
-D a -iu > ~/urchin_af/analysis/cmh.genes.txt

# convert cmh.genes.txt gene names to those compatible with echinobase 

#get id from 5th column 

cat ~/urchin_af/analysis/cmh.genes.txt | grep   "ID=" > ~/urchin_af/analysis/cmh.genes_name.txt

cat ~/urchin_af/analysis/cmh.genes_name.txt | grep -ohP "ID=.*?;" | sed 's/ID=//g' | sed 's/;//g' | sed 's/"//g'| \
	cut -f 1-2 -d '.' > ~/urchin_af/analysis/cmh.id 

# whl22.v1.0.tmap.gz

#convert based on whl22.v1.0.tmap.gz conversion table.

zcat ~/reference/whl22.v1.0.tmap.gz | awk 'BEGIN{FS="\t";OFS="\t"}FNR==NR { a[$4]=$2; next } $1 in a { $1=a[$1] }1' - ~/urchin_af/analysis/cmh.id  > ~/urchin_af/analysis/cmh.id.converted
 
paste ~/urchin_af/analysis/cmh.genes_name.txt  ~/urchin_af/analysis/cmh.id.converted > ~/urchin_af/analysis/cmh.converted_genes

```



```r

genes <- read.table("~/urchin_af/analysis/cmh.converted_genes", header=FALSE)

out <- (table(genes$V1)[which(table(genes$V1) > 4)])

genes.sub <- subset(genes, V1 %in% names(out))

grep("SPU_000825", genes$V15)

```

Acid-sensitive potassium channel protein TASK-2
	- "Scaffold542:34370" gta to ata, valine to Ile- not a big change
	- "Scaffold542:34445" aaa to taa, lysine to STOP- um. big change
	- "Scaffold542:34450" agc to aga, Serine to Arginine- polar to basic
	- "Scaffold542:34467" gtt to gct, Valine to alanine- not big
	- "Scaffold542:34469" ttt to gtt, Phe to valine - not big
	- SPU_018872
	- WHL22.580025.1

### use snpEFF to determine effect of variants

This is currently not working well. not pulling out transcripts in echinobase

it looks like the Strongylocentrotus_purpuratus genome is already in the database. good. 


I think that was old version. add to config file:
 - Strongylocentrotus_purpuratus genome, version SpBase3.1_build8
   SpBase3.1_build8.genome : Strongylocentrotus_purpuratus

```
~/bin/SNPdat_package_v1.0.5/SNPdat_v1.0.5.pl -i ~/urchin_af/variants/urchin_filtered.vcf -f /data/OASV2/Spur_3.1.LinearScaffold.fa -g ~/reference/whl22.v1.0.gtf.gz -o test.out

```

```bash

# convert gff to gtf

/data/programs/cufflinks-2.2.1.Linux_x86_64/gffread -E ~/reference/SpBase3.1_build8.gff3/GLEAN-UTR-3.1_annotated08282015.gff3 -T -o- test.out  > snpEff.stdout 2> snpEff.stderr



cp ~/reference/SpBase3.1_build8.gff3/GLEAN-UTR-3.1_annotated08282015.gff3 ~/bin/snpEff/data/SpBase3.1_build8/
mv ~/bin/snpEff/data/SpBase3.1_build8/GLEAN-UTR-3.1_annotated08282015.gff3 ~/bin/snpEff/data/SpBase3.1_build8/genes.gff

# put reference fasta in ~/bin/snpEff/data/genomes

cp /data/OASV2/Spur_3.1.LinearScaffold.fa ~/bin/snpEff/data/genomes/
mv ~/bin/snpEff/data/genomes/Spur_3.1.LinearScaffold.fa ~/bin/snpEff/data/genomes/SpBase3.1_build8.fa 

# add fasta to gff
echo "###"  >> ~/bin/snpEff/data/SpBase3.1_build8/genes.gff
echo "##FASTA"  >> ~/bin/snpEff/data/SpBase3.1_build8/genes.gff
cat ~/bin/snpEff/data/genomes/SpBase3.1_build8.fa  >> ~/bin/snpEff/data/SpBase3.1_build8/genes.gff

#mv ~/bin/snpEff/data/SpBase3.1_build8/genes.gff ~/bin/snpEff/data/genomes/


java -Xmx16g -jar ~/bin/snpEff/snpEff.jar build -gff3 -v SpBase3.1_build8   > snpEff.stdout 2> snpEff.stderr

```

```bash

#-v makes it verbose
# specify config file with -c
java -Xmx4g -jar ~/bin/snpEff/snpEff.jar -c ~/bin/snpEff/snpEff.config  -v -noGenome SpBase3.1_build8  ~/urchin_af/variants/urchin_filtered.vcf > ~/urchin_af/variants/urchin_ann.vcf

```



run with ensemble genome

java -Xmx4g -jar ~/bin/snpEff/snpEff.jar -c ~/bin/snpEff/snpEff.config  -v Strongylocentrotus_purpuratus  ~/urchin_af/variants/urchin_filtered.vcf > ~/urchin_af/variants/urchin_ann.vcf > snpEff.stdout 2> snpEff.stderr


where I left of fon fri nov 17:
	- trying to figure out if using snpEff with ensemble genome is ok. 
	- check with gff/gtf from ensemble: Strongylocentrotus_purpuratus.Spur_3.1.37.gff3.gz
	- there is a whl gtf. I could use that? in reference directory/.




## simulations still a work in progress



### simulate data to check power for AF changes

using fastsimcoal

run simulation: `fastsimcoal.sh`

convert to mimicree format: `fastsimcoal_to_hap.py`

need to space delimit 5th col:

```bash

cat ~/urchin_af/analysis/simulations/hap.out | sed 's/ \+ /\t/g' | cut -f 1-4 > ~/urchin_af/analysis/simulations/hap.out.1

cat ~/urchin_af/analysis/simulations/hap.out | sed 's/ \+ /\t/g' | cut -f 5- | sed 's/\t/ /g' > ~/urchin_af/analysis/simulations/hap.out.2

paste ~/urchin_af/analysis/simulations/hap.out.1 ~/urchin_af/analysis/simulations/hap.out.2 > ~/urchin_af/analysis/simulations/mimicree.hap

rm ~/urchin_af/analysis/simulations/hap.out.2
rm ~/urchin_af/analysis/simulations/hap.out.1

```

recombination rate file:

col1: genomic locus in the {{form chromosome:start..end}}
col2: recombination rate at the beginning of the window
col3: recombination rate in the middle of the window
col4: recombination rate in the end of the window



run mimicree: 

```bash


#pick 10 snps to be under selection:

# -s = selection
# -n number selected
# -e effect of het
# -m max frequency of selected allele

for i in {1..10}; do python ~/bin/MimicrEE/scripts/general/pick-random-addtive-snps.py -n 10 -s 0.1 -e 0.5 -m 0.8 --loci-count 72561 --input ~/urchin_af/analysis/simulations/mimicree.hap > ~/urchin_af/analysis/simulations/mimicree_out/snps/s01-n$i.txt; done

##########
#### generate recombination file
##########



for i in {1..10}; do nohup java -Xmx1g -jar  ~/bin/MimicrEE/MimicrEESummary.jar --haplotypes-g0 ~/urchin_af/analysis/simulations/mimicree.hap --recombination-rate 0.195 --output-mode 1 --replicate-runs 1 --output-format sync --threads 6 --additive ~/urchin_af/analysis/simulations/mimicree_out/snps/s01-n$i.txt --output-file ~/urchin_af/analysis/simulations/mimicree_out/sim/s01-n$i.sync; done

# run cmh

mkdir cmh
for i in {1..10}; do perl <path>/popoolation2/cmh-test.pl --min-count 1 --min-coverage 1 --max-coverage 100000 --min-logpvalue 0.0 --population 1-2,3-4,5-6 --input sim/s01-n$i.sync --output cmh/s01-n$i.cmh; done &

# create labs and predictions for ROCR
mkdir simres
python <path>/rocr-generate-labellist.py sim/s01-n1.sync snps/s01-n{1..10}.txt > simres/s01.labels

python <path>/rocr-generate-predictionlist.py cmh/s01-n{1..10}.cmh > simres/s01.predictions

#visualize results

R --vanilla --args simres/s01.labels simres/s01.predictions simres/s0025.labels simres/s0025.predictions roc < ../toysoft/toyroc.R


```

### forqs

see forgs.md



















