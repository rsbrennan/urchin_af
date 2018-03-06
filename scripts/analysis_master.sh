# driver bash script to exectute urchin selection experiment analysis
# the following executes data analysis following filtering

DATE=`date +%Y-%m-%d`
now=$(date +"%T")

echo "data analysis script run on ${DATE} at ${now}"

Rscript ~/urchin_af/scripts/AF_change.R 2> ~/urchin_af/log_out/AF_change.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/AF_change.stdout_$(date +"%F_%R").txt

wait $!

if [ $? -eq 0 ]
then
  echo "cmh analysis successful"
else
  echo "cmh analysis failed"
  exit 1
fi

##############
## snp eff
#############
# format data base
cd ~/bin/snpEff/data/Strongylocentrotus_purpuratus
zcat snpEffectPredictor.bin | grep -v 'EMSPUG' | gzip > snpEffectPredictor.bin1
mv snpEffectPredictor.bin oldannotation
mv snpEffectPredictor.bin1 snpEffectPredictor.bin

# execute snp eff

java -Xmx4g -jar ~/bin/snpEff/snpEff.jar -c ~/bin/snpEff/snpEff.config  -v Strongylocentrotus_purpuratus  ~/urchin_af/variants/urchin_final.vcf.gz  > ~/urchin_af/variants/urchin_ann.vcf > snpEff.stdout 2> snpEff.stderr

wait

# look at distribution of af, annotations, balancing selection with

Rscript ~/urchin_af/scripts/balancing_sel.R 2> ~/urchin_af/log_out/balancing_sel.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/balancing_sel.stdout_$(date +"%F_%R").txt

########################
#
# go assignment
#
########################

python  ~/urchin_af/scripts/snp_go_assign.py 2> ~/urchin_af/log_out/go_assign_stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/go_assign_stdout_$(date +"%F_%R").txt

# outputs ~/urchin_af/analysis/cmh.master.out

wait $!

if [ $? -eq 0 ]
then
  echo "go assign successful"
else
  echo "go assign failed"
  exit 1
fi


#######################
#
# go enrichment
#
#######################

#format 1
Rscript ~/urchin_af/scripts/go_format.R 2> ~/urchin_af/log_out/go_format_1.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/go_format_1.stdout_$(date +"%F_%R").txt

wait $!

# format 2

bash ~/urchin_af/scripts/go_format.sh 2> ~/urchin_af/log_out/go_format_2_stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/go_format_2.stdout_$(date +"%F_%R").txt

wait $!

# run enrichment:

Rscript ~/urchin_af/scripts/go_enrich.R 2> ~/urchin_af/log_out/go_enrich.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/go_enrich.stdout_$(date +"%F_%R").txt


#











# prep files for ldx
# isolate each group for comparisons. get bam out. these should be around 15g total.

bash ~/urchin_af/scripts/ld.sh 2> ~/urchin_af/log_out/ld_stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/ld_stdout_$(date +"%F_%R").txt

wait $!

if [ $? -eq 0 ]
then
  echo "ldx successful"
else
  echo "ldx failed"
  exit 1
fi
# plot ld results

Rscript ~/urchin_af/scripts/ld.R 2> ~/urchin_af/log_out/depth_filt.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/depth_filt.stdout_$(date +"%F_%R").txt

wait $!

if [ $? -eq 0 ]
then
  echo "ldx plotting successful"
else
  echo "ldx plotting failed"
  exit 1
fi
# variance in af

Rscript ~/urchin_af/scripts/af_variance.R 2> ~/urchin_af/log_out/af_variance.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/af_variance.stdout_$(date +"%F_%R").txt

wait $!

if [ $? -eq 0 ]
then
  echo "allele freq variance calc successful"
else
  echo "allele freq variance calc failed"
  exit 1
fi

### characterize prop of regulatory vs coding SNPs

bash ~/urchin_af/scripts/coding_vs_non.sh 2> ~/urchin_af/log_out/coding_vs_non.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/coding_vs_non.stdout_$(date +"%F_%R").txt

wait $!

if [ $? -eq 0 ]
then
 echo "coding vs non-coding count successful"
else
  echo "coding vs non-coding count failed"
  exit 1
fi

# overlap of sig loci with genes

wait $!

if [ $? -eq 0 ]
then
 echo "overlap of loci with genes successful"
else
  echo "overlap with loci with genes failed"
  exit 1
fi

bash ~/urchin_af/scripts/snp_gene_overlap.sh 2> ~/urchin_af/log_out/snp_gene_overlap.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/snp_gene_overlap.stdout_$(date +"%F_%R").txt


wait $!

if [ $? -eq 0 ]
then
     echo "coding vs non-coding count successful"
 else
       echo "coding vs non-coding count failed"
         exit 1
     fi

     bash ~/urchin_af/scripts/snp_gene_overlap.sh 2> ~/urchin_af/log_out/snp_gene_overlap.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/snp_gene_overlap.stdout_$(date +"%F_%R").txt

