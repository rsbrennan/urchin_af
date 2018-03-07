cd ~/bin/snpEff/data/Strongylocentrotus_purpuratus
zcat snpEffectPredictor.bin | grep -v 'EMSPUG' | gzip > snpEffectPredictor.bin1
mv snpEffectPredictor.bin oldannotation
mv snpEffectPredictor.bin1 snpEffectPredictor.bin

java -Xmx4g -jar ~/bin/snpEff/snpEff.jar -c ~/bin/snpEff/snpEff.config  -v Strongylocentrotus_purpuratus  ~/urchin_af/variants/urchin_final.vcf.gz  > ~/urchin_af/variants/urchin_ann.vcf > ~/urchin_af/log_out/snpEff.stdout 2> ~/urchin_af/log_out/snpEff.stderr

