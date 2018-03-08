#!/bin/bash


# this script will run the entire filtering process. simiply execute the entire script from the command line

# initial filtering.

cd ~/urchin_af/variants/

bash ~/urchin_af/scripts/01.2_vcf_filter.sh 2> ~/urchin_af/log_out/01.2_vcf_filter.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/01.2_vcf_filter.stdout_$(date +"%F_%R").txt

wait $!

# filter by depth

Rscript ~/urchin_af/scripts/01.3_filter.R 2> ~/urchin_af/log_out/01.3_filter.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/01.3_filter.stdout_$(date +"%F_%R").txt

wait $!

bash ~/urchin_af/scripts/01.4_filter_depth.sh 2> ~/urchin_af/log_out/01.4_filter_depth.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/01.4_filter_depth.stdout_$(date +"%F_%R").txt

wait $!


# filter to on target only

bash ~/urchin_af/scripts/01.5_probe_distance.sh 2> ~/urchin_af/log_out/01.5_probe_distance.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/01.5_probe_distance.stdout_$(date +"%F_%R").txt

wait $!

Rscript ~/urchin_af/scripts/01.6_probe_distance.R 2> ~/urchin_af/log_out/01.6_probe_distance.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/01.6_probe_distance.stdout_$(date +"%F_%R").txt

wait $!

# final filtering

bash ~/urchin_af/scripts/01.7_filter_final.sh 2> ~/urchin_af/log_out/01.7_filter_final.stderr_$(date +"%F_%R").txt 1> ~/urchin_af/log_out/01.7_filter_final.stdout_$(date +"%F_%R").txt
