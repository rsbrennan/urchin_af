# 07


## cmh.annotations.out contains spu assignments
# ~/reference/ensembl_goterms.txt were downloaded on 2/18/18

# want to link spu from cmh.annotation.out to spu in ensembl_goterms.txt to get GO terms
# secondly, can pull down gene description from wiki gene

# some of these are no longer necessary...
import pandas as pd
import numpy as np
import gzip
import csv
import itertools
import re
from collections import OrderedDict
import time

# assign GO terms
go_path = '/users/r/b/rbrennan/reference/ensembl_goterms.txt'

#make empty array to save output. This speeds things up compared to appending each iteration.
# appending requires resaving entire array
go_out = np.empty(shape=(75368,5), dtype = object)

i=0 # start counter. not necessary, but useful when trouble shooting

# note that in the go terms, there are multiple rows for each spu if mnultiple GO terms are assigned

with open('/users/r/b/rbrennan/urchin_af/analysis/cmh.annotations.out') as master_file:
    header_line = next(master_file) # skip header row
    #start_time = time.time()
    for idx, line in enumerate(master_file):
        tmp_snp = line.split("\t")[0] # pull out snp name
        tmp_spu = line.split("\t")[5] # pull out sput
        tmp_go = ""
        tmp_short = ""

        with open(go_path) as go_file:
            for go_line in go_file:
                # check if spu matches each row of go terms
                # but some spu are between two genes. want to use both
                if(len(tmp_spu.split("-")) == 2):
                    spu1 = tmp_spu.split("-")[0]
                    spu2 = tmp_spu.split("-")[1]
                    # match spu1
                    if spu1 == go_line.split("\t")[0]:
                        if len(tmp_go) == 0:  # if this is the first match, fill
                            tmp_go = go_line.split("\t")[6]
                        if len(tmp_short) == 0:
                            tmp_short = go_line.split("\t")[4]
                        if len(tmp_go) > 0: # if this isn't the first match, append.
                            tmp_go = tmp_go + ";" +  go_line.split("\t")[6]
                        if len(tmp_short) > 0:
                            tmp_short = go_line.split("\t")[4] + ";" + tmp_short

                    if spu2 == go_line.split("\t")[0]:
                        if len(tmp_go) == 0:
                            tmp_go = go_line.split("\t")[6]
                        if len(tmp_short) == 0:
                            tmp_short = go_line.split("\t")[4]
                        if len(tmp_go) > 0:
                            tmp_go = tmp_go + ";" +  go_line.split("\t")[6]
                        if len(tmp_short) > 0:
                            tmp_short = tmp_short + ";" +  go_line.split("\t")[4]
                else:
                    if tmp_spu == go_line.split("\t")[0]:
                    # if tmp go is empty, assign all values to it
                        if len(tmp_go) == 0:
                            tmp_go = go_line.split("\t")[6]
                            tmp_short= go_line.split("\t")[4]
                        # is not empty, already have names, etc. just add additional go terms
                        if len(tmp_go) > 0:
                            tmp_go = tmp_go + ";" +  go_line.split("\t")[6]
        # pull out class
        tmp_class = line.split("\t")[19].split("\n")[0]
        tmp_short1 = tmp_short.split(";")
        tmp_short2 = ";".join(list(OrderedDict.fromkeys(tmp_short1)))
        tmp_go1 = tmp_go.split(";")
        tmp_go2 = ";".join(list(OrderedDict.fromkeys(tmp_go1)))
        # if no GO term will be blank
        out_1 = tmp_snp + "\t" + tmp_spu + "\t" + tmp_short2 + "\t" + tmp_class + "\t" + tmp_go2
        go_out[idx]= out_1.split("\t")
        i=i+1
        if i % 10000 == 0: print(i)


# pull out common gene name ~/reference/annotation.build8/gene_info_table.txt
# in column 8 of this table
# spu in column 2

nm_out = np.empty(shape=(75368,2), dtype = object)
i=0 # start counter

for idx, go_line in enumerate(go_out):
    with open('/users/r/b/rbrennan/reference/annotation.build8/gene_info_table.txt') as nm_file:
        header_line = next(nm_file) # skip header row
        for nm_line in nm_file:
            tmp_spu = nm_line.split("\t")[1]
            if tmp_spu == go_line[1]:
                tmp_nm = nm_line.split("\t")[7]
    out1 = go_line[1] + "\t" + tmp_nm
    nm_out[idx] = out1.split("\t")
    i=i+1
    if i % 1000 == 0: print(i)


# also need to pull from ~/urchin_af/analysis/cmh.out.txt to get sig, etc.

sig_out = np.empty(shape=(75368,14), dtype = object)
i=0 # start counter

for idx, go_line in enumerate(go_out):
    tmp_chr = ""
    tmp_pos = ""
    tmp_pval = ""
    tmp_sig = ""
    tmp_snpnm = ""
    with open('/users/r/b/rbrennan/urchin_af/analysis/adaptive_allelefreq.txt') as sig_file:
        header_line = next(sig_file) # skip header row
        for sig_line in sig_file:
            tmp_snp = sig_line.split("\t")[0]+ ":" + sig_line.split("\t")[1]
            if tmp_snp == go_line[0]:
                tmp_chr = sig_line.split("\t")[0]
                tmp_pos = sig_line.split("\t")[1]
                tmp_ph7_pval = sig_line.split("\t")[51]
                tmp_ph8_pval = sig_line.split("\t")[50]
                tmp_sig_ph7 = sig_line.split("\t")[52]
                tmp_sig_ph8 = sig_line.split("\t")[53]
                tmp_af = sig_line.split("\t")[54:57]
                tmp_delta = sig_line.split("\t")[65:68]
                tmp_delta2 = sig_line.split("\t")[68].split("\n")[0]
                tmp_snpnm = tmp_chr + ":" + tmp_pos
    out_1 = tmp_chr + "\t" + tmp_pos + "\t" + tmp_snpnm + "\t" + tmp_ph7_pval + "\t" + tmp_ph8_pval + "\t" +  tmp_sig_ph7 + "\t" + tmp_sig_ph8 + "\t" + "\t".join(tmp_af) + "\t" + "\t".join(tmp_delta) + "\t" + tmp_delta2
    sig_out[idx]= out_1.split("\t")
    i=i+1
    if i % 5000 == 0: print(i)




# combine and save all
out_1 = np.column_stack((sig_out, nm_out, go_out))
head = "CHR" + "\t" + "POS" + "\t" + "SNP_1" + "\t" + "pval_pH75" + "\t" + "pval_pH80" + "\t" + "sig_pH75"+ "\t" + "sig_pH80" + "\t" + "D1_8_af" + "\t" + "D7_7_af" + "\t" + "D7_8_af" + "\t" + "mean_delta_75" + "\t" + "mean_delta_80" + "\t" + "sd_delta_75" + "\t" + "sd_delta_80" + "\t" + "SPU_1" + "\t" + "description" + "\t" + "SNP_2"  "\t" + "SPU_2"  + "\t" + "short_name" + "\t" + "class" + "\t" + "GO_term"
head = head.split("\t")
out_final = np.vstack((head, out_1))
np.savetxt('/users/r/b/rbrennan/urchin_af/analysis/cmh.master.out', out_final,fmt='%5s', delimiter='\t')
