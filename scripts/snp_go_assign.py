
## cmh.annotations.out contains spu assignments
# ~/reference/ensembl_goterms.txt were downloaded on 2/18/18

# want to link spu from cmh.annotation.out to spu in ensembl_goterms.txt to get GO terms
#### this should also let me
# secondly, can pull down gene description from wiki gene

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
#spu_path = '/users/r/b/rbrennan/reference/whl22.v1.0.tmap.gz'

#make empty array to save output
go_out = np.empty(shape=(41494,5), dtype = object)

i=0 # start counter
# 41495 rows

# note that in the go terms, there are multiple rows for each spu if mnultiple GO terms are assigned

with open('/users/r/b/rbrennan/urchin_af/analysis/cmh.annotations.out') as master_file:
    header_line = next(master_file) # skip header row
    #start_time = time.time()
    for idx, line in enumerate(master_file):
        tmp_snp = line.split("\t")[0]
        tmp_spu = line.split("\t")[5]
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
                        if len(tmp_go) == 0:
                            tmp_go = go_line.split("\t")[6]
                        if len(tmp_short) == 0:
                            tmp_short = go_line.split("\t")[4]
                        if len(tmp_go) > 0:
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
        tmp_class = line.split("\t")[18].split("\n")[0]
        tmp_short1 = tmp_short.split(";")
        tmp_short2 = ";".join(list(OrderedDict.fromkeys(tmp_short1)))
        tmp_go1 = tmp_go.split(";")
        tmp_go2 = ";".join(list(OrderedDict.fromkeys(tmp_go1)))
        # I think that if no GO term will be blank... but double check
        out_1 = tmp_snp + "\t" + tmp_spu + "\t" + tmp_short2 + "\t" + tmp_class + "\t" + tmp_go2
        #go_out = np.vstack((go_out, out_1.split("\t")))
        go_out[idx]= out_1.split("\t")
        i=i+1
        if i % 1000 == 0: print(i)
        #if i % 1000 == 0: end_time2 = time.time()
        #if i % 1000 == 0: print("total time taken for 1000 loops: ", end_time2 - start_time)


# pull out common gene name ~/reference/annotation.build8/gene_info_table.txt
# in column 8 of this table
# spu in column 2

nm_out = np.empty(shape=(41494,2), dtype = object)
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

sig_out = np.empty(shape=(41494,4), dtype = object)
i=0 # start counter

with open('/users/r/b/rbrennan/urchin_af/analysis/cmh.out.txt') as sig_file:
    header_line = next(sig_file) # skip header row
    for idx, sig_line in enumerate(sig_file):
        tmp_snp = sig_line.split("\t")[0] + ":" + sig_line.split("\t")[1]
        for go_line in go_out:
            if tmp_snp == go_line[0]:
                tmp_chr = sig_line.split("\t")[0]
                tmp_pos = sig_line.split("\t")[1]
                tmp_pval = sig_line.split("\t")[74]
                tmp_sig = sig_line.split("\t")[77].split("\n")[0]
        out_1 = tmp_chr + "\t" + tmp_pos + "\t" + tmp_pval + "\t" + tmp_sig
        sig_out[idx]= out_1.split("\t")
        i=i+1
        if i % 1000 == 0: print(i)

# combine and save all
out_1 = np.column_stack((sig_out, nm_out, go_out))
head = "CHR" + "\t" + "POS" + "\t" + "PVAL"  + "\t" + "sig" + "\t" + "SPU_1" + "\t" + "description" + "\t" + "SNP"  "\t" + "SPU_2"  + "\t" + "short_name" + "\t" + "class" + "\t" + "GO_term"
head = head.split("\t")
out_final = np.vstack((head, out_1))
np.savetxt('/users/r/b/rbrennan/urchin_af/analysis/cmh.master.out', out_final,fmt='%5s', delimiter='\t')
