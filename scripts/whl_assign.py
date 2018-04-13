import pandas as pd
import numpy as np
import gzip
import csv
import itertools
import re
from collections import OrderedDict
import time

# assign sig or not:
#with open('/users/r/b/rbrennan/urchin_af/analysis/cmh.out.sorted.txt') as inf:
#    reader = csv.reader(inf, delimiter="\t")
#    sig_col = list(zip(*reader))[78]

#sig_col = sig_col[1:]
# assign GO terms
filepath = '/users/r/b/rbrennan/reference/blast2go-whl.nospace.txt.gz'
spu_path = '/users/r/b/rbrennan/reference/whl22.v1.0.tmap.gz'
####
## all loci
####

#make empty array
go_out = np.empty(shape=(75368,2), dtype = object)
spu_out = np.empty(shape=(75368,2), dtype = object)

i=0 # start counter

# in col 3 of whl-annot, there are gene names. add these
with open('/users/r/b/rbrennan/urchin_af/analysis/cmh.all.id') as master_file:
        for idx, line in enumerate(master_file):
            tmp_go = line.split("\n")[0]
            tmp_spu = line.split("\n")[0]
            # assign spu values
            with gzip.open(spu_path) as spu_file:
                for spu_line in spu_file:
                    #make empty array to add results to
                    if line.split("\n")[0] == spu_line.split("\t")[3].split("\n")[0]:
                        if len(tmp_spu.split()) == 1:
                            #print(gene_line)
                            tmp_spu = tmp_spu + "\t" + spu_line.split("\t")[1].split("\n")[0]

                        if len(tmp_spu.split()) > 1:
                            tmp_spu = tmp_spu + ";" +  spu_line.split("\t")[1].split("\n")[0]

            if len(tmp_spu.split()) == 1:
                tmp_spu = tmp_spu + "\t" + "NA"
            tmp_spu1 = tmp_spu.split()
            tmp_spu2 = [s.replace('tr', '') for s in tmp_spu1]
            tmp_spu3 = tmp_spu2[1].split(";")
            seen = set()
            result = []
            for item in tmp_spu3:
                if item not in seen:
                    seen.add(item)
                    result.append(item)
            tmp_spu4 = ";".join(result)
            tmp_spu5 = tmp_spu.split("\t")[0] + "\t" + tmp_spu4.split()[0]

            spu_out[idx] = tmp_spu5.split()

            # match GO terms
            with gzip.open(filepath) as file:
                for gene_line in file:
                    #make empty array to add results to
                    if line.split("\n")[0] == gene_line.split("\t")[0]:
                        if len(tmp_go.split()) == 1:
                            #print(gene_line)
                            if len(gene_line.split("\t")) == 3:
                                tmp_go = tmp_go + "\t" + gene_line.split("\t")[2].split("\n")[0]
                                #print(gene_line)
            if len(tmp_go.split()) == 1:
                tmp_go = tmp_go + "\t" + "NA"

            go_out[idx] = tmp_go.split()
            i=i+1
            if i % 1000 == 0: print(i)

gene_name_out = [i[1] for i in go_out]

#np.savetxt('/users/r/b/rbrennan/urchin_af/analysis/cmh.all.GO', go_out,fmt='%5s', delimiter='\t')

#combine go and spu
go_spu_out = np.column_stack((spu_out,gene_name_out))


# pull out snp ID
i=0
snp_out = np.empty(shape=(75368,14), dtype = object)

with open('/users/r/b/rbrennan/urchin_af/analysis/cmh.master.sort.out') as master_file:
    head1 = next(master_file)
    for idx, line in enumerate(master_file):
        snp_out[idx] =  line.split("\t")
        tmp_go = snp_out[idx,13].split("\n")[0]
        snp_out[idx,13] = tmp_go
        i=i+1
        if i % 10000 == 0: print(i)

# then join together sig or not, snp name, etc

out_1 = np.column_stack((snp_out, go_spu_out))
head2 = head1.split("\n")[0]
head3 = head2 + "\t" + "WHL" + "\t" + "WHL_SPU" + "\t" + "WHL_name"
head_final = head3.split("\t")
head_final = [x.strip(' ') for x in head_final]
out_final = np.vstack((head_final, out_1))
np.savetxt('/users/r/b/rbrennan/urchin_af/analysis/cmh.master_goodnm.out', out_final,fmt='%5s', delimiter='\t')
