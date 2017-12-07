import pandas as pd
import numpy as np
import gzip
import csv

# assign sig or not:
with open('/users/r/b/rbrennan/urchin_af/analysis/cmh.out.sorted.txt') as inf:
    reader = csv.reader(inf, delimiter="\t")
    sig_col = list(zip(*reader))[78]

sig_col = sig_col[1:]
# assign GO terms
filepath = '/users/r/b/rbrennan/reference/blast2go-whl.annot.txt.gz'

####
## all loci
####

#make empty array
go_out = np.empty(shape=(0,2))
i=0
with open('/users/r/b/rbrennan/urchin_af/analysis/cmh.all.id') as master_file:
        for line in master_file:
            tmp_go = line.split("\n")[0]
            with gzip.open(filepath) as file:
                for gene_line in file:
                    #print(gene_line)
                    #make empty array to add results to
                    if line.split("\n")[0] == gene_line.split("\t")[0]:
                        if len(tmp_go.split()) == 1:
                            tmp_go = tmp_go + "\t" + gene_line.split("\t")[1].split("\n")[0]
                        if len(tmp_go.split()) > 1:
                            tmp_go = tmp_go + ";" +  gene_line.split("\t")[1].split("\n")[0]
            if len(tmp_go.split()) == 1:
                tmp_go = tmp_go + "\t" + "unknown"
            go_out = np.vstack((go_out, tmp_go.split()))
            i=i+1
            if i % 500 == 0: print(i)

np.savetxt('/users/r/b/rbrennan/urchin_af/analysis/cmh.all.GO', go_out,fmt='%5s', delimiter='\t')

# then join together sig or not, snp name, etc

out_1 = np.column_stack((go_out, sig_col))

# pull out snp ID

snp_out = np.empty(shape=(0,4))
with open('/users/r/b/rbrennan/urchin_af/analysis/cmh.out.sorted.txt') as master_file:
    next(master_file)
    for line in master_file:
        # make snp name
        snp_nm = line.split("\t")[0] + ":" + line.split("\t")[1]
        # chr and pos
        chr_out = line.split("\t")[0]
        pos_out = line.split("\t")[1]
        # p value
        p_val = line.split("\t")[75]
        # combine snp and pval
        out_2 = np.column_stack((snp_nm, chr_out,pos_out, p_val))
        snp_out = np.vstack((snp_out, out_2))

out_3 = np.column_stack((snp_out, out_1))
head = "SNP" + "\t" + "CHR" + "\t" + "POS" + "\t" + "PVAL" + "\t" "WHL" + "\t" + "GO" + "\t" + "sig"
head = head.split("\t")
out_final = np.vstack((head, out_3))
np.savetxt('/users/r/b/rbrennan/urchin_af/analysis/cmh.master.out', out_final,fmt='%5s', delimiter='\t')
