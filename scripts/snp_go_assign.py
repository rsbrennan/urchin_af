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
filepath = '/users/r/b/rbrennan/reference/blast2go-whl.nospace.txt.gz'
spu_path = '/users/r/b/rbrennan/reference/whl22.v1.0.tmap.gz'
####
## all loci
####

#make empty array
go_out = np.empty(shape=(0,3))
spu_out = np.empty(shape=(0,2))

i=0 # start counter

# in col 3 of whl-annot, there are gene names. add these
with open('/users/r/b/rbrennan/urchin_af/analysis/cmh.all.id') as master_file:
        for line in master_file:
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
            spu_out = np.vstack((spu_out, tmp_spu.split()))

            # match GO terms
            with gzip.open(filepath) as file:
                for gene_line in file:
                    #make empty array to add results to
                    if line.split("\n")[0] == gene_line.split("\t")[0]:
                        if len(tmp_go.split()) == 1:
                            #print(gene_line)
                            if len(gene_line.split("\t")) == 3:
                                tmp_go = tmp_go + "\t" + gene_line.split("\t")[2].split("\n")[0]
                                tmp_go = tmp_go + "\t" + gene_line.split("\t")[1].split("\n")[0]
                                #print(gene_line)
                            if len(gene_line.split("\t")) == 2:
                                #tmp_go = tmp_go + "\t" + "NA"
                                tmp_go = tmp_go + "\t" + gene_line.split("\t")[1].split("\n")[0]
                        if len(tmp_go.split()) > 1:
                            tmp_go = tmp_go + ";" +  gene_line.split("\t")[1].split("\n")[0]
            if len(tmp_go.split()) == 1:
                tmp_go = tmp_go + "\t" + "NA" + "\t" + "unknown"

            go_out = np.vstack((go_out, tmp_go.split()))
            i=i+1
            if i % 50 == 0: print(i)

go_out_go = [i[2] for i in go_out]
gene_name_out = [i[1] for i in go_out]

np.savetxt('/users/r/b/rbrennan/urchin_af/analysis/cmh.all.GO', go_out,fmt='%5s', delimiter='\t')

#combine go and spu
go_spu_out = np.column_stack((gene_name_out,spu_out, go_out_go))


# pull out snp ID
i=0
snp_out = np.empty(shape=(0,5))
with open('/users/r/b/rbrennan/urchin_af/analysis/cmh.out.sorted.txt') as master_file:
    next(master_file)
    for line in master_file:
        # make snp name
        snp_nm = line.split("\t")[0] + ":" + line.split("\t")[1]
        # chr and pos
        chr_out = line.split("\t")[0]
        pos_out = line.split("\t")[1]
        # p value: corresponds to col pH_selection_pval
        p_val = line.split("\t")[79]
        #test statistic corresponds to col pH_selection_stat
        test_stat = line.split("\t")[83]
        # combine snp and pval
        out_2 = np.column_stack((snp_nm, chr_out,pos_out, p_val,test_stat ))
        snp_out = np.vstack((snp_out, out_2))
        i=i+1
        if i % 10000 == 0: print(i)

# then join together sig or not, snp name, etc

out_1 = np.column_stack((snp_out, sig_col))

out_3 = np.column_stack((out_1, go_spu_out))
head = "SNP" + "\t" + "CHR" + "\t" + "POS" + "\t" + "PVAL" + "\t" "test_stat" + "\t" + "sig"  + "\t" "gene_name" + "\t" "WHL" + "\t" "SPU"  + "\t" + "GO"
head = head.split("\t")
out_final = np.vstack((head, out_3))
np.savetxt('/users/r/b/rbrennan/urchin_af/analysis/cmh.master.out', out_final,fmt='%5s', delimiter='\t')


