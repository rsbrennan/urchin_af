import pandas as pd
import numpy as np
import gzip
import csv
import itertools
import re
from collections import OrderedDict

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
with open('/users/r/b/rbrennan/urchin_af/analysis/cmh.all_dups_in.ID.txt') as master_file:
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
                    # check if whl matches
                    if line.split("\n")[0] == gene_line.split("\t")[0]:
                        # if tmp_go is 1, this means only whl is present
                        if len(tmp_go.split()) == 1:
                            # check if gene name is present
                            if len(gene_line.split("\t")) == 3:
                                # if gene name is present, add it
                                tmp_go = tmp_go + "\t" + gene_line.split("\t")[2].split("\n")[0]
                                # then add GO terms
                                tmp_go = tmp_go + "\t" + gene_line.split("\t")[1].split("\n")[0]
                            # if gene_line length is 2, there is no gene name on that line. just add go Term
                            if len(gene_line.split("\t")) == 2:
                                #tmp_go = tmp_go + "\t" + "NA"
                                tmp_go = tmp_go + "\t" + gene_line.split("\t")[1].split("\n")[0]
                        # if GO term has already been added, just append additional GO terms
                        if len(tmp_go.split()) > 1:
                            tmp_go = tmp_go + ";" +  gene_line.split("\t")[1].split("\n")[0]
            # if tmp_go length is 1, there are no matches. Add unknown, etc.
            if len(tmp_go.split()) == 1:
                tmp_go = tmp_go + "\t" + "NA" + "\t" + "unknown"
            # if length is 2, there is no gene name. fill it in with NA
            if len(tmp_go.split()) == 2:
                tmp_go = tmp_go.split("\t")[0] + "\t" + "NA" + "\t" + tmp_go.split("\t")[1]
            go_out = np.vstack((go_out, tmp_go.split()))
            i=i+1
            if i % 1000 == 0: print(i)

np.savetxt('/users/r/b/rbrennan/urchin_af/analysis/cmh.dupsincl.GO', go_out, fmt='%5s', delimiter='\t')

# only if not rerunning the whole thing
#go_out = open('/users/r/b/rbrennan/urchin_af/analysis/cmh.dupsincl.GO', 'r').read().split("\t")
#go_outn = np.asarray(go_out)

# after this, have 42840 rows in go_out.
# this is because some are duplicates. Go through and merge those

# ~/urchin_af/analysis/cmh.all_dups_in.genes.txt
# this file only has duplicates where multple gene assignments exist
# corresponds to go_out
# so merge based on SNP name in ~/urchin_af/analysis/cmh.all_dups_in.genes.txt

# first, read in file of all variants, including duplicates and make snp name
with open('/users/r/b/rbrennan/urchin_af/analysis/cmh.all_dups_in.sorted.genes.txt') as inf:
    reader = csv.reader(inf, delimiter="\t")
    snp_1, snp_2 = list(zip(*reader))[0:2]

snp_id = np.empty(shape=(0,1))

for f,b in itertools.izip(snp_1,snp_2):
    nm_tmp = f + ":" + b
    snp_id = np.vstack((snp_id, nm_tmp))

# parse down this snp_id list to unique snps only.
new_array = [tuple(row) for row in snp_id]

# the following thsould keep the order
snp_uniq = list(OrderedDict.fromkeys(new_array))

# pull out matches between unique snps, dup snps,
  # then make new array from go_out

go_uniq = np.empty(shape=(0,3))
spu_uniq = np.empty(shape=(0,2))

i=0 # start counter
dups=0
for line in snp_uniq:
    B = np.where(snp_id == line)[0]
    # if length of match == 1, just add the original line from go_out
    if len(B) == 1:
        #print("start")
        tmp_A = go_out[B]
        #print("tmpA done")
        tmp_B = spu_out[B]
        #print("tmpB")
    # if length of match is >1, need to combine entries
    elif len(B) > 1:
        dups=dups+1
        #print("more than 1!")
        tmp_whl = ['']
        tmp_nm = ['']
        tmp_gterm = ['']
        tmp_spu = ['']
        for un_len in range(len(B)):
            tmp_1 = go_out[B[un_len]]
            tmp_2 = spu_out[B[un_len]]
            if(len(tmp_1[0]) > 0):
                tmp_whl = tmp_whl[0] + tmp_1[0] + ";"
                tmp_whl = tmp_whl.split()
            if(len(tmp_1[1]) > 0):
                tmp_nm = tmp_nm[0] + tmp_1[1] + ";"
                tmp_nm = tmp_nm.split()
            if(len(tmp_1[2]) > 0):
                tmp_gterm = tmp_gterm[0] + tmp_1[2] + ";"
                tmp_gterm = tmp_gterm.split()
            if(len(tmp_2[1]) > 0):
                tmp_spu = tmp_spu[0] + tmp_2[1] + ";"
                tmp_spu = tmp_spu.split()
        # save go output
        tmp_A = tmp_whl[0] + "\t" + tmp_nm[0] + "\t" + tmp_gterm[0]
        # save spu output
        tmp_B = tmp_whl[0] + "\t" + tmp_spu[0]
    #for some reason, sometimes tmp_A is a string, sometimes a nparray. the following deals with it.
    if isinstance(tmp_A, basestring):
        go_uniq = np.vstack((go_uniq, tmp_A.split()))
    else:
        go_uniq = np.vstack((go_uniq, tmp_A.tolist()))
    if isinstance(tmp_B, basestring):
        spu_uniq = np.vstack((spu_uniq, tmp_B.split()))
    else:
        spu_uniq = np.vstack((spu_uniq, tmp_B.tolist()))
    i=i+1
    if i % 2500 == 0: print(i)

# save output
np.savetxt('/users/r/b/rbrennan/urchin_af/analysis/cmh.all.GO', go_uniq, fmt='%5s', delimiter='\t')

# format spu
# remove duplicate spus, get rid of tr at end
spu_new = np.empty(shape=(0,2))

i=0 # start counter
for line in spu_uniq:
    tmp_whl = re.sub(';$', '', line[0])
    tmp_spu = re.sub('tr', '', line[1])
    tmp_spu = re.sub(';$', '', tmp_spu)
    #remove duplicates
    tmp_uniq = list(set(tmp_spu.split(";"))) # note that this won't preserve order of spus
    tmp_uniq = ";".join(tmp_uniq)
    tmp_new = tmp_whl + "\t" + tmp_uniq
    spu_new = np.vstack((spu_new, tmp_new.split()))
    i=i+1
    if i % 2500 == 0: print(i)

# pull out go term, gene name
go_out_go = [i[2] for i in go_uniq]
gene_name_out = [i[1] for i in go_uniq]
#spu_out = [i[1] for i in spu_new]

#combine go and spu
go_spu_out = np.column_stack((gene_name_out,spu_new, go_out_go))

# pull out snp ID
i=0
snp_out = np.empty(shape=(0,4))
with open('/users/r/b/rbrennan/urchin_af/analysis/cmh.out.sorted.txt') as master_file:
    next(master_file)
    for line in master_file:
        # make snp name
        snp_nm = line.split("\t")[0] + ":" + line.split("\t")[1]
        # chr and pos
        chr_out = line.split("\t")[0]
        pos_out = line.split("\t")[1]
        # p value: corresponds to col pH_selection_pval
        p_val = line.split("\t")[74]
        # combine snp and pval
        out_2 = np.column_stack((snp_nm, chr_out,pos_out, p_val))
        snp_out = np.vstack((snp_out, out_2))
        i=i+1
        if i % 10000 == 0: print(i)


# assign sig or not:
with open('/users/r/b/rbrennan/urchin_af/analysis/cmh.out.sorted.txt') as inf:
    reader = csv.reader(inf, delimiter="\t")
    sig_col = list(zip(*reader))[77]

sig_col = sig_col[1:]


# then join together sig or not, snp name, etc

out_1 = np.column_stack((snp_out, sig_col))



out_3 = np.column_stack((out_1, go_spu_out))
head = "SNP" + "\t" + "CHR" + "\t" + "POS" + "\t" + "PVAL"  + "\t" + "sig"  + "\t" "gene_name" + "\t" "WHL" + "\t" "SPU"  + "\t" + "GO"
head = head.split("\t")
out_final = np.vstack((head, out_3))
np.savetxt('/users/r/b/rbrennan/urchin_af/analysis/cmh.master.out', out_final,fmt='%5s', delimiter='\t')
