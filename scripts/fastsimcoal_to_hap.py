mport pandas as pd
import numpy as np

###
#
# need the following format:

#chrom position ancestral_allele ancestral/derived haplotupes(in AA AC)


#read in file:

sim_file = open('/users/r/b/rbrennan/urchin_af/analysis/simulations/test/test_1_2.arp', "r")

print sim_file.readline()

#for line in sim_file:
#   print line

#first: pull out chromosomes and their positions. save in 2 cols

#in following format:

# 2 polymorphic positions on chromosome 1
#151, 533

#make empty matrix. row then col
result_a = np.empty(shape=(0,2))

filepath = '/users/r/b/rbrennan/urchin_af/analysis/simulations/test/test_1_2.arp'
with open(filepath) as fp:
    for line in fp:
        if ("polymorphic positions on chromosome") in line:
            chr_out = line.split()[5:]
            chr_out = '_'.join(chr_out)
            if line.split()[1] != '0': # only pull out polymorphic
                snp_num = (next(fp, '').strip()[1:]).split()
                snp_num = [s.replace(',', '') for s in snp_num]
                chr_new = [chr_out] * (len(snp_num))
                a = np.column_stack((chr_new, snp_num))
                result_a = np.append(result_a, a, axis=0)

#should have 72605

# fastsimcoal is generating haplotypes
# so need to assign 2 haplotypes to each individual

# first assign the haplotypes. need to keep consistent bc of recombination

haps =list(range(1000)) # could modify this to choose from larger range of samples
hap_list = np.empty(shape=(0,2)) #make array to hold assignments
hap_samp = np.random.choice(haps,50,replace=False)# randomly sample haplotypes

# assign haplotypes to indivs
for i in xrange(0,50,2):
    hap_ind = [hap_samp[i:i+2]]
    hap_list = np.append(hap_list, hap_ind, axis=0).astype(int)


# now assign actual values to indivs
geno_all = np.empty(shape=(0,1))

filepath = '/users/r/b/rbrennan/urchin_af/analysis/simulations/test/test_1_2.arp'
with open(filepath) as fp:
    for line in fp:
        if line.startswith('1_'):
            # need to make new object with genos
            geno = line.split()[2:]
            geno_all = np.append(geno_all,geno[0])

# make empty array to hold haplotypes
ind_hap = np.empty(shape=(result_a.shape[0],0))
#now cycle through each individual
for line in hap_list:
    b1 = geno_all[line[0]]
    b2 = geno_all[line[1]]
    gt_all = np.empty(shape=(0,0))
    for i in xrange(0,result_a.shape[0], 1): # want to paste two bases together to make a genotype
        gt = b1[i]+b2[i] # make specific genotype
        gt_all = np.append(gt_all, gt)
    #append to ind_hap
    ind_hap = np.column_stack((ind_hap, gt_all))



# now need ancestral_allele and possible variants (in A/A A/C)

anc_list = np.empty(shape=(0,2))

for i in xrange(0,ind_hap.shape[0],1):
    anc_base = list(set(''.join(ind_hap[i])))
    if len(anc_base) == 2:
        bases = '/'.join((anc_base[0], anc_base[1]))
        ac1 = anc_base[0]
        a_temp = np.column_stack((ac1, bases))
        anc_list = np.vstack((anc_list, a_temp))
    else:
        bases = '/'.join((anc_base[0], anc_base[0]))
        ac1 = anc_base[0]
        a_temp = np.column_stack((ac1, bases))
        anc_list = np.vstack((anc_list, a_temp))

###
# now joint everything together
###

# in this order: result_a anc_list ind_hap

final_out = np.column_stack((result_a, anc_list, ind_hap))

# remove invariant sites

final_new = np.empty(shape=(0,final_out.shape[1]))

for line in final_out:
    geno = line[3].split("/")
    if geno[0] != geno[1]:
        final_new = np.vstack((final_new, line))



np.savetxt('/users/r/b/rbrennan/urchin_af/analysis/simulations/hap.out', final_new,fmt='%5s')
# the format option is necessary. not exactly sure what it is doing though.







