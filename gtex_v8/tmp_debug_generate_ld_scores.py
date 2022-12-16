import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import os
import pdb
import numpy as np
from pandas_plink import read_plink1_bin
import pickle







genotype_stem = sys.argv[1]
hapmap3_rsid_file = sys.argv[2]




# Load in Reference Genotype data
G_obj = read_plink1_bin(genotype_stem + '.bed', genotype_stem + '.bim', genotype_stem + '.fam', verbose=False)
G = G_obj.values # Numpy 2d array of dimension num samples X num snps
ref_chrom = np.asarray(G_obj.chrom)
ref_pos = np.asarray(G_obj.pos)
# For our purposes, a0 is the effect allele
# For case of plink package, a0 is the first column in the plink bim file
ref_a0 = np.asarray(G_obj.a0)
ref_a1 = np.asarray(G_obj.a1)
# RSids
rsids = np.asarray(G_obj.snp)
# Centimorgan distances
cm = np.asarray(G_obj.cm)

# Genotype sample size
NN = G.shape[0]


cur_rs_id = 'rs2839376'

cur_snp_index = np.where(rsids==cur_rs_id)[0][0]

nearby_snp_indices = np.abs(cm - cm[cur_snp_index]) <= 1.0

nearby_snp_names = rsids[nearby_snp_indices]

r_squared = np.square(np.corrcoef(np.transpose(G[:,nearby_snp_indices])))

adj_r_squared = r_squared - ((1.0-r_squared)/(NN-2.0))

new_cur_snp_index = np.where(nearby_snp_names == cur_rs_id)[0][0]

score = np.sum(adj_r_squared[new_cur_snp_index,:])

pdb.set_trace()