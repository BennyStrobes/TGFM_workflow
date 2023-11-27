import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import pdb
import h5py
from sklearn.decomposition import PCA
import scanpy as sc
from anndata import AnnData
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.linear_model import LinearRegression
from sklearn.cluster import KMeans
import pandas as pd
from scipy.sparse import csr_matrix
import scipy.io


def extract_dictionary_list_of_rna_seq_individuals(input_h5py_file):
	adata = sc.read_h5ad(input_h5py_file)
	indi_arr = np.asarray(adata.obs['ind_cov'])
	all_indis = {}
	for indi in indi_arr:
		all_indis[indi] = 1
	return all_indis



def extract_dictionary_list_of_genotyped_individuals(pysam_file):
	f = open(pysam_file)
	all_indis = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		all_indis[data[0]] = 1
	f.close()
	return all_indis


#######################
# Command line args
#######################
updated_individual_info_file = sys.argv[1]
individual_info_file = sys.argv[2]
genotype_data_dir = sys.argv[3]
input_h5py_file = sys.argv[4]



# Extract dictionary list of genotyped individuals
genotyped_individuals = extract_dictionary_list_of_genotyped_individuals(genotype_data_dir + 'inds_v2_header.maf10.psam')

# Extract dictionary list of rna-seq individuals
rna_seq_individuals = extract_dictionary_list_of_rna_seq_individuals(input_h5py_file)



f = open(individual_info_file)
t = open(updated_individual_info_file,'w')
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	if data[1] != data[3]:
		print('assumption eroror')
		pdb.set_trace()
	ind_id = data[1]
	if data[2] != 'European':
		continue
	if ind_id in genotyped_individuals and ind_id in rna_seq_individuals:
		t.write(line + '\n')
	else:
		# SEEMS TO BE MISSING ENTIRELY FROM GENOTYPE
		# WILL PUSH AHEAD WITH 120 individuals
		if ind_id not in rna_seq_individuals:
			pdb.set_trace()
		if data[4] in genotyped_individuals:
			pdb.set_trace()
		if data[5] in genotyped_individuals:
			pdb.set_trace()
f.close()
t.close()







'''
output_file = processed_genotype_dir + 'sc_rna_seq_ea_individual_list.txt'

# Load in processed-SC Ann-Data file
input_h5py_file = processed_sc_expression_dir + 'scran_normalization_hvg_2000_regress_batch_True_2.h5ad'
adata = sc.read_h5ad(input_h5py_file)


indiz = np.asarray(adata.obs['ind_cov'])
popz = np.asarray(adata.obs['pop_cov'])

unique_individuals_eur = {}
unique_individuals_aa = {}
unique_individuals_asian = {}
unique_individuals_hispanic = {}


for i, indi in enumerate(indiz):
	if popz[i] == 'European':
		unique_individuals_eur[indi] = 1
	elif popz[i] == 'African American':
		unique_individuals_aa[indi] = 1
	elif popz[i] == 'Asian':
		unique_individuals_asian[indi] = 1
	elif popz[i] == 'Hispanic':
		unique_individuals_hispanic[indi] = 1
	else:
		print('assumption eororor')
		pdb.set_trace()

t = open(output_file,'w')
eur_donors = np.sort(np.asarray([*unique_individuals_eur]))

for eur_donor in eur_donors:
	t.write(eur_donor + '\n')
t.close()
'''