import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import pickle






def extract_tissue_names(gtex_pseudotissue_file, remove_testis=False):
	f = open(gtex_pseudotissue_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0] == 'Testis' and remove_testis:
			continue
		arr.append(data[0])
	f.close()
	return np.asarray(arr)








#####################
# Command line args
#####################
tgfm_input_summary_file = sys.argv[1]
gtex_pseudotissue_file = sys.argv[2]
output_file = sys.argv[3]


# Extract ordered tissue information
ordered_tissue_names = extract_tissue_names(gtex_pseudotissue_file, remove_testis=True)

# Create mapping to keep track of correlations
mapping = {}
num_tissues = len(ordered_tissue_names)
for ii in range(num_tissues):
	for jj in range(num_tissues):
		if ii != jj:
			mapping[ordered_tissue_names[ii] + ':' + ordered_tissue_names[jj]] = []





# Now loop through windows
# In each window run TGFM independently
# Loop through trait components
tgfm_input_data = np.loadtxt(tgfm_input_summary_file,dtype=str,delimiter='\t')
tgfm_input_data = tgfm_input_data[1:,:]


# Get n_windows on this run
n_windows = tgfm_input_data.shape[0]

for window_iter in range(n_windows):
	data = tgfm_input_data[window_iter, :]

	##############################
	# Extract relevent fields
	###############################
	window_name = data[0]
	print(window_name)

	ld_file = data[1]
	tgfm_input_pkl = data[2]
	tgfm_trait_input_pkl = data[3]

	##############################
	# Load in Data
	###############################
	# Load in tgfm trait input data
	g = open(tgfm_trait_input_pkl, "rb")
	tgfm_trait_data = pickle.load(g)
	g.close()


	# Load in LD
	ld_mat = np.load(ld_file)
	# Load in tgfm input data
	g = open(tgfm_input_pkl, "rb")
	tgfm_data = pickle.load(g)
	g.close()
	# Add ld to tgfm_data obj
	tgfm_data['reference_ld'] = ld_mat

	# Skip windows with no genes
	if len(tgfm_data['genes']) <= 1:
		continue

	# PMCES for this window limiting to genes in the middle of the window
	middle_gene_eqtl_pmces = tgfm_data['gene_eqtl_pmces'][tgfm_data['middle_gene_indices'], :]
	# Gene names limiting to genes in the middle of the window
	middle_gene_names = tgfm_data['genes'][tgfm_data['middle_gene_indices']]

	# tgfm_data['gene_eqtl_pmces'], tgfm_data['reference_ld']
	gene_gene_ld = np.dot(np.dot(middle_gene_eqtl_pmces, ld_mat), np.transpose(middle_gene_eqtl_pmces))

	n_genes = len(middle_gene_names)
	for ii in range(n_genes):
		for jj in range(n_genes):
			if ii != jj:
				gene_ii = middle_gene_names[ii].split('_')[0]
				gene_jj = middle_gene_names[jj].split('_')[0]
				if gene_ii == gene_jj:
					tissue_ii = '_'.join(middle_gene_names[ii].split('_')[1:])
					tissue_jj = '_'.join(middle_gene_names[jj].split('_')[1:])
					corry = gene_gene_ld[ii,jj]
					mapping[tissue_ii + ':' + tissue_jj].append(corry)


t = open(output_file,'w')
t.write('tissue_ii\ttissue_jj\taverage_correlation\tn_genes\n')
num_tissues = len(ordered_tissue_names)
for ii in range(num_tissues):
	for jj in range(num_tissues):
		if ii == jj:
			t.write(ordered_tissue_names[ii] + '\t' + ordered_tissue_names[jj] + '\t' + '1.0' + '\t' + 'NA' + '\n')
		if ii != jj:
			arr = mapping[ordered_tissue_names[ii] + ':' + ordered_tissue_names[jj]]
			if len(arr) > 0:
				t.write(ordered_tissue_names[ii] + '\t' + ordered_tissue_names[jj] + '\t' + str(np.mean(arr)) + '\t' + str(len(arr)) + '\n')
			else:
				t.write(ordered_tissue_names[ii] + '\t' + ordered_tissue_names[jj] + '\t' + 'NA' + '\t' + '0' + '\n')
t.close()

print(output_file)
















