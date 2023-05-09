import numpy as np 
import os
import sys
import pdb
import pickle



def extract_iterative_variant_gene_tissue_log_priors(mapping, variants, genes):
	n_var = len(variants)
	n_genes = len(genes)
	n_bs = len(mapping['variant'])
	var_probs = np.zeros((n_var, n_bs))
	gene_probs = np.zeros((n_genes, n_bs))

	for var_iter in range(n_var):
		var_probs[var_iter, :] = mapping['variant']
	for gene_iter, gene_name in enumerate(genes):
		tissue_name = gene_name.split('_')[1]
		gene_probs[gene_iter, :] = mapping[tissue_name]
	
	# Normalize rows
	normalizers = np.sum(gene_probs,axis=0) + np.sum(var_probs,axis=0)
	norm_var_probs = var_probs/normalizers
	norm_gene_probs = gene_probs/normalizers

	e_ln_pi_var = np.mean(np.log(norm_var_probs),axis=1)
	e_ln_pi_gene = np.mean(np.log(norm_gene_probs),axis=1)
	
	return e_ln_pi_var, e_ln_pi_gene







tgfm_input_summary_file = sys.argv[1]
iterative_prior_summary_file = sys.argv[2]
version = sys.argv[3]

# Create mapping from element name to bs-probs
mapping = {}
f = open(iterative_prior_summary_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	ele_name = data[0]
	bs_probs = np.asarray(data[3].split(';')).astype(float)
	mapping[ele_name] = bs_probs
f.close()



f = open(tgfm_input_summary_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	window_name = data[0]
	tgfm_input_pkl = data[2]
	ln_pi_stem = data[3]

	# Load in tgfm input data
	g = open(tgfm_input_pkl, "rb")
	tgfm_data = pickle.load(g)
	g.close()
	print(window_name)

	# extract prior for this window
	variant_log_prior, gene_log_prior = extract_iterative_variant_gene_tissue_log_priors(mapping, tgfm_data['variants'], tgfm_data['genes'])

	# print to output
	output_file = ln_pi_stem + '_iterative_variant_gene_tissue_' + version + '.txt'
	t = open(output_file,'w')
	t.write('element_name\tln_pi\n')
	for variant_iter, var_name in enumerate(tgfm_data['variants']):
		t.write(var_name + '\t' + str(variant_log_prior[variant_iter]) + '\n')
	for gene_iter, gene_name in enumerate(tgfm_data['genes']):
		t.write(gene_name + '\t' + str(gene_log_prior[gene_iter]) + '\n')
	t.close()
f.close()
