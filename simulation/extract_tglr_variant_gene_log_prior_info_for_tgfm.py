import numpy as np 
import os
import sys
import pdb
import pickle




def extract_tglr_variant_gene_log_priors(per_variant_h2, per_gene_h2, variants, genes):
	# Deal with case where estimated per element h2 is < 0
	if per_variant_h2 <= 1e-15:
		per_variant_h2 = 1e-15
	if per_gene_h2 <= 1e-15:
		per_gene_h2 = 1e-15
	# Get n genes and n-variants
	n_genes = len(genes)
	n_variants = len(variants)

	# Fill in genetic elements h2s
	gene_h2s = np.zeros(n_genes) + per_gene_h2
	variant_h2s = np.zeros(n_variants) + per_variant_h2

	# Calculate h2 from region (sum of element h2s)
	region_h2 = np.sum(gene_h2s) + np.sum(variant_h2s)


	# Normalize to get probabilities
	gene_probs = gene_h2s/region_h2
	variant_probs = variant_h2s/region_h2	

	return np.log(variant_probs), np.log(gene_probs)



tgfm_input_summary_file = sys.argv[1]
eqtl_sample_size = sys.argv[2]
simulation_name_string = sys.argv[3]
simulated_sldsc_results_dir = sys.argv[4]
identifier_string = sys.argv[5]


# First extract taus and h2s
tglr_tau_file = simulated_sldsc_results_dir + simulation_name_string + '_eqtl_ss_' + eqtl_sample_size + '_susie_pmces_sldsc_results_organized_res.txt'
tglr_h2_file = simulated_sldsc_results_dir + simulation_name_string + '_eqtl_ss_' + eqtl_sample_size + '_susie_pmces_sldsc_results_organized_mediated_h2.txt'
tglr_tau_data = np.loadtxt(tglr_tau_file, dtype=str,delimiter='\t')
tglr_h2_data = np.loadtxt(tglr_h2_file, dtype=str, delimiter='\t')

# Extract taus and heritabilities
anno_names = tglr_tau_data[1:,0]
taus = tglr_tau_data[1:,1].astype(float)
h2s = tglr_h2_data[1:,1].astype(float)
# Extract m-vec
m_vec = h2s/taus
# Eqtl start index
eqtl_start_index = np.where(anno_names=='tissue0')[0][0]

# Get number of snps
n_snps = m_vec[0]
n_genes = np.sum(m_vec[eqtl_start_index:])

# Extract variant and gene h2s
jacknifed_geno_h2 = np.sum(h2s[:eqtl_start_index])
jacknifed_expr_h2 = np.sum(h2s[eqtl_start_index:])


# Per snp and per-gene h2 across jacknifed samples
per_snp_h2 = jacknifed_geno_h2/n_snps
per_gene_h2 = jacknifed_expr_h2/n_genes


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

	# extract prior for this window
	variant_log_prior, gene_log_prior = extract_tglr_variant_gene_log_priors(per_snp_h2, per_gene_h2, tgfm_data['variants'], tgfm_data['genes'])

	# print to output
	output_file = ln_pi_stem + '_tglr_variant_gene_' + identifier_string + '.txt'
	t = open(output_file,'w')
	t.write('element_name\tln_pi\n')
	for variant_iter, var_name in enumerate(tgfm_data['variants']):
		t.write(var_name + '\t' + str(variant_log_prior[variant_iter]) + '\n')
	for gene_iter, gene_name in enumerate(tgfm_data['genes']):
		t.write(gene_name + '\t' + str(gene_log_prior[gene_iter]) + '\n')
	t.close()
f.close()


