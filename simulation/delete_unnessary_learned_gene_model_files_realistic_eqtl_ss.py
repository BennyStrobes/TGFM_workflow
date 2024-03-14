import numpy as np
import os
import sys
import pdb











simulated_gene_expression_dir = sys.argv[1]
simulated_learned_gene_models_dir = sys.argv[2]
simulation_name_string = sys.argv[3]

gene_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'


n_tiss = 10

# Loop through genes (note: not gene tissue pairs)
f = open(gene_summary_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split()
	# Skip header
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Extract relevent fields from line
	ensamble_id = data[0]
	gene_tss = int(data[2])
	cis_snp_indices_file = data[5]
	total_n_genome_snps = int(data[6])

	fitted_gene_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_gene_model_pmces.npy'
	gene_model_mat = np.load(fitted_gene_file)
				
	# Relevent info
	n_tiss = gene_model_mat.shape[0]
	n_cis_snps = gene_model_mat.shape[1]

	# Remove marginal beta and beta var files
	marginal_beta_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_marginal_beta.npy'
	marginal_betavar_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_marginal_beta_var.npy'
	os.system('rm ' + marginal_beta_file)
	os.system('rm ' + marginal_betavar_file)

	# Loop through tissues
	for tiss_iter in range(n_tiss):
		susie_stem = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_tissue_' + str(tiss_iter) + '_gene_model_susie_'
		susie_alpha_file = susie_stem + 'alpha.npy'
		susie_mu_file = susie_stem + 'mu.npy'
		susie_mu_var_file = susie_stem + 'mu_var.npy'

		# Delete files
		os.system('rm ' + susie_alpha_file)
		os.system('rm ' + susie_mu_file)
		os.system('rm ' + susie_mu_var_file)


f.close()