import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import sys
import pdb
import sumstats
import rpy2
import math
import rpy2.robjects.numpy2ri as numpy2ri
import rpy2.robjects as ro
import scipy.stats
ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
#colocR_pkg = importr('coloc')


def extract_gwas_pvalues(merged_gwas_summary_stat_file):
	f = open(merged_gwas_summary_stat_file)
	beta_vec = []
	beta_var_vec = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		beta = float(data[6])
		beta_var = float(data[7])

		beta_vec.append(beta)
		beta_var_vec.append(beta_var)
	f.close()
	return np.asarray(beta_vec), np.asarray(beta_var_vec)

def get_approx_log_bf_estimates(beta_vec, beta_var_vec, sd_prior=0.15):
	z = beta_vec/np.sqrt(beta_var_vec)
	r = np.square(sd_prior)/(np.square(sd_prior) + beta_var_vec)
	lABF = 0.5*(np.log(1.0 - r) + (r*np.square(z)))
	return lABF

def run_coloc(trait1_lnbfs,trait2_lnbfs, prior1=1e-4,prior2 = 1e-4, prior12= 1e-5):
	log_numerators = (
		0,
		math.log(prior1) + sumstats.log_sum(trait1_lnbfs),
		math.log(prior2) + sumstats.log_sum(trait2_lnbfs),
		math.log(prior1) + math.log(prior2) + sumstats.log_sum(
			trait1_lnbf + trait2_lnbf
			for i, trait1_lnbf in enumerate(trait1_lnbfs)
			for j, trait2_lnbf in enumerate(trait2_lnbfs)
			if i != j
		),
		math.log(prior12) + sumstats.log_sum(
			trait1_lnbf + trait2_lnbf
			for trait1_lnbf, trait2_lnbf in zip(trait1_lnbfs, trait2_lnbfs)
		)
	)

	coloc_probs = np.exp(np.asarray(log_numerators) - sumstats.log_sum(log_numerators))
	return coloc_probs


merged_gwas_summary_stat_file = sys.argv[1]
simulation_number = sys.argv[2]
chrom_num = sys.argv[3]
simulation_name_string = sys.argv[4]
n_gwas_individuals = sys.argv[5]
simulated_gene_expression_dir = sys.argv[6]
simulated_learned_gene_models_dir = sys.argv[7]
simulated_coloc_results_dir = sys.argv[8]

#################
# Extract gwas pvalues across all snps on this chromosome
gwas_beta, gwas_beta_var = extract_gwas_pvalues(merged_gwas_summary_stat_file)

##################
# Loop through gene summary file
gene_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'

# open outpout file handle
output_file = simulated_coloc_results_dir + simulation_name_string + '_coloc_results.txt'
print(output_file)
t = open(output_file,'w')
t.write('gene_name\tcoloc_pph4\n')

# Open gene summary file-handle
f = open(gene_summary_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Extract relevent fields
	ensamble_id = data[0]
	gene_snp_indices_file = data[5]
	total_n_genome_snps = int(data[6])

	gene_snp_indices_raw = np.load(gene_snp_indices_file)
	gene_snp_indices = np.asarray([False]*total_n_genome_snps)
	gene_snp_indices[gene_snp_indices_raw] = True



	# Quick error check
	if len(gwas_beta) != len(gene_snp_indices):
		print('assumption eroorr')
		pdb.set_trace()

	# Extract gene window gwas pvalues
	gene_window_gwas_beta = gwas_beta[gene_snp_indices]
	gene_window_gwas_beta_var = gwas_beta_var[gene_snp_indices]

	# Get ABF for gwas
	gene_window_gwas_lbf = get_approx_log_bf_estimates(gene_window_gwas_beta, gene_window_gwas_beta_var)

	# Load in eqtl pvalues for this gene across tissues
	eqtl_betas_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_marginal_beta.npy'
	eqtl_betas = np.load(eqtl_betas_file)
	eqtl_beta_var_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_marginal_beta_var.npy'
	eqtl_beta_vars = np.load(eqtl_beta_var_file)


	# Now loop through tissues
	for tissue_iter in range(10):
		# Name of gene-tissue pair
		gene_tissue_name = ensamble_id + '_' + 'tissue' + str(tissue_iter)

		# Now extract eqtl pvalues for this gene tissue pair
		gene_tissue_eqtl_beta = eqtl_betas[tissue_iter,:]
		gene_tissue_eqtl_beta_var = eqtl_beta_vars[tissue_iter,:]

		gene_tissue_eqtl_lbf = get_approx_log_bf_estimates(gene_tissue_eqtl_beta, gene_tissue_eqtl_beta_var)
		coloc_prob = run_coloc(gene_window_gwas_lbf, gene_tissue_eqtl_lbf)
		coloc_pph4 = coloc_prob[-1]
		t.write(gene_tissue_name + '\t' + str(coloc_pph4) + '\n')

f.close()
t.close()





