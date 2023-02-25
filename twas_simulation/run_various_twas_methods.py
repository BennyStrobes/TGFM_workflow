import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os
import pdb
from pandas_plink import read_plink1_bin
import scipy.stats




def extract_list_of_gwas_sum_stats(gwas_results_file):
	f = open(gwas_results_file)
	rsids = []
	z_scores = []
	beta = []
	beta_se = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsids.append(data[0])
		beta.append(float(data[1]))
		beta_se.append(float(data[2]))
		z_scores.append(float(data[3]))
	f.close()
	return np.asarray(beta), np.asarray(beta_se), np.asarray(z_scores), np.asarray(rsids)

def extract_list_of_simulated_causal_gene_effect_sizes(simulated_causal_gene_file):
	causal_effect_sizes = []
	sim_gene_names = []
	f = open(simulated_causal_gene_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		causal_effect_sizes.append(float(data[1]))
		sim_gene_names.append(data[0])
	f.close()
	return np.asarray(causal_effect_sizes), np.asarray(sim_gene_names)

def load_in_ordered_array_of_rsids_from_bim_file(genotype_bim):
	f = open(genotype_bim)
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		arr.append(data[1])
	f.close()
	arr = np.asarray(arr)

	# Quick error check 
	if len(np.unique(arr)) != len(arr):
		print('assumption eroror')
		pdb.set_trace()

	return arr


def compute_ld_on_window(gwas_plink_stem, window_rsids):
	# Load in ordered array of rsids
	genotype_bim = gwas_plink_stem + '.bim'
	ordered_rsids = load_in_ordered_array_of_rsids_from_bim_file(genotype_bim)

	# Load in variant names in plink format
	window_variant_names = []
	for snp_iter, snp_rsid in enumerate(ordered_rsids):
		# Skip snps that are not hapmap3 snps
		if snp_rsid not in window_rsids:
			continue
		window_variant_names.append('variant' + str(snp_iter))
	window_variant_names = np.asarray(window_variant_names)

	# Load in genotype object
	genotype_obj = read_plink1_bin(gwas_plink_stem + '.bed', gwas_plink_stem + '.bim', gwas_plink_stem + '.fam', verbose=False)


	# Extract genotype data for this batch of snps
	window_variant_genotype = np.asarray(genotype_obj.sel(variant=window_variant_names))

	# Standardize and mean impute genotype
	nvary = window_variant_genotype.shape[1]
	for var_iter in range(nvary):
		variant_genotype = window_variant_genotype[:,var_iter]
		std_variant_genotype = np.copy(variant_genotype)
		nan_indices = np.isnan(variant_genotype)
		non_nan_mean = np.mean(variant_genotype[nan_indices==False])
		std_variant_genotype[nan_indices] = non_nan_mean
		std_variant_genotype = (std_variant_genotype - np.mean(std_variant_genotype))/np.std(std_variant_genotype)
		window_variant_genotype[:,var_iter] = std_variant_genotype


	# Compute LD
	LD_mat = np.corrcoef(np.transpose(window_variant_genotype))

	return LD_mat, window_variant_genotype

def calculate_gene_variance_according_to_susie_distribution(susie_mu, susie_alpha, susie_mu_sd, ld):
	gene_var = 0.0

	# Component level eqtl effect sizes for this gene		
	gene_component_effect_sizes = (susie_mu)*susie_alpha

	# eQTL effect sizes for this gene
	gene_eqtl_effect_sizes = np.sum(gene_component_effect_sizes,axis=0)


	num_susie_components = susie_mu.shape[0]
	for k_index in range(num_susie_components):
		gene_var = gene_var + np.sum((np.square(susie_mu[k_index,:]) + np.square(susie_mu_sd[k_index,:]))*np.diag(ld)*susie_alpha[k_index,:])
		eqtl_component_pmces = (susie_mu[k_index,:])*(susie_alpha[k_index,:])
		gene_var = gene_var - np.dot(np.dot(eqtl_component_pmces,ld), eqtl_component_pmces)
	gene_var = gene_var + np.dot(np.dot(gene_eqtl_effect_sizes,ld), gene_eqtl_effect_sizes)
				
	return gene_var



def susie_distr_twas(ww, susie_alpha, susie_mu, susie_mu_var, gwas_gene_z_scores, ld_mat, gene_gwas_s_vec):
	gene_variance = calculate_gene_variance_according_to_susie_distribution(susie_mu, susie_alpha, np.sqrt(susie_mu_var), ld_mat)
	ww_stand = ww/np.sqrt(gene_variance)
	
	twas_z = np.dot(ww_stand, gwas_gene_z_scores)
	twas_coef_se = np.sqrt(np.square(np.mean(gene_gwas_s_vec))/1.0)
	twas_coef = twas_z*twas_coef_se
	twas_p = scipy.stats.norm.sf(abs(twas_z))*2

	variance_ratio = np.dot(np.dot(ww, ld_mat),ww)/gene_variance

	return twas_coef, twas_z, twas_p, variance_ratio


def pmces_twas(ww, zz, ld_mat, gene_gwas_s_vec):
	# Standardize predicted gene expression
	gene_variance = np.dot(np.dot(ww,ld_mat),ww)
	ww_stand = ww/np.sqrt(gene_variance)

	# Compute standardard error of twas effect size
	#pdb.set_trace()
	twas_coef_se = np.sqrt(1.0/np.dot(np.dot(ww_stand/gene_gwas_s_vec, ld_mat), ww_stand/gene_gwas_s_vec))


	#np.sqrt(np.square(np.mean(gene_gwas_s_vec))/np.dot(np.dot(ww_stand, ld_mat), ww_stand))
	#np.sqrt(1.0/(np.square(np.mean(1.0/gene_gwas_s_vec))*np.dot(np.dot(ww_stand, ld_mat), ww_stand)))
	#S_inv_mat = np.diag(1.0/gene_gwas_s_vec)
	#D_mat = np.multiply(np.multiply(np.diag(S_inv_mat)[:, None], ld_mat), np.diag(S_inv_mat))
	twas_z = np.dot(ww_stand, zz)
	twas_coef = twas_z*twas_coef_se
	twas_p = scipy.stats.norm.sf(abs(twas_z))*2
	return twas_coef, twas_z, twas_p


def lava_style_pmces_twas(marginal_eqtl_beta, marginal_eqtl_beta_se, eqtl_sample_size, pruned_svd_Q, pruned_svd_lambda, ld_mat, gwas_gene_z_scores, gene_gwas_s_vec, sample_size_fraction=.75):
	# Potentially filter more pcs
	n_pcs = len(pruned_svd_lambda)
	max_num_pcs = int(np.floor(eqtl_sample_size*sample_size_fraction))
	if n_pcs >= max_num_pcs:
		pc_filter = np.asarray([False]*n_pcs)
		pc_filter[:max_num_pcs] = True
		pruned_svd_lambda = pruned_svd_lambda[pc_filter]
		pruned_svd_Q = pruned_svd_Q[:, pc_filter]

	# Compute predicted causal eqtl effects (alpha)
	S_inv_mat = np.dot(np.dot(pruned_svd_Q, np.diag(1.0/pruned_svd_lambda)), np.transpose(pruned_svd_Q))
	alpha = np.dot(S_inv_mat, marginal_eqtl_beta)

	# Compute sparse-pc predicted causal effects (delta)
	delta = np.dot(np.dot(np.diag(np.sqrt(pruned_svd_lambda)), np.transpose(pruned_svd_Q)), alpha)
	
	# compute variances
	r_sq = np.dot(alpha, marginal_eqtl_beta)
	sigma = (1.0-r_sq)/(eqtl_sample_size - n_pcs - 1.0)

	if sigma < 0.0:
		sigma = 0.0


	# Observed h2
	obs_h2 = 1.0 - (1.0 - np.dot(delta, delta))*((eqtl_sample_size-1.0)/(eqtl_sample_size-n_pcs-1.0))

	# Partition genetic gene variance
	genetic_gene_pmces_var = np.dot(alpha, marginal_eqtl_beta) # This is equivelent to running alpha*S*alpha 
	genetic_gene_noise_var = sigma*n_pcs

	unbiased_genetic_gene_var_est = genetic_gene_pmces_var - genetic_gene_noise_var
	if unbiased_genetic_gene_var_est < .01:
		unbiased_genetic_gene_var_est = .01

	expected_total_genetic_var = genetic_gene_pmces_var + genetic_gene_noise_var



	# RUN TWAS
	alpha_stand = alpha/np.sqrt(genetic_gene_pmces_var)
	expected_total_genetic_var_stand = genetic_gene_pmces_var/genetic_gene_pmces_var
	
	twas_z = np.dot(alpha_stand, gwas_gene_z_scores)/np.sqrt(expected_total_genetic_var_stand)
	twas_coef_se = np.sqrt(np.square(np.mean(gene_gwas_s_vec))/expected_total_genetic_var_stand)
	twas_coef = twas_z*twas_coef_se
	twas_p = scipy.stats.norm.sf(abs(twas_z))*2

	return twas_coef, twas_z, twas_p



def lava_style_distribution_twas(marginal_eqtl_beta, marginal_eqtl_beta_se, eqtl_sample_size, pruned_svd_Q, pruned_svd_lambda, ld_mat, gwas_gene_z_scores, gene_gwas_s_vec, sample_size_fraction=.75):
	# Potentially filter more pcs
	n_pcs = len(pruned_svd_lambda)
	max_num_pcs = int(np.floor(eqtl_sample_size*sample_size_fraction))
	if n_pcs >= max_num_pcs:
		pc_filter = np.asarray([False]*n_pcs)
		pc_filter[:max_num_pcs] = True
		pruned_svd_lambda = pruned_svd_lambda[pc_filter]
		pruned_svd_Q = pruned_svd_Q[:, pc_filter]

	# Compute predicted causal eqtl effects (alpha)
	S_inv_mat = np.dot(np.dot(pruned_svd_Q, np.diag(1.0/pruned_svd_lambda)), np.transpose(pruned_svd_Q))
	alpha = np.dot(S_inv_mat, marginal_eqtl_beta)

	# Compute sparse-pc predicted causal effects (delta)
	delta = np.dot(np.dot(np.diag(np.sqrt(pruned_svd_lambda)), np.transpose(pruned_svd_Q)), alpha)
	
	# compute variances
	r_sq = np.dot(alpha, marginal_eqtl_beta)
	sigma = (1.0-r_sq)/(eqtl_sample_size - n_pcs - 1.0)

	if sigma < 0.0:
		sigma = 0.0

	# Observed h2
	obs_h2 = 1.0 - (1.0 - np.dot(delta, delta))*((eqtl_sample_size-1.0)/(eqtl_sample_size-n_pcs-1.0))

	# Partition genetic gene variance
	genetic_gene_pmces_var = np.dot(alpha, marginal_eqtl_beta) # This is equivelent to running alpha*S*alpha 
	genetic_gene_noise_var = sigma*n_pcs

	unbiased_genetic_gene_var_est = genetic_gene_pmces_var - genetic_gene_noise_var
	if unbiased_genetic_gene_var_est < .01:
		unbiased_genetic_gene_var_est = .01

	expected_total_genetic_var = genetic_gene_pmces_var + genetic_gene_noise_var

	# RUN TWAS
	att_lambda = unbiased_genetic_gene_var_est/(genetic_gene_pmces_var+genetic_gene_noise_var)
	print(att_lambda)
	twas_z = np.dot(alpha, gwas_gene_z_scores)/np.sqrt(expected_total_genetic_var)
	twas_coef_se = np.sqrt(np.square(np.mean(gene_gwas_s_vec))/((expected_total_genetic_var/unbiased_genetic_gene_var_est)))/att_lambda
	twas_coef = twas_z*twas_coef_se
	twas_p = scipy.stats.norm.sf(abs(twas_z))*2

	# Attenuation bias correction
	variance_ratio = genetic_gene_pmces_var/expected_total_genetic_var

	return twas_coef, twas_z, twas_p, variance_ratio



def run_twas_shell(gene_name, simulation_name_string, twas_method_name, eqtl_ss, ld_mat, gwas_gene_z_scores, gene_sim_eqtl_effects_file, simulated_learned_gene_models_dir, gene_gwas_s_vec, pruned_svd_Q, pruned_svd_lambda):
	if twas_method_name == 'susie_pmces':
		gene_eqtl_pmces_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_susie_pmces_gene_model.npy'
		gene_eqtl_pmces = np.load(gene_eqtl_pmces_file)[0,:]
		if np.var(gene_eqtl_pmces) == 0.0:
			twas_coef = 'nan'
			twas_z = 'nan'
			twas_p = 'nan'
			variance_ratio = 'nan'
		else:
			twas_coef, twas_z, twas_p = pmces_twas(gene_eqtl_pmces, gwas_gene_z_scores, ld_mat, gene_gwas_s_vec)
			variance_ratio = 1.0
	elif twas_method_name == 'fusion_lasso_pmces':
		gene_eqtl_pmces_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_fusion_lasso_pmces_gene_model.npy'
		gene_eqtl_pmces = np.load(gene_eqtl_pmces_file)[0,:]
		if np.var(gene_eqtl_pmces) == 0.0:
			twas_coef = 'nan'
			twas_z = 'nan'
			twas_p = 'nan'
			variance_ratio = 'nan'
		else:
			twas_coef, twas_z, twas_p = pmces_twas(gene_eqtl_pmces, gwas_gene_z_scores, ld_mat, gene_gwas_s_vec)
			variance_ratio = 1.0

	elif twas_method_name == 'susie_distr':
		gene_eqtl_pmces_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_susie_pmces_gene_model.npy'
		gene_eqtl_pmces = np.load(gene_eqtl_pmces_file)[0,:]
		if np.var(gene_eqtl_pmces) == 0.0:
			twas_coef = 'nan'
			twas_z = 'nan'
			twas_p = 'nan'
			variance_ratio = 'nan'
		else:
			susie_alpha_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_susie_alpha_gene_model.npy'
			susie_alpha = np.load(susie_alpha_file)
			susie_mu_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_susie_mu_gene_model.npy'
			susie_mu = np.load(susie_mu_file)
			susie_mu_var_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_susie_mu_var_gene_model.npy'
			susie_mu_var = np.load(susie_mu_var_file)
			twas_coef, twas_z, twas_p, variance_ratio = susie_distr_twas(gene_eqtl_pmces, susie_alpha, susie_mu, susie_mu_var, gwas_gene_z_scores, ld_mat, gene_gwas_s_vec)
	elif twas_method_name == 'true_causal_effects':
		gene_eqtl_pmces = np.load(gene_sim_eqtl_effects_file)[:,0]
		twas_coef, twas_z, twas_p = pmces_twas(gene_eqtl_pmces, gwas_gene_z_scores, ld_mat, gene_gwas_s_vec)
		variance_ratio = 1.0
	elif twas_method_name == 'marginal_distr':
		gene_eqtl_pmces = np.load(gene_sim_eqtl_effects_file)[:,0]
		# Load in data
		marginal_eqtl_beta_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_marginal_effects_gene_model.npy'
		marginal_eqtl_beta = np.load(marginal_eqtl_beta_file)[0,:]
		marginal_eqtl_beta_se_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_marginal_effects_se_gene_model.npy'
		marginal_eqtl_beta_se = np.load(marginal_eqtl_beta_se_file)[0,:]
		twas_coef, twas_z, twas_p, variance_ratio = lava_style_distribution_twas(marginal_eqtl_beta, marginal_eqtl_beta_se, int(eqtl_ss), pruned_svd_Q, pruned_svd_lambda, ld_mat, gwas_gene_z_scores, gene_gwas_s_vec)
		print(twas_coef)
	elif twas_method_name == 'marginal_pmces':
		# Load in data
		marginal_eqtl_beta_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_marginal_effects_gene_model.npy'
		marginal_eqtl_beta = np.load(marginal_eqtl_beta_file)[0,:]
		marginal_eqtl_beta_se_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_marginal_effects_se_gene_model.npy'
		marginal_eqtl_beta_se = np.load(marginal_eqtl_beta_se_file)[0,:]
		twas_coef, twas_z, twas_p = lava_style_pmces_twas(marginal_eqtl_beta, marginal_eqtl_beta_se, int(eqtl_ss), pruned_svd_Q, pruned_svd_lambda, ld_mat, gwas_gene_z_scores, gene_gwas_s_vec)		
		variance_ratio = 1.0
	else:
		print('twas method ' + twas_method_name + ' not currently implemented')
		pdb.set_trace()
	return twas_coef, twas_z, twas_p, variance_ratio


def extract_indices_of_pcs_that_explain_specified_fraction_of_variance(svd_lambda, fraction):
	if np.sum(svd_lambda < 0.0) != 0:
		print('negative lambda assumption error')
		pdb.set_trace()
	pve = svd_lambda/np.sum(svd_lambda)
	cum_sum = []
	cum = 0.0
	for pve_val in pve:
		cum = cum + pve_val
		cum_sum.append(cum)
	cum_sum = np.asarray(cum_sum)

	valid_pcs = cum_sum <= fraction
	if np.sum(valid_pcs) < 2:
		print('assumption eororor')
		pdb.set_trace()

	return valid_pcs



###########################
# Command line args
##########################
simulation_number = sys.argv[1]
chrom_num = sys.argv[2]
simulation_name_string = sys.argv[3]
processed_genotype_data_dir = sys.argv[4]
simulated_gwas_dir = sys.argv[5]
simulated_gene_expression_dir = sys.argv[6]
simulated_learned_gene_models_dir = sys.argv[7]
simulated_trait_dir = sys.argv[8]
simulated_twas_dir = sys.argv[9]


# Get list of methods and eqtl sample sizes
method_names = []
eqtl_sample_sizes = []
#for eqtl_ss in [100,200,300,500,1000]:
for eqtl_ss in [100]:	
	for method in ['marginal_distr']:
	#for method in ['susie_pmces', 'susie_distr', 'marginal_pmces', 'marginal_distr', 'fusion_lasso_pmces']:
		method_names.append(method)
		eqtl_sample_sizes.append(eqtl_ss)
method_names.append('true_causal_effects')
eqtl_sample_sizes.append('inf')
method_names = np.asarray(method_names)
eqtl_sample_sizes = np.asarray(eqtl_sample_sizes)


# Extract list of gene causal effect sizes
simulated_causal_gene_file = simulated_trait_dir + simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
simulated_causal_gene_effect_sizes, sim_gene_names = extract_list_of_simulated_causal_gene_effect_sizes(simulated_causal_gene_file)



# Open output file handle and print header
twas_output_file = simulated_twas_dir + simulation_name_string + '_simualated_twas_results.txt'
t = open(twas_output_file,'w')
t.write('gene_name\teqtl_ss\ttwas_method\tcausal_gene\tsimulated_gene_effect_size\ttwas_effect_size\ttwas_z_score\ttwas_p_value\tvariance_ratio\n')

# Loop through genes
gene_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'
f = open(gene_summary_file)
head_count = 0
gene_counter = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue

	# Parse line
	gene_name = data[0]

	print(gene_name)
	sim_gene_name = sim_gene_names[gene_counter]
	gene_causal_effect_size = simulated_causal_gene_effect_sizes[gene_counter]
	print(gene_causal_effect_size)

	causal_gene = 'True'
	if gene_causal_effect_size == 0.0:
		causal_gene = 'False'
	cis_snp_id_file = data[4]
	cis_snp_index_file = data[5]
	gene_sim_eqtl_effects_file = data[3]
	gene_counter = gene_counter + 1

	# Quick error check
	if gene_name != sim_gene_name:
		print('assumption error')
		pdb.set_trace()
	if np.abs(gene_causal_effect_size) == 0.0:
		continue

	# Extract data for this line
	cis_snp_rsids = np.load(cis_snp_id_file, allow_pickle=True)[:,1]
	cis_snp_indices = np.load(cis_snp_index_file)
	n_cis_snps = np.sum(cis_snp_indices)


	# Extract list of gwas z-scores
	gene_gwas_results_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results_' + gene_name + '.txt'
	gene_gwas_beta, gene_gwas_beta_se, gwas_gene_z_scores, sumstat_rsids = extract_list_of_gwas_sum_stats(gene_gwas_results_file)
	gene_gwas_s_squared_vec = np.square(gene_gwas_beta_se) + (np.square(gene_gwas_beta)/100000.0)
	gwas_gene_s_vec = np.sqrt(gene_gwas_s_squared_vec)

	# QUick error check
	if np.array_equal(sumstat_rsids, cis_snp_rsids) == False:
		print('assumption eroroor')
		pdb.set_trace()

	# Compute gene LD 
	gwas_plink_stem = processed_genotype_data_dir + 'simulated_gwas_data_' + str(chrom_num)  # Genotype files
	ld_mat, window_variant_genotype = compute_ld_on_window(gwas_plink_stem, cis_snp_rsids)

	# Reduce dimensionality of window variant genotype with svd
	n_ref = window_variant_genotype.shape[0]
	uu, ss, vh = np.linalg.svd(window_variant_genotype, full_matrices=False)
	svd_lambda = ss*ss/(n_ref-1)
	svd_Q = np.transpose(vh)
	# Prune to PCs that explain 99% of variance
	pc_indices = extract_indices_of_pcs_that_explain_specified_fraction_of_variance(svd_lambda, .99)
	pruned_svd_lambda = svd_lambda[pc_indices]
	pruned_svd_Q = svd_Q[:, pc_indices]

	# Loop through methods and 
	for method_iter in range(len(method_names)):
		# Name of method and eqtl sample size
		twas_method_name = method_names[method_iter]
		eqtl_ss = eqtl_sample_sizes[method_iter]

		# Run twas
		twas_coef,twas_z, twas_p, variance_ratio = run_twas_shell(gene_name, simulation_name_string, twas_method_name, eqtl_ss, ld_mat, gwas_gene_z_scores, gene_sim_eqtl_effects_file, simulated_learned_gene_models_dir, gwas_gene_s_vec, pruned_svd_Q, pruned_svd_lambda)
		t.write(gene_name + '\t' + str(eqtl_ss) + '\t' + twas_method_name + '\t' + causal_gene + '\t' + str(gene_causal_effect_size) + '\t' + str(twas_coef) + '\t' + str(twas_z) + '\t' + str(twas_p) + '\t' + str(variance_ratio) + '\n')
	t.flush()
f.close()
t.close()


