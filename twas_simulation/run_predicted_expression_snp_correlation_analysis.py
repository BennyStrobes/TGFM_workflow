import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os
import pdb
from pandas_plink import read_plink1_bin
import scipy.stats









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

	# Quick error check
	if np.sum(np.isnan(window_variant_genotype)) != 0.0:
		print('assumption eroror')
		pdb.set_trace()

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

def pmces_correlation(gene_eqtl_pmces, ld_mat):
	# Standardize simulated eqtl effects
	stand_eqtl_effects = gene_eqtl_pmces/np.sqrt(np.dot(np.dot(gene_eqtl_pmces, ld_mat), gene_eqtl_pmces))

	# Extract simulated correlations
	corrz = np.dot(stand_eqtl_effects, ld_mat)

	return corrz

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




def susie_distr_correlation(ww, susie_alpha, susie_mu, susie_mu_var, ld_mat):
	gene_variance = calculate_gene_variance_according_to_susie_distribution(susie_mu, susie_alpha, np.sqrt(susie_mu_var), ld_mat)
	ww_stand = ww/np.sqrt(gene_variance)

	variance_ratio = np.dot(np.dot(ww, ld_mat),ww)/gene_variance

	corrz = np.dot(ww_stand, ld_mat)
	return corrz, variance_ratio

def unbiased_r_squared(marginal_effects, bootstrapped_marginal_effects, pruned_svd_ld_inv_mat, pruned_svd_lambda, pruned_svd_Q, eqtl_sample_size):
	n_pcs = len(pruned_svd_lambda)
	n_bootstraps = bootstrapped_marginal_effects.shape[0]

	alpha_mat = np.dot(pruned_svd_ld_inv_mat, np.transpose(bootstrapped_marginal_effects))

	bootstrapped_gene_variances = np.sum(alpha_mat*np.transpose(bootstrapped_marginal_effects),axis=0)

	bootstrapped_standardized_r_squared = np.transpose(np.square(bootstrapped_marginal_effects))/bootstrapped_gene_variances
	expected_bootstrapped_standardized_r_squared = np.mean(bootstrapped_standardized_r_squared,axis=1)

	#expected_pmces_standardized_r_squared = np.square(np.mean(bootstrapped_marginal_effects,axis=0))/np.dot(np.mean(alpha_mat,axis=1), np.mean(bootstrapped_marginal_effects,axis=0))
	expected_pmces_standardized_r_squared = np.square(marginal_effects)/np.dot(np.dot(pruned_svd_ld_inv_mat, np.transpose(marginal_effects)), marginal_effects)

	noise_terms = expected_bootstrapped_standardized_r_squared - expected_pmces_standardized_r_squared

	unbiased_standardized_r_squared = expected_pmces_standardized_r_squared - noise_terms

	return unbiased_standardized_r_squared


def lava_style_distribution_correlation(marginal_eqtl_beta, marginal_eqtl_beta_se, eqtl_sample_size, pruned_svd_Q, pruned_svd_lambda, ld_mat, bootstrapped_marginal_effects, sample_size_fraction=.75):
	# Potentially filter more pcs
	n_pcs = len(pruned_svd_lambda)
	max_num_pcs = int(np.floor(eqtl_sample_size*sample_size_fraction))
	if n_pcs >= max_num_pcs:
		pc_filter = np.asarray([False]*n_pcs)
		pc_filter[:max_num_pcs] = True
		pruned_svd_lambda = pruned_svd_lambda[pc_filter]
		pruned_svd_Q = pruned_svd_Q[:, pc_filter]
		n_pcs = len(pruned_svd_lambda)

	# Compute predicted causal eqtl effects (alpha)
	S_inv_mat = np.dot(np.dot(pruned_svd_Q, np.diag(1.0/pruned_svd_lambda)), np.transpose(pruned_svd_Q))
	alpha = np.dot(S_inv_mat, marginal_eqtl_beta)

	unbiased_r_squared_est = unbiased_r_squared(marginal_eqtl_beta, bootstrapped_marginal_effects, S_inv_mat, pruned_svd_lambda, pruned_svd_Q, eqtl_sample_size)

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

	att_lambda = unbiased_genetic_gene_var_est/(genetic_gene_pmces_var+genetic_gene_noise_var)

	# Compute variance of predictor
	standardized_predictor_variance = expected_total_genetic_var/unbiased_genetic_gene_var_est
	# Now correct variance for attenuation bias
	desired_variance = standardized_predictor_variance*np.square(att_lambda)

	# Scale alpha
	#scaling_ratio = desired_variance/expected_total_genetic_var
	#scaling_ratio = desired_variance/genetic_gene_pmces_var
	alpha_stand = alpha/np.sqrt(unbiased_genetic_gene_var_est)

	# Compute correlation
	corrz = np.dot(alpha_stand, ld_mat)

	# Attenuation bias correction
	variance_ratio = genetic_gene_pmces_var/expected_total_genetic_var
	return corrz, variance_ratio, unbiased_genetic_gene_var_est, unbiased_r_squared_est

def lava_style_pmces_correlation(marginal_eqtl_beta, marginal_eqtl_beta_se, eqtl_sample_size, pruned_svd_Q, pruned_svd_lambda, ld_mat, sample_size_fraction=.75):
	# Potentially filter more pcs
	n_pcs = len(pruned_svd_lambda)
	max_num_pcs = int(np.floor(eqtl_sample_size*sample_size_fraction))
	if n_pcs >= max_num_pcs:
		pc_filter = np.asarray([False]*n_pcs)
		pc_filter[:max_num_pcs] = True
		pruned_svd_lambda = pruned_svd_lambda[pc_filter]
		pruned_svd_Q = pruned_svd_Q[:, pc_filter]
		n_pcs = len(pruned_svd_lambda)

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

	corrz = np.dot(alpha_stand, ld_mat)

	return corrz


def run_correlation_shell(gene_name, simulation_name_string, method_name, eqtl_ss, ld_mat, simulated_learned_gene_models_dir, pruned_svd_Q, pruned_svd_lambda):
	n_snps = ld_mat.shape[0]
	if method_name == 'susie_pmces':
		gene_eqtl_pmces_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_susie_pmces_gene_model.npy'
		gene_eqtl_pmces = np.load(gene_eqtl_pmces_file)[0,:]
		if np.var(gene_eqtl_pmces) == 0.0:
			corrz = np.asarray(['nan']*n_snps)
			corrz_squared = np.asarray(['nan']*n_snps)
			variance_ratio = 'nan'
		else:
			corrz = pmces_correlation(gene_eqtl_pmces, ld_mat)
			corrz_squared = np.square(corrz)
			variance_ratio = 1.0
	elif method_name == 'fusion_lasso_pmces':
		gene_eqtl_pmces_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_fusion_lasso_pmces_gene_model.npy'
		gene_eqtl_pmces = np.load(gene_eqtl_pmces_file)[0,:]
		if np.var(gene_eqtl_pmces) == 0.0:
			corrz = np.asarray(['nan']*n_snps)
			corrz_squared = np.asarray(['nan']*n_snps)
			variance_ratio = 'nan'
		else:
			corrz = pmces_correlation(gene_eqtl_pmces, ld_mat)
			corrz_squared = np.square(corrz)
			variance_ratio = 1.0
	elif method_name == 'susie_distr':
		gene_eqtl_pmces_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_susie_pmces_gene_model.npy'
		gene_eqtl_pmces = np.load(gene_eqtl_pmces_file)[0,:]
		if np.var(gene_eqtl_pmces) == 0.0:
			corrz = np.asarray(['nan']*n_snps)
			corrz_squared = np.asarray(['nan']*n_snps)
			variance_ratio = 'nan'
		else:
			susie_alpha_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_susie_alpha_gene_model.npy'
			susie_alpha = np.load(susie_alpha_file)
			susie_mu_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_susie_mu_gene_model.npy'
			susie_mu = np.load(susie_mu_file)
			susie_mu_var_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_susie_mu_var_gene_model.npy'
			susie_mu_var = np.load(susie_mu_var_file)
			corrz, variance_ratio = susie_distr_correlation(gene_eqtl_pmces, susie_alpha, susie_mu, susie_mu_var, ld_mat)
			corrz_squared = np.square(corrz)
	elif method_name == 'marginal_distr':
		# Load in data
		marginal_eqtl_beta_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_marginal_effects_gene_model.npy'
		marginal_eqtl_beta = np.load(marginal_eqtl_beta_file)[0,:]
		marginal_eqtl_beta_se_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_marginal_effects_se_gene_model.npy'
		marginal_eqtl_beta_se = np.load(marginal_eqtl_beta_se_file)[0,:]
		bootstrapped_marginal_effects_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_marginal_effects_bootstrapped.npy'
		bootstrapped_marginal_effects = np.load(bootstrapped_marginal_effects_file)
		#unbiased_r_squared_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_unbiased_standardized_r_squared.npy'
		#unbiased_r_squared = np.load(unbiased_r_squared_file)
		corrz, variance_ratio, unbiased_genetic_gene_var_est, corrz_squared = lava_style_distribution_correlation(marginal_eqtl_beta, marginal_eqtl_beta_se, int(eqtl_ss), pruned_svd_Q, pruned_svd_lambda, ld_mat, bootstrapped_marginal_effects)
	elif method_name == 'marginal_eqtl_distr':
		# Load in data
		marginal_eqtl_beta_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_marginal_effects_gene_model.npy'
		marginal_eqtl_beta = np.load(marginal_eqtl_beta_file)[0,:]
		marginal_eqtl_beta_se_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_marginal_effects_se_gene_model.npy'
		marginal_eqtl_beta_se = np.load(marginal_eqtl_beta_se_file)[0,:]
		bootstrapped_marginal_effects_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_marginal_effects_bootstrapped.npy'
		bootstrapped_marginal_effects = np.load(bootstrapped_marginal_effects_file)
		unbiased_r_squared_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_unbiased_standardized_r_squared.npy'
		unbiased_r_squared = np.load(unbiased_r_squared_file)
		corrz, variance_ratio, unbiased_genetic_gene_var_est, corrz_squared = lava_style_distribution_correlation(marginal_eqtl_beta, marginal_eqtl_beta_se, int(eqtl_ss), pruned_svd_Q, pruned_svd_lambda, ld_mat, bootstrapped_marginal_effects)
		corrz_squared = unbiased_r_squared[0,:]
	elif method_name == 'marginal_pmces':
		# Load in data
		marginal_eqtl_beta_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_marginal_effects_gene_model.npy'
		marginal_eqtl_beta = np.load(marginal_eqtl_beta_file)[0,:]
		marginal_eqtl_beta_se_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_ss) + '_marginal_effects_se_gene_model.npy'
		marginal_eqtl_beta_se = np.load(marginal_eqtl_beta_se_file)[0,:]
		corrz = lava_style_pmces_correlation(marginal_eqtl_beta, marginal_eqtl_beta_se, int(eqtl_ss), pruned_svd_Q, pruned_svd_lambda, ld_mat)
		variance_ratio = 1.0
		corrz_squared = np.square(corrz)
	else:
		print('twas method ' + method_name + ' not currently implemented')
		pdb.set_trace()
	return corrz, corrz_squared, variance_ratio

######################
# Command line args
######################
simulation_number = sys.argv[1]
chrom_num = sys.argv[2]
simulation_name_string = sys.argv[3]
processed_genotype_data_dir = sys.argv[4]
simulated_gwas_dir = sys.argv[5]
simulated_gene_expression_dir = sys.argv[6]
simulated_learned_gene_models_dir = sys.argv[7]
simulated_expr_snp_corr_dir = sys.argv[8]



# Get list of methods and eqtl sample sizes
method_names = []
eqtl_sample_sizes = []
for eqtl_ss in [100,300,500,1000]:
	for method in ['susie_pmces', 'susie_distr', 'marginal_pmces', 'marginal_distr', 'marginal_eqtl_distr', 'fusion_lasso_pmces']:
		method_names.append(method)
		eqtl_sample_sizes.append(eqtl_ss)
#method_names.append('true_causal_effects')
#eqtl_sample_sizes.append('inf')
method_names = np.asarray(method_names)
eqtl_sample_sizes = np.asarray(eqtl_sample_sizes)




# Open output file handle and print header
twas_output_file = simulated_expr_snp_corr_dir + simulation_name_string + '_simualated_expr_snp_corr_results.txt'
t = open(twas_output_file,'w')
t.write('gene_name\tsnp_name\teqtl_ss\teqtl_method\tknown_correlation\testimated_correlation\tknow_squared_correlation\testimated_squared_correlation\tvariance_ratio\n')

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

	# Extract info from this line
	cis_snp_id_file = data[4]
	cis_snp_index_file = data[5]
	gene_sim_eqtl_effects_file = data[3]
	gene_counter = gene_counter + 1

	# Extract data for this line
	cis_snp_rsids = np.load(cis_snp_id_file, allow_pickle=True)[:,1]
	cis_snp_indices = np.load(cis_snp_index_file)
	n_cis_snps = np.sum(cis_snp_indices)

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

	# Extract simulated eQTL effects
	sim_eqtl_effects = np.load(gene_sim_eqtl_effects_file)[:,0]
	# Standardize simulated eqtl effects
	sim_stand_eqtl_effects = sim_eqtl_effects/np.sqrt(np.dot(np.dot(sim_eqtl_effects, ld_mat), sim_eqtl_effects))

	# Extract simulated correlations
	simulated_corrz = np.dot(sim_stand_eqtl_effects, ld_mat)

	# Loop through methods and 
	for method_iter in range(len(method_names)):
		# Name of method and eqtl sample size
		method_name = method_names[method_iter]
		eqtl_ss = eqtl_sample_sizes[method_iter]

		# Run twas
		method_corrz, method_corr_squared, variance_ratio = run_correlation_shell(gene_name, simulation_name_string, method_name, eqtl_ss, ld_mat, simulated_learned_gene_models_dir, pruned_svd_Q, pruned_svd_lambda)

		# Quick error check
		if len(method_corrz) != len(cis_snp_rsids) or len(cis_snp_rsids) != len(simulated_corrz):
			print('assumption eroror')
			pdb.set_trace()

		# Print results
		for snp_iter in range(len(cis_snp_rsids)):
			snp_rsid = cis_snp_rsids[snp_iter]
			simulated_corr = simulated_corrz[snp_iter]
			method_corr = method_corrz[snp_iter]
			simulated_squared_corr = np.square(simulated_corrz[snp_iter])
			method_squared_corr = method_corr_squared[snp_iter]
			t.write(gene_name + '\t' + snp_rsid + '\t' + str(eqtl_ss) + '\t' + method_name + '\t' + str(simulated_corr) + '\t' + str(method_corr) + '\t' + str(simulated_squared_corr) + '\t' + str(method_squared_corr) + '\t' + str(variance_ratio) + '\n')
	t.flush()
f.close()
t.close()


















