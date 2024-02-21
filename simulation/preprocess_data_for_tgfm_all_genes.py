import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os
import pdb
import pickle
import time
import gzip
from pandas_plink import read_plink1_bin






def create_dictionary_mapping_from_rsid_to_genomic_annotation_vector(annotation_file):
	dicti = {}
	rs_id_vec = []
	variant_position_vec = []
	f = open(annotation_file)
	head_count=0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[2]
		annotation = np.asarray(data[4:]).astype(float)
		variant_position_vec.append(int(data[1]))
		rs_id_vec.append(rsid)
		# Quick error check
		if rsid in dicti:
			print('assumption eroror')
			pdb.set_trace()

		dicti[rsid] = annotation

	f.close()

	return dicti, np.asarray(variant_position_vec), np.asarray(rs_id_vec)

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

def mean_impute_and_standardize_genotype(G_obj_geno):
	# Fill in missing values
	G_obj_geno_stand = np.copy(G_obj_geno)
	ncol = G_obj_geno_stand.shape[1]
	n_missing = []
	for col_iter in range(ncol):
		nan_indices = np.isnan(G_obj_geno[:,col_iter])
		non_nan_mean = np.mean(G_obj_geno[nan_indices==False, col_iter])
		G_obj_geno_stand[nan_indices, col_iter] = non_nan_mean
		n_missing.append(np.sum(nan_indices))
	n_missing = np.asarray(n_missing)

	# Quick error check
	if np.sum(np.std(G_obj_geno_stand,axis=0) == 0) > 0:
		print('no variance genotype assumption error')
		pdb.set_trace()

	# Standardize genotype
	G_obj_geno_stand = (G_obj_geno_stand -np.mean(G_obj_geno_stand,axis=0))/np.std(G_obj_geno_stand,axis=0)

	return G_obj_geno_stand

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

def extract_sparse_genotype_matrix_for_this_gene_from_ld(gene_ld_mat, eqtl_sample_size):
	svd_lambda, svd_Q = np.linalg.eig(gene_ld_mat)
	# Filter out complex eigenvalues
	real_valued_eigs = np.iscomplex(svd_lambda) == False
	svd_lambda = svd_lambda[real_valued_eigs].astype(float)
	if np.sum(np.iscomplex(svd_Q[:, real_valued_eigs])) != 0.0:
		print('asssumption eroror')
		pdb.set_trace()
	svd_Q = svd_Q[:, real_valued_eigs].astype(float)
	# Filter out negative eigen values
	positive_eigs = svd_lambda > 0.0
	svd_lambda = svd_lambda[positive_eigs]
	svd_Q = svd_Q[:, positive_eigs]


	# Prune to PCs that explain 99% of variance
	pc_indices = extract_indices_of_pcs_that_explain_specified_fraction_of_variance(svd_lambda, .99)
	pruned_svd_lambda = svd_lambda[pc_indices]
	pruned_svd_Q = svd_Q[:, pc_indices]

	n_pcs = len(pruned_svd_lambda)
	max_num_pcs = int(np.floor(eqtl_sample_size*.75))
	if n_pcs >= max_num_pcs:
		pc_filter = np.asarray([False]*n_pcs)
		pc_filter[:max_num_pcs] = True
		pruned_svd_lambda = pruned_svd_lambda[pc_filter]
		pruned_svd_Q = pruned_svd_Q[:, pc_filter]
		n_pcs = len(pruned_svd_lambda)

	pruned_svd_ld_inv_mat = np.dot(np.dot(pruned_svd_Q, np.diag(1.0/pruned_svd_lambda)), np.transpose(pruned_svd_Q))

	return pruned_svd_lambda, pruned_svd_Q, pruned_svd_ld_inv_mat

def extract_sparse_genotype_matrix_for_this_gene(genotype_obj, cis_snp_indices, eqtl_sample_size):
	total_snps = len(cis_snp_indices)
	gene_snps = np.arange(total_snps)[cis_snp_indices]
	gene_snp_names = []
	for gene_snp in gene_snps:
		gene_snp_names.append('variant' + str(gene_snp))
	gene_snp_names = np.asarray(gene_snp_names)
	gene_variant_genotype = np.asarray(genotype_obj.sel(variant=gene_snp_names))
	gene_variant_genotype_stand = mean_impute_and_standardize_genotype(gene_variant_genotype)

	pruned_svd_lambda, pruned_svd_Q = reduce_dimensionality_of_genotype_data_with_pca(gene_variant_genotype_stand, eqtl_sample_size)
	pruned_svd_ld_inv_mat = np.dot(np.dot(pruned_svd_Q, np.diag(1.0/pruned_svd_lambda)), np.transpose(pruned_svd_Q))
	return pruned_svd_lambda, pruned_svd_Q, pruned_svd_ld_inv_mat

def reduce_dimensionality_of_genotype_data_with_pca(window_variant_genotype, eqtl_ss):
	# Reduce dimensionality of window variant genotype with svd
	n_ref = window_variant_genotype.shape[0]
	uu, ss, vh = np.linalg.svd(window_variant_genotype, full_matrices=False)
	svd_lambda = ss*ss/(n_ref-1)
	svd_Q = np.transpose(vh)
	# Prune to PCs that explain 99% of variance
	pc_indices = extract_indices_of_pcs_that_explain_specified_fraction_of_variance(svd_lambda, .99)
	pruned_svd_lambda = svd_lambda[pc_indices]
	pruned_svd_Q = svd_Q[:, pc_indices]

	n_pcs = len(pruned_svd_lambda)
	max_num_pcs = int(np.floor(eqtl_ss*.75))
	if n_pcs >= max_num_pcs:
		pc_filter = np.asarray([False]*n_pcs)
		pc_filter[:max_num_pcs] = True
		pruned_svd_lambda = pruned_svd_lambda[pc_filter]
		pruned_svd_Q = pruned_svd_Q[:, pc_filter]
		n_pcs = len(pruned_svd_lambda)

	return pruned_svd_lambda, pruned_svd_Q

def calculate_gene_variance_according_to_marginal_distribution(marginal_effects, pruned_svd_ld_inv_mat, pruned_svd_lambda, pruned_svd_Q, eqtl_sample_size):
	n_pcs = len(pruned_svd_lambda)

	alpha = np.dot(pruned_svd_ld_inv_mat, marginal_effects)

	# Compute sparse-pc predicted causal effects (delta)
	delta = np.dot(np.dot(np.diag(np.sqrt(pruned_svd_lambda)), np.transpose(pruned_svd_Q)), alpha)

	# compute variances
	r_sq = np.dot(alpha, marginal_effects)
	sigma = (1.0-r_sq)/(eqtl_sample_size - n_pcs - 1.0)

	# Observed h2
	obs_h2 = 1.0 - (1.0 - np.dot(delta, delta))*((eqtl_sample_size-1.0)/(eqtl_sample_size-n_pcs-1.0))

	genetic_gene_pmces_var = np.dot(alpha, marginal_effects) # This is equivelent to running alpha*S*alpha 
	genetic_gene_noise_var = sigma*n_pcs

	total_var = genetic_gene_pmces_var + genetic_gene_noise_var
	#stat = (np.sum(np.square(delta))/n_pcs)/sigma

	#pval = 1.0 - scipy.stats.f.cdf(stat, n_pcs, eqtl_sample_size-n_pcs-1)

	return total_var

def sample_eqtl_effects_from_susie_distribution(gene_susie_mu, gene_susie_alpha, gene_susie_mu_var):
	n_components = gene_susie_mu.shape[0]
	n_snps = gene_susie_mu.shape[1]

	sampled_eqtl_effects = np.zeros(n_snps)

	for component_iter in range(n_components):
		# Randomly draw snp index for component
		random_snp_index = np.random.choice(np.arange(n_snps).astype(int), replace=False, p=gene_susie_alpha[component_iter,:])
		
		effect_size_mean = gene_susie_mu[component_iter,random_snp_index]
		effect_size_var = gene_susie_mu_var[component_iter, random_snp_index]

		random_effect_size = np.random.normal(loc=effect_size_mean, scale=np.sqrt(effect_size_var))

		sampled_eqtl_effects[random_snp_index] = sampled_eqtl_effects[random_snp_index] + random_effect_size
	return sampled_eqtl_effects


def extract_gene_tissue_pairs_and_associated_gene_models_in_window(window_start, window_end, gene_summary_file, simulated_learned_gene_models_dir, simulation_name_string, eqtl_sample_size, window_indices, simulated_gene_expression_dir, ld_mat, eqtl_type, n_bs):
	# Initialize output vectors
	gene_tissue_pairs = []
	weight_vectors = []
	gene_tss_arr = []
	gene_variances = []
	full_gene_variances = []
	pmces_weights = []

	susie_mus = []
	susie_vars = []
	susie_alphas = []
	susie_indices = []

	bs_arr = []

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


		if gene_tss >= (window_start + 100000) and gene_tss <= (window_end - 100000):
			# Get indices corresponding to cis snps for this gene
			cis_snp_indices_raw = np.load(cis_snp_indices_file)
			cis_snp_indices = np.asarray([False]*total_n_genome_snps)
			cis_snp_indices[cis_snp_indices_raw] = True

			# Quick error check
			if len(cis_snp_indices) != len(window_indices):
				print('assumption eroror')
				pdb.set_trace()


			# Gene is in cis with respect to window

			# Fitted gene file
			if eqtl_type == 'susie':
				fitted_gene_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_all_genes_gene_model_pmces.npy'
				gene_model_mat = np.load(fitted_gene_file)
				
				#sim_gene_file = simulated_gene_expression_dir + simulation_name_string + '_' + ensamble_id + '_causal_eqtl_effects.npy'
				#sim_gene_model_mat = np.transpose(np.load(sim_gene_file))

			# Relevent info
			n_tiss = gene_model_mat.shape[0]
			n_cis_snps = gene_model_mat.shape[1]

			# QUick error check
			if np.sum(cis_snp_indices) != n_cis_snps:
				print('assumption eroror')
				pdb.set_trace()

			# Loop through tissues
			for tiss_iter in range(n_tiss):
				# Skip gene, tissue pairs with no gene models
				if np.array_equal(gene_model_mat[tiss_iter,:], np.zeros(n_cis_snps)):
					continue

				# Compute variance of gene
				if eqtl_type == 'susie':
					gene_variance = np.dot(np.dot(gene_model_mat[tiss_iter,:], ld_mat[cis_snp_indices[window_indices],:][:,cis_snp_indices[window_indices]]), gene_model_mat[tiss_iter,:])
					gene_susie_mu_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_tissue_' + str(tiss_iter) + '_all_genes_gene_model_susie_mu.npy'
					gene_susie_alpha_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_tissue_' + str(tiss_iter) + '_all_genes_gene_model_susie_alpha.npy'
					gene_susie_mu_var_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_tissue_' + str(tiss_iter) + '_all_genes_gene_model_susie_mu_var.npy'
					gene_susie_mu = np.load(gene_susie_mu_file)
					gene_susie_alpha = np.load(gene_susie_alpha_file)
					gene_susie_mu_var = np.load(gene_susie_mu_var_file)


					window_level_eqtl_effect_sizes = np.zeros(np.sum(window_indices))
					window_level_eqtl_effect_sizes[cis_snp_indices[window_indices]] = gene_model_mat[tiss_iter,:]/np.sqrt(gene_variance)

					#full_gene_variance = calculate_gene_variance_according_to_susie_distribution(gene_susie_mu, gene_susie_alpha, np.sqrt(gene_susie_mu_var), ld_mat[cis_snp_indices[window_indices],:][:,cis_snp_indices[window_indices]])

					n_snps = gene_susie_mu.shape[1]
					bs_alpha_effects = np.zeros((n_snps, n_bs))

					for bs_iter in range(n_bs):
						susie_sampled_eqtl_effects = sample_eqtl_effects_from_susie_distribution(gene_susie_mu, gene_susie_alpha, gene_susie_mu_var)
						bs_alpha_effects[:, bs_iter] = susie_sampled_eqtl_effects
					
					bs_gene_variances = np.diag(np.dot(np.dot(np.transpose(bs_alpha_effects), ld_mat[cis_snp_indices[window_indices],:][:,cis_snp_indices[window_indices]]), bs_alpha_effects))
					bs_std_alpha_effects = bs_alpha_effects/np.sqrt(bs_gene_variances)
					bs_arr.append((bs_std_alpha_effects, cis_snp_indices[window_indices]))
					susie_mus.append(gene_susie_mu)
					susie_vars.append(gene_susie_mu_var)
					susie_alphas.append(gene_susie_alpha)
					susie_indices.append(cis_snp_indices[window_indices])

				else:
					print('assumption erooror')
					pdb.set_trace()

				# Add info to arrays
				gene_tss_arr.append(gene_tss)
				gene_tissue_pairs.append(ensamble_id + '_' + 'tissue' + str(tiss_iter))
				gene_variances.append(gene_variance)
				#full_gene_variances.append(full_gene_variance)
				pmces_weights.append(window_level_eqtl_effect_sizes)

	f.close()

	return np.asarray(gene_tissue_pairs), bs_arr, np.asarray(gene_tss_arr), np.asarray(pmces_weights), np.asarray(gene_variances), susie_mus, susie_vars, susie_alphas, susie_indices


def create_anno_matrix_for_set_of_rsids(rsid_to_genomic_annotation, window_rsids):
	anno_mat_arr = []
	for window_rsid in window_rsids:
		anno_mat_arr.append(rsid_to_genomic_annotation[window_rsid])
	return np.asarray(anno_mat_arr)


def compute_log_expected_probability(full_anno_mat, sldsc_tau_mean, threshold):
	expected_per_ele_h2 = np.dot(full_anno_mat, sldsc_tau_mean)
	expected_per_ele_h2[expected_per_ele_h2 <= threshold] = threshold

	prob = expected_per_ele_h2/np.sum(expected_per_ele_h2)

	return np.log(prob)

def compute_log_expected_probability_with_ratio_to_max(full_anno_mat, sldsc_tau_mean, ratio_to_max, max_variant_value, max_tissue_value):
	expected_per_ele_h2 = np.dot(full_anno_mat, sldsc_tau_mean)
	variant_indices = full_anno_mat[:,0] == 1
	gene_indices = full_anno_mat[:,0] == 0

	variant_h2 = np.copy(expected_per_ele_h2[variant_indices])
	variant_h2[variant_h2 <= ratio_to_max*max_variant_value] = ratio_to_max*max_variant_value
	
	gene_h2 = np.copy(expected_per_ele_h2[gene_indices])
	gene_h2[gene_h2 <= ratio_to_max*max_tissue_value] = ratio_to_max*max_tissue_value

	expected_per_ele_h2[variant_indices] = variant_h2
	expected_per_ele_h2[gene_indices] = gene_h2

	prob = expected_per_ele_h2/np.sum(expected_per_ele_h2)

	return np.log(prob)

def compute_shared_variant_log_expected_probability_with_ratio_to_max(full_anno_mat, sldsc_tau_mean, ratio_to_max, max_variant_value, max_tissue_value):
	expected_per_ele_h2 = np.dot(full_anno_mat, sldsc_tau_mean)
	variant_indices = full_anno_mat[:,0] == 1
	gene_indices = full_anno_mat[:,0] == 0

	variant_h2 = np.copy(expected_per_ele_h2[variant_indices])
	per_variant_h2 = np.sum(variant_h2)/len(variant_h2)
	if per_variant_h2 < 1e-8:
		per_variant_h2 = 1e-8
	shared_variant_h2 = variant_h2*0.0 + per_variant_h2
	gene_h2 = np.copy(expected_per_ele_h2[gene_indices])
	gene_h2[gene_h2 <= max_tissue_value*ratio_to_max] = max_tissue_value*ratio_to_max

	expected_per_ele_h2[variant_indices] = shared_variant_h2
	expected_per_ele_h2[gene_indices] = gene_h2

	prob = expected_per_ele_h2/np.sum(expected_per_ele_h2)
	return np.log(prob)

def compute_shared_variant_expected_log_probability_with_ratio_to_max(full_anno_mat, sldsc_tau_mean, sldsc_tau_cov, ratio_to_max, max_variant_value, max_tissue_value, n_samples=10000):
	# Sample a whole bunch of taus
	sampled_taus = np.random.multivariate_normal(mean=sldsc_tau_mean, cov=sldsc_tau_cov, size=n_samples)

	sampled_expected_per_ele_h2 = np.dot(full_anno_mat, np.transpose(sampled_taus))

	variant_indices = full_anno_mat[:,0] == 1
	gene_indices = full_anno_mat[:,0] == 0

	sampled_expected_per_variant_h2 = sampled_expected_per_ele_h2[variant_indices,:]
	n_var = np.sum(variant_indices)
	sampled_shared_variant_h2 = np.sum(sampled_expected_per_variant_h2,axis=0)/n_var

def compute_shared_variant_expected_log_probability(full_anno_mat, sldsc_tau_mean,sldsc_tau_cov, threshold, n_samples=10000):
	# Sample a whole bunch of taus
	sampled_taus = np.random.multivariate_normal(mean=sldsc_tau_mean, cov=sldsc_tau_cov, size=n_samples)

	sampled_expected_per_ele_h2 = np.dot(full_anno_mat, np.transpose(sampled_taus))

	variant_indices = full_anno_mat[:,0] == 1
	gene_indices = full_anno_mat[:,0] == 0

	sampled_expected_per_variant_h2 = sampled_expected_per_ele_h2[variant_indices,:]
	n_var = np.sum(variant_indices)
	sampled_shared_variant_h2 = np.sum(sampled_expected_per_variant_h2,axis=0)/n_var

	prob = np.copy(sampled_expected_per_ele_h2)
	for sample_iter in range(n_samples):
		prob[variant_indices, sample_iter] = prob[variant_indices, sample_iter]*0.0 + sampled_shared_variant_h2[sample_iter]
		neg_indices = prob[:, sample_iter] <= threshold
		prob[neg_indices, sample_iter] = threshold
		prob[:, sample_iter] = prob[:, sample_iter]/np.sum(prob[:, sample_iter])

	log_prob = np.log(prob)

	return np.mean(log_prob,axis=1)


def compute_shared_variant_log_expected_probability(full_anno_mat, sldsc_tau_mean, threshold):
	expected_per_ele_h2 = np.dot(full_anno_mat, sldsc_tau_mean)
	variant_indices = full_anno_mat[:,0] == 1
	gene_indices = full_anno_mat[:,0] == 0

	variant_h2 = np.copy(expected_per_ele_h2[variant_indices])
	per_variant_h2 = np.sum(variant_h2)/len(variant_h2)

	shared_variant_h2 = variant_h2*0.0 + per_variant_h2
	gene_h2 = np.copy(expected_per_ele_h2[gene_indices])

	expected_per_ele_h2[variant_indices] = shared_variant_h2
	expected_per_ele_h2[gene_indices] = gene_h2

	expected_per_ele_h2[expected_per_ele_h2 <= threshold] = threshold

	prob = expected_per_ele_h2/np.sum(expected_per_ele_h2)
	return np.log(prob)

def compute_expected_log_probability_with_ratio_to_max(full_anno_mat, sldsc_tau_mean, sldsc_tau_cov, ratio_to_max, max_variant_value, max_tissue_value, n_samples=10000):
	# Sample a whole bunch of taus
	sampled_taus = np.random.multivariate_normal(mean=sldsc_tau_mean, cov=sldsc_tau_cov, size=n_samples)

	sampled_expected_per_ele_h2 = np.dot(full_anno_mat, np.transpose(sampled_taus))

	variant_indices = full_anno_mat[:,0] == 1
	gene_indices = full_anno_mat[:,0] == 0

	sampled_expected_per_variant_h2 = sampled_expected_per_ele_h2[variant_indices,:]
	sampled_expected_per_variant_h2[sampled_expected_per_variant_h2 <= ratio_to_max*max_variant_value] = ratio_to_max*max_variant_value

	sampled_expected_per_gene_h2 = sampled_expected_per_ele_h2[gene_indices,:]
	sampled_expected_per_gene_h2[sampled_expected_per_gene_h2 <= ratio_to_max*max_tissue_value] = ratio_to_max*max_variant_value

	sampled_expected_per_ele_h2[variant_indices,:] = sampled_expected_per_variant_h2
	sampled_expected_per_ele_h2[gene_indices,:] = sampled_expected_per_gene_h2

	prob = np.copy(sampled_expected_per_ele_h2)
	for sample_iter in range(n_samples):
		prob[:, sample_iter] = prob[:, sample_iter]/np.sum(prob[:, sample_iter])

	log_prob = np.log(prob)

	return np.mean(log_prob,axis=1)

def compute_expected_log_probability(full_anno_mat, sldsc_tau_mean, sldsc_tau_cov, threshold, n_samples=10000):
	# Sample a whole bunch of taus
	sampled_taus = np.random.multivariate_normal(mean=sldsc_tau_mean, cov=sldsc_tau_cov, size=n_samples)

	sampled_expected_per_ele_h2 = np.dot(full_anno_mat, np.transpose(sampled_taus))

	sampled_expected_per_ele_h2[sampled_expected_per_ele_h2 <= threshold] = threshold

	prob = np.copy(sampled_expected_per_ele_h2)
	for sample_iter in range(n_samples):
		prob[:, sample_iter] = prob[:, sample_iter]/np.sum(prob[:, sample_iter])

	log_prob = np.log(prob)

	return np.mean(log_prob,axis=1)

def compute_log_expected_probability_variant_v_gene_only(full_anno_mat, sldsc_tau_mean, threshold):
	expected_per_ele_h2 = np.dot(full_anno_mat, sldsc_tau_mean)
	variant_elements = full_anno_mat[:,0] == 1
	gene_elements = full_anno_mat[:,0] != 1
	# Quick error check
	if np.sum(variant_elements) + np.sum(gene_elements) != full_anno_mat.shape[0]:
		print('assumption eroror')
		pdb.set_trace()
	window_snp_h2 = np.sum(expected_per_ele_h2[variant_elements])
	window_gene_h2 = np.sum(expected_per_ele_h2[gene_elements])
	if window_snp_h2 < threshold:
		window_snp_h2 = threshold
	if window_gene_h2 < threshold:
		window_gene_h2 = threshold

	avg_per_ele_h2 = np.ones(len(expected_per_ele_h2))
	avg_per_ele_h2[variant_elements] = window_snp_h2/np.sum(variant_elements)
	avg_per_ele_h2[gene_elements] = window_gene_h2/np.sum(gene_elements)

	prob = avg_per_ele_h2/np.sum(avg_per_ele_h2)

	return np.log(prob)


def compute_various_versions_of_log_prior_probabilities_with_ratio_to_max(window_rsids, window_snp_anno_mat, gene_tissue_pairs, sldsc_tau_mean, sldsc_tau_cov, sparse_sldsc_tau, ratio_to_max):
	# Merge together window_snp_anno_mat with tissue_anno_mat
	n_snps = len(window_rsids)
	n_genes = len(gene_tissue_pairs)
	tissue_anno_mat = np.zeros((n_genes, 10))
	for gene_iter, gene_tissue_pair in enumerate(gene_tissue_pairs):
		tissue_index = int(gene_tissue_pair.split('_')[1].split('issue')[1])
		tissue_anno_mat[gene_iter, tissue_index] = 1

	full_anno_mat_top = np.hstack((window_snp_anno_mat, np.zeros((n_snps,10))))
	full_anno_mat_bottom = np.hstack((np.zeros((n_genes, window_snp_anno_mat.shape[1])), tissue_anno_mat))
	full_anno_mat = np.vstack((full_anno_mat_top, full_anno_mat_bottom))

	# Quick error check
	if full_anno_mat.shape[1] != len(sldsc_tau_mean):
		print('assumption error')
		pdb.set_trace()

	max_variant_value = np.max(np.dot(window_snp_anno_mat, sldsc_tau_mean[0:-10]))
	if max_variant_value < 1e-8:
		max_variant_value = 1e-8
	max_tissue_value = np.max(sldsc_tau_mean[-10:])
	if max_tissue_value < 1e-8:
		max_tissue_value = 1e-8

	#point_estimate_ln_pi = compute_log_expected_probability_with_ratio_to_max(full_anno_mat, sldsc_tau_mean, ratio_to_max, max_variant_value, max_tissue_value)
	#distribution_estimate_ln_pi = compute_expected_log_probability_with_ratio_to_max(full_anno_mat, sldsc_tau_mean, sldsc_tau_cov, ratio_to_max, max_variant_value, max_tissue_value)
	#shared_variant_point_estimate_ln_pi = compute_shared_variant_log_expected_probability_with_ratio_to_max(full_anno_mat, sldsc_tau_mean, ratio_to_max, max_variant_value, max_tissue_value)
	shared_variant_distribution_estimate_ln_pi = compute_shared_variant_expected_log_probability_with_ratio_to_max(full_anno_mat, sldsc_tau_mean, sldsc_tau_cov, ratio_to_max, max_variant_value, max_tissue_value)


def compute_various_versions_of_log_prior_probabilities(window_rsids, window_snp_anno_mat, gene_tissue_pairs, sldsc_tau_mean, sldsc_tau_cov, sparse_sldsc_tau, threshold=1e-30):
	# Merge together window_snp_anno_mat with tissue_anno_mat
	n_snps = len(window_rsids)
	n_genes = len(gene_tissue_pairs)
	tissue_anno_mat = np.zeros((n_genes, 10))
	for gene_iter, gene_tissue_pair in enumerate(gene_tissue_pairs):
		tissue_index = int(gene_tissue_pair.split('_')[1].split('issue')[1])
		tissue_anno_mat[gene_iter, tissue_index] = 1

	full_anno_mat_top = np.hstack((window_snp_anno_mat, np.zeros((n_snps,10))))
	full_anno_mat_bottom = np.hstack((np.zeros((n_genes, window_snp_anno_mat.shape[1])), tissue_anno_mat))
	full_anno_mat = np.vstack((full_anno_mat_top, full_anno_mat_bottom))

	# Quick error check
	if full_anno_mat.shape[1] != len(sldsc_tau_mean):
		print('assumption error')
		pdb.set_trace()

	shared_variant_point_estimate_ln_pi = compute_shared_variant_log_expected_probability(full_anno_mat, sldsc_tau_mean, threshold)
	shared_variant_sparse_estimate_ln_pi = compute_shared_variant_log_expected_probability(full_anno_mat, sparse_sldsc_tau, threshold)
	shared_variant_distribution_estimate_ln_pi = compute_shared_variant_expected_log_probability(full_anno_mat, sldsc_tau_mean,sldsc_tau_cov, threshold)
	variant_v_gene_only_ln_pi = compute_log_expected_probability_variant_v_gene_only(full_anno_mat, sldsc_tau_mean, threshold)	
	point_estimate_ln_pi = compute_log_expected_probability(full_anno_mat, sldsc_tau_mean, threshold)
	sparse_estimate_ln_pi = compute_log_expected_probability(full_anno_mat, sparse_sldsc_tau, threshold)
	distribution_estimate_ln_pi = compute_expected_log_probability(full_anno_mat, sldsc_tau_mean, sldsc_tau_cov, threshold)

	return point_estimate_ln_pi, sparse_estimate_ln_pi, distribution_estimate_ln_pi, variant_v_gene_only_ln_pi, shared_variant_point_estimate_ln_pi, shared_variant_distribution_estimate_ln_pi, shared_variant_sparse_estimate_ln_pi

def load_in_window_gwas_betas_and_ses(window_gwas_summary_file, window_rsids):
	f = open(window_gwas_summary_file)
	beta_vec =[]
	beta_se_vec = []
	rsid_vec = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[0]
		beta = float(data[1])
		beta_se = float(data[2])
		beta_vec.append(beta)
		beta_se_vec.append(beta_se)
		rsid_vec.append(rsid)

	f.close()

	if np.array_equal(np.asarray(rsid_vec), window_rsids) == False:
		print('assumption eroror')
		pdb.set_trace()

	return np.asarray(beta_vec), np.asarray(beta_se_vec)

# Convert gwas summary statistics to *STANDARDIZED* effect sizes
# Following SuSiE code found in these two places:
########1. https://github.com/stephenslab/susieR/blob/master/R/susie_rss.R  (LINES 277-279)
########2. https://github.com/stephenslab/susieR/blob/master/R/susie_ss.R (LINES 148-156 AND 203-205)
def convert_to_standardized_summary_statistics(gwas_beta_raw, gwas_beta_se_raw, gwas_sample_size, R, sigma2=1.0):
	gwas_z_raw = gwas_beta_raw/gwas_beta_se_raw

	XtX = (gwas_sample_size-1)*R
	Xty = np.sqrt(gwas_sample_size-1)*gwas_z_raw
	var_y = 1

	dXtX = np.diag(XtX)
	csd = np.sqrt(dXtX/(gwas_sample_size-1))
	csd[csd == 0] = 1

	XtX = (np.transpose((1/csd) * XtX) / csd)
	Xty = Xty / csd

	dXtX2 = np.diag(XtX)

	beta_scaled = (1/dXtX2)*Xty
	beta_se_scaled = np.sqrt(sigma2/dXtX2)

	return beta_scaled, beta_se_scaled, XtX


def save_ln_pi_output_file(ln_pi, output_file, window_rsids, gene_tissue_pairs):
	t2 = open(output_file,'w')
	t2.write('element_name\tln_pi\n')

	genetic_element_names = np.hstack((window_rsids, gene_tissue_pairs))

	# Quick error check
	if len(genetic_element_names) != len(ln_pi):
		print('assumption error')
		pdb.set_trace()

	for itera, genetic_element_name in enumerate(genetic_element_names):
		t2.write(genetic_element_name + '\t' + str(ln_pi[itera]) + '\n')

	t2.close()

	return


#####################
# Command line args
#####################
simulation_number = int(sys.argv[1])
chrom_num = sys.argv[2]
simulation_name_string = sys.argv[3]
n_gwas_individuals = int(sys.argv[4])
eqtl_sample_size = sys.argv[5]
simulation_window_list_file = sys.argv[6]
annotation_file = sys.argv[7]
simulated_gwas_dir = sys.argv[8]
simulated_gene_expression_dir = sys.argv[9]
simulated_learned_gene_models_dir = sys.argv[10]
simulated_tgfm_input_data_dir = sys.argv[11]
eqtl_type = sys.argv[12]
processed_genotype_data_dir = sys.argv[13]


n_bs = 100


# Create dictionary mapping from rsid to genomic annotation vector
rsid_to_genomic_annotation, variant_position_vec, rsids = create_dictionary_mapping_from_rsid_to_genomic_annotation_vector(annotation_file)


# Open outputful summarizing TGFM input (one line for each window)
tgfm_input_data_summary_file = simulated_tgfm_input_data_dir + simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + eqtl_type + '_all_genes_bootstrapped_tgfm_input_data_summary.txt'
t = open(tgfm_input_data_summary_file,'w')
# Write header
t.write('window_name\tLD_npy_file\tTGFM_input_pkl\n')


# Now loop through windows
f = open(simulation_window_list_file)
head_count = 0

for line in f:
	line = line.rstrip()
	data = line.split('\t')
	# Skip header
	if head_count == 0:
		head_count = head_count + 1
		continue

	# Extract relevent info for this window
	window_name = data[0]
	print(window_name)

	window_start = int(data[1])
	window_middle_start = int(data[2])
	window_middle_end = int(data[3])
	window_end = int(data[4])

	# Extract indices of variants in this window
	window_indices = (variant_position_vec >= window_start) & (variant_position_vec < window_end)
	window_rsids = rsids[window_indices]
	window_variant_position_vec = variant_position_vec[window_indices]

	if len(window_rsids) < 20:
		print('skipped window: ' + window_name)
		continue

	# Create annotation matrix for window rsids
	#window_anno_mat = create_anno_matrix_for_set_of_rsids(rsid_to_genomic_annotation, window_rsids)

	# Extract LD
	ld_mat_file = processed_genotype_data_dir + window_name + '_in_sample_ld.npy'
	ld_mat = np.load(ld_mat_file)

	# Extract gene-tissue pairs and fitted models in this window
	gene_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'

	gene_tissue_pairs, bs_eqtl_arr, gene_tissue_pairs_tss, pmces_weights, gene_variances, gene_susie_mus, gene_susie_mu_vars, gene_susie_alphas, gene_susie_indices = extract_gene_tissue_pairs_and_associated_gene_models_in_window(window_start, window_end, gene_summary_file, simulated_learned_gene_models_dir, simulation_name_string,  eqtl_sample_size, window_indices, simulated_gene_expression_dir,ld_mat, eqtl_type, n_bs)

	# Option to save some memory
	# Change to current version used in real data
	#bs_gene_eqtl_files = []
	sparse_sampled_gene_eqtl_pmces = []
	for bs_iter in range(n_bs):
		eqtl_mat = []
		for gene_itera, bs_eqtl_tuple in enumerate(bs_eqtl_arr):
			eqtl_gene_window = bs_eqtl_tuple[0][:, bs_iter]
			boolean_indices = bs_eqtl_tuple[1]
			eqtl_indices = np.arange(len(boolean_indices))[boolean_indices]
			for ii, eqtl_effect in enumerate(eqtl_gene_window):
				if eqtl_effect == 0.0:
					continue
				eqtl_index = eqtl_indices[ii]
				eqtl_mat.append(np.asarray([gene_itera, eqtl_index, eqtl_effect]))
		eqtl_mat = np.asarray(eqtl_mat)
		sparse_sampled_gene_eqtl_pmces.append(eqtl_mat)
		#bs_eqtl_output_file = simulated_tgfm_input_data_dir + simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + eqtl_type + '_' + window_name + '_' + str(bs_iter) + '.npy'
		#np.save(bs_eqtl_output_file, eqtl_mat)
		#bs_gene_eqtl_files.append(bs_eqtl_output_file)
	#bs_gene_eqtl_files = np.asarray(bs_gene_eqtl_files)


	# Get middle variant indices and middle gene indices
	middle_variant_indices = np.where((window_variant_position_vec >= window_middle_start) & (window_variant_position_vec < window_middle_end))[0]
	middle_gene_indices = np.where((gene_tissue_pairs_tss >= window_middle_start) & (gene_tissue_pairs_tss < window_middle_end))[0]

	# Load in GWAS betas and standard errors
	window_gwas_summary_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results_window_' + window_name + '.txt'
	gwas_beta, gwas_beta_se = load_in_window_gwas_betas_and_ses(window_gwas_summary_file, window_rsids)
	# Standardize gwas beta and se
	beta_scaled, beta_se_scaled, XtX = convert_to_standardized_summary_statistics(gwas_beta, gwas_beta_se, n_gwas_individuals, ld_mat)


	# Compute various ln(pi) and save to output # and save those results to output files
	#ln_pi_output_stem = simulated_tgfm_input_data_dir + simulation_name_string + '_' + window_name + '_eqtl_ss_' + str(eqtl_sample_size)+ '_' + eqtl_type + '_ln_pi'
	# Uniform prior
	#n_window_elements = len(window_rsids) + len(gene_tissue_pairs)
	#uniform_pi = np.ones(n_window_elements)*(1.0/n_window_elements)
	#uniform_ln_pi = np.log(uniform_pi)
	#save_ln_pi_output_file(uniform_ln_pi, ln_pi_output_stem + '_uniform.txt', window_rsids, gene_tissue_pairs)

	# Organize TGFM data into nice data structure
	tgfm_data = {}
	tgfm_data['genes'] = gene_tissue_pairs
	tgfm_data['variants'] = window_rsids
	tgfm_data['gwas_beta'] = beta_scaled
	tgfm_data['gwas_beta_se'] = beta_se_scaled
	tgfm_data['gwas_sample_size'] = n_gwas_individuals
	tgfm_data['sparse_sampled_gene_eqtl_pmces'] = sparse_sampled_gene_eqtl_pmces
	tgfm_data['middle_gene_indices'] = middle_gene_indices
	tgfm_data['middle_variant_indices'] = middle_variant_indices
	tgfm_data['gene_eqtl_pmces'] = pmces_weights
	tgfm_data['gene_variances'] = gene_variances
	#tgfm_data['full_gene_variances'] = full_gene_variances
	#tgfm_data['annotation'] = window_anno_mat
	tgfm_data['tss'] = gene_tissue_pairs_tss
	tgfm_data['variant_positions'] = window_variant_position_vec
	tgfm_data['gene_susie_mu'] = gene_susie_mus
	tgfm_data['gene_susie_mu_var'] = gene_susie_mu_vars
	tgfm_data['gene_susie_alpha'] = gene_susie_alphas
	tgfm_data['gene_susie_indices'] = gene_susie_indices

	# Save TGFM output data to pickle
	window_pickle_output_file = simulated_tgfm_input_data_dir + simulation_name_string + '_' + window_name + '_eqtl_ss_' + str(eqtl_sample_size)+ '_' + eqtl_type + '_all_genes_tgfm_input_data.pkl'
	g = open(window_pickle_output_file, "wb")
	pickle.dump(tgfm_data, g)
	g.close()

	# Write to output summary file
	t.write(window_name + '\t' + ld_mat_file + '\t' + window_pickle_output_file + '\n')

t.close()
f.close()





'''
	# Compute various ln(pi) and save to output # and save those results to output files
	ln_pi_output_stem = simulated_tgfm_input_data_dir + simulation_name_string + '_' + window_name + '_eqtl_ss_' + str(eqtl_sample_size)+ '_' + eqtl_type + '_ln_pi'
	# Uniform prior
	n_window_elements = len(window_rsids) + len(gene_tissue_pairs)
	uniform_pi = np.ones(n_window_elements)*(1.0/n_window_elements)
	uniform_ln_pi = np.log(uniform_pi)
	save_ln_pi_output_file(uniform_ln_pi, ln_pi_output_stem + '_uniform.txt', window_rsids, gene_tissue_pairs)
	## V2
	#ratio_to_maxs = [0.01, 0.001]
	#for ratio_to_max in ratio_to_maxs:
	#	point_estimate_ln_pi, distribution_estimate_ln_pi, shared_variant_point_estimate_ln_pi, shared_variant_distribution_estimate_ln_pi = compute_various_versions_of_log_prior_probabilities_with_ratio_to_max(window_rsids, window_anno_mat, gene_tissue_pairs, sldsc_tau_mean, sldsc_tau_cov, sparse_sldsc_tau, ratio_to_max)
	# v1
	thresholds = [1e-8,1e-10,1e-30]
	thresholds = [1e-8]

	for threshold in thresholds:
		point_estimate_ln_pi, sparse_estimate_ln_pi, distribution_estimate_ln_pi, variant_v_gene_only_ln_pi,shared_variant_point_estimate_ln_pi, shared_variant_distribution_estimate_ln_pi, shared_variant_sparse_estimate_ln_pi  = compute_various_versions_of_log_prior_probabilities(window_rsids, window_anno_mat, gene_tissue_pairs, sldsc_tau_mean, sldsc_tau_cov, sparse_sldsc_tau, threshold=threshold)
		save_ln_pi_output_file(point_estimate_ln_pi, ln_pi_output_stem + '_point_estimate_' + str(threshold) + '.txt', window_rsids, gene_tissue_pairs)
		save_ln_pi_output_file(sparse_estimate_ln_pi, ln_pi_output_stem + '_sparse_estimate_' + str(threshold) + '.txt', window_rsids, gene_tissue_pairs)
		save_ln_pi_output_file(distribution_estimate_ln_pi, ln_pi_output_stem + '_distribution_estimate_' + str(threshold) + '.txt', window_rsids, gene_tissue_pairs)
		save_ln_pi_output_file(variant_v_gene_only_ln_pi, ln_pi_output_stem + '_variant_v_gene_only_' + str(threshold) + '.txt', window_rsids, gene_tissue_pairs)
		save_ln_pi_output_file(shared_variant_point_estimate_ln_pi, ln_pi_output_stem + '_shared_variant_point_estimate_' + str(threshold) + '.txt', window_rsids, gene_tissue_pairs)
		save_ln_pi_output_file(shared_variant_distribution_estimate_ln_pi, ln_pi_output_stem + '_shared_variant_distribution_estimate_' + str(threshold) + '.txt', window_rsids, gene_tissue_pairs)
		save_ln_pi_output_file(shared_variant_sparse_estimate_ln_pi, ln_pi_output_stem + '_shared_variant_sparse_estimate_' + str(threshold) + '.txt', window_rsids, gene_tissue_pairs)

'''



