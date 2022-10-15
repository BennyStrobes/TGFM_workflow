import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pyreadr
import numpy as np 
import os
import sys
import pdb
import pickle
import time




def get_pseudotissue_names(gtex_pseudotissue_file):
	f = open(gtex_pseudotissue_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0])
	f.close()
	return np.asarray(arr)

def get_dictionary_list_of_hapmap3_snpids(hapmap3_snpid_file):
	f = open(hapmap3_snpid_file)
	dicti = {}
	for line in f:
		line = line.rstrip()
		dicti[line] = 1
	f.close()
	return dicti

def create_tissue_to_gene_model_df(pseudotissues, gtex_susie_gene_models_dir):
	mapping = {}
	for pseudotissue in pseudotissues:
		gene_model_summary_file = gtex_susie_gene_models_dir + pseudotissue + '/' + pseudotissue + '_cis_heritable_gene_pos_file.txt'
		raw_table = np.loadtxt(gene_model_summary_file, dtype=str, delimiter='\t')[1:]
		mapping[pseudotissue] = raw_table
	return mapping

def get_list_of_gene_tissue_pairs_and_weight_files_for_window(window_chrom_num, window_start, window_end, pseudotissues, tissue_to_gene_model_df):
	# Initialize output vectors
	gene_tissue_pairs = []
	weight_files = []
	gene_tss_arr = []
	# Loop through tissues
	for pseudotissue in pseudotissues:
		# Identify which genes are in cis of window
		pseudotissue_gene_model_df = tissue_to_gene_model_df[pseudotissue]
		df_chrom_num = (pseudotissue_gene_model_df[:,2]).astype(int)
		df_gene_name = pseudotissue_gene_model_df[:,1]
		df_weight_files = pseudotissue_gene_model_df[:,0]
		df_tss = (pseudotissue_gene_model_df[:,3]).astype(int)
		# Loop through genes
		num_genes = pseudotissue_gene_model_df.shape[0]
		for gene_num in range(num_genes):
			# Check if gene in cis with window
			if df_chrom_num[gene_num] == window_chrom_num:
				gene_tss = df_tss[gene_num]
				if gene_tss >= (window_start + 500000) and gene_tss <= (window_end - 500000):
					# Gene is in cis with respect to window
					gene_name = df_gene_name[gene_num]
					wgt_file = df_weight_files[gene_num]
					gene_tissue_pairs.append(gene_name + '_' + pseudotissue)
					weight_files.append(wgt_file)
					gene_tss_arr.append(gene_tss)
	return np.asarray(gene_tissue_pairs), np.asarray(weight_files), np.asarray(gene_tss_arr)

def create_mapping_from_variant_name_to_global_variant_arr_pos(global_variant_arr):
	dicti = {}
	for i, val in enumerate(global_variant_arr):
		dicti[val] = i
	return dicti

def create_gene_local_to_global_mapping(variant_name_to_global_variant_arr_pos, gene_variants):
	local_to_global = []
	flips = []

	for gene_variant in gene_variants:
		gene_variant_info = gene_variant.split('_')
		gene_variant_norm = gene_variant_info[0] + '_' + gene_variant_info[1] + '_' + gene_variant_info[2] + '_' + gene_variant_info[3]
		gene_variant_alt = gene_variant_info[0] + '_' + gene_variant_info[1] + '_' + gene_variant_info[3] + '_' + gene_variant_info[2]
		if gene_variant_norm in variant_name_to_global_variant_arr_pos:
			local_to_global.append(variant_name_to_global_variant_arr_pos[gene_variant_norm])
			flips.append(1.0)
		elif gene_variant_alt in variant_name_to_global_variant_arr_pos:
			local_to_global.append(variant_name_to_global_variant_arr_pos[gene_variant_alt])
			flips.append(-1.0)
		else:
			print('gene variant mapping assumption erooror')
			pdb.set_trace()
	return np.asarray(local_to_global), np.asarray(flips)

def create_global_susie_no_flips(gene_local_to_global_mapping, local_susie_data, num_global_variants):
	global_susie_data = np.zeros((local_susie_data.shape[0], num_global_variants))
	global_susie_data[:, gene_local_to_global_mapping] = local_susie_data
	return global_susie_data

def create_global_susie_w_flips(gene_local_to_global_mapping, gene_flips, local_susie_data, num_global_variants):
	for var_index, flip_value in enumerate(gene_flips):
		if flip_value == -1.0:
			local_susie_data[:, var_index] = local_susie_data[:, var_index]*-1.0
	global_susie_data = create_global_susie_no_flips(gene_local_to_global_mapping, local_susie_data, num_global_variants)
	return global_susie_data

def organize_tgfm_trait_agnostic_data_for_single_window(window_name, window_chrom_num, window_start, window_end, global_variant_arr, LD, pseudotissues, tissue_to_gene_model_df):
	# Get list of gene_tissue pairs and weight files for this window
	gene_tissue_pairs, gene_tissue_pair_weight_files, gene_tss_arr = get_list_of_gene_tissue_pairs_and_weight_files_for_window(window_chrom_num, window_start, window_end, pseudotissues, tissue_to_gene_model_df)

	# Get number of total variants
	num_global_variants = len(global_variant_arr)

	# Create dictionary mapping from variant name to global-variant arr position
	variant_name_to_global_variant_arr_pos = create_mapping_from_variant_name_to_global_variant_arr_pos(global_variant_arr)
	

	# Initialize arrays to keep track of eqtl-susie data
	global_susie_mu_data = []
	global_susie_alpha_data = []
	global_susie_mu_sd_data = []

	# Loop through gene-tissue pairs
	for gene_tissue_index, gene_tissue_name in enumerate(gene_tissue_pairs):
		# Get gene-tissue weight file
		gene_tissue_weight_file = gene_tissue_pair_weight_files[gene_tissue_index]

		# Load gene-tissue model from gene-tissue weight file
		gene_tissue_model = pyreadr.read_r(gene_tissue_weight_file)

		# Create mapping from gene coordinates to global coordinates
		gene_variants = np.asarray(gene_tissue_model['variant_names'])[:,0]
		gene_local_to_global_mapping, gene_flips = create_gene_local_to_global_mapping(variant_name_to_global_variant_arr_pos, gene_variants)

		# Create global susie alpha
		global_susie_alpha = create_global_susie_no_flips(gene_local_to_global_mapping, np.asarray(gene_tissue_model['susie_alpha']), num_global_variants)
		global_susie_alpha_data.append(global_susie_alpha)
		# Create global susie mu_sd
		local_susie_mu_sd = np.sqrt(np.asarray(gene_tissue_model['susie_mu2']) - np.square(np.asarray(gene_tissue_model['susie_mu'])))
		global_susie_mu_sd = create_global_susie_no_flips(gene_local_to_global_mapping, local_susie_mu_sd, num_global_variants)
		global_susie_mu_sd_data.append(global_susie_mu_sd)

		# Create global mu
		global_susie_mu = create_global_susie_w_flips(gene_local_to_global_mapping, gene_flips, np.asarray(gene_tissue_model['susie_mu']), num_global_variants)
		global_susie_mu_data.append(global_susie_mu)

	global_dictionary = {'genes': gene_tissue_pairs,'tss':gene_tss_arr, 'variants': global_variant_arr, 'reference_ld':LD, 'susie_mu':global_susie_mu_data, 'susie_alpha': global_susie_alpha_data, 'susie_mu_sd': global_susie_mu_sd_data}

	return global_dictionary

def extract_regression_snp_indices(window_variant_names, hapmap3_snps, window_start, window_end):
	middle_start = window_start + 1000000.0
	middle_end = window_end + 2000000.0
	regression_snp_indices = []
	hapmap3_snp_indices = []

	for snp_index, variant_name in enumerate(window_variant_names):
		if variant_name not in hapmap3_snps:
			continue
		hapmap3_snp_indices.append(snp_index)
		variant_pos = int(variant_name.split('_')[1])
		if variant_pos >= middle_end or variant_pos < middle_start:
			continue
		regression_snp_indices.append(snp_index)
	return np.asarray(regression_snp_indices), np.asarray(hapmap3_snp_indices)

def get_tissue_names_from_gene_tissue_names_arr(gene_tissue_names):
	tissue_names = []

	for ele in gene_tissue_names:
		tissue_names.append('_'.join(ele.split('_')[1:]))
	return np.asarray(tissue_names)


def compute_gene_variance(susie_mu, susie_mu_sd, susie_alpha, ld):
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

def extract_ld_annotations_for_this_region_old(ld, susie_mu, susie_mu_sd, susie_alpha, ordered_tissue_names, region_tissue_names, variant_indices):
	# LD scores
	ld_scores = np.sum(np.square(ld),axis=0)
	# Number of variants
	num_var = len(ld_scores)

	# Create eqtl weighted scores for each gene, independently
	per_gene_eqtl_weighted_ld_scores_arr = []
	per_gene_eqtl_cross_terms = []
	per_gene_variance = []
	for gene_iter, region_tissue_name in enumerate(region_tissue_names):

		# Component level eqtl effect sizes for this gene		
		gene_component_effect_sizes = (susie_mu[gene_iter])*susie_alpha[gene_iter]

		# eQTL effect sizes for this gene
		gene_eqtl_effect_sizes = np.sum(gene_component_effect_sizes,axis=0)

		# Compute squared eqtl effect sizes for this gene
		gene_squared_eqtl_effect_sizes = np.sum((np.square(susie_mu[gene_iter]) + np.square(susie_mu_sd[gene_iter]))*susie_alpha[gene_iter],axis=0) + gene_eqtl_effect_sizes*gene_eqtl_effect_sizes - np.sum(gene_component_effect_sizes*gene_component_effect_sizes,axis=0)

		# Variance w/o modeling distribution
		#print(np.dot(np.dot(gene_eqtl_effect_sizes, ld), gene_eqtl_effect_sizes))

		# E[beta_k*beta_j]
		cross_terms = np.dot(np.reshape(gene_eqtl_effect_sizes, (num_var,1)), np.reshape(gene_eqtl_effect_sizes, (1,num_var))) - np.dot(np.transpose(gene_component_effect_sizes), gene_component_effect_sizes)


		gene_variance = compute_gene_variance(susie_mu[gene_iter], susie_mu_sd[gene_iter], susie_alpha[gene_iter], ld)
		per_gene_variance.append(gene_variance)

		gene_eqtl_weighted_ld_scores = np.sum((np.square(ld)*gene_squared_eqtl_effect_sizes),axis=1)

		per_gene_eqtl_cross_terms.append(cross_terms)
		per_gene_eqtl_weighted_ld_scores_arr.append(gene_eqtl_weighted_ld_scores[variant_indices])


	for regression_var_num, global_var_num in enumerate(variant_indices):
		directional_ld_scores = np.dot(np.reshape(ld[global_var_num,:], (num_var,1)), np.reshape(ld[global_var_num,:], (1,num_var)))
		
		for gene_iter, region_tissue_name in enumerate(region_tissue_names):
			anno_weighted_directional_ld_scores = directional_ld_scores*per_gene_eqtl_cross_terms[gene_iter]
			per_gene_eqtl_weighted_ld_scores_arr[gene_iter][regression_var_num] = per_gene_eqtl_weighted_ld_scores_arr[gene_iter][regression_var_num] + np.sum(anno_weighted_directional_ld_scores) - np.sum(np.diag(anno_weighted_directional_ld_scores))


	filtered_ld_scores = ld_scores[variant_indices]

	tissue_eqtl_ld_scores = []
	standardized_tissue_eqtl_ld_scores = []


	for tissue_name in ordered_tissue_names:
		eqtl_ld_score = np.zeros(len(filtered_ld_scores))
		standardized_eqtl_ld_score = np.zeros(len(filtered_ld_scores))
		# Get gene indices corresponding to the tissue
		gene_indices = np.where(region_tissue_names==tissue_name)[0]
		for gene_index in gene_indices:
			eqtl_ld_score = eqtl_ld_score + per_gene_eqtl_weighted_ld_scores_arr[gene_index]
			standardized_eqtl_ld_score = standardized_eqtl_ld_score + per_gene_eqtl_weighted_ld_scores_arr[gene_index]/per_gene_variance[gene_index]

		tissue_eqtl_ld_scores.append(eqtl_ld_score)
		standardized_tissue_eqtl_ld_scores.append(standardized_eqtl_ld_score)

	regression_weights = np.sum(np.square(ld)[variant_indices,:][:,variant_indices],axis=0)

	return filtered_ld_scores, np.transpose(np.asarray(tissue_eqtl_ld_scores)), np.transpose(np.asarray(standardized_tissue_eqtl_ld_scores)), regression_weights


def extract_ld_annotations_for_this_region(ld, susie_mu, susie_mu_sd, susie_alpha, ordered_tissue_names, region_tissue_names, variant_indices, hapmap3_snp_indices):
	# LD scores
	ld_scores = np.sum(np.square(ld),axis=0)
	# Number of variants
	num_var = len(ld_scores)

	# Create eqtl weighted scores for each gene, independently
	per_gene_eqtl_weighted_ld_scores_arr = []
	per_gene_variance = []
	for gene_iter, region_tissue_name in enumerate(region_tissue_names):

		# Component level eqtl effect sizes for this gene		
		gene_component_effect_sizes = (susie_mu[gene_iter])*susie_alpha[gene_iter]

		# eQTL effect sizes for this gene
		gene_eqtl_effect_sizes = np.sum(gene_component_effect_sizes,axis=0)

		# Compute squared eqtl effect sizes for this gene
		gene_squared_eqtl_effect_sizes = np.sum((np.square(susie_mu[gene_iter]) + np.square(susie_mu_sd[gene_iter]))*susie_alpha[gene_iter],axis=0) + gene_eqtl_effect_sizes*gene_eqtl_effect_sizes - np.sum(gene_component_effect_sizes*gene_component_effect_sizes,axis=0)

		# E[beta_k*beta_j]
		cross_terms = np.dot(np.reshape(gene_eqtl_effect_sizes, (num_var,1)), np.reshape(gene_eqtl_effect_sizes, (1,num_var))) - np.dot(np.transpose(gene_component_effect_sizes), gene_component_effect_sizes)


		gene_variance = compute_gene_variance(susie_mu[gene_iter], susie_mu_sd[gene_iter], susie_alpha[gene_iter], ld)
		per_gene_variance.append(gene_variance)

		gene_eqtl_weighted_ld_scores = np.sum((np.square(ld)*gene_squared_eqtl_effect_sizes),axis=1)

		# Temp
		cross_terms = cross_terms - np.diag(np.diag(cross_terms))
		cross_effects = np.sum(np.dot(cross_terms, ld[:, variant_indices])*ld[:,variant_indices],axis=0)

		per_gene_eqtl_weighted_ld_scores_arr.append(gene_eqtl_weighted_ld_scores[variant_indices] + cross_effects)


	filtered_ld_scores = ld_scores[variant_indices]

	tissue_eqtl_ld_scores = []
	standardized_tissue_eqtl_ld_scores = []


	for tissue_name in ordered_tissue_names:
		eqtl_ld_score = np.zeros(len(filtered_ld_scores))
		standardized_eqtl_ld_score = np.zeros(len(filtered_ld_scores))
		# Get gene indices corresponding to the tissue
		if len(region_tissue_names) > 0:
			gene_indices = np.where(region_tissue_names==tissue_name)[0]
			for gene_index in gene_indices:
				eqtl_ld_score = eqtl_ld_score + per_gene_eqtl_weighted_ld_scores_arr[gene_index]
				standardized_eqtl_ld_score = standardized_eqtl_ld_score + per_gene_eqtl_weighted_ld_scores_arr[gene_index]/per_gene_variance[gene_index]

		tissue_eqtl_ld_scores.append(eqtl_ld_score)
		standardized_tissue_eqtl_ld_scores.append(standardized_eqtl_ld_score)


	regression_weights = np.sum(np.square(ld)[variant_indices,:][:,hapmap3_snp_indices],axis=1)

	return filtered_ld_scores, np.transpose(np.asarray(tissue_eqtl_ld_scores)), np.transpose(np.asarray(standardized_tissue_eqtl_ld_scores)), regression_weights


def extract_tgfm_ld_score_annotation_file(window_ld_score_annotation_file, window_unstandardized_ld_score_annotation_file, tgfm_trait_agnostic_obj, regression_snp_indices, hapmap3_snp_indices, ordered_tissue_names):
	t = open(window_ld_score_annotation_file,'w')
	t.write('variant_name\tregression_weight\tld_score')
	for tissue_name in ordered_tissue_names:
		t.write('\t' + tissue_name + '_eqtl_ld_score')
	t.write('\n')

	t2 = open(window_unstandardized_ld_score_annotation_file,'w')
	t2.write('variant_name\tregression_weight\tld_score')
	for tissue_name in ordered_tissue_names:
		t2.write('\t' + tissue_name + '_eqtl_ld_score')
	t2.write('\n')


	region_tissue_names = get_tissue_names_from_gene_tissue_names_arr(tgfm_trait_agnostic_obj['genes'])

	ld_scores, eqtl_ld_scores, standardized_eqtl_ld_scores, regression_weights = extract_ld_annotations_for_this_region(tgfm_trait_agnostic_obj['reference_ld'], tgfm_trait_agnostic_obj['susie_mu'], tgfm_trait_agnostic_obj['susie_mu_sd'], tgfm_trait_agnostic_obj['susie_alpha'], ordered_tissue_names, region_tissue_names, regression_snp_indices, hapmap3_snp_indices)	

	for variant_index, variant_id in enumerate(tgfm_trait_agnostic_obj['variants'][regression_snp_indices]):
		t.write(variant_id + '\t' + str(regression_weights[variant_index]) + '\t' + str(ld_scores[variant_index]) + '\t' + '\t'.join(standardized_eqtl_ld_scores[variant_index,:].astype(str)) + '\n')
		t2.write(variant_id + '\t' + str(regression_weights[variant_index]) + '\t' + str(ld_scores[variant_index]) + '\t' + '\t'.join(eqtl_ld_scores[variant_index,:].astype(str)) + '\n')
	
	t.close()
	t2.close()


def print_study_chi_sq_vec(study_chi_sq_window_file, chi_sq_vec, window_variant_names, regression_snp_indices, study_sample_size):
	t = open(study_chi_sq_window_file, 'w')
	t.write('variant_name\tchi_square_stat\tsample_size\n')

	regression_chi_sq_vec = chi_sq_vec[regression_snp_indices]

	for variant_index, variant_id in enumerate(window_variant_names[regression_snp_indices]):
		t.write(variant_id + '\t' + str(regression_chi_sq_vec[variant_index]) + '\t' + study_sample_size + '\n')

	t.close()

def extract_rss_likelihood_trait_agnostic_data(twas_data_obj, window_name):
	# Total number of genes
	num_genes = len(twas_data_obj['genes'])
	
	# Get eqtl PMCES for each gene
	gene_eqtl_pmces = []
	for g_index in range(num_genes):
		eqtl_pmces = np.sum((twas_data_obj['susie_mu'][g_index])*(twas_data_obj['susie_alpha'][g_index]),axis=0)
		gene_eqtl_pmces.append(eqtl_pmces)

	# Get gene variances
	gene_variances = []
	for g_index in range(num_genes):
		gene_variance = compute_gene_variance(twas_data_obj['susie_mu'][g_index], twas_data_obj['susie_mu_sd'][g_index], twas_data_obj['susie_alpha'][g_index], twas_data_obj['reference_ld'])
		gene_variances.append(gene_variance)
	gene_variances = np.asarray(gene_variances)

	# Get indices of genes and variants in the middle of window
	middle_window_start = int(window_name.split(':')[1]) + 1000000
	middle_window_end = int(window_name.split(':')[1]) + 2000000
	middle_gene_indices = np.where((twas_data_obj['tss'] >= middle_window_start) & (twas_data_obj['tss'] < middle_window_end))[0]
	middle_variant_indices = []
	for ii, variant in enumerate(twas_data_obj['variants']):
		variant_pos = int(variant.split('_')[1])
		if variant_pos >= middle_window_start and variant_pos < middle_window_end:
			middle_variant_indices.append(ii)
	middle_variant_indices = np.asarray(middle_variant_indices)

	# Save to storage dictionary
	rss_trait_agnostic_data = {'gene_variances': gene_variances, 'middle_gene_indices': middle_gene_indices, 'middle_variant_indices': middle_variant_indices, 'gene_eqtl_pmces': gene_eqtl_pmces}
	return rss_trait_agnostic_data

def extract_rss_likelihood_data_for_single_trait(study_pickle_output_file, gwas_beta, gwas_beta_se, gwas_sample_size, twas_data_obj, rss_likelihood_trait_agnostic_data):
	# Initialize data obj
	data_obj = {}
	# Number of genes
	data_obj['G'] = len(twas_data_obj['genes'])
	# Number of variants
	data_obj['K'] = len(twas_data_obj['variants'])
	# Gene names
	data_obj['genes'] = np.copy(twas_data_obj['genes'])
	# Variant names
	data_obj['variants'] = np.copy(twas_data_obj['variants'])

	# Add info from trait agnostic data
	data_obj['gene_variances'] = rss_likelihood_trait_agnostic_data['gene_variances']
	data_obj['middle_gene_indices'] = rss_likelihood_trait_agnostic_data['middle_gene_indices']
	data_obj['middle_variant_indices'] = rss_likelihood_trait_agnostic_data['middle_variant_indices']
	data_obj['gene_eqtl_pmces'] = rss_likelihood_trait_agnostic_data['gene_eqtl_pmces']

	# Save gwas beta, se and sample size
	data_obj['gwas_beta'] = gwas_beta
	data_obj['gwas_beta_se'] = gwas_beta_se
	data_obj['gwas_sample_size'] = gwas_sample_size

	# Generate S matrix
	s_squared_vec = np.square(gwas_beta_se) + (np.square(gwas_beta)/gwas_sample_size)
	s_vec = np.sqrt(s_squared_vec)
	S_mat = np.diag(s_vec)
	S_inv_mat = np.diag(1.0/s_vec)
	S_inv_2_mat = np.diag(1.0/np.square(s_vec))
	# Compute (S^-1)R(S^-1) taking advantage of fact that S^-1 is a diagonal matrix
	D_mat = np.multiply(np.multiply(np.diag(S_inv_mat)[:, None], twas_data_obj['reference_ld']), np.diag(S_inv_mat))
	# Compute (S)R(S^-1) taking advantage of fact that S and S^-1 is a diagonal matrix
	srs_inv_mat = np.multiply(np.multiply(np.diag(S_mat)[:, None], twas_data_obj['reference_ld']), np.diag(S_inv_mat))
	# Generate data object containing statistics that are precomputed
	data_obj['srs_inv'] = srs_inv_mat
	data_obj['s_inv_2_diag'] = np.diag(S_inv_2_mat)
	data_obj['D_diag'] = np.diag(D_mat)

	# Precompute gene's a terms (where a term is a term in variational updates)
	precomputed_a_terms = np.zeros(data_obj['G'])
	for g_index in range(data_obj['G']):
		num_susie_components = twas_data_obj['susie_mu'][g_index].shape[0]
		for k_index in range(num_susie_components):
			precomputed_a_terms[g_index] = precomputed_a_terms[g_index] - np.sum(.5*(np.square(twas_data_obj['susie_mu'][g_index][k_index,:]) + np.square(twas_data_obj['susie_mu_sd'][g_index][k_index,:]))*np.diag(D_mat)*twas_data_obj['susie_alpha'][g_index][k_index,:])
			eqtl_component_pmces = (twas_data_obj['susie_mu'][g_index][k_index,:])*(twas_data_obj['susie_alpha'][g_index][k_index,:])
			precomputed_a_terms[g_index] = precomputed_a_terms[g_index] + .5*np.dot(np.dot(eqtl_component_pmces,D_mat), eqtl_component_pmces)
		precomputed_a_terms[g_index] = precomputed_a_terms[g_index] - .5*np.dot(np.dot(data_obj['gene_eqtl_pmces'][g_index],D_mat), data_obj['gene_eqtl_pmces'][g_index])
	data_obj['precomputed_a_terms'] = precomputed_a_terms

	# Write pickle file
	g = open(study_pickle_output_file, "wb")
	pickle.dump(data_obj, g)
	g.close()

	return


def extract_rss_likelihood_data_for_multiple_traits(window_name, tgfm_trait_agnostic_obj, gwas_beta, gwas_beta_se, study_names, study_sample_sizes, rss_likelihood_data_output_root):
	# First extract trait agnostic rss likelihood data
	rss_likelihood_trait_agnostic_data = extract_rss_likelihood_trait_agnostic_data(tgfm_trait_agnostic_obj, window_name)

	valid_studies = {'biochemistry_Cholesterol': 1, 'blood_WHITE_COUNT': 1, 'body_WHRadjBMIz':1, 'bp_DIASTOLICadjMEDz':1, 'body_BMIz':1, 'blood_RED_COUNT':1}
	# Now loop through traits and save data object containing trait specific rss likelihood data
	for study_index, study_name in enumerate(study_names):
		if study_name not in valid_studies:
			continue
		# Output file to save study-specific info
		study_pickle_output_file = rss_likelihood_data_output_root + study_name + '_data.pkl'
		# Extract and save data
		extract_rss_likelihood_data_for_single_trait(study_pickle_output_file, gwas_beta[study_index,:], gwas_beta_se[study_index,:], float(study_sample_sizes[study_index]), tgfm_trait_agnostic_obj, rss_likelihood_trait_agnostic_data)


ukkbb_window_summary_file = sys.argv[1]
hapmap3_snpid_file = sys.argv[2]
gtex_pseudotissue_file = sys.argv[3]
gtex_susie_gene_models_dir = sys.argv[4]
preprocessed_tgfm_data_dir = sys.argv[5] # Output dir
job_number = int(sys.argv[6])  # For parallelization purposes
num_jobs = int(sys.argv[7])  # For parallelization purposes

# Load in UKBB Genome wide windows file
ukbb_windows = np.loadtxt(ukkbb_window_summary_file, dtype=str,delimiter='\t')[1:,:]
# Subset to just windows in this parallel run
ukbb_windows_parr = np.array_split(ukbb_windows, num_jobs)[job_number]

# Get array of pseudotissue names
pseudotissues = get_pseudotissue_names(gtex_pseudotissue_file)

# Create dictionary from pseudotissue to array of gene models for that tissue
tissue_to_gene_model_df = create_tissue_to_gene_model_df(pseudotissues, gtex_susie_gene_models_dir)

# Get dictionary list of hapmap3 snpids
hapmap3_snps = get_dictionary_list_of_hapmap3_snpids(hapmap3_snpid_file)

start_time = time.time()
# Loop through windows
num_windows = ukbb_windows_parr.shape[0]
for window_iter in range(num_windows):
	# Extract relevent info for this window
	window_chrom_num = int(ukbb_windows_parr[window_iter,0])
	window_start = int(ukbb_windows_parr[window_iter,1])
	window_end = int(ukbb_windows_parr[window_iter,2])
	window_gwas_beta_file = ukbb_windows_parr[window_iter,5]
	window_gwas_beta_se_file = ukbb_windows_parr[window_iter,6]
	window_variant_id_file = ukbb_windows_parr[window_iter,7]
	window_study_name_file = ukbb_windows_parr[window_iter,8]
	window_variant_in_sample_ld_file = ukbb_windows_parr[window_iter,10]
	window_study_sample_size_file = ukbb_windows_parr[window_iter,11]

	# Name of window
	window_name = str(window_chrom_num) + ':' + str(window_start) + ':' + str(window_end)


	print(window_name)
	end_time = time.time()
	print((end_time-start_time)/60.0)
	start_time = end_time

	# Load in LD
	LD = np.loadtxt(window_variant_in_sample_ld_file)

	# Load variant names
	window_variant_names = np.loadtxt(window_variant_id_file, dtype=str)

	# GET HAPMAP3 SNPS and middle snps
	regression_snp_indices, hapmap3_snp_indices = extract_regression_snp_indices(window_variant_names, hapmap3_snps, window_start, window_end)
	if len(regression_snp_indices) == 0:
		continue

	# study names
	study_names = np.loadtxt(window_study_name_file, dtype=str)
	study_sample_sizes = np.loadtxt(window_study_sample_size_file,dtype=str)

	# Gwas beta and beta_se
	gwas_beta = np.loadtxt(window_gwas_beta_file)
	gwas_beta_se = np.loadtxt(window_gwas_beta_se_file)

	# Save chi-squared stats of regression snps for each study
	for study_index, study_name in enumerate(study_names):
		study_chi_sq_window_file = preprocessed_tgfm_data_dir + window_name + '_' + study_name + '_tgfm_ldscore_chi_squared_stats.txt'
		chi_sq_vec = np.square(gwas_beta[study_index,:]/gwas_beta_se[study_index,:])
		print_study_chi_sq_vec(study_chi_sq_window_file, chi_sq_vec, window_variant_names, regression_snp_indices, study_sample_sizes[study_index])

	# Organize TGFM trait agnostic input for single window
	tgfm_trait_agnostic_obj = organize_tgfm_trait_agnostic_data_for_single_window(window_name, window_chrom_num, window_start, window_end, window_variant_names, LD, pseudotissues, tissue_to_gene_model_df)


	# Save TGFM trait agnostic input data
	pkl_file = preprocessed_tgfm_data_dir + window_name + '_tgfm_trait_agnostic_data_obj.pkl'
	g = open(pkl_file, "wb")
	pickle.dump(tgfm_trait_agnostic_obj, g)
	g.close()
	'''
	pkl_file = preprocessed_tgfm_data_dir + window_name + '_tgfm_trait_agnostic_data_obj.pkl'
	f = open(pkl_file, "rb")
	tgfm_trait_agnostic_obj = pickle.load(f)
	f.close()
	'''

	# Extract RSS-likelihood formatted data
	rss_likelihood_data_output_root = preprocessed_tgfm_data_dir + window_name + '_rss_likelihood_'
	extract_rss_likelihood_data_for_multiple_traits(window_name, tgfm_trait_agnostic_obj, gwas_beta, gwas_beta_se, study_names, study_sample_sizes, rss_likelihood_data_output_root)


	'''
	# Extract ld score annotation file
	window_ld_score_annotation_file = preprocessed_tgfm_data_dir + window_name + '_tgfm_ldscore_annotation_file.txt'
	window_unstandardized_ld_score_annotation_file = preprocessed_tgfm_data_dir + window_name + '_tgfm_ldscore_not_standardized_annotation_file.txt'
	extract_tgfm_ld_score_annotation_file(window_ld_score_annotation_file, window_unstandardized_ld_score_annotation_file, tgfm_trait_agnostic_obj, regression_snp_indices, hapmap3_snp_indices, pseudotissues)
	'''