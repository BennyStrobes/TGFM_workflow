import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pyreadr
import numpy as np 
import os
import sys
import pdb
import pickle
import time
import gzip




def get_pseudotissue_names(gtex_pseudotissue_file, remove_testis=False):
	f = open(gtex_pseudotissue_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if remove_testis and data[0] == 'Testis':
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

def create_tissue_to_gene_model_df(pseudotissues, gtex_susie_gene_models_dir, gene_type):
	mapping = {}
	for pseudotissue in pseudotissues:
		gene_model_summary_file = gtex_susie_gene_models_dir + pseudotissue + '/' + pseudotissue + '_' + gene_type + '_pos_file.txt'
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
	return global_susie_data, local_susie_data

def create_global_susie_w_flips(gene_local_to_global_mapping, gene_flips, local_susie_data, num_global_variants):
	for var_index, flip_value in enumerate(gene_flips):
		if flip_value == -1.0:
			local_susie_data[:, var_index] = local_susie_data[:, var_index]*-1.0
	global_susie_data, local_susie_data = create_global_susie_no_flips(gene_local_to_global_mapping, local_susie_data, num_global_variants)
	return global_susie_data, local_susie_data

def extract_annotation_mat_for_window_variants(short_variant_to_annotation_vector, global_variant_arr):
	annotation_arr = []
	for full_variant in global_variant_arr:
		variant_info = full_variant.split('_')
		short_variant_name = variant_info[0] + ':' + variant_info[1]
		annotation_vec = short_variant_to_annotation_vector[short_variant_name]
		annotation_arr.append(annotation_vec)
	annotation_arr = np.asarray(annotation_arr)
	return annotation_arr

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


def extract_variant_positions_from_variant_array(global_variant_arr):
	variant_positions = []
	for variant_name in global_variant_arr:
		variant_positions.append(int(variant_name.split('_')[1]))
	return np.asarray(variant_positions)

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

def organize_tgfm_trait_agnostic_data_for_single_window(window_name, window_chrom_num, window_start, window_end, global_variant_arr, LD, pseudotissues, tissue_to_gene_model_df, n_bs):
	# Get list of gene_tissue pairs and weight files for this window
	gene_tissue_pairs, gene_tissue_pair_weight_files, gene_tss_arr = get_list_of_gene_tissue_pairs_and_weight_files_for_window(window_chrom_num, window_start, window_end, pseudotissues, tissue_to_gene_model_df)

	# Get number of total variants
	num_global_variants = len(global_variant_arr)

	# Create dictionary mapping from variant name to global-variant arr position
	variant_name_to_global_variant_arr_pos = create_mapping_from_variant_name_to_global_variant_arr_pos(global_variant_arr)

	# Extract variant positions from variant array
	variant_positions = extract_variant_positions_from_variant_array(global_variant_arr)
	

	# Initialize arrays to keep track of eqtl-susie data
	global_weight_vectors = []
	gene_variances = []
	gene_full_variances = []
	global_components = []

	bs_arr = []

	# Loop through gene-tissue pairs
	print(len(gene_tissue_pairs))
	for gene_tissue_index, gene_tissue_name in enumerate(gene_tissue_pairs):
		#print(gene_tissue_index)
		# Get gene-tissue weight file
		gene_tissue_weight_file = gene_tissue_pair_weight_files[gene_tissue_index]

		# Load gene-tissue model from gene-tissue weight file
		gene_tissue_model = pyreadr.read_r(gene_tissue_weight_file)

		component_bool = np.asarray(gene_tissue_model['component_bool'])[0,0]
		if component_bool == False:
			global_components.append(np.asarray([]))
		else:
			gene_components = np.asarray(gene_tissue_model['susie_cs'])[:,0] - 1
			global_components.append(gene_components)

		# Create mapping from gene coordinates to global coordinates
		gene_variants = np.asarray(gene_tissue_model['variant_names'])[:,0]
		gene_local_to_global_mapping, gene_flips = create_gene_local_to_global_mapping(variant_name_to_global_variant_arr_pos, gene_variants)

		# Create global susie alpha
		global_susie_alpha, local_susie_alpha = create_global_susie_no_flips(gene_local_to_global_mapping, np.asarray(gene_tissue_model['susie_alpha']), num_global_variants)
		#global_susie_alpha_data.append(global_susie_alpha)
		# Create global susie mu_sd
		local_susie_mu_sd_raw = np.sqrt(np.asarray(gene_tissue_model['susie_mu2']) - np.square(np.asarray(gene_tissue_model['susie_mu'])))
		global_susie_mu_sd, local_susie_mu_sd = create_global_susie_no_flips(gene_local_to_global_mapping, local_susie_mu_sd_raw, num_global_variants)
		#global_susie_mu_sd_data.append(global_susie_mu_sd)

		# Create global mu
		global_susie_mu, local_susie_mu = create_global_susie_w_flips(gene_local_to_global_mapping, gene_flips, np.asarray(gene_tissue_model['susie_mu']), num_global_variants)
		#global_susie_mu_data.append(global_susie_mu)

		# Create susie PMCES
		local_susie_pmces = np.sum(local_susie_mu*local_susie_alpha,axis=0)
		global_susie_pmces = np.sum(global_susie_mu*global_susie_alpha,axis=0)
		#global_susie_pmces2 = np.copy(global_susie_pmces)*0.0
		#global_susie_pmces2[gene_local_to_global_mapping] = local_susie_pmces

		# Construct local LD
		local_LD = LD[gene_local_to_global_mapping,:][:, gene_local_to_global_mapping]

		# Compute gene variance
		#gene_variance_alt = np.dot(np.dot(global_susie_pmces, LD), global_susie_pmces)
		gene_variance  = np.dot(np.dot(local_susie_pmces, local_LD), local_susie_pmces)

		# Calculate full gene variance (variance according to susie posterior distribution)
		full_gene_variance = calculate_gene_variance_according_to_susie_distribution(local_susie_mu, local_susie_alpha, local_susie_mu_sd, local_LD)

		# Standardize susie PMCES
		global_susie_pmces_standardized = global_susie_pmces/np.sqrt(gene_variance)

		# Store data
		global_weight_vectors.append(global_susie_pmces_standardized)
		gene_variances.append(gene_variance)
		gene_full_variances.append(full_gene_variance)


		# Run sampling analysis
		n_snps = local_susie_mu.shape[1]
		bs_alpha_effects = np.zeros((n_snps, n_bs))

		# For number of samplese (n_bs), sample
		for bs_iter in range(n_bs):
			susie_sampled_eqtl_effects = sample_eqtl_effects_from_susie_distribution(local_susie_mu, local_susie_alpha, np.square(local_susie_mu_sd))
			bs_alpha_effects[:, bs_iter] = susie_sampled_eqtl_effects
		# Standardize and store data		
		bs_gene_variances = np.diag(np.dot(np.dot(np.transpose(bs_alpha_effects), local_LD), bs_alpha_effects))
		bs_std_alpha_effects = bs_alpha_effects/np.sqrt(bs_gene_variances)
		bs_arr.append((bs_std_alpha_effects, gene_local_to_global_mapping))

	# Organize sample analysis
	sparse_sampled_gene_eqtl_pmces = []
	for bs_iter in range(n_bs):
		eqtl_mat = []
		for gene_itera, bs_eqtl_tuple in enumerate(bs_arr):
			eqtl_gene_window = bs_eqtl_tuple[0][:, bs_iter]
			eqtl_indices = bs_eqtl_tuple[1]
			#eqtl_indices = np.arange(len(boolean_indices))[boolean_indices]
			for ii, eqtl_effect in enumerate(eqtl_gene_window):
				if eqtl_effect == 0.0:
					continue
				eqtl_index = eqtl_indices[ii]
				eqtl_mat.append(np.asarray([gene_itera, eqtl_index, eqtl_effect]))
		eqtl_mat = np.asarray(eqtl_mat)
		sparse_sampled_gene_eqtl_pmces.append(eqtl_mat)



	# Extract annotation mat for these variants
	#annotation_mat = extract_annotation_mat_for_window_variants(short_variant_to_annotation_vector, global_variant_arr)

	# Get middle indices
	window_middle_start = window_start + 1000000
	window_middle_end = window_end - 1000000
	middle_gene_indices = np.where((gene_tss_arr >= window_middle_start) & (gene_tss_arr < window_middle_end))[0]
	middle_variant_indices = np.where((variant_positions >= window_middle_start) & (variant_positions < window_middle_end))[0]

	# Save to global object
	global_dictionary = {'genes': gene_tissue_pairs,'tss':gene_tss_arr, 'valid_susie_components':global_components, 'variants': global_variant_arr, 'varriant_positions': variant_positions, 'gene_eqtl_pmces': np.asarray(global_weight_vectors), 'gene_variances': np.asarray(gene_variances), 'full_gene_variances': np.asarray(gene_full_variances), 'sparse_sampled_gene_eqtl_pmces':sparse_sampled_gene_eqtl_pmces, 'middle_gene_indices': middle_gene_indices, 'middle_variant_indices': middle_variant_indices}

	return global_dictionary

def extract_regression_snp_indices(window_variant_names, hapmap3_snps, window_start, window_end):
	middle_start = window_start + 1000000.0
	middle_end = window_start + 2000000.0
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

def extract_ld_annotations_for_this_region_v2(ld, susie_mu, susie_mu_sd, susie_alpha, ordered_tissue_names, region_tissue_names, variant_indices, hapmap3_snp_indices):
	# LD scores
	ld_scores = np.sum(np.square(ld),axis=0)
	# Number of variants
	num_var = len(ld_scores)

	# Create eqtl weighted scores for each gene, independently
	per_gene_eqtl_weighted_ld_scores_arr = []
	per_gene_variance = []
	for gene_iter, region_tissue_name in enumerate(region_tissue_names):

		gene_variance_init = compute_gene_variance(susie_mu[gene_iter], susie_mu_sd[gene_iter], susie_alpha[gene_iter], ld)
		susie_mu[gene_iter] = susie_mu[gene_iter]/np.sqrt(gene_variance_init)
		susie_mu_sd[gene_iter] = susie_mu_sd[gene_iter]/np.sqrt(gene_variance_init)

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


def extract_ld_annotations_for_this_region_point_estimate_gene_models(ld, susie_mu, susie_alpha, ordered_tissue_names, region_tissue_names, variant_indices, hapmap3_snp_indices):
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
		#gene_squared_eqtl_effect_sizes = np.sum((np.square(susie_mu[gene_iter]) + np.square(susie_mu_sd[gene_iter]))*susie_alpha[gene_iter],axis=0) + gene_eqtl_effect_sizes*gene_eqtl_effect_sizes - np.sum(gene_component_effect_sizes*gene_component_effect_sizes,axis=0)
		gene_squared_eqtl_effect_sizes = np.square(gene_eqtl_effect_sizes)

		# E[beta_k*beta_j]
		#cross_terms = np.dot(np.reshape(gene_eqtl_effect_sizes, (num_var,1)), np.reshape(gene_eqtl_effect_sizes, (1,num_var))) - np.dot(np.transpose(gene_component_effect_sizes), gene_component_effect_sizes)
		cross_terms = np.dot(np.reshape(gene_eqtl_effect_sizes, (num_var,1)), np.reshape(gene_eqtl_effect_sizes, (1,num_var)))

		#gene_variance = compute_gene_variance(susie_mu[gene_iter], susie_mu_sd[gene_iter], susie_alpha[gene_iter], ld)
		gene_variance = np.dot(np.dot(gene_eqtl_effect_sizes, ld), gene_eqtl_effect_sizes)
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


def extract_tgfm_ld_score_annotation_file_point_estimate_gene_model(window_ld_score_annotation_file, window_unstandardized_ld_score_annotation_file, tgfm_trait_agnostic_obj, regression_snp_indices, hapmap3_snp_indices, ordered_tissue_names):
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

	ld_scores, eqtl_ld_scores, standardized_eqtl_ld_scores, regression_weights = extract_ld_annotations_for_this_region_point_estimate_gene_models(tgfm_trait_agnostic_obj['reference_ld'], tgfm_trait_agnostic_obj['susie_mu'], tgfm_trait_agnostic_obj['susie_alpha'], ordered_tissue_names, region_tissue_names, regression_snp_indices, hapmap3_snp_indices)

	for variant_index, variant_id in enumerate(tgfm_trait_agnostic_obj['variants'][regression_snp_indices]):
		t.write(variant_id + '\t' + str(regression_weights[variant_index]) + '\t' + str(ld_scores[variant_index]) + '\t' + '\t'.join(standardized_eqtl_ld_scores[variant_index,:].astype(str)) + '\n')
		t2.write(variant_id + '\t' + str(regression_weights[variant_index]) + '\t' + str(ld_scores[variant_index]) + '\t' + '\t'.join(eqtl_ld_scores[variant_index,:].astype(str)) + '\n')
	
	t.close()
	t2.close()


def extract_indices_of_genes_with_components(global_components):
	indices = []
	for ii, temp_arr in enumerate(global_components):
		if len(temp_arr) > 0:
			indices.append(ii)
	return np.asarray(indices)

def filter_list_to_indices(listy, indices):
	new_listy = []
	for index in indices:
		new_listy.append(listy[index])
	return new_listy

def filter_array_list_to_indices(listy, row_indices, indices):
	new_listy = []
	for index in indices:
		new_listy.append(listy[index][row_indices[index],:])
	return new_listy


def extract_tgfm_ld_score_annotation_file_component_only_gene_model(window_ld_score_annotation_file, window_unstandardized_ld_score_annotation_file, tgfm_trait_agnostic_obj, regression_snp_indices, hapmap3_snp_indices, ordered_tissue_names):
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

	indices_of_genes_with_components = extract_indices_of_genes_with_components(tgfm_trait_agnostic_obj['valid_susie_components'])

	if len(indices_of_genes_with_components) > 0:
		region_tissue_names = get_tissue_names_from_gene_tissue_names_arr(tgfm_trait_agnostic_obj['genes'][indices_of_genes_with_components])
	else:
		region_tissue_names = np.asarray([])

	# Filter gene models to just predicted expression in a component
	filtered_susie_mu = filter_array_list_to_indices(tgfm_trait_agnostic_obj['susie_mu'], tgfm_trait_agnostic_obj['valid_susie_components'], indices_of_genes_with_components)
	filtered_susie_mu_sd = filter_array_list_to_indices(tgfm_trait_agnostic_obj['susie_mu_sd'], tgfm_trait_agnostic_obj['valid_susie_components'], indices_of_genes_with_components)
	filtered_susie_alpha = filter_array_list_to_indices(tgfm_trait_agnostic_obj['susie_alpha'], tgfm_trait_agnostic_obj['valid_susie_components'], indices_of_genes_with_components)

	# Compute LD-scores
	ld_scores, eqtl_ld_scores, standardized_eqtl_ld_scores, regression_weights = extract_ld_annotations_for_this_region(tgfm_trait_agnostic_obj['reference_ld'], filtered_susie_mu, filtered_susie_mu_sd, filtered_susie_alpha, ordered_tissue_names, region_tissue_names, regression_snp_indices, hapmap3_snp_indices)	

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

def extract_rss_likelihood_trait_agnostic_data(twas_data_obj, window_name, shared_pickle_output_file):
	# Total number of genes
	num_genes = len(twas_data_obj['genes'])
	
	# Get eqtl PMCES for each gene
	gene_eqtl_pmces = []
	for g_index in range(num_genes):
		eqtl_pmces = np.sum((twas_data_obj['susie_mu'][g_index])*(twas_data_obj['susie_alpha'][g_index]),axis=0)
		gene_eqtl_pmces.append(eqtl_pmces)

	# Get gene variances
	if num_genes != 0:
		expression_covariance = np.dot(np.dot(gene_eqtl_pmces, twas_data_obj['reference_ld']), np.transpose(gene_eqtl_pmces))
		gene_variances = []
		for g_index in range(num_genes):
			gene_variance = compute_gene_variance(twas_data_obj['susie_mu'][g_index], twas_data_obj['susie_mu_sd'][g_index], twas_data_obj['susie_alpha'][g_index], twas_data_obj['reference_ld'])
			gene_variances.append(gene_variance)
			expression_covariance[g_index, g_index] = gene_variance
		gene_variances = np.asarray(gene_variances)
		# normalize covairance
		dd = np.diag(1.0/np.sqrt(np.diag(expression_covariance)))
		ge_ld = np.dot(np.dot(dd, expression_covariance),dd)
	else:
		gene_variances = np.asarray([])
		ge_ld = np.zeros((0,0))


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

	# Number of genes
	rss_trait_agnostic_data['G'] = len(twas_data_obj['genes'])
	# Number of variants
	rss_trait_agnostic_data['K'] = len(twas_data_obj['variants'])
	# Gene names
	rss_trait_agnostic_data['genes'] = np.copy(twas_data_obj['genes'])
	# Variant names
	rss_trait_agnostic_data['variants'] = np.copy(twas_data_obj['variants'])
	# LD
	rss_trait_agnostic_data['reference_ld'] = np.copy(twas_data_obj['reference_ld'])
	# Gene expression LD
	rss_trait_agnostic_data['ge_ld'] = ge_ld
	# gene eqtl pmces in np form
	rss_trait_agnostic_data['gene_eqtl_pmces_np'] = np.asarray(rss_trait_agnostic_data['gene_eqtl_pmces'])

	# Write pickle file
	g = open(shared_pickle_output_file, "wb")
	pickle.dump(rss_trait_agnostic_data, g)
	g.close()

	return rss_trait_agnostic_data

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

	return beta_scaled, beta_se_scaled

def extract_rss_likelihood_data_for_single_trait(study_pickle_output_file, gwas_beta, gwas_beta_se, gwas_sample_size, twas_data_obj, rss_likelihood_trait_agnostic_data):
	# Initialize data obj
	data_obj = {}

	# Save gwas beta, se and sample size
	data_obj['gwas_beta'] = gwas_beta
	data_obj['gwas_beta_se'] = gwas_beta_se
	data_obj['gwas_sample_size'] = gwas_sample_size

	# Standardize betas and beta_se
	beta_scaled, beta_se_scaled, XtX = convert_to_standardized_summary_statistics(gwas_beta, gwas_beta_se, gwas_sample_size, twas_data_obj['reference_ld'])
	data_obj['standardized_gwas_beta'] = beta_scaled
	data_obj['standardized_gwas_beta_se'] = beta_se_scaled
	
	# Generate S matrix
	#s_squared_vec = np.square(gwas_beta_se) + (np.square(gwas_beta)/gwas_sample_size)
	s_squared_vec = np.square(beta_se_scaled) + (np.square(beta_scaled)/gwas_sample_size)
	s_vec = np.sqrt(s_squared_vec)
	S_mat = np.diag(s_vec)
	S_inv_mat = np.diag(1.0/s_vec)
	S_inv_2_mat = np.diag(1.0/np.square(s_vec))
	# Compute (S^-1)R(S^-1) taking advantage of fact that S^-1 is a diagonal matrix
	D_mat = np.multiply(np.multiply(np.diag(S_inv_mat)[:, None], twas_data_obj['reference_ld']), np.diag(S_inv_mat))
	# Compute (S)R(S^-1) taking advantage of fact that S and S^-1 is a diagonal matrix
	#srs_inv_mat = np.multiply(np.multiply(np.diag(S_mat)[:, None], twas_data_obj['reference_ld']), np.diag(S_inv_mat))
	# Generate data object containing statistics that are precomputed
	#data_obj['srs_inv'] = srs_inv_mat
	data_obj['s_diag'] = np.diag(S_mat)
	data_obj['s_inv_2_diag'] = np.diag(S_inv_2_mat)
	data_obj['D_diag'] = np.diag(D_mat)

	# Precompute gene-gene interaction terms (where gene-gene interaction terms are used in variational updates)
	# Note that diagonal elements are currently wrong as they do not account for probability distribution governing susie eqtl effect sizes
	# Next section (precomputing a_terms will fix this issue)
	if rss_likelihood_trait_agnostic_data['G'] > 0:
		precomputed_gene_gene_terms = -np.dot(np.dot(rss_likelihood_trait_agnostic_data['gene_eqtl_pmces_np'], D_mat),np.transpose(rss_likelihood_trait_agnostic_data['gene_eqtl_pmces_np']))
	else:
		precomputed_gene_gene_terms = np.zeros((rss_likelihood_trait_agnostic_data['G'], rss_likelihood_trait_agnostic_data['G']))

	# Precompute gene's a terms (where a term is a term in variational updates)
	precomputed_a_terms = np.zeros(rss_likelihood_trait_agnostic_data['G'])
	for g_index in range(rss_likelihood_trait_agnostic_data['G']):
		precomputed_a_terms[g_index] = (precomputed_gene_gene_terms[g_index,g_index]/2.0) 
		precomputed_a_terms[g_index] = precomputed_a_terms[g_index] - (np.sum(np.sum((np.square(twas_data_obj['susie_mu'][g_index]) + np.square(twas_data_obj['susie_mu_sd'][g_index]))*twas_data_obj['susie_alpha'][g_index],axis=0)*np.diag(D_mat))/2.0)
		eqtl_component_pmces = twas_data_obj['susie_mu'][g_index]*twas_data_obj['susie_alpha'][g_index]
		precomputed_a_terms[g_index] = precomputed_a_terms[g_index] + np.sum(np.diag(np.dot(np.dot(eqtl_component_pmces,D_mat), np.transpose(eqtl_component_pmces))))/2.0
		precomputed_gene_gene_terms[g_index, g_index] = 2.0*precomputed_a_terms[g_index]

	'''
	precomputed_a_terms2 = np.zeros(rss_likelihood_trait_agnostic_data['G'])
	precomputed_gene_gene_terms2 = np.zeros((rss_likelihood_trait_agnostic_data['G'], rss_likelihood_trait_agnostic_data['G']))
	for g_index in range(rss_likelihood_trait_agnostic_data['G']):
		print(g_index)
		num_susie_components = twas_data_obj['susie_mu'][g_index].shape[0]
		for k_index in range(num_susie_components):
			precomputed_a_terms2[g_index] = precomputed_a_terms2[g_index] - np.sum(.5*(np.square(twas_data_obj['susie_mu'][g_index][k_index,:]) + np.square(twas_data_obj['susie_mu_sd'][g_index][k_index,:]))*np.diag(D_mat)*twas_data_obj['susie_alpha'][g_index][k_index,:])
			eqtl_component_pmces = (twas_data_obj['susie_mu'][g_index][k_index,:])*(twas_data_obj['susie_alpha'][g_index][k_index,:])
			precomputed_a_terms2[g_index] = precomputed_a_terms2[g_index] + .5*np.dot(np.dot(eqtl_component_pmces,D_mat), eqtl_component_pmces)
		precomputed_a_terms2[g_index] = precomputed_a_terms2[g_index] - .5*np.dot(np.dot(rss_likelihood_trait_agnostic_data['gene_eqtl_pmces'][g_index],D_mat), rss_likelihood_trait_agnostic_data['gene_eqtl_pmces'][g_index])
			

		# Cross terms
		precomputed_gene_gene_terms2[g_index, g_index] = 2.0*precomputed_a_terms2[g_index]
		for g_prime_index in range((g_index+1), rss_likelihood_trait_agnostic_data['G']):
			temp_val = -np.dot(np.dot(rss_likelihood_trait_agnostic_data['gene_eqtl_pmces'][g_prime_index],D_mat), rss_likelihood_trait_agnostic_data['gene_eqtl_pmces'][g_index])
			precomputed_gene_gene_terms2[g_index, g_prime_index] = temp_val
			precomputed_gene_gene_terms2[g_prime_index, g_index] = temp_val
	'''

	# Add to data obj
	data_obj['precomputed_a_terms'] = precomputed_a_terms
	data_obj['precomputed_gene_gene_terms'] = precomputed_gene_gene_terms

	# Write pickle file
	g = open(study_pickle_output_file, "wb")
	pickle.dump(data_obj, g)
	g.close()

	return


def extract_rss_likelihood_data_for_multiple_traits(window_name, tgfm_trait_agnostic_obj, gwas_beta, gwas_beta_se, study_names, study_sample_sizes, rss_likelihood_data_output_root, standardize):
	# Standardize eQTL effects
	if standardize == True:
		for g_index in range(len(tgfm_trait_agnostic_obj['susie_mu'])):
			gene_variance = compute_gene_variance(tgfm_trait_agnostic_obj['susie_mu'][g_index], tgfm_trait_agnostic_obj['susie_mu_sd'][g_index], tgfm_trait_agnostic_obj['susie_alpha'][g_index], tgfm_trait_agnostic_obj['reference_ld'])
			tgfm_trait_agnostic_obj['susie_mu'][g_index] = tgfm_trait_agnostic_obj['susie_mu'][g_index]/np.sqrt(gene_variance)
			tgfm_trait_agnostic_obj['susie_mu_sd'][g_index] = tgfm_trait_agnostic_obj['susie_mu_sd'][g_index]/np.sqrt(gene_variance)

	# First extract trait agnostic rss likelihood data
	if standardize == False:
		shared_pickle_output_file = rss_likelihood_data_output_root + 'shared_data.pkl'
	else:
		shared_pickle_output_file = rss_likelihood_data_output_root + 'shared_standardized_data.pkl'
	rss_likelihood_trait_agnostic_data = extract_rss_likelihood_trait_agnostic_data(tgfm_trait_agnostic_obj, window_name, shared_pickle_output_file)

	# Now loop through traits and save data object containing trait specific rss likelihood data
	for study_index, study_name in enumerate(study_names):
		# Output file to save study-specific info
		if standardize == False:
			study_pickle_output_file = rss_likelihood_data_output_root + study_name + '_data.pkl'
		else:
			study_pickle_output_file = rss_likelihood_data_output_root + study_name + '_standardized_data.pkl'
		# Extract and save data
		extract_rss_likelihood_data_for_single_trait(study_pickle_output_file, gwas_beta[study_index,:], gwas_beta_se[study_index,:], float(study_sample_sizes[study_index]), tgfm_trait_agnostic_obj, rss_likelihood_trait_agnostic_data)


def get_dictionary_list_of_all_regression_snps(hapmap3_snps, ukbb_windows):
	all_regression_snps = {}
	num_windows = ukbb_windows.shape[0]
	for window_iter in range(num_windows):
		# Extract relevent info for this window
		window_start = int(ukbb_windows[window_iter,1])
		window_end = int(ukbb_windows[window_iter,2])
		window_variant_id_file = ukbb_windows[window_iter,7]
		# Load variant names
		window_variant_names = np.loadtxt(window_variant_id_file, dtype=str)
		regression_snp_indices, hapmap3_snp_indices = extract_regression_snp_indices(window_variant_names, hapmap3_snps, window_start, window_end)
		if len(regression_snp_indices) == 0:
			continue
		regression_snps = window_variant_names[regression_snp_indices]
		for regression_snp in regression_snps:
			regression_snp_info = regression_snp.split('_')
			regression_snp_alt = regression_snp_info[0] + '_' + regression_snp_info[1] + '_' + regression_snp_info[3] + '_' + regression_snp_info[2]
			all_regression_snps[regression_snp] = 1
			all_regression_snps[regression_snp_alt] = 1
	return all_regression_snps

def get_chromosome_names_used_in_this_analysis(ukbb_windows_parr):
	num_windows = ukbb_windows_parr.shape[0]
	chrom_nums = []
	for window_iter in range(num_windows):
		chrom_num = ukbb_windows_parr[window_iter,0]
		chrom_nums.append(chrom_num)
	return np.sort(np.unique(np.asarray(chrom_nums)))

def create_mapping_from_short_variant_name_to_annotation_vec(used_chrom_arr, annotation_dir):
	dicti_mapping = {}
	for chrom_string in used_chrom_arr:
		chrom_annot_file = annotation_dir + 'baselineLD_no_qtl.' + chrom_string + '.annot'
		f = open(chrom_annot_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			short_variant_name = 'chr' + data[0] + ':' + data[1]
			if short_variant_name in dicti_mapping:
				print('assumption error')
			annotation_vec = np.asarray(data[4:]).astype(float)
			dicti_mapping[short_variant_name] = annotation_vec
	return dicti_mapping

def extract_uniform_log_prior_probabilities(gene_names, variant_names, study_names):
	# General info
	n_genes = len(gene_names)
	n_variants = len(variant_names)
	n_studies = len(study_names)

	# Initialize prior mats
	gene_prior = np.zeros((n_studies, n_genes))
	variant_prior = np.zeros((n_studies, n_variants))

	# Calculate prior prob for each element under uniform prior
	prior_prob = 1.0/(n_genes + n_variants)

	gene_prior = gene_prior + prior_prob
	variant_prior = variant_prior + prior_prob

	return np.log(gene_prior), np.log(variant_prior)

def extract_sparse_variant_gene_tissue_log_prior_probabilities(gene_names, variant_names, study_names, tgfm_sldsc_results_dir):
	# General info
	n_genes = len(gene_names)
	n_variants = len(variant_names)
	n_studies = len(study_names)

	# Initialize prior mats
	gene_prior = np.zeros((n_studies, n_genes))
	variant_prior = np.zeros((n_studies, n_variants))

	for study_iter, study_name in enumerate(study_names):
		# Extract per gene and per variant h2
		tglr_results_file = tgfm_sldsc_results_dir + study_name + '_baseline_no_qtl_component_gene_no_testis_pmces_gene_adj_ld_scores_0.5_sparse_ard_eqtl_coefficients_mv_update_avg_per_snp_and_gene_tissue_h2.txt'
		tmp = np.loadtxt(tglr_results_file,dtype=str,delimiter='\t')

		per_variant_h2 = float(tmp[1,1])
		per_gene_tissue_h2 = tmp[2:,1].astype(float)
		tissue_names = tmp[2:,0]
		tissue_to_gene_h2 = {}

		# Deal with case where estimated per element h2 is < 0
		if per_variant_h2 <= 1e-15:
			per_variant_h2 = 1e-15
		per_gene_tissue_h2[per_gene_tissue_h2 <= 1e-15] = 1e-15

		# create mapping from tissue to gene
		for tissue_iter, tissue_name in enumerate(tissue_names):
			tissue_to_gene_h2[tissue_name] = per_gene_tissue_h2[tissue_iter]


		# Fill in genetic elements h2s
		gene_h2s = np.zeros(n_genes)
		variant_h2s = np.zeros(n_variants) + per_variant_h2

		for gene_iter, gene_name in enumerate(gene_names):
			tissue_name = '_'.join(gene_name.split('_')[1:])
			gene_h2s[gene_iter] = tissue_to_gene_h2[tissue_name]

		# Calculate h2 from region (sum of element h2s)
		region_h2 = np.sum(gene_h2s) + np.sum(variant_h2s)


		# Normalize to get probabilities
		gene_probs = gene_h2s/region_h2
		variant_probs = variant_h2s/region_h2

		# Add to matrix
		gene_prior[study_iter,:] = gene_probs
		variant_prior[study_iter,:] = variant_probs
	return np.log(gene_prior), np.log(variant_prior)


def extract_variant_gene_log_prior_probabilities(gene_names, variant_names, study_names, tgfm_sldsc_results_dir):
	# General info
	n_genes = len(gene_names)
	n_variants = len(variant_names)
	n_studies = len(study_names)

	# Initialize prior mats
	gene_prior = np.zeros((n_studies, n_genes))
	variant_prior = np.zeros((n_studies, n_variants))

	for study_iter, study_name in enumerate(study_names):
		# Extract per gene and per variant h2
		tglr_results_file = tgfm_sldsc_results_dir + study_name + '_baseline_no_qtl_component_gene_no_testis_pmces_gene_adj_ld_scores_avg_per_snp_and_gene_h2.txt'
		tmp = np.loadtxt(tglr_results_file,dtype=str,delimiter='\t')
		per_variant_h2 = float(tmp[1,1])
		per_gene_h2 = float(tmp[2,1])
		# Deal with case where estimated per element h2 is < 0
		if per_variant_h2 <= 1e-15:
			per_variant_h2 = 1e-15
		if per_gene_h2 <= 1e-15:
			per_gene_h2 = 1e-15

		# Fill in genetic elements h2s
		gene_h2s = np.zeros(n_genes) + per_gene_h2
		variant_h2s = np.zeros(n_variants) + per_variant_h2

		# Calculate h2 from region (sum of element h2s)
		region_h2 = np.sum(gene_h2s) + np.sum(variant_h2s)


		# Normalize to get probabilities
		gene_probs = gene_h2s/region_h2
		variant_probs = variant_h2s/region_h2

		# Add to matrix
		gene_prior[study_iter,:] = gene_probs
		variant_prior[study_iter,:] = variant_probs

	return np.log(gene_prior), np.log(variant_prior)



ukkbb_window_summary_file = sys.argv[1]
gtex_pseudotissue_file = sys.argv[2]
gtex_susie_gene_models_dir = sys.argv[3]
preprocessed_tgfm_data_dir = sys.argv[4] # Output dir
job_number = int(sys.argv[5])  # For parallelization purposes
num_jobs = int(sys.argv[6])  # For parallelization purposes
gene_type = sys.argv[7]

# Number of bootstrap samples
n_bs = 100

# Append gene type to output root
preprocessed_tgfm_data_dir = preprocessed_tgfm_data_dir + gene_type + '_'

# Load in UKBB Genome wide windows file
ukbb_windows = np.loadtxt(ukkbb_window_summary_file, dtype=str,delimiter='\t')[1:,:]
# Subset to just windows in this parallel run
ukbb_windows_parr = np.array_split(ukbb_windows, num_jobs)[job_number]

print(ukbb_windows_parr.shape)

# Get array of pseudotissue names
pseudotissues = get_pseudotissue_names(gtex_pseudotissue_file, remove_testis=False)


# Create dictionary from pseudotissue to array of gene models for that tissue
tissue_to_gene_model_df = create_tissue_to_gene_model_df(pseudotissues, gtex_susie_gene_models_dir, gene_type)

# Get dictionary list of hapmap3 snpids
#hapmap3_snps = get_dictionary_list_of_hapmap3_snpids(hapmap3_snpid_file)
# This is all hapmap3 snps that fall in the middle of a window that is tested
#all_regression_snps = get_dictionary_list_of_all_regression_snps(hapmap3_snps, ukbb_windows)

# Get chromosome nums in this parallel window
used_chrom_arr = get_chromosome_names_used_in_this_analysis(ukbb_windows_parr)
print(used_chrom_arr)
#short_variant_name_to_annotation_vec = create_mapping_from_short_variant_name_to_annotation_vec(used_chrom_arr, annotation_dir)



# Open outputful summarizing TGFM input (one line for each window)
tgfm_input_data_summary_file = preprocessed_tgfm_data_dir + 'tgfm_input_data_summary_' + str(job_number) + '_' + str(num_jobs) + '.txt'
t = open(tgfm_input_data_summary_file,'w')
# Write header
t.write('window_name\tLD_npy_file\tTGFM_input_pkl\tTGFM_trait_input_pkl\n')



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
	window_1kg_genotype_file = ukbb_windows_parr[window_iter,9]
	window_variant_in_sample_ld_file = ukbb_windows_parr[window_iter,10]
	window_study_sample_size_file = ukbb_windows_parr[window_iter,11]

	# Name of window
	window_name = str(window_chrom_num) + ':' + str(window_start) + ':' + str(window_end)

	print(window_name)
	end_time = time.time()
	print((end_time-start_time)/60.0)
	start_time = end_time

	# Load in LD
	LD = np.load(window_variant_in_sample_ld_file)

	# Load variant names
	window_variant_names = np.loadtxt(window_variant_id_file, dtype=str)

	# study names
	study_names = np.loadtxt(window_study_name_file, dtype=str)
	study_sample_sizes = np.loadtxt(window_study_sample_size_file,dtype=str)

	# Gwas beta and beta_se
	gwas_beta = np.loadtxt(window_gwas_beta_file)
	gwas_beta_se = np.loadtxt(window_gwas_beta_se_file)

	# Organize TGFM trait agnostic input for single window
	tgfm_data_obj = organize_tgfm_trait_agnostic_data_for_single_window(window_name, window_chrom_num, window_start, window_end, window_variant_names, LD, pseudotissues, tissue_to_gene_model_df, n_bs)

	# Save TGFM input data object
	pkl_file = preprocessed_tgfm_data_dir + window_name + '_tgfm_trait_agnostic_input_data_obj.pkl'
	g = open(pkl_file, "wb")
	pickle.dump(tgfm_data_obj, g)
	g.close()

	##***********## TEMP #####****###
	#g = open(pkl_file, "rb")
	#tgfm_data_obj = pickle.load(g)
	#g.close()
	##***********## TEMP #####****###


	# Extract log-prior probabilities of each element being causal
	uniform_ln_prior_gene, uniform_ln_prior_variant = extract_uniform_log_prior_probabilities(tgfm_data_obj['genes'], tgfm_data_obj['variants'], study_names)
	#variant_gene_ln_prior_gene, variant_gene_ln_prior_variant = extract_variant_gene_log_prior_probabilities(tgfm_data_obj['genes'], tgfm_data_obj['variants'], study_names, tgfm_sldsc_results_dir)
	#sparse_variant_gene_tissue_ln_prior_gene, sparse_variant_gene_tissue_ln_prior_variant = extract_sparse_variant_gene_tissue_log_prior_probabilities(tgfm_data_obj['genes'], tgfm_data_obj['variants'], study_names, tgfm_sldsc_results_dir)

	# Add trait data
	tgfm_trait_data = {}
	tgfm_trait_data['gwas_beta'] = gwas_beta
	tgfm_trait_data['gwas_beta_se'] = gwas_beta_se
	tgfm_trait_data['gwas_sample_size'] = study_sample_sizes
	tgfm_trait_data['gwas_study_names'] = study_names
	tgfm_trait_data['uniform_ln_prior_gene'] = uniform_ln_prior_gene
	tgfm_trait_data['uniform_ln_prior_variant'] = uniform_ln_prior_variant
	#tgfm_trait_data['variant_gene_ln_prior_gene'] = variant_gene_ln_prior_gene
	#tgfm_trait_data['variant_gene_ln_prior_variant'] = variant_gene_ln_prior_variant
	#tgfm_trait_data['sparse_variant_gene_tissue_ln_prior_gene'] = sparse_variant_gene_tissue_ln_prior_gene
	#tgfm_trait_data['sparse_variant_gene_tissue_ln_prior_variant'] = sparse_variant_gene_tissue_ln_prior_variant

	# Save TGFM input data object
	trait_pkl_file = preprocessed_tgfm_data_dir + window_name + '_tgfm_ukbb_data_obj.pkl'
	g = open(trait_pkl_file, "wb")
	pickle.dump(tgfm_trait_data, g)
	g.close()

	t.write(window_name + '\t' + window_variant_in_sample_ld_file + '\t' + pkl_file + '\t' + trait_pkl_file + '\n')
	t.flush()


t.close()

