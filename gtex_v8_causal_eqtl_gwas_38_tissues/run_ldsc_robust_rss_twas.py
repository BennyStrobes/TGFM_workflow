import numpy as np 
import os
import sys
import pdb
import pickle

def extract_tissue_names(gtex_pseudotissue_file):
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

def extract_ld_annotations_for_this_region(ld, susie_mu, susie_mu_sd, susie_alpha, ordered_tissue_names, region_tissue_names, variant_indices):
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

def get_tissue_names_from_gene_tissue_names_arr(gene_tissue_names):
	tissue_names = []

	for ele in gene_tissue_names:
		tissue_names.append('_'.join(ele.split('_')[1:]))
	return np.asarray(tissue_names)

def extract_middle_independent_variant_indices(region_name, bim_mat, ld_mat, ld_sq_thresh):
	region_name_info = region_name.split(':')
	region_start = float(region_name_info[1])
	region_end = float(region_name_info[2])

	if region_end - region_start != 3000000.0:
		print('assumption eroror')
		pdb.set_trace()

	variant_positions = bim_mat[:,3].astype(float)

	middle_start = region_start + 1000000.0
	middle_end = region_start + 2000000.0

	middle_indices = np.where((variant_positions >= middle_start) & (variant_positions < middle_end))[0]

	valid_snps = np.asarray([True]*len(middle_indices))
	picked_snps = np.asarray([False]*len(middle_indices))
	middle_ld_sq = np.square(ld_mat)[middle_indices, :][:, middle_indices]



	if np.sum(valid_snps) == 0:
		keep_going = False
		return middle_indices
	else:
		keep_going = True

	while keep_going:
		snp_index = np.random.choice(np.arange(len(valid_snps))[valid_snps])
		picked_snps[snp_index] = True

		new_valid_snps = middle_ld_sq[snp_index,:] <= ld_sq_thresh
		valid_snps = valid_snps*new_valid_snps

		if np.sum(valid_snps) == 0:
			keep_going = False


	independent_middle_indices = middle_indices[picked_snps]

	return independent_middle_indices


def extract_middle_variant_indices(region_name, bim_mat):
	region_name_info = region_name.split(':')
	region_start = float(region_name_info[1])
	region_end = float(region_name_info[2])

	if region_end - region_start != 3000000.0:
		print('assumption eroror')
		pdb.set_trace()

	variant_positions = bim_mat[:,3].astype(float)

	middle_start = region_start + 1000000.0
	middle_end = region_start + 2000000.0

	middle_indices = np.where((variant_positions >= middle_start) & (variant_positions < middle_end))[0]

	if len(middle_indices) > 300:
		middle_indices = np.random.choice(middle_indices, 300, replace=False)
		if len(np.unique(middle_indices)) != 300:
			print('assumption erorroorro')
			pdb.set_trace()

	return middle_indices


def extract_tgfm_ldscore_annotations(ordered_tissue_names, chrom_component_data_file, output_file_anno, output_file_std_anno, output_file_chi_sq):
	used_regions = {}
	f = open(chrom_component_data_file)

	t = open(output_file_anno,'w')
	t.write('variant_name\tregression_weight\tld_score')
	for tissue_name in ordered_tissue_names:
		t.write('\t' + tissue_name + '_eqtl_ld_score')
	t.write('\n')
	t2 = open(output_file_chi_sq,'w')
	t2.write('variant_name\tchi_square_stat\tsample_size\n')

	t3 = open(output_file_std_anno,'w')
	t3.write('variant_name\tregression_weight\tld_score')
	for tissue_name in ordered_tissue_names:
		t3.write('\t' + tissue_name + '_eqtl_ld_score')
	t3.write('\n')


	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent fields from line corresponding to trait component
		component_name = data[1]
		num_genes = int(data[2])
		pickled_data_file = data[3]

		region_name = data[1].split('_')[-1]
		if region_name in used_regions:
			continue
		used_regions[region_name] = 1


		if pickled_data_file == 'NA':
			continue
		else:
			# Load in pickled data for this component
			f = open(pickled_data_file, "rb")
			twas_data = pickle.load(f)
			f.close()
			region_tissue_names = get_tissue_names_from_gene_tissue_names_arr(twas_data['genes'])

			#middle_variant_indices = extract_middle_variant_indices(region_name, twas_data['bim'])
			middle_variant_indices = extract_middle_independent_variant_indices(region_name, twas_data['bim'], twas_data['reference_ld'], .1)


			ld_scores, eqtl_ld_scores, standardized_eqtl_ld_scores, regression_weights = extract_ld_annotations_for_this_region(twas_data['reference_ld'], twas_data['susie_mu'], twas_data['susie_mu_sd'], twas_data['susie_alpha'], ordered_tissue_names, region_tissue_names, middle_variant_indices)

			chi_squared_stats = np.square(twas_data['gwas_beta']/twas_data['gwas_beta_se'])
			chi_squared_stats = chi_squared_stats[middle_variant_indices]


			for variant_index, variant_id in enumerate(twas_data['variants'][middle_variant_indices]):
				t.write(variant_id + '\t' + str(regression_weights[variant_index]) + '\t' + str(ld_scores[variant_index]) + '\t' + '\t'.join(eqtl_ld_scores[variant_index,:].astype(str)) + '\n')
				t2.write(variant_id + '\t' + str(chi_squared_stats[variant_index]) + '\t' + str(twas_data['gwas_sample_size']) + '\n')
				t3.write(variant_id + '\t' + str(regression_weights[variant_index]) + '\t' + str(ld_scores[variant_index]) + '\t' + '\t'.join(standardized_eqtl_ld_scores[variant_index,:].astype(str)) + '\n')

	f.close()
	t.close()
	t2.close()
	t3.close()


trait_name = sys.argv[1]
gtex_pseudotissue_file = sys.argv[2]
pseudotissue_gtex_rss_multivariate_twas_data_dir = sys.argv[3]
ukbb_genome_wide_susie_organized_results_dir = sys.argv[4]
pseudotissue_gtex_rss_multivariate_twas_dir = sys.argv[5]
gene_version = sys.argv[6]
chrom_num = int(sys.argv[7])


# Get gtex tissue names
ordered_tissue_names = extract_tissue_names(gtex_pseudotissue_file)


# Extract LD-score annotations on each chromosome seperately
chrom_component_data_file=pseudotissue_gtex_rss_multivariate_twas_data_dir + trait_name + '_' + gene_version + '_' + str(chrom_num) + '_component_rss_multivariate_twas_data_organized.txt'
output_file_anno = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + gene_version + '_' + str(chrom_num) + '_tgfm_ldsc_annotations.txt'
output_file_std_anno = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + gene_version + '_' + str(chrom_num) + '_tgfm_ldsc_standardized_annotations.txt'
output_file_chi_sq = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + gene_version + '_' + str(chrom_num) + '_tgfm_ldsc_chi_sq.txt'

extract_tgfm_ldscore_annotations(ordered_tissue_names, chrom_component_data_file, output_file_anno, output_file_std_anno, output_file_chi_sq)
