import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os
import pdb
import pickle



def extract_data_from_tgfm_twas_pickle_file(tgfm_tiss_prior_twas_pickle_file):
	f = open(tgfm_tiss_prior_twas_pickle_file, "rb")
	tgfm_twas_data = pickle.load(f)
	f.close()

	# Extract relevent fields from twas object
	gene_tissue_pairs = tgfm_twas_data.genes
	twas_fusion_nominal_z_scores = tgfm_twas_data.nominal_twas_z
	twas_tgfm_nominal_z_scores = tgfm_twas_data.nominal_twas_rss_z
	twas_tgfm_multivariate_z_scores = tgfm_twas_data.alpha_mu/np.sqrt(tgfm_twas_data.alpha_var)
	
	# Get gene names and tissue names
	gene_names = []
	tissue_names = []
	for gene_tissue_pair in gene_tissue_pairs:
		gene_name = gene_tissue_pair.split('_')[0]
		tissue_name = '_'.join(gene_tissue_pair.split('_')[1:])
		gene_names.append(gene_name)
		tissue_names.append(tissue_name)
	gene_names = np.asarray(gene_names)
	tissue_names = np.asarray(tissue_names)
	
	return gene_names, tissue_names, twas_fusion_nominal_z_scores, twas_tgfm_multivariate_z_scores, tgfm_twas_data.alpha_mu, tgfm_twas_data.alpha_var

def extract_data_from_tgfm_fine_mapping_table(tgfm_tiss_prior_fine_mapping_file):
	data = np.loadtxt(tgfm_tiss_prior_fine_mapping_file, dtype=str, delimiter='\t')[1:,:]

	# Create mapping from tissue-gene to sum_posterior prob
	tg_to_sum_posterior_prob = {}
	tg_to_max_posterior_prob = {}
	nrows = data.shape[0]

	for row_num in range(nrows):
		tissue_gene_name = data[row_num, 0] + '_' + data[row_num, 1]
		posterior_prob = float(data[row_num, 3])
		if tissue_gene_name not in tg_to_sum_posterior_prob:
			tg_to_sum_posterior_prob[tissue_gene_name] = posterior_prob
		else:
			tg_to_sum_posterior_prob[tissue_gene_name] = tg_to_sum_posterior_prob[tissue_gene_name] + posterior_prob
		if tissue_gene_name not in tg_to_max_posterior_prob:
			tg_to_max_posterior_prob[tissue_gene_name] = posterior_prob
		else:
			tg_to_max_posterior_prob[tissue_gene_name] = np.max([posterior_prob, tg_to_max_posterior_prob[tissue_gene_name]])
	return tg_to_sum_posterior_prob

def get_var_dicti(dir_name, trait_name):
	file_name = dir_name + trait_name + '_component_organized_multivariate_twas_overlaps.txt'
	dicti = {}
	if os.path.exists(file_name) == False:
		return dicti
	f = open(file_name)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissue_gene_name = data[1] + '_' + data[2]
		sd = float(data[5])
		dicti[tissue_gene_name] = sd
	f.close()
	return dicti

def extract_gene_variances(twas_data_file):
	f = open(twas_data_file, "rb")
	twas_data = pickle.load(f)
	f.close()

	n_genes = len(twas_data['susie_mu'])

	varz = []
	for g_index in range(n_genes):
		var = np.sum((np.square(twas_data['susie_mu_sd'][g_index]) + np.square(twas_data['susie_mu'][g_index]))*twas_data['susie_alpha'][g_index])
		varz.append(var)
	return np.asarray(varz)

def get_coloc_prob_dicti_for_set_of_genes(coloc_results_dir, unique_genes, trait_name, adaptive=False):
	dicti = {}
	for gene_name in unique_genes:
		if adaptive:
			file_name = coloc_results_dir + trait_name + '_' + gene_name + '_adaptive_coloc_posterior_probabilities.txt'
		else:
			file_name = coloc_results_dir + trait_name + '_' + gene_name + '_coloc_posterior_probabilities.txt'
		f = open(file_name)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			tissue_name = data[0]
			posterior_prob = (data[-1])
			tissue_gene_name = tissue_name + '_' + gene_name
			dicti[tissue_gene_name] = posterior_prob
		f.close()
	return dicti

####################
# Command line args
####################

trait_name = sys.argv[1]
gtex_pseudotissue_file = sys.argv[2]
pseudotissue_gtex_rss_multivariate_twas_dir = sys.argv[3]
gene_version = sys.argv[4]
coloc_results_dir = sys.argv[5]


# Version with tissue specific prior (NEED TO FIX ISSUE STILLLLLLL)
output_file = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + gene_version + '_component_organized_tgfm_results.txt'
t = open(output_file,'w')
t.write('trait_component\ttissue\tgene\tnominal_fusion_twas_z_score\ttgfm_const_prior_twas_alpha\ttgfm_tissue_prior_twas_alpha\trobust_tgfm_const_prior_twas_alpha\trobust_tgfm_tissue_prior_twas_alpha\ttgfm_const_prior_twas_alpha_var\ttgfm_tissue_prior_twas_alpha_var\trobust_tgfm_const_prior_twas_alpha_var\trobust_tgfm_tissue_prior_twas_alpha_var\ttgfm_const_prior_twas_z_score\ttgfm_tissue_prior_twas_z_score\trobust_tgfm_const_prior_twas_z_score\trobust_tgfm_tissue_prior_twas_z_score\ttgfm_rss_regression_const_prior_posterior_prob\ttgfm_rss_regression_tissue_prior_posterior_prob\trobust_tgfm_rss_regression_const_prior_posterior_prob\trobust_tgfm_rss_regression_tissue_prior_posterior_prob\tcoloc_posterior_prob\tadaptive_prior_coloc_posterior_prob\n')

# Loop through chromosomes
for chrom_num in range(1,23):
	print(chrom_num)
	chromosome_summary_file = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + str(chrom_num) + '_summary_tgfm_results_' + gene_version + '_const_1e-5_prior.txt'
	# Stream chromosome summary file
	f = open(chromosome_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent data
		component_name = data[0]
		num_genes = data[1]

		# Several different versions of analysis
		# Version 1: TGFM with constant prior
		tgfm_const_prior_twas_pickle_file = data[2]
		#tgfm_const_prior_fine_mapping_file = data[3]
		tgfm_const_prior_rss_regression_fine_mapping_file = data[5]
		#tgfm_const_prior_rss_prediction_fine_mapping_file = data[7]
		#tgfm_const_prior_rss_regression_nominal_fine_mapping_file = data[9]
		if tgfm_const_prior_twas_pickle_file == 'NA':
			continue


		# Version 2: TGFM with tissue prior
		tgfm_tiss_prior_twas_pickle_file = tgfm_const_prior_twas_pickle_file.split('tgfm_twas_results.pkl')[0] + 'tgfm_twas_tissue_specific_prior_results.pkl'
		#tgfm_tiss_prior_fine_mapping_file = tgfm_const_prior_fine_mapping_file.split('table.txt')[0] + 'tissue_specific_prior_table.txt'
		tgfm_tiss_prior_rss_regression_fine_mapping_file = tgfm_const_prior_rss_regression_fine_mapping_file.split('table.txt')[0] + 'tissue_specific_prior_table.txt'
		#tgfm_tiss_prior_rss_prediction_fine_mapping_file = tgfm_const_prior_rss_prediction_fine_mapping_file.split('table.txt')[0] + 'tissue_specific_prior_table.txt'

		# Version 3: Robust TGFM with constant prior
		robust_tgfm_const_prior_twas_pickle_file = tgfm_const_prior_twas_pickle_file.split('tgfm_twas_results.pkl')[0] + 'tgfm_robust_twas_pleiotropic_only_prior_results.pkl'
		robust_tgfm_const_prior_rss_regression_fine_mapping_file = tgfm_const_prior_rss_regression_fine_mapping_file.split('tgfm_rss_regression_fine_mapping_table.txt')[0] + 'tgfm_robust_rss_regression_fine_mapping_pleiotropic_only_prior_table.txt'

		# Version 4: Robust TGFM with tissue prior
		robust_tgfm_tissue_prior_twas_pickle_file = tgfm_const_prior_twas_pickle_file.split('tgfm_twas_results.pkl')[0] + 'tgfm_robust_twas_tissue_specific_prior_results.pkl'
		robust_tgfm_tissue_prior_rss_regression_fine_mapping_file = tgfm_const_prior_rss_regression_fine_mapping_file.split('tgfm_rss_regression_fine_mapping_table.txt')[0] + 'tgfm_robust_rss_regression_fine_mapping_tissue_specific_prior_table.txt'


		# Extract TGFM TWAS data for this disease component
		tp_gene_names, tp_tissue_names, tp_twas_fusion_nominal_z_scores, tp_twas_tgfm_multivariate_z_scores, tp_twas_tgfm_multivariate_alpha_mu, tp_twas_tgfm_multivariate_alpha_var  = extract_data_from_tgfm_twas_pickle_file(tgfm_tiss_prior_twas_pickle_file)
		cp_gene_names, cp_tissue_names, cp_twas_fusion_nominal_z_scores, cp_twas_tgfm_multivariate_z_scores, cp_twas_tgfm_multivariate_alpha_mu, cp_twas_tgfm_multivariate_alpha_var = extract_data_from_tgfm_twas_pickle_file(tgfm_const_prior_twas_pickle_file)
		r_cp_gene_names, r_cp_tissue_names, r_cp_twas_fusion_nominal_z_scores, r_cp_twas_tgfm_multivariate_z_scores, r_cp_twas_tgfm_multivariate_alpha_mu, r_cp_twas_tgfm_multivariate_alpha_var  = extract_data_from_tgfm_twas_pickle_file(robust_tgfm_const_prior_twas_pickle_file)
		r_tp_gene_names, r_tp_tissue_names, r_tp_twas_fusion_nominal_z_scores, r_tp_twas_tgfm_multivariate_z_scores, r_tp_twas_tgfm_multivariate_alpha_mu, r_tp_twas_tgfm_multivariate_alpha_var  = extract_data_from_tgfm_twas_pickle_file(robust_tgfm_tissue_prior_twas_pickle_file)



		# Quick error checking
		if np.array_equal(tp_gene_names, cp_gene_names) == False or np.array_equal(tp_tissue_names, cp_tissue_names) == False:
			print('assumption error')
			pdb.set_trace()

		# Extract TGFM fine mapping data for this disease component
		#tp_sum_post_prob_dicti = extract_data_from_tgfm_fine_mapping_table(tgfm_tiss_prior_fine_mapping_file)
		tp_sum_post_prob_rss_regression_dicti = extract_data_from_tgfm_fine_mapping_table(tgfm_tiss_prior_rss_regression_fine_mapping_file)
		#tp_sum_post_prob_rss_prediction_dicti = extract_data_from_tgfm_fine_mapping_table(tgfm_tiss_prior_rss_prediction_fine_mapping_file)


		#cp_sum_post_prob_dicti = extract_data_from_tgfm_fine_mapping_table(tgfm_const_prior_fine_mapping_file)
		cp_sum_post_prob_rss_regression_dicti = extract_data_from_tgfm_fine_mapping_table(tgfm_const_prior_rss_regression_fine_mapping_file)
		#cp_sum_post_prob_rss_prediction_dicti = extract_data_from_tgfm_fine_mapping_table(tgfm_const_prior_rss_prediction_fine_mapping_file)
		#cp_sum_post_prob_rss_regression_nom_dicti = extract_data_from_tgfm_fine_mapping_table(tgfm_const_prior_rss_regression_nominal_fine_mapping_file)

		r_tp_sum_post_prob_rss_regression_dicti = extract_data_from_tgfm_fine_mapping_table(robust_tgfm_tissue_prior_rss_regression_fine_mapping_file)

		r_cp_sum_post_prob_rss_regression_dicti = extract_data_from_tgfm_fine_mapping_table(robust_tgfm_const_prior_rss_regression_fine_mapping_file)


		# Extract Unique genes
		unique_genes = np.unique(cp_gene_names)
		adapive_coloc_prob_dicti = get_coloc_prob_dicti_for_set_of_genes(coloc_results_dir, unique_genes, trait_name, adaptive=True)
		coloc_prob_dicti = get_coloc_prob_dicti_for_set_of_genes(coloc_results_dir, unique_genes, trait_name, adaptive=False)


		for g_index in range(len(tp_gene_names)):
			tissue_gene_name = cp_tissue_names[g_index] + '_' + cp_gene_names[g_index]
			t.write(component_name + '\t' + cp_tissue_names[g_index] + '\t' + cp_gene_names[g_index] + '\t')
			t.write(str(cp_twas_fusion_nominal_z_scores[g_index]) + '\t')

			t.write(str(cp_twas_tgfm_multivariate_alpha_mu[g_index]) + '\t')
			t.write(str(tp_twas_tgfm_multivariate_alpha_mu[g_index]) + '\t')
			t.write(str(r_cp_twas_tgfm_multivariate_alpha_mu[g_index]) + '\t')
			t.write(str(r_tp_twas_tgfm_multivariate_alpha_mu[g_index]) + '\t')

			t.write(str(cp_twas_tgfm_multivariate_alpha_var[g_index]) + '\t')
			t.write(str(tp_twas_tgfm_multivariate_alpha_var[g_index]) + '\t')
			t.write(str(r_cp_twas_tgfm_multivariate_alpha_var[g_index]) + '\t')
			t.write(str(r_tp_twas_tgfm_multivariate_alpha_var[g_index]) + '\t')

			t.write(str(cp_twas_tgfm_multivariate_z_scores[g_index]) + '\t')
			t.write(str(tp_twas_tgfm_multivariate_z_scores[g_index]) + '\t')
			t.write(str(r_cp_twas_tgfm_multivariate_z_scores[g_index]) + '\t')
			t.write(str(r_tp_twas_tgfm_multivariate_z_scores[g_index]) + '\t')


			t.write(str(cp_sum_post_prob_rss_regression_dicti[tissue_gene_name]) + '\t')
			t.write(str(tp_sum_post_prob_rss_regression_dicti[tissue_gene_name]) + '\t')
			t.write(str(r_cp_sum_post_prob_rss_regression_dicti[tissue_gene_name]) + '\t')
			t.write(str(r_tp_sum_post_prob_rss_regression_dicti[tissue_gene_name]) + '\t')

			t.write(coloc_prob_dicti[tissue_gene_name] + '\t')
			t.write(adapive_coloc_prob_dicti[tissue_gene_name] + '\n')

	f.close()
t.close()
