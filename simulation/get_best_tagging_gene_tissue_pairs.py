import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os
import pdb
from pandas_plink import read_plink1_bin
import pickle


def extract_gene_gene_ld(standardized_eqtl_effects, variant_ld):
	expression_covariance = np.dot(np.dot(standardized_eqtl_effects, variant_ld), np.transpose(standardized_eqtl_effects))
	np.fill_diagonal(expression_covariance, 1.0)
	dd = np.diag(1.0/np.sqrt(np.diag(expression_covariance)))
	ge_ld = np.dot(np.dot(dd, expression_covariance),dd)
	return ge_ld



def get_best_tagging_gt_pairs_for_missing_tissue0(tgfm_input_summary_file, best_tagging_gt_output_file):
	# Now loop through n_windows
	# In each window run TGFM independently
	# Loop through trait components
	tgfm_input_data = np.loadtxt(tgfm_input_summary_file,dtype=str,delimiter='\t')
	tgfm_input_data = tgfm_input_data[1:,:]

	# Open output file
	t = open(best_tagging_gt_output_file,'w')
	t.write('window\tgene_tissue_pair\tbest_tagging_gene_tissue_pair\tbest_abs_corr\tall_tagging_gene_tissue_pairs\tall_abs_corrs\n')

	# Get n_windows on this run
	n_windows = tgfm_input_data.shape[0]

	for window_iter in range(n_windows):
		data = tgfm_input_data[window_iter, :]

		##############################
		# Extract relevent fields
		###############################
		window_name = data[0]
		print(window_name)

		ld_file = data[1]
		tgfm_input_pkl = data[2]

		##############################
		# Load in Data
		###############################
		# Load in tgfm input data
		g = open(tgfm_input_pkl, "rb")
		tgfm_data = pickle.load(g)
		g.close()

		# Skip windows with no genes
		if len(tgfm_data['genes']) == 0:
			print('skipped because of no genes')
			continue

		# Load in Variant LD
		ld_mat = np.load(ld_file)
		# Add ld to tgfm_data obj
		tgfm_data['reference_ld'] = ld_mat


		# Get gene-gene LD
		gene_gene_ld = extract_gene_gene_ld(tgfm_data['gene_eqtl_pmces'], tgfm_data['reference_ld'])

		# Loop through genes
		for gt_index, gt_pair in enumerate(tgfm_data['genes']):
			# Extract tissue name from gene-tissue pair
			tissue_name = gt_pair.split('_')[1]

			# Ignore gene if not from tissue
			if tissue_name != 'tissue0':
				continue

			# Get gene name
			gene_name = gt_pair.split('_')[0]

			# Now get best tagging gene-tissue pair
			best_tagging_gt_pair = 'none'
			best_abs_corr = 0.0

			all_tagging_gt_pairs = []
			all_tagging_gt_pair_abs_corrs = []

			# Loop through genes again!
			# Get best tagging tissue
			for gt_prime_index, gt_prime_pair in enumerate(tgfm_data['genes']):
				# Get tissue name of gt_prime
				gt_prime_tissue_name = gt_prime_pair.split('_')[1]
				# Get gene name of gt_prime
				gt_prime_gene_name = gt_prime_pair.split('_')[0]

				# Ignore gene if from tissue0
				if gt_prime_tissue_name == 'tissue0':
					continue
				# Ignore gene if not from gene_name
				if gt_prime_gene_name != gene_name:
					continue

				# Get abs correlation
				abs_corr = np.abs(gene_gene_ld[gt_index, gt_prime_index])

				all_tagging_gt_pairs.append(gt_prime_pair)
				all_tagging_gt_pair_abs_corrs.append(abs_corr)

				if abs_corr > best_abs_corr:
					best_abs_corr = abs_corr
					best_tagging_gt_pair = gt_prime_pair
			if best_tagging_gt_pair == 'none':
				all_tagging_gt_pairs.append('none')
				all_tagging_gt_pair_abs_corrs.append('none')
			all_tagging_gt_pairs = np.asarray(all_tagging_gt_pairs)
			all_tagging_gt_pair_abs_corrs = np.asarray(all_tagging_gt_pair_abs_corrs)
			t.write(window_name + '\t' + gt_pair + '\t' + best_tagging_gt_pair + '\t' + str(best_abs_corr) + '\t' + ';'.join(all_tagging_gt_pairs) + '\t' + ';'.join(all_tagging_gt_pair_abs_corrs.astype(str)) + '\n')
	t.close()
	return

def mean_impute_and_standardize_genotype(G_obj_geno):
	G_obj_geno_stand = np.copy(G_obj_geno)
	ncol = G_obj_geno_stand.shape[1]
	#n_missing = []
	for col_iter in range(ncol):
		nan_indices = np.isnan(G_obj_geno[:,col_iter])
		non_nan_mean = np.mean(G_obj_geno[nan_indices==False, col_iter])
		G_obj_geno_stand[nan_indices, col_iter] = non_nan_mean
		G_obj_geno_stand[:, col_iter] = (G_obj_geno_stand[:, col_iter] - np.mean(G_obj_geno_stand[:, col_iter]))/np.std(G_obj_geno_stand[:, col_iter])
		#n_missing.append(np.sum(nan_indices))
	#n_missing = np.asarray(n_missing)

	#G_obj_geno_stand = (G_obj_geno_stand -np.mean(G_obj_geno_stand,axis=0))/np.std(G_obj_geno_stand,axis=0)

	return G_obj_geno_stand



def get_best_simulated_causal_effect_tagging_gt_pairs_for_missing_tissue0(sim_gene_expression_summary, output_file, gwas_plink_stem, sim_gene_trait_effect_size_file):
	valid_genes = {}
	f = open(sim_gene_trait_effect_size_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		gene_id = data[0]
		tissue0_effect_size = float(data[1])

		if tissue0_effect_size == 0:
			continue
		valid_genes[gene_id] = 1
	f.close()


	# Load in genotype object
	genotype_obj = read_plink1_bin(gwas_plink_stem + '.bed', gwas_plink_stem + '.bim', gwas_plink_stem + '.fam', verbose=False)

	# Open output file
	t = open(output_file,'w')
	t.write('gene_tissue_pair\tbest_tagging_gene_tissue_pair\tbest_abs_corr\tall_tagging_gene_tissue_pairs\tall_abs_corrs\n')


	# Loop through genes
	f = open(sim_gene_expression_summary)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent fields from line
		gene_id = data[0]
		if gene_id not in valid_genes:
			continue
		print(gene_id)
		causal_eqtl_effect_summary_file = data[3]
		eqtl_indices_file = data[5]
		total_n_genome_snps = int(data[6])

		# Load in causal eqtl effects
		causal_eqtl_effect_sizes = np.load(causal_eqtl_effect_summary_file)

		# Load in "global" indices of snps defining causal eqtl effect sizes
		eqtl_snp_indices_raw = np.load(eqtl_indices_file)
		eqtl_snp_indices = np.asarray([False]*total_n_genome_snps)
		eqtl_snp_indices[eqtl_snp_indices_raw] = True		

		# Total number of snps
		global_n_snps = len(eqtl_snp_indices)

		# Extract eqtl variant names
		eqtl_index_positions = np.arange(global_n_snps)[eqtl_snp_indices].astype(int).astype(str)
		eqtl_index_names = []
		for eqtl_index_position in eqtl_index_positions:
			eqtl_index_names.append('variant' + eqtl_index_position)
		eqtl_index_names = np.asarray(eqtl_index_names)

		# Extract genotype matrix for these snps
		eqtl_genotype = np.asarray(genotype_obj.sel(variant=eqtl_index_names))
		stand_eqtl_genotype = mean_impute_and_standardize_genotype(eqtl_genotype)
		
		# Get genetic ge
		genetically_predicted_gene_expression = np.dot(stand_eqtl_genotype, causal_eqtl_effect_sizes)

		corr_mat = np.corrcoef(np.transpose(genetically_predicted_gene_expression))

		corrz = corr_mat[0,:]

		if np.isnan(corrz[0]):
			continue

		gt_name = gene_id + '_' + 'tissue0'

		gt_prime_names = []
		gt_prime_abs_corrz = []

		for ii, corry in enumerate(corrz[1:]):
			if np.isnan(corry):
				continue
			tissue_name = 'tissue' + str(ii+1)
			gene_tissue_name = gene_id + '_' + tissue_name
			gt_prime_names.append(gene_tissue_name)
			gt_prime_abs_corrz.append(np.abs(corry))

		if len(gt_prime_names) == 0:
			t.write(gt_name + '\t' + 'none' + '\t' + '0.0' + '\t' + 'none\tnone\n')
		else:
			gt_prime_names = np.asarray(gt_prime_names)
			gt_prime_abs_corrz = np.asarray(gt_prime_abs_corrz)
			best_gt_prime_index = np.argmax(gt_prime_abs_corrz)
			t.write(gt_name + '\t' + gt_prime_names[best_gt_prime_index] + '\t' + str(gt_prime_abs_corrz[best_gt_prime_index]) + '\t' + ';'.join(gt_prime_names) + '\t' + ';'.join(gt_prime_abs_corrz.astype(str)) + '\n')
	f.close()
	t.close()
	return






#################
# Command line args
tgfm_input_summary_file = sys.argv[1]
tgfm_tissues = sys.argv[2]
best_tagging_gt_output_stem = sys.argv[3]
sim_gene_expression_summary = sys.argv[4]
processed_genotype_data_dir = sys.argv[5]
sim_gene_trait_effect_size_file = sys.argv[6]
eqtl_sample_size = sys.argv[7]



gwas_plink_stem = processed_genotype_data_dir + 'simulated_gwas_data_' + str(1)  # Genotype directory
eqtl_plink_stem = processed_genotype_data_dir + 'simulated_eqtl_' + str(eqtl_sample_size) + '_data_' + str(1)  # Genotype directory


if tgfm_tissues == 'no_t0':
	get_best_tagging_gt_pairs_for_missing_tissue0(tgfm_input_summary_file, best_tagging_gt_output_stem + '.txt')
	get_best_simulated_causal_effect_tagging_gt_pairs_for_missing_tissue0(sim_gene_expression_summary, best_tagging_gt_output_stem + '_from_sim_causal_effect.txt', gwas_plink_stem, sim_gene_trait_effect_size_file)
	get_best_simulated_causal_effect_tagging_gt_pairs_for_missing_tissue0(sim_gene_expression_summary, best_tagging_gt_output_stem + '_from_sim_causal_eqtl_effect.txt', eqtl_plink_stem, sim_gene_trait_effect_size_file)

print(best_tagging_gt_output_stem + '_from_sim_causal_eqtl_effect.txt')

