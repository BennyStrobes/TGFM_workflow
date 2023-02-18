import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os
import pdb
from pandas_plink import read_plink1_bin




def extract_ordered_regression_snps_names(simulation_gwas_file):
	arr = []
	f = open(simulation_gwas_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0])
	f.close()
	arr = np.asarray(arr)
	# Quick error check
	if len(arr) != len(np.unique(arr)):
		print('assumptino eroror')
		pdb.set_trace()
	return arr

def load_in_genotype_data(genotype_stem):
	# Load in Reference Genotype data
	G_obj = read_plink1_bin(genotype_stem + '.bed', genotype_stem + '.bim', genotype_stem + '.fam', verbose=False)

	G_obj_geno = G_obj.values # Numpy 2d array of dimension num samples X num snps
	G_obj_chrom = np.asarray(G_obj.chrom)
	G_obj_pos = np.asarray(G_obj.pos)
	# For our purposes, a0 is the effect allele
	# For case of plink package, a0 is the first column in the plink bim file
	G_obj_a0 = np.asarray(G_obj.a0)
	G_obj_a1 = np.asarray(G_obj.a1)
	# RSids
	G_obj_rsids = np.asarray(G_obj.snp)
	# Centimorgan distances
	G_obj_cm = np.asarray(G_obj.cm)

	# Put geno into organized dictionary
	genotype_obj = {'G': G_obj_geno, 'rsid': G_obj_rsids, 'position': G_obj_pos, 'cm': G_obj_cm}

	return genotype_obj

def create_regression_snp_names_dictionary(regression_snp_names):
	dicti = {}
	for snp_name in regression_snp_names:
		if snp_name in dicti:
			print('assumption eroror')
		dicti[snp_name] = 1
	return dicti

def get_ordered_list_of_gene_names_and_gene_indices_files_from_gene_summary_file(gene_summary_file):
	f = open(gene_summary_file)
	gene_names = []
	gene_tss = []
	gene_indices_files = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[0]
		tss = int(data[2])
		snp_index_file = data[5]
		gene_names.append(ensamble_id)
		gene_tss.append(tss)
		gene_indices_files.append(snp_index_file)
	f.close()
	return np.asarray(gene_names), np.asarray(gene_tss), np.asarray(gene_indices_files)


def extract_regression_snp_indices(global_rsids, regression_snp_names_set):
	boolean_vector = []
	for global_rsid in global_rsids:
		if global_rsid in regression_snp_names_set:
			boolean_vector.append(True)
		else:
			boolean_vector.append(False)
	return np.asarray(boolean_vector)

def extract_gene_window_snp_indices(global_snp_positions, global_snp_cms, eqtl_snp_indices):
	# Eqtl min cm position
	eqtl_min_cm = np.min(global_snp_cms[eqtl_snp_indices])
	eqtl_max_cm = np.max(global_snp_cms[eqtl_snp_indices])

	# Gene window cm lb and ub
	gene_window_cm_lb = eqtl_min_cm - 1.0
	gene_window_cm_ub = eqtl_max_cm + 1.0

	gene_window_snp_indices = (global_snp_cms >= gene_window_cm_lb) & (global_snp_cms < gene_window_cm_ub)

	return gene_window_snp_indices


def mean_impute_genotype(geno_mat):
	geno_mat_stand = np.copy(geno_mat)
	ncol = geno_mat_stand.shape[1]
	#n_missing = []
	for col_iter in range(ncol):
		nan_indices = np.isnan(geno_mat[:,col_iter])
		non_nan_mean = np.mean(geno_mat[nan_indices==False, col_iter])
		geno_mat_stand[nan_indices, col_iter] = non_nan_mean

	return geno_mat_stand



def extract_gene_window_ld(global_genotype_matrix, gene_window_snp_indices):
	# Subset global genotype matrix to just variants in gene window
	gene_window_genotype_matrix = global_genotype_matrix[:, gene_window_snp_indices]

	# Mean imput genotype matrix
	mean_imputed_window_genotype_matrix = mean_impute_genotype(gene_window_genotype_matrix)

	# Compute ld across snps in window
	gene_window_ld = np.corrcoef(np.transpose(mean_imputed_window_genotype_matrix))

	return gene_window_ld


def extract_ld_annotations_for_this_gene_region_with_eqtl_susie_distribution(ld, gene_eqtl_effect_sizes, gene_variance, variant_indices, n_ref_panel_samples):
	# Get standardized eqtl effect sizes (A little hacky because some tissues are all 0)
	standardized_gene_eqtl_effect_sizes = np.copy(gene_eqtl_effect_sizes)
	n_tiss = standardized_gene_eqtl_effect_sizes.shape[1]
	for tiss_iter in range(n_tiss):
		if gene_variance[tiss_iter] > 0.0:
			standardized_gene_eqtl_effect_sizes[:,tiss_iter] = standardized_gene_eqtl_effect_sizes[:,tiss_iter]/np.sqrt(gene_variance[tiss_iter])
	
	ld_scores = np.square(np.dot(ld[variant_indices,:], standardized_gene_eqtl_effect_sizes))

	# Get adjusted ld scores
	adj_ld_scores = np.copy(ld_scores)
	for tiss_iter in range(n_tiss):
		if gene_variance[tiss_iter] > 0.0:
			adj_ld_scores[:,tiss_iter] = ld_scores[:,tiss_iter] - ((1.0-ld_scores[:,tiss_iter])/(n_ref_panel_samples-2.0))
	return adj_ld_scores

def extract_ld_annotations_for_this_gene_region_with_eqtl_point_estimate(ld, gene_eqtl_effect_sizes, variant_indices, n_ref_panel_samples):
	gene_variance = np.diag(np.dot(np.dot(np.transpose(gene_eqtl_effect_sizes), ld), gene_eqtl_effect_sizes))

	# Get standardized eqtl effect sizes (A little hacky because some tissues are all 0)
	standardized_gene_eqtl_effect_sizes = np.copy(gene_eqtl_effect_sizes)
	n_tiss = standardized_gene_eqtl_effect_sizes.shape[1]
	for tiss_iter in range(n_tiss):
		if gene_variance[tiss_iter] > 0.0:
			standardized_gene_eqtl_effect_sizes[:,tiss_iter] = standardized_gene_eqtl_effect_sizes[:,tiss_iter]/np.sqrt(gene_variance[tiss_iter])
	
	ld_scores = np.square(np.dot(ld[variant_indices,:], standardized_gene_eqtl_effect_sizes))

	# Get adjusted ld scores
	adj_ld_scores = np.copy(ld_scores)
	for tiss_iter in range(n_tiss):
		if gene_variance[tiss_iter] > 0.0:
			adj_ld_scores[:,tiss_iter] = ld_scores[:,tiss_iter] - ((1.0-ld_scores[:,tiss_iter])/(n_ref_panel_samples-2.0))

	return adj_ld_scores

def extract_boolean_on_whether_gene_model_exists(gene_eqtl_pmces):
	n_tiss = gene_eqtl_pmces.shape[0]
	n_snps = gene_eqtl_pmces.shape[1]
	boolean_arr = []
	for tiss_iter in range(n_tiss):
		if np.array_equal(gene_eqtl_pmces[tiss_iter,:], np.zeros(n_snps)):
			boolean_arr.append(0.0)
		else:
			boolean_arr.append(1.0)
	return np.asarray(boolean_arr)

def print_gene_weighted_ld_scores_to_output(gene_weighted_ld_scores, regression_snp_names, gene_weighted_ld_score_output_file):
	t = open(gene_weighted_ld_score_output_file,'w')

	# Quick error check
	if gene_weighted_ld_scores.shape[0] != len(regression_snp_names):
		print('assumption error')

	for snp_iter, snp_name in enumerate(regression_snp_names):
		t.write(snp_name + '\t' + '\t'.join(gene_weighted_ld_scores[snp_iter, :].astype(str)) + '\n')
	t.close()

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




def generate_susie_distr_gene_weighted_ld_scores(regression_snp_names, genotype_obj, eqtl_sample_sizes, gene_summary_file, simulated_learned_gene_models_dir, simulation_name_string, simulated_gene_expression_dir, gene_weighted_ld_score_output_root):
	# Create dictionary of regression_snp_names
	regression_snp_names_set = create_regression_snp_names_dictionary(regression_snp_names)

	# Extract boolean vector of regression snps
	regression_snp_indices = extract_regression_snp_indices(genotype_obj['rsid'], regression_snp_names_set)

	# Initialize gene_weighted_ld_scores matrix and num_genes_counters for various eqtl data sets
	n_regression_snps = len(regression_snp_names)
	gene_weighted_ld_scores_across_eqtl_data_sets = []
	num_genes_across_eqtl_data_sets = []
	# Loop through eqtl data sets to initialize matrix for each data set
	for eqtl_sample_size in eqtl_sample_sizes:
		gene_weighted_ld_scores_across_eqtl_data_sets.append(np.zeros((n_regression_snps, 10)))
		num_genes_across_eqtl_data_sets.append(np.zeros(10))

	# Get ordered list of gene names and gene-index files
	gene_names,gene_tss_positions,gene_indices_files = get_ordered_list_of_gene_names_and_gene_indices_files_from_gene_summary_file(gene_summary_file)

	# Number of reference panel samples
	n_ref_panel_samples = genotype_obj['G'].shape[0]

	# Loop through genes
	for gene_iter, gene_name in enumerate(gene_names):
		# Extract other relevent information for this gene
		gene_tss = gene_tss_positions[gene_iter]
		gene_snp_index_file = gene_indices_files[gene_iter]

		# Extract indices of snps nearby gene used for eQTL calling
		eqtl_snp_indices = np.load(gene_snp_index_file)

		# Extract indices of snps within 1CM of gene start and end
		gene_window_snp_indices = extract_gene_window_snp_indices(genotype_obj['position'], genotype_obj['cm'], eqtl_snp_indices)
		n_gene_window_snps = np.sum(gene_window_snp_indices)

		# Extract gene window LD
		gene_window_ld = extract_gene_window_ld(genotype_obj['G'], gene_window_snp_indices)

		# Get regression snp indices corresponding to this gene window
		gene_window_regression_snp_indices = regression_snp_indices[gene_window_snp_indices]
		# Skip windows with no regression snps 
		if np.sum(gene_window_regression_snp_indices) == 0:
			continue

		# Subloop through eqtl sample sizes
		for eqtl_sample_size_iter, eqtl_sample_size in enumerate(eqtl_sample_sizes):
			# Extract estimated eQTL PMCES for this gene
			# Use learned eqtl effect sizes
			gene_eqtl_pmces_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_sample_size) + '_gene_model_pmces.npy'
			gene_eqtl_pmces = np.load(gene_eqtl_pmces_file)

			# Extract boolean whether gene model exists for each tissue
			gene_model_boolean_cross_tissues = extract_boolean_on_whether_gene_model_exists(gene_eqtl_pmces)
			num_genes_across_eqtl_data_sets[eqtl_sample_size_iter] = num_genes_across_eqtl_data_sets[eqtl_sample_size_iter] + gene_model_boolean_cross_tissues

			# Get eQTL PMCES FOR GENE WINDOW (now just gene window space, global space and regression snp space)
			gene_window_eqtl_pmces = np.zeros((n_gene_window_snps, 10))
			gene_window_eqtl_pmces[eqtl_snp_indices[gene_window_snp_indices], :] = np.transpose(gene_eqtl_pmces)

			# Extract variance for each gene
			gene_variances = []
			for tiss_iter, gene_model_boolean in enumerate(gene_model_boolean_cross_tissues):
				# Gene-tissue pairs with no gene model get no variance
				if gene_model_boolean == 0.0:
					gene_variances.append(0.0)
				else:
					# Load in susie files for this gene
					gene_susie_mu_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_sample_size) + '_tissue_' + str(tiss_iter) + '_gene_model_susie_mu.npy'
					gene_susie_alpha_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_sample_size) + '_tissue_' + str(tiss_iter) + '_gene_model_susie_alpha.npy'
					gene_susie_mu_var_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_sample_size) + '_tissue_' + str(tiss_iter) + '_gene_model_susie_mu_var.npy'
					gene_susie_mu = np.load(gene_susie_mu_file)
					gene_susie_alpha = np.load(gene_susie_alpha_file)
					gene_susie_mu_var = np.load(gene_susie_mu_var_file)

					# quick error checking
					if np.array_equal(np.sum(gene_susie_mu*gene_susie_alpha,axis=0), gene_eqtl_pmces[tiss_iter,:]) == False:
						print('PMCES mismatch assumption error')
						pdb.set_trace()

					# Load in susie data for gene window
					window_susie_mu = np.zeros((gene_susie_mu.shape[0], n_gene_window_snps))
					window_susie_mu[:, eqtl_snp_indices[gene_window_snp_indices]] = gene_susie_mu
					window_susie_alpha = np.zeros((gene_susie_alpha.shape[0], n_gene_window_snps))
					window_susie_alpha[:, eqtl_snp_indices[gene_window_snp_indices]] = gene_susie_alpha
					window_susie_mu_var = np.zeros((gene_susie_mu_var.shape[0], n_gene_window_snps))
					window_susie_mu_var[:, eqtl_snp_indices[gene_window_snp_indices]] = gene_susie_mu_var

					# Compute gene variance
					gene_variance = calculate_gene_variance_according_to_susie_distribution(window_susie_mu, window_susie_alpha, np.sqrt(window_susie_mu_var), gene_window_ld)
					gene_variances.append(gene_variance)
			gene_variances = np.asarray(gene_variances)

			# Extract gene window gene-ld scores
			gene_window_ld_scores_susie_distr = extract_ld_annotations_for_this_gene_region_with_eqtl_susie_distribution(gene_window_ld, gene_window_eqtl_pmces, gene_variances, gene_window_regression_snp_indices, n_ref_panel_samples)
			
			# Update global ld scores
			gene_weighted_ld_scores_across_eqtl_data_sets[eqtl_sample_size_iter][gene_window_snp_indices[regression_snp_indices], :] = gene_weighted_ld_scores_across_eqtl_data_sets[eqtl_sample_size_iter][gene_window_snp_indices[regression_snp_indices], :] + gene_window_ld_scores_susie_distr


	# Save everything to output
	for eqtl_sample_size_iter, eqtl_sample_size in enumerate(eqtl_sample_sizes):
		# Extract relevent info for this data set
		gene_weighted_ld_scores = gene_weighted_ld_scores_across_eqtl_data_sets[eqtl_sample_size_iter]
		num_genes = num_genes_across_eqtl_data_sets[eqtl_sample_size_iter]

		# Print
		gene_weighted_ld_score_output_file = gene_weighted_ld_score_output_root + '_eqtlss_' + str(eqtl_sample_size) + '_ld_scores.txt'
		print_gene_weighted_ld_scores_to_output(gene_weighted_ld_scores, regression_snp_names, gene_weighted_ld_score_output_file)
		num_genes_output_file = gene_weighted_ld_score_output_root + '_eqtlss_' + str(eqtl_sample_size) + '_M.txt'
		np.savetxt(num_genes_output_file, num_genes, fmt="%s", delimiter='\t')

	return

def generate_gene_weighted_ld_scores(regression_snp_names, genotype_obj, eqtl_sample_sizes, gene_summary_file, simulated_learned_gene_models_dir, simulation_name_string, simulated_gene_expression_dir, gene_weighted_ld_score_output_root):
	# Create dictionary of regression_snp_names
	regression_snp_names_set = create_regression_snp_names_dictionary(regression_snp_names)

	# Extract boolean vector of regression snps
	regression_snp_indices = extract_regression_snp_indices(genotype_obj['rsid'], regression_snp_names_set)

	# Initialize gene_weighted_ld_scores matrix and num_genes_counters for various eqtl data sets
	n_regression_snps = len(regression_snp_names)
	gene_weighted_ld_scores_across_eqtl_data_sets = []
	num_genes_across_eqtl_data_sets = []
	# Loop through eqtl data sets to initialize matrix for each data set
	for eqtl_sample_size in eqtl_sample_sizes:
		gene_weighted_ld_scores_across_eqtl_data_sets.append(np.zeros((n_regression_snps, 10)))
		num_genes_across_eqtl_data_sets.append(np.zeros(10))

	# Get ordered list of gene names and gene-index files
	gene_names,gene_tss_positions,gene_indices_files = get_ordered_list_of_gene_names_and_gene_indices_files_from_gene_summary_file(gene_summary_file)

	# Number of reference panel samples
	n_ref_panel_samples = genotype_obj['G'].shape[0]

	# Loop through genes
	for gene_iter, gene_name in enumerate(gene_names):
		# Extract other relevent information for this gene
		gene_tss = gene_tss_positions[gene_iter]
		gene_snp_index_file = gene_indices_files[gene_iter]

		# Extract indices of snps nearby gene used for eQTL calling
		eqtl_snp_indices = np.load(gene_snp_index_file)

		# Extract indices of snps within 1CM of gene start and end
		gene_window_snp_indices = extract_gene_window_snp_indices(genotype_obj['position'], genotype_obj['cm'], eqtl_snp_indices)
		n_gene_window_snps = np.sum(gene_window_snp_indices)

		# Extract gene window LD
		gene_window_ld = extract_gene_window_ld(genotype_obj['G'], gene_window_snp_indices)

		# Get regression snp indices corresponding to this gene window
		gene_window_regression_snp_indices = regression_snp_indices[gene_window_snp_indices]
		# Skip windows with no regression snps 
		if np.sum(gene_window_regression_snp_indices) == 0:
			continue

		# Subloop through eqtl sample sizes
		for eqtl_sample_size_iter, eqtl_sample_size in enumerate(eqtl_sample_sizes):
			# Extract estimated eQTL PMCES for this gene
			if eqtl_sample_size == 'inf':
				# Use known (simulated eqtl effect sizes)
				gene_eqtl_ce_file = simulated_gene_expression_dir + simulation_name_string + '_' + gene_name + '_causal_eqtl_effects.npy'
				gene_eqtl_pmces = np.transpose(np.load(gene_eqtl_ce_file))
			else:
				# Use learned eqtl effect sizes
				gene_eqtl_pmces_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_sample_size) + '_gene_model_pmces.npy'
				gene_eqtl_pmces = np.load(gene_eqtl_pmces_file)

			# Extract boolean whether gene model exists for each tissue
			gene_model_boolean_cross_tissues = extract_boolean_on_whether_gene_model_exists(gene_eqtl_pmces)
			num_genes_across_eqtl_data_sets[eqtl_sample_size_iter] = num_genes_across_eqtl_data_sets[eqtl_sample_size_iter] + gene_model_boolean_cross_tissues


			# Get eQTL PMCES FOR GENE WINDOW (now just gene window space, global space and regression snp space)
			gene_window_eqtl_pmces = np.zeros((n_gene_window_snps, 10))
			gene_window_eqtl_pmces[eqtl_snp_indices[gene_window_snp_indices], :] = np.transpose(gene_eqtl_pmces)

			# Extract gene window gene-ld scores
			gene_window_adj_ld_scores_pe = extract_ld_annotations_for_this_gene_region_with_eqtl_point_estimate(gene_window_ld, gene_window_eqtl_pmces, gene_window_regression_snp_indices, n_ref_panel_samples)
			
			# Update global ld scores
			gene_weighted_ld_scores_across_eqtl_data_sets[eqtl_sample_size_iter][gene_window_snp_indices[regression_snp_indices], :] = gene_weighted_ld_scores_across_eqtl_data_sets[eqtl_sample_size_iter][gene_window_snp_indices[regression_snp_indices], :] + gene_window_adj_ld_scores_pe



	# Save everything to output
	for eqtl_sample_size_iter, eqtl_sample_size in enumerate(eqtl_sample_sizes):
		# Extract relevent info for this data set
		gene_weighted_ld_scores = gene_weighted_ld_scores_across_eqtl_data_sets[eqtl_sample_size_iter]
		num_genes = num_genes_across_eqtl_data_sets[eqtl_sample_size_iter]

		# Print
		gene_weighted_ld_score_output_file = gene_weighted_ld_score_output_root + '_eqtlss_' + str(eqtl_sample_size) + '_ld_scores.txt'
		print_gene_weighted_ld_scores_to_output(gene_weighted_ld_scores, regression_snp_names, gene_weighted_ld_score_output_file)
		num_genes_output_file = gene_weighted_ld_score_output_root + '_eqtlss_' + str(eqtl_sample_size) + '_M.txt'
		np.savetxt(num_genes_output_file, num_genes, fmt="%s", delimiter='\t')

	return


##############################
# Command line argumemnts
##############################
simulation_number = int(sys.argv[1])
chrom_num = sys.argv[2]
simulation_name_string = sys.argv[3]
processed_genotype_data_dir = sys.argv[4]
simulated_gwas_dir = sys.argv[5]
simulated_gene_expression_dir = sys.argv[6]
simulated_learned_gene_models_dir = sys.argv[7]
simulated_ld_scores_dir = sys.argv[8]


####################################################
# Extract names and number of regression snps
####################################################
# Need to know number of regression snps
simulation_gwas_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results_hm3_noHMC_snps_only.txt'
regression_snp_names = extract_ordered_regression_snps_names(simulation_gwas_file)
n_regression_snps = len(regression_snp_names)

# Save regression snps to output
np.savetxt(simulated_ld_scores_dir + simulation_name_string + '_regression_snp_ids.txt', regression_snp_names, fmt="%s", delimiter="\n")


####################################################
# Load in genotype data for whole chromosome in EUR ancestry 1KG individuals
####################################################
genotype_stem = processed_genotype_data_dir + '100G.EUR.QC.filtered.' + chrom_num
genotype_obj = load_in_genotype_data(genotype_stem)


####################################################
# Generate gene-weighted ld scores
####################################################
gene_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'  # names of genes
eqtl_sample_sizes = np.asarray([100,200,300,500,1000,'inf'])  # Various eqtl data sets

gene_weighted_ld_score_output_root = simulated_ld_scores_dir + simulation_name_string + '_gene_weighted_ld_scores'  # output root
generate_gene_weighted_ld_scores(regression_snp_names, genotype_obj, eqtl_sample_sizes, gene_summary_file, simulated_learned_gene_models_dir, simulation_name_string,simulated_gene_expression_dir, gene_weighted_ld_score_output_root)


eqtl_sample_sizes = np.asarray([100, 200, 300, 500, 1000])  # Various eqtl data sets
susie_distr_gene_weighted_ld_score_output_root = simulated_ld_scores_dir + simulation_name_string + '_susie_distr_gene_weighted_ld_scores'  # output root
generate_susie_distr_gene_weighted_ld_scores(regression_snp_names, genotype_obj, eqtl_sample_sizes, gene_summary_file, simulated_learned_gene_models_dir, simulation_name_string,simulated_gene_expression_dir, susie_distr_gene_weighted_ld_score_output_root)


