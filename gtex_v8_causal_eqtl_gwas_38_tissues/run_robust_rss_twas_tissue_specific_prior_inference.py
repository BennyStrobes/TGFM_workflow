import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import pickle


def get_tissue_names_from_gene_tissue_string_array(genes):
	tissues = []
	for gene in genes:
		info = gene.split('_')
		tissue = '_'.join(info[1:])
		tissues.append(tissue)
	return np.asarray(tissues)

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

def extract_twas_pickle_file_names(trait_name, pseudotissue_gtex_rss_multivariate_twas_dir, gene_version, fusion_weights):
	file_names = []
	for chrom_num in range(1,23):
		file_name = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + str(chrom_num) + '_summary_tgfm_robust_results_' + gene_version + '_const_1e-5_prior.txt'
		f = open(file_name)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			if fusion_weights == "False":
				file_name = data[2]
			elif fusion_weights == "True":
				file_name = data[3]
			else:
				print('assumption error, unknown fusion_weights variable: ' + fusion_weights)
			if file_name == 'NA':
				continue
			file_names.append(file_name)
		f.close()
	return np.asarray(file_names)

def create_alpha_mu_and_alpha_var_objects(twas_pickle_file_names, tissue_to_position_mapping, gene_count_method, init_version):
	twas_obj_file_names = []
	alpha_mu_object = []
	alpha_var_object = []
	beta_mu_object = []
	beta_var_object = []
	tissue_object = []
	valid_gene_object = []
	used_genes = {}
	print(len(twas_pickle_file_names))
	for itera, twas_pickle_file_name in enumerate(twas_pickle_file_names):
		#print(itera)
		f = open(twas_pickle_file_name, "rb")
		twas_obj = pickle.load(f)
		f.close()

		# Throw out components where expected squared effect sizes are huge (currently believe it is indication of poor model fitting)
		# Need to debug more.
		max_expected_alpha_squared = np.max(np.square(twas_obj.alpha_mu) + twas_obj.alpha_var)
		#if max_expected_alpha_squared >= 1e-3:
			#continue

		twas_obj_file_names.append(twas_pickle_file_name)
		if init_version == 'trained_init':
			alpha_mu_object.append(np.copy(twas_obj.alpha_mu))
			alpha_var_object.append(np.copy(twas_obj.alpha_var))
			beta_mu_object.append(np.copy(twas_obj.beta_mu))
			beta_var_object.append(np.copy(twas_obj.beta_var))
		elif init_version == 'null_init':
			alpha_mu_object.append(np.zeros(len(twas_obj.alpha_mu)))
			alpha_var_object.append(np.ones(len(twas_obj.alpha_var)))
			beta_mu_object.append(np.zeros(len(twas_obj.beta_mu)))
			beta_var_object.append(np.ones(len(twas_obj.beta_var)))
		# get tissue names
		tissue_names = get_tissue_names_from_gene_tissue_string_array(twas_obj.genes)
		tissue_positions = []
		valid_genes = []
		for i,tissue_name in enumerate(tissue_names):
			tissue_positions.append(tissue_to_position_mapping[tissue_name])
			gene_name = twas_obj.genes[i]
			if gene_count_method == 'count_genes_once':
				if gene_name in used_genes:
					#valid_genes.append(False)
					valid_genes.append(True)
				else:
					valid_genes.append(False)
					used_genes[gene_name] = 1
			elif gene_count_method == 'count_all_genes':
				valid_genes.append(True)
		tissue_object.append(np.asarray(tissue_positions))
		valid_gene_object.append(np.asarray(valid_genes))
	return alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object, tissue_object, valid_gene_object, np.asarray(twas_obj_file_names)

def create_tissue_to_position_mapping(ordered_tissue_names):
	mapping = {}
	for i, tissue in enumerate(ordered_tissue_names):
		mapping[tissue] = i
	return mapping

def update_gamma_alpha(alpha_mu_object, alpha_var_object, tissue_object, valid_gene_object, num_tiss, ard_a_prior, ard_b_prior):
	#alpha_squared_expected_val = np.square(self.alpha_mu) + self.alpha_var
	# VI updates
	#self.gamma_alpha_a = self.ard_a_prior + (self.G/2.0)
	#self.gamma_alpha_b = self.ard_b_prior + (np.sum(alpha_squared_expected_val)/2.0)
	
	# Number of components to loop over
	num_components = len(alpha_mu_object)
	# Initialize vector to keep track of number of tissue specific genes
	num_genes_vec = np.zeros(num_tiss)
	# Initialize vector to keep track of sum of alpha-squared expected val
	alpha_squared_vec = np.zeros(num_tiss)
	for nn in range(num_components):
		num_genes = len(tissue_object[nn])
		for kk in range(num_genes):
			tissue_index = tissue_object[nn][kk]
			valid_gene_bool = valid_gene_object[nn][kk]
			# Use this so we don't count genes twice
			if valid_gene_bool:
				num_genes_vec[tissue_index] = num_genes_vec[tissue_index] + 1
				alpha_squared_vec[tissue_index] = alpha_squared_vec[tissue_index] + np.square(alpha_mu_object[nn][kk]) + alpha_var_object[nn][kk]
	gamma_alpha_a = ard_a_prior + (num_genes_vec/2.0)
	gamma_alpha_b = ard_b_prior + (alpha_squared_vec/2.0)
	return gamma_alpha_a, gamma_alpha_b

def update_gamma_beta(beta_mu_object, beta_var_object, ard_a_prior, ard_b_prior):
	# Number of components to loop over
	num_components = len(beta_mu_object)
	# Initialize floats to keep track of beta-squared and total number of variants
	beta_squared = 0.0
	num_variants = 0.0

	# loop over components
	for nn in range(num_components):
		beta_squared_vec = np.square(beta_mu_object[nn]) + beta_var_object[nn]
		beta_squared = beta_squared + np.sum(np.square(beta_mu_object[nn]) + beta_var_object[nn])
		num_variants = num_variants + len(beta_mu_object[nn])
	
	gamma_beta_a = ard_a_prior + (num_variants/2.0)
	gamma_beta_b = ard_b_prior + (beta_squared/2.0)
	return gamma_beta_a, gamma_beta_b	


def update_alphas_and_betas(alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object, tissue_object, twas_obj_file_names, expected_gamma_alpha, expected_gamma_beta):
	# Number of components to loop over
	num_components = len(alpha_mu_object)


	for nn in range(num_components):

		# First create vector of length number of genes where element is tissue-specific prior precision corresponding to that gene
		component_gamma_alpha = np.zeros(len(tissue_object[nn]))
		for pos in range(len(tissue_object[nn])):
			component_gamma_alpha[pos] = expected_gamma_alpha[tissue_object[nn][pos]]
		# Load in twas obj
		f = open(twas_obj_file_names[nn], "rb")
		twas_obj = pickle.load(f)
		f.close()


		# Set twas_obj.residual to new values per new values of alpha_mu
		temp_gene_pred = np.zeros(len(twas_obj.residual))
		for g_index in range(twas_obj.G):
			# Remove effect of the gene corresponding to g_index from the residaul and add new effect
			gene_trait_pred_old = twas_obj.gene_eqtl_pmces[g_index]*twas_obj.alpha_mu[g_index]
			gene_trait_pred_new = twas_obj.gene_eqtl_pmces[g_index]*alpha_mu_object[nn][g_index]
			#twas_obj.residual = twas_obj.residual + np.dot(twas_obj.srs_inv, gene_trait_pred_old) - np.dot(twas_obj.srs_inv, gene_trait_pred_new)
			#twas_obj.residual = twas_obj.residual + np.dot(twas_obj.srs_inv, (gene_trait_pred_old - gene_trait_pred_new))
			temp_gene_pred = temp_gene_pred + (gene_trait_pred_old - gene_trait_pred_new)
		twas_obj.residual = twas_obj.residual + np.dot(twas_obj.srs_inv, temp_gene_pred)

		# Set twas_obj.residual to new values per new values of beta_mu
		twas_obj.residual = twas_obj.residual + np.dot(twas_obj.srs_inv, twas_obj.beta_mu - beta_mu_object[nn])

		# Set twas_obj alpha_mu and alpha_var to current estimates of alpha_mu and alpha_var
		twas_obj.alpha_mu = np.copy(alpha_mu_object[nn])
		twas_obj.alpha_var = np.copy(alpha_var_object[nn])

		# Set twas_obj beta_mu and beta_var to current estimates of alpha_mu and alpha_var
		twas_obj.beta_mu = np.copy(beta_mu_object[nn])
		twas_obj.beta_var = np.copy(beta_var_object[nn])

		# Perform VI UPDATES on alpha_mu and alpha_var
		twas_obj.update_alpha(component_gamma_alpha)

		# Now set alpha_mu_object and alpha_var_object to current estimates of alpha_mu and alpha_var
		alpha_mu_object[nn] = np.copy(twas_obj.alpha_mu)
		alpha_var_object[nn] = np.copy(twas_obj.alpha_var)

		# Perform VI UPDATES on beta_mu and beta_var
		twas_obj.update_beta(expected_gamma_beta)

		# Now set beta_mu_object and beta_var_object to current estimates of beta_mu and beta_var
		beta_mu_object[nn] = np.copy(twas_obj.beta_mu)
		beta_var_object[nn] = np.copy(twas_obj.beta_var)


	return alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object

def infer_rss_twas_tissue_specific_priors(ordered_tissue_names, twas_pickle_file_names, temp_output_file, temp_output_file2, gene_count_method, init_version, max_iter=600):
	# Number of tissues
	num_tiss = len(ordered_tissue_names)
	# Create mapping from tissue name to position
	tissue_to_position_mapping = create_tissue_to_position_mapping(ordered_tissue_names)
	# First create initial alpha mu object
	alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object, tissue_object, valid_gene_object, twas_obj_file_names = create_alpha_mu_and_alpha_var_objects(twas_pickle_file_names, tissue_to_position_mapping, gene_count_method, init_version)
	
	# INITIALIZE GAMMA
	if init_version == 'trained_init':
		# Update gamma
		gamma_alpha_a, gamma_alpha_b = update_gamma_alpha(alpha_mu_object, alpha_var_object, tissue_object, valid_gene_object, num_tiss, 1e-16, 1e-16)
		gamma_beta_a, gamma_beta_b = update_gamma_beta(beta_mu_object, beta_var_object, 1e-16, 1e-16)

	elif init_version == 'null_init':
		gamma_alpha_a = np.ones(num_tiss)
		gamma_alpha_b = np.ones(num_tiss)
		gamma_beta_a = 1.0
		gamma_beta_b = 1.0
	# Start variational updates
	for itera in range(max_iter):
		print('Variational iteration ' + str(itera))

		# Update alpha and beta
		expected_gamma_alpha = gamma_alpha_a/gamma_alpha_b
		expected_gamma_beta = gamma_beta_a/gamma_beta_b
		alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object = update_alphas_and_betas(alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object, tissue_object, twas_obj_file_names, expected_gamma_alpha, expected_gamma_beta)
		

		# Update gamma_alpha
		gamma_alpha_a, gamma_alpha_b = update_gamma_alpha(alpha_mu_object, alpha_var_object, tissue_object, valid_gene_object, num_tiss, 1e-16, 1e-16)

		# Update gamma_beta
		gamma_beta_a, gamma_beta_b = update_gamma_beta(beta_mu_object, beta_var_object, 1e-16, 1e-16)

		print(gamma_alpha_a/gamma_alpha_b)	
		print(gamma_beta_a/gamma_beta_b)	

		if np.mod(itera, 5) == 0:
			df = pd.DataFrame(data={'tissue': ordered_tissue_names, 'expected_precision': (gamma_alpha_a/gamma_alpha_b), 'gamma_alpha_a': gamma_alpha_a, 'gamma_alpha_b': gamma_alpha_b, 'iteration': np.ones(len(gamma_alpha_b))*itera})
			df.to_csv(temp_output_file, sep='\t', index=False)

			gamma_beta_output_mat = np.vstack((np.asarray(['expected_precision', 'gamma_beta_a', 'gamma_beta_b']), np.asarray([gamma_beta_a/gamma_beta_b, gamma_beta_a, gamma_beta_b]).astype(str)))
			np.savetxt(temp_output_file2, gamma_beta_output_mat, fmt="%s", delimiter='\t')


	return gamma_alpha_a/gamma_alpha_b, gamma_alpha_a, gamma_alpha_b, gamma_beta_a, gamma_beta_b


trait_name = sys.argv[1]
gtex_pseudotissue_file = sys.argv[2]
pseudotissue_gtex_rss_multivariate_twas_dir = sys.argv[3]
gene_version = sys.argv[4]
gene_count_method = sys.argv[5]
init_version = sys.argv[6]
fusion_weights = sys.argv[7]

ordered_tissue_names = extract_tissue_names(gtex_pseudotissue_file)

twas_pickle_file_names = extract_twas_pickle_file_names(trait_name, pseudotissue_gtex_rss_multivariate_twas_dir, gene_version, fusion_weights)

temp_output_file = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + gene_version + '_' + gene_count_method + '_' + init_version + '_fusion_weights_' + fusion_weights  + '_robust_tissue_specific_prior_precision_temp.txt'
temp_output_file2 = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + gene_version + '_' + gene_count_method + '_' + init_version + '_fusion_weights_' + fusion_weights  + '_robust_pleiotropic_prior_precision_temp.txt'

expected_gamma_alpha, gamma_alpha_a, gamma_alpha_b, gamma_beta_a, gamma_beta_b = infer_rss_twas_tissue_specific_priors(ordered_tissue_names, twas_pickle_file_names, temp_output_file, temp_output_file2, gene_count_method, init_version)


output_file = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + gene_version + '_' + gene_count_method + '_' + init_version + '_fusion_weights_' + fusion_weights  + '_robust_tissue_specific_prior_precision.txt'
df = pd.DataFrame(data={'tissue': ordered_tissue_names, 'expected_precision': expected_gamma_alpha, 'gamma_alpha_a': gamma_alpha_a, 'gamma_alpha_b': gamma_alpha_b})
df.to_csv(output_file, sep='\t', index=False)

output_file2 = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + gene_version + '_' + gene_count_method + '_' + init_version + '_fusion_weights_' + fusion_weights  + '_robust_pleiotropic_prior_precision.txt'
gamma_beta_output_mat = np.vstack((np.asarray(['expected_precision', 'gamma_beta_a', 'gamma_beta_b']), np.asarray([gamma_beta_a/gamma_beta_b, gamma_beta_a, gamma_beta_b]).astype(str)))
np.savetxt(output_file2, gamma_beta_output_mat, fmt="%s", delimiter='\t')


