import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import pickle
import rss_likelihood_updates as rss





def get_tissue_names(tissue_name_file):
	aa = np.loadtxt(tissue_name_file,dtype=str,delimiter='\t')
	return aa[1:,0]



def get_window_names(ukkbb_window_summary_file):
	f = open(ukkbb_window_summary_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0] + ':' + data[1] + ':' + data[2]
		arr.append(window_name)
	f.close()
	return np.asarray(arr)

def extract_window_pickle_file_names(trait_name, preprocessed_tgfm_data_dir, window_names, standardize_expression_boolean):
	file_names_arr = []
	count = 0
	for window_name in window_names:
		count = count + 1
		if np.mod(count, 100) == 0:
			print(count)
		if standardize_expression_boolean == "False":
			window_file = preprocessed_tgfm_data_dir + window_name + '_rss_likelihood_' + trait_name + '_data.pkl'
		elif standardize_expression_boolean == "True":
			window_file = preprocessed_tgfm_data_dir + window_name + '_rss_likelihood_' + trait_name + '_standardized_data.pkl'

		# Make sure file exists
		if os.path.isfile(window_file) == False:
			continue

		# load in pickled data
		f = open(window_file, "rb")
		twas_obj = pickle.load(f)
		f.close()

		# Ignore file if there are not many variants
		if len(twas_obj['middle_variant_indices']) < 50:
			continue

		# Ignore incompletely finished window files
		if 'standardized_gwas_beta' not in twas_obj:
			continue

		file_names_arr.append(window_file)

	return np.asarray(file_names_arr)

def create_tissue_to_position_mapping(ordered_tissue_names):
	mapping = {}
	for i, tissue in enumerate(ordered_tissue_names):
		mapping[tissue] = i
	return mapping


def get_tissue_names_from_gene_tissue_string_array(genes):
	tissues = []
	for gene in genes:
		info = gene.split('_')
		tissue = '_'.join(info[1:])
		tissues.append(tissue)
	return np.asarray(tissues)

def create_alpha_mu_and_alpha_var_objects(twas_pickle_file_names, tissue_to_position_mapping, num_tiss):
	twas_obj_file_names = []
	alpha_mu_object = []
	alpha_var_object = []
	beta_mu_object = []
	beta_var_object = []
	tissue_object = []
	valid_gene_object = []
	valid_variant_object = []
	used_genes = {}
	total_genes = np.zeros(num_tiss)
	total_variants = 0.0
	for itera, twas_pickle_file_name in enumerate(twas_pickle_file_names):
		#print(itera)
		f = open(twas_pickle_file_name, "rb")
		twas_obj = pickle.load(f)
		f.close()

		twas_obj_file_names.append(twas_pickle_file_name)

		alpha_mu_object.append(np.zeros(twas_obj['G']))
		alpha_var_object.append(np.ones(twas_obj['G']))
		beta_mu_object.append(np.zeros(twas_obj['K']))
		beta_var_object.append(np.ones(twas_obj['K']))

		total_variants = total_variants + len(twas_obj['middle_variant_indices'])
		
		# get tissue names
		tissue_names = get_tissue_names_from_gene_tissue_string_array(twas_obj['genes'])
		tissue_positions = []
		for i,tissue_name in enumerate(tissue_names):
			tissue_positions.append(tissue_to_position_mapping[tissue_name])
			gene_name = twas_obj['genes'][i]

		for gene_index in twas_obj['middle_gene_indices']:
			total_genes[tissue_to_position_mapping[tissue_names[gene_index]]] = total_genes[tissue_to_position_mapping[tissue_names[gene_index]]] + 1.0

		tissue_object.append(np.asarray(tissue_positions))
		valid_gene_object.append(np.copy(twas_obj['middle_gene_indices']))
		valid_variant_object.append(np.copy(twas_obj['middle_variant_indices']))
	return alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object, tissue_object, valid_gene_object, valid_variant_object, np.asarray(twas_obj_file_names), total_variants, total_genes


def update_alphas_and_betas(alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object, tissue_object, twas_obj_file_names, expected_gamma_alpha, expected_gamma_beta, iteration_svi_window_indices):
	# Number of windows to loop over
	num_windows = len(alpha_mu_object)

	# loop through windows
	for nn in iteration_svi_window_indices:

		# First create vector of length number of genes where element is tissue-specific prior precision corresponding to that gene
		component_gamma_alpha = np.zeros(len(tissue_object[nn]))
		for pos in range(len(tissue_object[nn])):
			component_gamma_alpha[pos] = expected_gamma_alpha[tissue_object[nn][pos]]
		# Load in twas obj
		f = open(twas_obj_file_names[nn], "rb")
		twas_obj = pickle.load(f)
		f.close()

		##########################
		# Compute residual vector
		#########################
		# Initialize residual to beta vectors
		residual = np.copy(twas_obj['standardized_gwas_beta'])
		# Remove predicted effect from each gene (alphas)
		temp_gene_pred = np.zeros(len(residual))
		for g_index in range(twas_obj['G']):
			gene_trait_pred_new = twas_obj['gene_eqtl_pmces'][g_index]*alpha_mu_object[nn][g_index]
			temp_gene_pred = temp_gene_pred - gene_trait_pred_new
		residual = residual + np.dot(twas_obj['srs_inv'], (temp_gene_pred-beta_mu_object[nn]))  # Beta_mu remove predicted effect of genotype

		# Set twas_obj alpha_mu and alpha_var to current estimates of alpha_mu and alpha_var
		alpha_mu = np.copy(alpha_mu_object[nn])
		alpha_var = np.copy(alpha_var_object[nn])

		# Set twas_obj beta_mu and beta_var to current estimates of alpha_mu and alpha_var
		beta_mu = np.copy(beta_mu_object[nn])
		beta_var = np.copy(beta_var_object[nn])

		##########################
		# Update alpha and beta in this window
		#########################
		# Perform VI UPDATES on alpha_mu and alpha_var
		alpha_mu, alpha_var, residual = rss.update_alpha(alpha_mu, alpha_var, residual, twas_obj['gene_eqtl_pmces'], twas_obj['srs_inv'], twas_obj['G'], twas_obj['s_inv_2_diag'], twas_obj['precomputed_a_terms'], component_gamma_alpha)

		# Now set alpha_mu_object and alpha_var_object to current estimates of alpha_mu and alpha_var
		alpha_mu_object[nn] = np.copy(alpha_mu)
		alpha_var_object[nn] = np.copy(alpha_var)

		# Perform VI UPDATES on beta_mu and beta_var
		beta_mu, beta_var, residual = rss.update_beta(beta_mu, beta_var, residual, twas_obj['srs_inv'], twas_obj['K'], twas_obj['D_diag'], twas_obj['s_inv_2_diag'], expected_gamma_beta)

		# Now set beta_mu_object and beta_var_object to current estimates of beta_mu and beta_var
		beta_mu_object[nn] = np.copy(beta_mu)
		beta_var_object[nn] = np.copy(beta_var)

	return alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object

def update_gamma_alpha(gamma_alpha_a, gamma_alpha_b, alpha_mu_object, alpha_var_object, tissue_object, valid_gene_object, num_tiss, ard_a_prior, ard_b_prior, iteration_svi_window_indices, svi_learning_rate, total_num_genes):
	# Number of components to loop over
	num_components = len(alpha_mu_object)
	# Initialize vector to keep track of number of tissue specific genes
	num_genes_vec = np.zeros(num_tiss)
	# Initialize vector to keep track of sum of alpha-squared expected val
	alpha_squared_vec = np.zeros(num_tiss)
	for nn in iteration_svi_window_indices:
		for kk in valid_gene_object[nn]:
			tissue_index = tissue_object[nn][kk]
			num_genes_vec[tissue_index] = num_genes_vec[tissue_index] + 1
			alpha_squared_vec[tissue_index] = alpha_squared_vec[tissue_index] + np.square(alpha_mu_object[nn][kk]) + alpha_var_object[nn][kk]

	gamma_alpha_a_update = ard_a_prior + (num_genes_vec/2.0)*(total_num_genes/num_genes_vec)
	gamma_alpha_b_update = ard_b_prior + (alpha_squared_vec/2.0)*(total_num_genes/num_genes_vec)

	gamma_alpha_a = (svi_learning_rate*gamma_alpha_a_update) + ((1.0 - svi_learning_rate)*gamma_alpha_a)
	gamma_alpha_b = (svi_learning_rate*gamma_alpha_b_update) + ((1.0 - svi_learning_rate)*gamma_alpha_b)

	return gamma_alpha_a, gamma_alpha_b

def update_gamma_beta(gamma_beta_a, gamma_beta_b, beta_mu_object, beta_var_object, valid_variant_object, ard_a_prior, ard_b_prior, iteration_svi_window_indices, svi_learning_rate, total_num_variants):
	# Number of components to loop over
	num_components = len(beta_mu_object)
	# Initialize floats to keep track of beta-squared and total number of variants
	beta_squared = 0.0
	num_variants = 0.0

	# loop over components
	for nn in iteration_svi_window_indices:
		beta_squared = beta_squared + np.sum(np.square(beta_mu_object[nn][valid_variant_object[nn]]) + beta_var_object[nn][valid_variant_object[nn]])
		num_variants = num_variants + len(valid_variant_object[nn])

	# Updates
	gamma_beta_a_update = ard_a_prior + (num_variants/2.0)*(total_num_variants/num_variants)
	gamma_beta_b_update = ard_b_prior + (beta_squared/2.0)*(total_num_variants/num_variants)

	gamma_beta_a = (svi_learning_rate*gamma_beta_a_update) + ((1.0 - svi_learning_rate)*gamma_beta_a)
	gamma_beta_b = (svi_learning_rate*gamma_beta_b_update) + ((1.0 - svi_learning_rate)*gamma_beta_b)


	return gamma_beta_a, gamma_beta_b	


def infer_rss_likelihood_genome_wide_heritabilities(ordered_tissue_names, twas_pickle_file_names, temp_alpha_output_file, temp_beta_output_file, svi_params, max_iter=600):
	# Number of tissues
	num_tiss = len(ordered_tissue_names)
	# Create mapping from tissue name to position
	tissue_to_position_mapping = create_tissue_to_position_mapping(ordered_tissue_names)
	# First create initial alpha mu object
	alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object, tissue_object, valid_gene_object, valid_variant_object, twas_obj_file_names, total_num_variants, total_num_genes = create_alpha_mu_and_alpha_var_objects(twas_pickle_file_names, tissue_to_position_mapping, num_tiss)

	# INITIALIZE Genome-wide parameters (Gammas)
	gamma_alpha_a = np.ones(num_tiss)
	gamma_alpha_b = np.ones(num_tiss)
	gamma_beta_a = 1.0
	gamma_beta_b = 1.0

	# Start variational updates
	for itera in range(max_iter):
		# starting new iteration
		print('Variational iteration ' + str(itera))

		# Compute SVI learning rate for this iteration
		svi_learning_rate = svi_params['tau']/np.power((1.0 + (svi_params['kappa']*itera)), .75)
		# Get array of windows for this SVI iteration
		iteration_svi_window_indices = np.random.choice(np.arange(len(alpha_mu_object)), size=svi_params['batch_size'], replace=False)

		# Take expectations of gamma alpha and gamma beta
		expected_gamma_alpha = gamma_alpha_a/gamma_alpha_b
		expected_gamma_beta = gamma_beta_a/gamma_beta_b

		# Update alpha and beta distributions
		alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object = update_alphas_and_betas(alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object, tissue_object, twas_obj_file_names, expected_gamma_alpha, expected_gamma_beta, iteration_svi_window_indices)

		# Update gamma_alpha
		gamma_alpha_a, gamma_alpha_b = update_gamma_alpha(gamma_alpha_a, gamma_alpha_b, alpha_mu_object, alpha_var_object, tissue_object, valid_gene_object, num_tiss, 1e-16, 1e-16, iteration_svi_window_indices, svi_learning_rate, total_num_genes)

		# Update gamma_beta
		gamma_beta_a, gamma_beta_b = update_gamma_beta(gamma_beta_a, gamma_beta_b, beta_mu_object, beta_var_object, valid_variant_object, 1e-16, 1e-16, iteration_svi_window_indices, svi_learning_rate, total_num_variants)

		# Save temporary results to output
		# Alpha
		df = pd.DataFrame(data={'tissue': ordered_tissue_names, 'expected_precision': (gamma_alpha_a/gamma_alpha_b), 'gamma_alpha_a': gamma_alpha_a, 'gamma_alpha_b': gamma_alpha_b, 'iteration': np.ones(len(gamma_alpha_b))*itera})
		df.to_csv(temp_alpha_output_file, sep='\t', index=False)
		#Beta
		gamma_beta_output_mat = np.vstack((np.asarray(['expected_precision', 'gamma_beta_a', 'gamma_beta_b', 'iteration']), np.asarray([gamma_beta_a/gamma_beta_b, gamma_beta_a, gamma_beta_b, float(itera)]).astype(str)))
		np.savetxt(temp_beta_output_file, gamma_beta_output_mat, fmt="%s", delimiter='\t')



##############################
# Command line args
##############################
trait_name = sys.argv[1]
ukkbb_window_summary_file = sys.argv[2]
tissue_name_file = sys.argv[3]
preprocessed_tgfm_data_dir = sys.argv[4]
standardize_expression_boolean = sys.argv[5]
output_root = sys.argv[6]


# SVI parameters
svi_params = {'batch_size':400, 'tau': .75, 'kappa':.5}

# Set seed
np.random.seed(1)

# Get names of tissues
ordered_tissue_names = get_tissue_names(tissue_name_file)

# Get array of names of windows
window_names = get_window_names(ukkbb_window_summary_file)

# Extract window pickle file names
window_pickle_file_names = extract_window_pickle_file_names(trait_name, preprocessed_tgfm_data_dir, window_names, standardize_expression_boolean)


# Temporary output files used to save intermediate results
temp_alpha_output_file = output_root + 'robust_tissue_specific_prior_precision_temp.txt'
temp_beta_output_file =output_root + 'robust_pleiotropic_prior_precision_temp.txt'

# Run iterative variational algorithm to get heritability estimates
expected_gamma_alpha, gamma_alpha_a, gamma_alpha_b, gamma_beta_a, gamma_beta_b = infer_rss_likelihood_genome_wide_heritabilities(ordered_tissue_names, window_pickle_file_names, temp_alpha_output_file, temp_beta_output_file, svi_params)



