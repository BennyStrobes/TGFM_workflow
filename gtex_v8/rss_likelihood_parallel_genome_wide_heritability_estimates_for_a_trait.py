import sys
sys.path.remove('/n/app/python/3.6.0/lib/python3.6/site-packages')
import numpy as np 
import pandas as pd
import os
import pdb
import tensorflow as tf
import gzip
import time
from numba import njit, prange, jit
from joblib import Parallel, delayed
from scipy.stats import gamma
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
from joblib.externals.loky import set_loky_pickler
from joblib import wrap_non_picklable_objects
import pickle
import rss_likelihood_updates as rss





def get_tissue_names(tissue_name_file):
	aa = np.loadtxt(tissue_name_file,dtype=str,delimiter='\t')
	return aa[1:,0]



def get_window_names(ukkbb_window_summary_file, preprocessed_tgfm_data_dir, standardize_expression_boolean):
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
		
		# Check if window_file_name was created
		temp_trait_name = 'repro_MENARCHE_AGE'
		if standardize_expression_boolean == "False":
			window_file_name = preprocessed_tgfm_data_dir + window_name + '_rss_likelihood_' + temp_trait_name + '_data.pkl'
		elif standardize_expression_boolean == "True":
			window_file_name = preprocessed_tgfm_data_dir + window_name + '_rss_likelihood_' + temp_trait_name + '_standardized_data.pkl'

		if os.path.isfile(window_file_name) == False:
			continue

		arr.append(window_name)
	f.close()
	return np.asarray(arr)

def extract_window_pickle_file_names(trait_name, preprocessed_tgfm_data_dir, window_names, standardize_expression_boolean):
	file_names_arr = []
	shared_file_names_arr = []
	count = 0
	new_window_names = []
	for window_name in window_names:
		count = count + 1
		if np.mod(count, 100) == 0:
			print(count)
		if standardize_expression_boolean == "False":
			window_file = preprocessed_tgfm_data_dir + window_name + '_rss_likelihood_' + trait_name + '_data.pkl'
			shared_window_file = preprocessed_tgfm_data_dir + window_name + '_rss_likelihood_shared_data.pkl'

		elif standardize_expression_boolean == "True":
			window_file = preprocessed_tgfm_data_dir + window_name + '_rss_likelihood_' + trait_name + '_standardized_data.pkl'
			shared_window_file = preprocessed_tgfm_data_dir + window_name + '_rss_likelihood_shared_standardized_data.pkl'
		# Make sure file exists
		if os.path.isfile(window_file) == False:
			continue

		# load in pickled data
		f = open(window_file, "rb")
		twas_obj = pickle.load(f)
		f.close()

		f = open(shared_window_file, "rb")
		shared_twas_obj = pickle.load(f)
		f.close()

		# Ignore file if there are not many variants
		if len(shared_twas_obj['middle_variant_indices']) < 5:
			continue

		# Ignore incompletely finished window files
		if 'standardized_gwas_beta' not in twas_obj:
			continue

		new_window_names.append(window_name)
		file_names_arr.append(window_file)
		shared_file_names_arr.append(shared_window_file)

	return np.asarray(new_window_names), np.asarray(file_names_arr), np.asarray(shared_file_names_arr)

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

def create_alpha_mu_and_alpha_var_objects(twas_pickle_file_names, shared_twas_pickle_file_names, tissue_to_position_mapping):
	twas_obj_file_names = []
	alpha_mu_object = []
	alpha_var_object = []
	beta_mu_object = []
	beta_var_object = []
	tissue_object = []
	valid_gene_object = []
	valid_variant_object = []
	used_genes = {}
	for itera, twas_pickle_file_name in enumerate(twas_pickle_file_names):
		#print(itera)
		f = open(twas_pickle_file_name, "rb")
		twas_obj = pickle.load(f)
		f.close()

		f = open(shared_twas_pickle_file_names[itera], "rb")
		shared_twas_obj = pickle.load(f)
		f.close()		

		alpha_mu_object.append(np.zeros(shared_twas_obj['G']))
		alpha_var_object.append(np.ones(shared_twas_obj['G']))
		beta_mu_object.append(np.zeros(shared_twas_obj['K']))
		beta_var_object.append(np.ones(shared_twas_obj['K']))
		
		# get tissue names
		tissue_names = get_tissue_names_from_gene_tissue_string_array(shared_twas_obj['genes'])
		tissue_positions = []
		for i,tissue_name in enumerate(tissue_names):
			tissue_positions.append(tissue_to_position_mapping[tissue_name])
			gene_name = shared_twas_obj['genes'][i]
		tissue_object.append(np.asarray(tissue_positions))
		valid_gene_object.append(np.copy(shared_twas_obj['middle_gene_indices']))
		valid_variant_object.append(np.copy(shared_twas_obj['middle_variant_indices']))
	return alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object, tissue_object, valid_gene_object, valid_variant_object

#def update_alpha_and_beta_in_a_single_window(alpha_mu_object[nn], alpha_var_object[nn], beta_mu_object[nn], beta_var_object[nn], tissue_object[nn], twas_obj_file_names[nn], expected_gamma_alpha, expected_gamma_beta):
def update_alpha_and_beta_in_a_single_window(alpha_mu_w, alpha_var_w, beta_mu_w, beta_var_w, tissue_object_w, twas_obj_file_name_w, shared_twas_obj_file_name_w, expected_gamma_alpha, expected_gamma_beta):
	# First create vector of length number of genes where element is tissue-specific prior precision corresponding to that gene
	component_gamma_alpha = np.zeros(len(tissue_object_w))
	for pos in range(len(tissue_object_w)):
		component_gamma_alpha[pos] = expected_gamma_alpha[tissue_object_w[pos]]
	# Load in twas obj
	f = open(twas_obj_file_name_w, "rb")
	twas_obj = pickle.load(f)
	f.close()
	# Load in shared twas obj
	f = open(shared_twas_obj_file_name_w, "rb")
	shared_twas_obj = pickle.load(f)
	f.close()

	# Compute SRS^-1/2
	srs_inv_mat = np.multiply(np.multiply(twas_obj['s_diag'][:, None], shared_twas_obj['reference_ld']), (1.0/(twas_obj['s_diag'])))

	##########################
	# Compute residual vector
	#########################
	# Initialize residual to beta vectors
	residual = np.copy(twas_obj['standardized_gwas_beta'])
	# Remove predicted effect from each gene (alphas)
	temp_gene_pred = np.zeros(len(residual))
	for g_index in range(shared_twas_obj['G']):
		gene_trait_pred_new = shared_twas_obj['gene_eqtl_pmces'][g_index]*alpha_mu_w[g_index]
		temp_gene_pred = temp_gene_pred - gene_trait_pred_new
	residual = residual + np.dot(srs_inv_mat, (temp_gene_pred-beta_mu_w))  # Beta_mu remove predicted effect of genotype

	##########################
	# Update alpha and beta in this window
	#########################
	# Perform VI UPDATES on alpha_mu and alpha_var
	alpha_mu_w, alpha_var_w, residual = rss.update_alpha(alpha_mu_w, alpha_var_w, residual, shared_twas_obj['gene_eqtl_pmces'], srs_inv_mat, shared_twas_obj['G'], twas_obj['s_inv_2_diag'], twas_obj['precomputed_a_terms'], component_gamma_alpha)

	# Perform VI UPDATES on beta_mu and beta_var
	beta_mu_w, beta_var_w, residual = rss.update_beta(beta_mu_w, beta_var_w, residual, srs_inv_mat, shared_twas_obj['K'], twas_obj['D_diag'], twas_obj['s_inv_2_diag'], expected_gamma_beta)

	return (alpha_mu_w, alpha_var_w, beta_mu_w, beta_var_w)

def update_alphas_and_betas(alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object, tissue_object, twas_obj_file_names, shared_twas_obj_file_names, expected_gamma_alpha, expected_gamma_beta):
	# Number of windows to loop over
	num_windows = len(alpha_mu_object)

	parr = False
	parr = True

	print('start')

	if parr == False:
		window_update_data = []
		# loop through windows
		for nn in range(num_windows):
			window_update_data.append(update_alpha_and_beta_in_a_single_window(alpha_mu_object[nn], alpha_var_object[nn], beta_mu_object[nn], beta_var_object[nn], tissue_object[nn], twas_obj_file_names[nn], shared_twas_obj_file_names[nn], expected_gamma_alpha, expected_gamma_beta))
	elif parr == True:
		window_update_data = Parallel(n_jobs=20)(delayed(update_alpha_and_beta_in_a_single_window)(alpha_mu_object[nn], alpha_var_object[nn], beta_mu_object[nn], beta_var_object[nn], tissue_object[nn], twas_obj_file_names[nn], shared_twas_obj_file_names[nn], expected_gamma_alpha, expected_gamma_beta) for nn in range(num_windows))


	for nn in range(num_windows):
		alpha_mu_object[nn] = window_update_data[nn][0]
		alpha_var_object[nn] = window_update_data[nn][1]
		beta_mu_object[nn] = window_update_data[nn][2]
		beta_var_object[nn] = window_update_data[nn][3]

	return alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object

def update_gamma_alpha(alpha_mu_object, alpha_var_object, tissue_object, valid_gene_object, num_tiss, ard_a_prior, ard_b_prior):
	# Number of components to loop over
	num_components = len(alpha_mu_object)
	# Initialize vector to keep track of number of tissue specific genes
	num_genes_vec = np.zeros(num_tiss)
	# Initialize vector to keep track of sum of alpha-squared expected val
	alpha_squared_vec = np.zeros(num_tiss)
	for nn in range(num_components):
		for kk in valid_gene_object[nn]:
			tissue_index = tissue_object[nn][kk]
			num_genes_vec[tissue_index] = num_genes_vec[tissue_index] + 1
			alpha_squared_vec[tissue_index] = alpha_squared_vec[tissue_index] + np.square(alpha_mu_object[nn][kk]) + alpha_var_object[nn][kk]
	gamma_alpha_a = ard_a_prior + (num_genes_vec/2.0)
	gamma_alpha_b = ard_b_prior + (alpha_squared_vec/2.0)
	return gamma_alpha_a, gamma_alpha_b

def update_gamma_beta(beta_mu_object, beta_var_object, valid_variant_object, ard_a_prior, ard_b_prior):
	# Number of components to loop over
	num_components = len(beta_mu_object)
	# Initialize floats to keep track of beta-squared and total number of variants
	beta_squared = 0.0
	num_variants = 0.0

	# loop over components
	for nn in range(num_components):
		beta_squared = beta_squared + np.sum(np.square(beta_mu_object[nn][valid_variant_object[nn]]) + beta_var_object[nn][valid_variant_object[nn]])
		num_variants = num_variants + len(valid_variant_object[nn])
	
	gamma_beta_a = ard_a_prior + (num_variants/2.0)
	gamma_beta_b = ard_b_prior + (beta_squared/2.0)
	return gamma_beta_a, gamma_beta_b	


def infer_rss_likelihood_genome_wide_heritabilities(ordered_tissue_names, twas_pickle_file_names, shared_twas_pickle_file_names, temp_alpha_output_file, temp_beta_output_file, max_iter=1000):
	# Number of tissues
	num_tiss = len(ordered_tissue_names)
	# Create mapping from tissue name to position
	tissue_to_position_mapping = create_tissue_to_position_mapping(ordered_tissue_names)
	# First create initial alpha mu object
	alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object, tissue_object, valid_gene_object, valid_variant_object = create_alpha_mu_and_alpha_var_objects(twas_pickle_file_names, shared_twas_pickle_file_names, tissue_to_position_mapping)
	
	# INITIALIZE Genome-wide parameters (Gammas)
	gamma_alpha_a = np.ones(num_tiss)
	gamma_alpha_b = np.ones(num_tiss)*1e-6
	#gamma_alpha_b = np.ones(num_tiss)*1e-5
	gamma_beta_a = 1.0
	gamma_beta_b = 5e-8

	# Start variational updates
	for itera in range(max_iter):
		# starting new iteration
		print('Variational iteration ' + str(itera))
		print(time.time())

		# Take expectations of gamma alpha and gamma beta
		expected_gamma_alpha = gamma_alpha_a/gamma_alpha_b
		expected_gamma_beta = gamma_beta_a/gamma_beta_b

		# Update alpha and beta distributions
		alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object = update_alphas_and_betas(alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object, tissue_object, twas_pickle_file_names, shared_twas_pickle_file_names, expected_gamma_alpha, expected_gamma_beta)

		if itera > 5:
			# Update gamma_alpha
			gamma_alpha_a, gamma_alpha_b = update_gamma_alpha(alpha_mu_object, alpha_var_object, tissue_object, valid_gene_object, num_tiss, 1e-16, 1e-16)

			# Update gamma_beta
			gamma_beta_a, gamma_beta_b = update_gamma_beta(beta_mu_object, beta_var_object, valid_variant_object, 1e-16, 1e-16)

		# Save temporary results to output
		# Alpha
		df = pd.DataFrame(data={'tissue': ordered_tissue_names, 'expected_precision': (gamma_alpha_a/gamma_alpha_b), 'gamma_alpha_a': gamma_alpha_a, 'gamma_alpha_b': gamma_alpha_b, 'iteration': np.ones(len(gamma_alpha_b))*itera})
		df.to_csv(temp_alpha_output_file, sep='\t', index=False)
		#Beta
		gamma_beta_output_mat = np.vstack((np.asarray(['expected_precision', 'gamma_beta_a', 'gamma_beta_b', 'iteration']), np.asarray([gamma_beta_a/gamma_beta_b, gamma_beta_a, gamma_beta_b, float(itera)]).astype(str)))
		np.savetxt(temp_beta_output_file, gamma_beta_output_mat, fmt="%s", delimiter='\t')


def extract_top_n_windows(window_names, trait_name, preprocessed_tgfm_data_dir, standardize_expression_boolean, num_top):
	# Keep track of maximum z score per window
	window_max_z_score_arr = []
	# Loop through windows
	for window_name in window_names:

		if standardize_expression_boolean == "False":
			window_file = preprocessed_tgfm_data_dir + window_name + '_rss_likelihood_' + trait_name + '_data.pkl'
		elif standardize_expression_boolean == "True":
			window_file = preprocessed_tgfm_data_dir + window_name + '_rss_likelihood_' + trait_name + '_standardized_data.pkl'
		# Make sure file exists
		if os.path.isfile(window_file) == False:
			print('assumption eroror')
			pdb.set_trace()
		# load in pickled data
		f = open(window_file, "rb")
		twas_obj = pickle.load(f)
		f.close()

		# Get window z scores
		window_z_scores = twas_obj['gwas_beta']/twas_obj['gwas_beta_se']

		# keep track of abs max z score
		window_max_z_score_arr.append(np.max(np.abs(window_z_scores)))
	window_max_z_score_arr = np.asarray(window_max_z_score_arr)

	# z score thresh
	z_score_thresh = -np.sort(-window_max_z_score_arr)[num_top]

	filtered_window_names = []
	for window_iter, window_name in enumerate(window_names):
		if window_max_z_score_arr[window_iter] > z_score_thresh:
			filtered_window_names.append(window_name)
	filtered_window_names = np.asarray(filtered_window_names)

	
	return filtered_window_names


##############################
# Command line args
##############################
trait_name = sys.argv[1]
ukkbb_window_summary_file = sys.argv[2]
tissue_name_file = sys.argv[3]
preprocessed_tgfm_data_dir = sys.argv[4]
standardize_expression_boolean = sys.argv[5]
output_root = sys.argv[6]



# Get names of tissues
ordered_tissue_names = get_tissue_names(tissue_name_file)

# Get array of names of windows
window_names = get_window_names(ukkbb_window_summary_file, preprocessed_tgfm_data_dir, standardize_expression_boolean)
window_names = extract_top_n_windows(window_names, trait_name, preprocessed_tgfm_data_dir, standardize_expression_boolean, 500)

# Extract window pickle file names
window_names, window_pickle_file_names, shared_window_pickle_file_names = extract_window_pickle_file_names(trait_name, preprocessed_tgfm_data_dir, window_names, standardize_expression_boolean)

# Temporary output files used to save intermediate results
temp_alpha_output_file = output_root + 'robust_tissue_specific_prior_precision_temp.txt'
temp_beta_output_file =output_root + '_robust_pleiotropic_prior_precision_temp.txt'

# Run iterative variational algorithm to get heritability estimates
expected_gamma_alpha, gamma_alpha_a, gamma_alpha_b, gamma_beta_a, gamma_beta_b = infer_rss_likelihood_genome_wide_heritabilities(ordered_tissue_names, window_pickle_file_names, shared_window_pickle_file_names, temp_alpha_output_file, temp_beta_output_file)



