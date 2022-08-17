import sys
sys.path.remove('/n/app/python/3.6.0/lib/python3.6/site-packages')
import pandas as pd
import numpy as np 
import os
import tensorflow as tf
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


def get_gene_names_from_gene_tissue_string_array(genes):
	tissues = []
	gene_arr = []
	for gene in genes:
		info = gene.split('_')
		gene_name = info[0]
		gene_arr.append(gene_name)
	return np.asarray(gene_arr)

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


def extract_twas_pickle_file_names_and_save_to_new_dir(trait_name, pseudotissue_gtex_rss_multivariate_twas_dir, gene_version, fusion_weights, new_pickle_directory):
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

			# Load pickle
			f = open(file_name, "rb")
			twas_obj = pickle.load(f)
			f.close()

			# Re-save pickle
			new_file_name = new_pickle_directory + file_name.split('/')[-1]
			file_pi = open(new_file_name, 'wb') 
			pickle.dump(twas_obj, file_pi)
			file_pi.close()

			file_names.append(new_file_name)
		f.close()
	return np.asarray(file_names)

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

def create_alpha_mu_and_alpha_var_objects(twas_pickle_file_names, tissue_to_position_mapping, gene_count_method, init_version, gene_to_annotation_vector, num_tiss):
	twas_obj_file_names = []
	alpha_mu_object = []
	alpha_var_object = []
	beta_mu_object = []
	beta_var_object = []
	gene_anno_object = []
	tissue_object = []
	tissue_anno_object = []
	valid_gene_object = []
	used_genes = {}
	num_anno = len(gene_to_annotation_vector[[*gene_to_annotation_vector][0]])
	null_anno = np.zeros(num_anno)
	null_anno[0] = 1.0
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
		# Get Gene info
		gene_names = get_gene_names_from_gene_tissue_string_array(twas_obj.genes)
		gene_anno = []
		for gene_name in gene_names:
			if gene_name.split('.')[0] not in gene_to_annotation_vector:
				gene_anno.append(np.copy(null_anno))
			else:
				gene_anno.append(gene_to_annotation_vector[gene_name.split('.')[0]])
		gene_anno_object.append(np.asarray(gene_anno))
		tissue_names = get_tissue_names_from_gene_tissue_string_array(twas_obj.genes)
		# Get tissue info
		tissue_positions = []
		valid_genes = []
		tissue_anno = []
		for i,tissue_name in enumerate(tissue_names):
			pos = tissue_to_position_mapping[tissue_name] + 1
			temp_tissue_anno = np.zeros(num_tiss+1)
			temp_tissue_anno[0] = 1.0
			temp_tissue_anno[pos] = 1.0
			tissue_anno.append(temp_tissue_anno)

			tissue_positions.append(tissue_to_position_mapping[tissue_name])
			gene_name = twas_obj.genes[i]
			if gene_count_method == 'count_genes_once':
				if gene_name in used_genes:
					valid_genes.append(False)
				else:
					valid_genes.append(True)
					used_genes[gene_name] = 1
			elif gene_count_method == 'count_all_genes':
				valid_genes.append(True)
		tissue_object.append(np.asarray(tissue_positions))
		tissue_anno_object.append(np.asarray(tissue_anno))
		valid_gene_object.append(np.asarray(valid_genes))
	return alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object, gene_anno_object, tissue_anno_object, tissue_object, valid_gene_object, np.asarray(twas_obj_file_names)

def create_tissue_to_position_mapping(ordered_tissue_names):
	mapping = {}
	for i, tissue in enumerate(ordered_tissue_names):
		mapping[tissue] = i
	return mapping

def update_non_linear_gamma_alpha(alpha_mu_object, alpha_var_object, tissue_anno_object, gene_anno_object, valid_gene_object, num_tiss):
	# initialize tensorflow model
	non_linear_gamma_alpha_model = init_nonlinear_gamma_alpha_model(tissue_anno_object[0].shape[1], gene_anno_object[0].shape[1])
	# Initialize vector to keep track of sum of alpha-squared expected val
	alpha_squared_vec = []
	# Initialize concatentated tissue and gene annotation matrix
	tiss_anno_mat = []
	gene_anno_mat = []
	# Number of components to loop over
	num_components = len(alpha_mu_object)

	for nn in range(num_components):
		num_genes = tissue_anno_object[nn].shape[0]
		for kk in range(num_genes):
			valid_gene_bool = valid_gene_object[nn][kk]
			if valid_gene_bool:
				gene_alpha_squared = np.square(alpha_mu_object[nn][kk]) + alpha_var_object[nn][kk]
				alpha_squared_vec.append(gene_alpha_squared)

				tiss_anno = tissue_anno_object[nn][kk,:]
				gene_anno = gene_anno_object[nn][kk,:]

				tiss_anno_mat.append(tiss_anno)
				gene_anno_mat.append(gene_anno)


	# Put in matrix form
	tiss_anno_mat = np.asarray(tiss_anno_mat)
	gene_anno_mat = np.asarray(gene_anno_mat)
	alpha_squared_vec = np.asarray(alpha_squared_vec)
	anno_mat = [tiss_anno_mat, gene_anno_mat]

	non_linear_gamma_alpha_model.fit(anno_mat, alpha_squared_vec, epochs=50)

	return non_linear_gamma_alpha_model


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


def update_alphas_and_betas(alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object, tissue_object, tissue_anno_object, gene_anno_object, twas_obj_file_names, non_linear_gamma_alpha_model, expected_gamma_beta):
	# Number of components to loop over
	num_components = len(alpha_mu_object)


	for nn in range(num_components):
		if non_linear_gamma_alpha_model == None:
			component_gamma_alpha = np.ones(len(tissue_object[nn]))
		else:
			component_var_alpha = non_linear_gamma_alpha_model.predict([tissue_anno_object[nn],gene_anno_object[nn]])[:,0]
			component_gamma_alpha = 1.0/component_var_alpha
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

def create_gene_to_annotation_vector_mapping(gene_set_annotation_file):
	gene_sets = ['Universe', 'MGI_essential', 'Haploinsufficient', 'highPLI_Exac', 'highShet_Cassa', 'Kinase_Uniprot', 'LoF_twohits', 'OMIM_autosomaldominant', 'Olfactory']
	indices = []
	f = open(gene_set_annotation_file)
	head_count = 0
	mapping = {}
	for line in f:
		line = line.rstrip()
		data = line.split(',')
		if len(data) != 165:
			print('assumption error')
		if head_count == 0:
			head_count = head_count + 1
			for gene_set in gene_sets:
				indices.append(np.where(np.asarray(data)==gene_set)[0][0])
			indices = np.asarray(indices)
			continue
		for gene_set_index, gene_set in enumerate(gene_sets):
			col_index = indices[gene_set_index]
			gene_name = data[col_index]
			if gene_name == '':
				continue
			if gene_name not in mapping:
				mapping[gene_name] = np.zeros(len(gene_sets))
			mapping[gene_name][gene_set_index] = 1.0
	f.close()
	return mapping, gene_sets

def gaussian_neg_log_likelihood_tf_padded_loss(beta_squared_true, gamma_pred):
	epsilon = 1e-30
	nll = tf.math.log(gamma_pred + epsilon) + tf.math.divide(beta_squared_true, (gamma_pred+epsilon))
	return nll

def init_nonlinear_gamma_alpha_model(num_tissue_anno, num_gene_anno):
	input_tissue = tf.keras.layers.Input(shape=(num_tissue_anno,))
	input_gene = tf.keras.layers.Input(shape=(num_gene_anno,))
	# the first branch operates on the first input
	x = tf.keras.layers.Dense(1, activation="linear")(input_tissue)
	x = tf.keras.models.Model(inputs=input_tissue, outputs=x)

	y = tf.keras.layers.Dense(1, activation="linear")(input_gene)
	y = tf.keras.models.Model(inputs=input_gene, outputs=y)

	combined = tf.keras.layers.concatenate([x.output, y.output])

	z = tf.keras.layers.Dense(1, activation="softplus")(combined)

	model = tf.keras.models.Model(inputs=[x.input, y.input], outputs=z)

	model.compile(loss=gaussian_neg_log_likelihood_tf_padded_loss, optimizer='adam')

	return model

def extract_expected_gene_tissue_precision_df(ordered_tissue_names, gene_set_names, non_linear_gamma_alpha_model, itera):
	tiss_arr = []
	gene_set_arr = []
	precision_arr = []

	num_tiss = len(ordered_tissue_names)
	num_gene_set = len(gene_set_names)

	gene_anno = []
	tissue_anno = []

	for ii, tissue_name in enumerate(ordered_tissue_names):
		tissue_vec = np.zeros(num_tiss +1)
		tissue_vec[0] = 1.0
		tissue_vec[(ii+1)] = 1.0

		for jj, gene_set_name in enumerate(gene_set_names):
			gene_vec = np.zeros(num_gene_set)
			gene_vec[0] = 1.0
			gene_vec[jj] = 1.0
			#non_linear_gamma_alpha_model.predict([tissue_vec, gene_vec])
			#non_linear_gamma_alpha_model.predict([np.reshape(tissue_vec, (1,len(tissue_vec))), np.reshape(gene_vec, (1,len(gene_vec)))])
			tiss_arr.append(tissue_name)
			gene_set_arr.append(gene_set_name)

			tissue_anno.append(tissue_vec)
			gene_anno.append(gene_vec)

	pred_var = non_linear_gamma_alpha_model.predict([np.asarray(tissue_anno), np.asarray(gene_anno)])[:,0]
	df = pd.DataFrame(data={'tissue': np.asarray(tiss_arr), 'gene_set': np.asarray(gene_set_arr), 'expected_precision': (1.0/pred_var),  'iteration': np.ones(len(pred_var))*itera})
	return df


def infer_rss_twas_tissue_specific_priors(ordered_tissue_names, twas_pickle_file_names, temp_output_file, temp_output_file2, gene_count_method, init_version, gene_set_annotation_file, max_iter=600):
	# Create gene to annotation mapping
	gene_to_annotation_vector, gene_set_names = create_gene_to_annotation_vector_mapping(gene_set_annotation_file)

	# Number of tissues
	num_tiss = len(ordered_tissue_names)
	# Create mapping from tissue name to position
	tissue_to_position_mapping = create_tissue_to_position_mapping(ordered_tissue_names)
	# First create initial alpha mu object
	alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object, gene_anno_object, tissue_anno_object, tissue_object, valid_gene_object, twas_obj_file_names = create_alpha_mu_and_alpha_var_objects(twas_pickle_file_names, tissue_to_position_mapping, gene_count_method, init_version, gene_to_annotation_vector, num_tiss)
	

	# INITIALIZE GAMMA
	if init_version == 'trained_init':
		# Update gamma
		print('not currently implemented')
		pdb.set_trace()

	elif init_version == 'null_init':
		#non_linear_gamma_alpha_model = init_nonlinear_gamma_alpha_model(tissue_anno_object[0].shape[1], gene_anno_object[0].shape[1])
		non_linear_gamma_alpha_model = None
		gamma_beta_a = 1.0
		gamma_beta_b = 1e-8
	# Start variational updates
	for itera in range(max_iter):
		print('Variational iteration ' + str(itera))

		# Update alpha and beta
		#expected_gamma_alpha = gamma_alpha_a/gamma_alpha_b
		expected_gamma_beta = gamma_beta_a/gamma_beta_b
		alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object = update_alphas_and_betas(alpha_mu_object, alpha_var_object, beta_mu_object, beta_var_object, tissue_object, tissue_anno_object, gene_anno_object, twas_obj_file_names, non_linear_gamma_alpha_model, expected_gamma_beta)
		

		# Update gamma_alpha
		#gamma_alpha_a, gamma_alpha_b = update_gamma_alpha(alpha_mu_object, alpha_var_object, tissue_object, valid_gene_object, num_tiss, 1e-16, 1e-16)
		non_linear_gamma_alpha_model = update_non_linear_gamma_alpha(alpha_mu_object, alpha_var_object, tissue_anno_object, gene_anno_object, valid_gene_object, num_tiss)

		# Update gamma_beta
		gamma_beta_a, gamma_beta_b = update_gamma_beta(beta_mu_object, beta_var_object, 1e-16, 1e-16)

		#print(gamma_alpha_a/gamma_alpha_b)	
		print(gamma_beta_a/gamma_beta_b)	

		if np.mod(itera, 5) == 0:
			expected_precision_df = extract_expected_gene_tissue_precision_df(ordered_tissue_names, gene_set_names, non_linear_gamma_alpha_model, itera)
			expected_precision_df.to_csv(temp_output_file, sep='\t', index=False)

			gamma_beta_output_mat = np.vstack((np.asarray(['expected_precision', 'gamma_beta_a', 'gamma_beta_b']), np.asarray([gamma_beta_a/gamma_beta_b, gamma_beta_a, gamma_beta_b]).astype(str)))
			np.savetxt(temp_output_file2, gamma_beta_output_mat, fmt="%s", delimiter='\t')

			# Save model
			non_linear_gamma_alpha_model.save(temp_output_file + '_tf_model')


	return gamma_alpha_a/gamma_alpha_b, gamma_alpha_a, gamma_alpha_b, gamma_beta_a, gamma_beta_b


trait_name = sys.argv[1]
gtex_pseudotissue_file = sys.argv[2]
pseudotissue_gtex_rss_multivariate_twas_dir = sys.argv[3]
gene_version = sys.argv[4]
gene_count_method = sys.argv[5]
init_version = sys.argv[6]
fusion_weights = sys.argv[7]
gene_set_annotation_file = sys.argv[8]


ordered_tissue_names = extract_tissue_names(gtex_pseudotissue_file)


new_pickle_directory = '/n/scratch3/users/b/bes710/causal_eqtl_gwas/component_data_scratch2/'
twas_pickle_file_names = extract_twas_pickle_file_names_and_save_to_new_dir(trait_name, pseudotissue_gtex_rss_multivariate_twas_dir, gene_version, fusion_weights, new_pickle_directory)

#twas_pickle_file_names = extract_twas_pickle_file_names(trait_name, pseudotissue_gtex_rss_multivariate_twas_dir, gene_version, fusion_weights)



temp_output_file = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + gene_version + '_' + gene_count_method + '_' + init_version + '_fusion_weights_' + fusion_weights  + '_robust_tissue_and_gene_specific_prior_precision_temp2.txt'
temp_output_file2 = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + gene_version + '_' + gene_count_method + '_' + init_version + '_fusion_weights_' + fusion_weights  + '_robust_tissue_and_gene_pleiotropic_prior_precision_temp2.txt'

expected_gamma_alpha, gamma_alpha_a, gamma_alpha_b, gamma_beta_a, gamma_beta_b = infer_rss_twas_tissue_specific_priors(ordered_tissue_names, twas_pickle_file_names, temp_output_file, temp_output_file2, gene_count_method, init_version, gene_set_annotation_file)


output_file = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + gene_version + '_' + gene_count_method + '_' + init_version + '_fusion_weights_' + fusion_weights  + '_robust_tissue_and_gene_specific_prior_precision2.txt'
df = pd.DataFrame(data={'tissue': ordered_tissue_names, 'expected_precision': expected_gamma_alpha, 'gamma_alpha_a': gamma_alpha_a, 'gamma_alpha_b': gamma_alpha_b})
df.to_csv(output_file, sep='\t', index=False)

output_file2 = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + gene_version + '_' + gene_count_method + '_' + init_version + '_fusion_weights_' + fusion_weights  + '_robust_tissue_and_gene_pleiotropic_prior_precision2.txt'
gamma_beta_output_mat = np.vstack((np.asarray(['expected_precision', 'gamma_beta_a', 'gamma_beta_b']), np.asarray([gamma_beta_a/gamma_beta_b, gamma_beta_a, gamma_beta_b]).astype(str)))
np.savetxt(output_file2, gamma_beta_output_mat, fmt="%s", delimiter='\t')


