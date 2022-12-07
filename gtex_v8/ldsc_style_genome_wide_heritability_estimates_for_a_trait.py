import sys
sys.path.remove('/n/app/python/3.6.0/lib/python3.6/site-packages')
import numpy as np 
import pandas as pd
import os
import pdb
import tensorflow as tf
import tensorflow_recommenders as tfrs
import gzip
import time
import scipy.optimize
import scipy.stats


def get_window_names(ukkbb_window_summary_file, preprocessed_tgfm_data_dir):
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
		window_file_name = preprocessed_tgfm_data_dir + gene_type + '_' + window_name + '_rss_likelihood_repro_MENARCHE_AGE_standardized_data.pkl'
		if os.path.isfile(window_file_name) == False:
			continue
		arr.append(window_name)
	f.close()
	return np.asarray(arr)

def load_in_ldsc_style_data(preprocessed_tgfm_data_dir, trait_name, window_names, gene_model_suffix, window_chi_sq_lb=-1, window_chi_sq_ub=100000000000000000.0):
	# Initialize output
	X = []
	y = []
	weights = []
	num_windows = 0

	# Loop through windows
	for window_name in window_names:
		# Input files for this window
		ldscore_file = preprocessed_tgfm_data_dir + gene_type + '_' + window_name + '_tgfm_ldscore_annotation_file' + gene_model_suffix + '.txt'
		#ldscore_file = preprocessed_tgfm_data_dir + window_name + '_tgfm_ldscore_not_standardized_annotation_file.txt'

		chi_sq_file = preprocessed_tgfm_data_dir + gene_type + '_' + window_name + '_' + trait_name + '_tgfm_ldscore_chi_squared_stats.txt'

		# Ignore windows where files don't exist (will need to be changed)
		if os.path.isfile(ldscore_file) == False or os.path.isfile(chi_sq_file) == False:
			continue

		# Extract info from input files
		anno_raw = np.loadtxt(ldscore_file,dtype=str,delimiter='\t')
		chi_sq_raw = np.loadtxt(chi_sq_file,dtype=str,delimiter='\t')	
		samp_size = float(chi_sq_raw[1,2])
		if anno_raw.shape[0] == 0:
			continue
		# Filter window based on chi-squared stats
		chi_sq_stats = chi_sq_raw[1:,1].astype(float)
		if np.max(chi_sq_stats) < window_chi_sq_lb:
			print('eroror')
			continue
		if np.max(chi_sq_stats) > window_chi_sq_ub:
			print('skip')
			continue

		num_windows = num_windows+1
		# Add window info to array
		X.append(samp_size*anno_raw[1:,2:].astype(float))
		y.append(chi_sq_stats)
		weights.append(anno_raw[1:,1].astype(float))
	print(num_windows)

	# Save to nicely formatted numpy matrices
	X = np.vstack(X)
	y=np.hstack(y)
	weights = np.hstack(weights)


	return X,y,weights

def ldsc_tf_loss_fxn(chi_sq, ldsc_model,X, snp_weights, intercept_variable):
	predy = tf.linalg.matmul(X, tf.math.softplus(ldsc_model))
	pred_chi_sq = (predy) + (tf.math.softplus(intercept_variable))


	log_like = (-.5)*tf.math.log(chi_sq) - tf.math.divide(chi_sq, 2.0*pred_chi_sq) - (.5*tf.math.log(2.0*pred_chi_sq))

	weighted_middle_indices_log_like = tf.math.divide(log_like, snp_weights)

	return -tf.math.reduce_sum(weighted_middle_indices_log_like) #+ tf.math.reduce_sum(tf.math.softplus(ldsc_model))

def ldsc_np_gradient_fxn(x, chi_sq, ld_scores, snp_weights):
	pred_chi_sq = np.dot(ld_scores, softplus_np(x))

	gradient_term_a = (np.exp(x)/(np.exp(x)+1.0))*np.dot((0.5*chi_sq/(np.square(pred_chi_sq)*snp_weights)), ld_scores)
	gradient_term_b = -(np.exp(x)/(np.exp(x)+1.0))*np.dot((0.5)/(pred_chi_sq*snp_weights), ld_scores)
	gradient = gradient_term_a + gradient_term_b
	return -gradient


def ldsc_np_loss_fxn(x, chi_sq,ld_scores, snp_weights):
	pred_chi_sq = np.dot(ld_scores, softplus_np(x))

	log_like = (-.5)*np.log(chi_sq) - np.divide(chi_sq, 2.0*pred_chi_sq) - (.5*np.log(2.0*pred_chi_sq))

	weighted_middle_indices_log_like = np.divide(log_like, snp_weights)

	return -np.sum(weighted_middle_indices_log_like)

def ldsc_np_loss_fxn_regularized(x, chi_sq, ld_scores, snp_weights, regularization_weight):
	pred_chi_sq = np.dot(ld_scores, softplus_np(x))


	log_like = (-.5)*np.log(chi_sq) - np.divide(chi_sq, 2.0*pred_chi_sq) - (.5*np.log(2.0*pred_chi_sq))	
	weighted_log_like = np.divide(log_like, snp_weights)
	loss = -np.sum(weighted_log_like) + np.sum(softplus_np(x[1:]))*regularization_weight
	return loss


def ldsc_np_loss_and_gradient_fxn_regularized(x, chi_sq, ld_scores, snp_weights, regularization_weight):
	pred_chi_sq = np.dot(ld_scores, softplus_np(x))

	gradient_term_a = (np.exp(x)/(np.exp(x)+1.0))*np.dot((0.5*chi_sq/(np.square(pred_chi_sq)*snp_weights)), ld_scores)

	gradient_term_b = -(np.exp(x)/(np.exp(x)+1.0))*np.dot((0.5)/(pred_chi_sq*snp_weights), ld_scores)

	gradient = -gradient_term_a - gradient_term_b 
	gradient[1:] = gradient[1:] + regularization_weight*np.exp(x[1:])/(1.0 + np.exp(x[1:]))

	log_like = (-.5)*np.log(chi_sq) - np.divide(chi_sq, 2.0*pred_chi_sq) - (.5*np.log(2.0*pred_chi_sq))	
	weighted_log_like = np.divide(log_like, snp_weights)
	loss = -np.sum(weighted_log_like) + np.sum(softplus_np(x[1:]))*regularization_weight
	return loss, gradient

def ldsc_np_loss_and_gradient_fxn(x, chi_sq, ld_scores, snp_weights):
	pred_chi_sq = np.dot(ld_scores, softplus_np(x))

	gradient_term_a = (np.exp(x)/(np.exp(x)+1.0))*np.dot((0.5*chi_sq/(np.square(pred_chi_sq)*snp_weights)), ld_scores)

	gradient_term_b = -(np.exp(x)/(np.exp(x)+1.0))*np.dot((0.5)/(pred_chi_sq*snp_weights), ld_scores)

	gradient = -gradient_term_a - gradient_term_b

	log_like = (-.5)*np.log(chi_sq) - np.divide(chi_sq, 2.0*pred_chi_sq) - (.5*np.log(2.0*pred_chi_sq))	
	weighted_log_like = np.divide(log_like, snp_weights)
	loss = -np.sum(weighted_log_like)
	return loss, gradient

def ldsc_np_loss_and_gradient_fxn_exp_link(x, chi_sq, ld_scores, snp_weights):
	pred_chi_sq = np.dot(ld_scores, np.exp(x))
	gradient_term_a = np.exp(x)*np.dot((0.5*chi_sq/(np.square(pred_chi_sq)*snp_weights)), ld_scores)
	gradient_term_b = -np.exp(x)*np.dot((0.5)/(pred_chi_sq*snp_weights), ld_scores)
	gradient = -gradient_term_a - gradient_term_b

	log_like = (-.5)*np.log(chi_sq) - np.divide(chi_sq, 2.0*pred_chi_sq) - (.5*np.log(2.0*pred_chi_sq))	
	weighted_log_like = np.divide(log_like, snp_weights)
	loss = -np.sum(weighted_log_like)
	return loss, gradient

def ldsc_np_loss_no_transform(x, chi_sq, ld_scores, snp_weights):
	pred_chi_sq = np.dot(ld_scores, np.exp(x))

	log_like = (-.5)*np.log(chi_sq) - np.divide(chi_sq, 2.0*pred_chi_sq) - (.5*np.log(2.0*pred_chi_sq))	
	weighted_log_like = np.divide(log_like, snp_weights)
	loss = -np.sum(weighted_log_like)
	return loss

def ldsc_np_loss_and_gradient_fxn_no_intercept(x, chi_sq, ld_scores, snp_weights):
	pred_chi_sq = np.dot(ld_scores, softplus_np(x)) + 1.0

	gradient_term_a = (np.exp(x)/(np.exp(x)+1.0))*np.dot((0.5*chi_sq/(np.square(pred_chi_sq)*snp_weights)), ld_scores)
	gradient_term_b = -(np.exp(x)/(np.exp(x)+1.0))*np.dot((0.5)/(pred_chi_sq*snp_weights), ld_scores)
	gradient = -gradient_term_a - gradient_term_b

	log_like = (-.5)*np.log(chi_sq) - np.divide(chi_sq, 2.0*pred_chi_sq) - (.5*np.log(2.0*pred_chi_sq))	
	weighted_log_like = np.divide(log_like, snp_weights)
	loss = -np.sum(weighted_log_like)
	return loss, gradient


def ldsc_np_loss_fxn_old(x, chi_sq,ld_scores, snp_weights):
	predy = np.dot(ld_scores, softplus_np(x[1:]))
	pred_chi_sq = (predy) + (softplus_np(x[0]))


	log_like = (-.5)*np.log(chi_sq) - np.divide(chi_sq, 2.0*pred_chi_sq) - (.5*np.log(2.0*pred_chi_sq))

	weighted_middle_indices_log_like = np.divide(log_like, snp_weights)

	return -np.sum(weighted_middle_indices_log_like) #+ tf.math.reduce_sum(tf.math.softplus(ldsc_model))

def ldsc_np_loss_fxn_no_intercept(x, chi_sq,ld_scores, snp_weights):
	predy = np.dot(ld_scores, softplus_np(x))
	pred_chi_sq = (predy) + 1.0


	log_like = (-.5)*np.log(chi_sq) - np.divide(chi_sq, 2.0*pred_chi_sq) - (.5*np.log(2.0*pred_chi_sq))

	weighted_middle_indices_log_like = np.divide(log_like, snp_weights)

	return -np.sum(weighted_middle_indices_log_like) #+ tf.math.reduce_sum(tf.math.softplus(ldsc_model))



#def softplus_np(x): return np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0)

def softplus_np(x): return np.log(np.exp(x) + 1)

def get_tissue_names(tissue_name_file):
	aa = np.loadtxt(tissue_name_file,dtype=str,delimiter='\t')
	return aa[1:,0]

def compute_genome_wide_heritability_estimates(X, chi_sq, snp_weights, learn_intercept='learn_intercept', max_epochs=20000):
	chi_sq = tf.convert_to_tensor(chi_sq.reshape(len(chi_sq),1), dtype=tf.float32)
	snp_weights = snp_weights.reshape(len(snp_weights),1)
	X = tf.convert_to_tensor(X, dtype=tf.float32)

	# Number of annotations
	annotation_data_dimension = X.shape[1]

	# Initialize mapping from annotations to per snp heritability
	ldsc_model = tf.Variable(np.ones((annotation_data_dimension,1))*-16, trainable=True,dtype=tf.float32)
	optimizer = tf.keras.optimizers.Adam()

	# Whether or not to learn intercept in LDSC
	# Initial value is np.log(np.exp(1)-1.0) [which equals 1 when put through softplus activation function]
	if learn_intercept == 'learn_intercept':
		log_intercept_variable = tf.Variable(initial_value=0.541324854612918,trainable=True, name='intercept')
	elif learn_intercept == 'fixed_intercept':
		log_intercept_variable = tf.Variable(initial_value=0.541324854612918,trainable=False, name='intercept')
	else:
		print('assumption error: intercept model called ' + learn_intercept + ' not currently implemented')
		pdb.set_trace()



	# Lopp through windows
	for epoch_iter in range(max_epochs):

		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = ldsc_tf_loss_fxn(chi_sq, ldsc_model, X, snp_weights, log_intercept_variable)

		# Define trainable variables
		trainable_variables = [ldsc_model]
		if learn_intercept == 'learn_intercept':
			trainable_variables.append(log_intercept_variable)
	
		# Compute and apply gradients
		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))

		if np.mod(epoch_iter, 100) == 0:
			print('###################################')
			print('epoch iter ' + str(epoch_iter))
			print('###################################')

			print(tf.math.softplus(ldsc_model)[:,0])

	h2 = np.asarray(tf.math.softplus(ldsc_model)[:,0])
	gene_h2 = h2[1:]
	return h2, gene_h2


def compute_genome_wide_heritability_estimates_lbfgs_exp_link(ld_scores, chi_sq, snp_weights, learn_intercept='learn_intercept', factr=10000000.0):
	if learn_intercept == 'learn_intercept':
		x0 = np.ones(ld_scores.shape[1] +1)*-15
		x0[0] = 0.0
		ld_scores_plus_intercept = np.hstack((np.ones((ld_scores.shape[0],1)), ld_scores))


		# scipy.optimize.approx_fprime(x0, ldsc_np_loss_no_transform, x0,chi_sq,ld_scores_plus_intercept,snp_weights)
		# ldsc_np_loss_and_gradient_fxn_no_transform(x0, chi_sq, ld_scores_plus_intercept,snp_weights)
		opti=scipy.optimize.fmin_l_bfgs_b(ldsc_np_loss_and_gradient_fxn_exp_link,x0, args=(chi_sq, ld_scores_plus_intercept, snp_weights), approx_grad=False)
		#opti=scipy.optimize.fmin_l_bfgs_b(ldsc_np_loss_no_transform,x0, args=(chi_sq, ld_scores_plus_intercept, snp_weights), approx_grad=True,  bounds=bounds)

		print('Warnflag: ' + str(opti[2]['warnflag']))
		opt_val = opti[0]
		t_opt_val = softplus_np(opt_val)
		opt_intercept = t_opt_val[0]
		opt_h2 = t_opt_val[1:]
		opt_h2_gene = opt_h2[1:]
	elif learn_intercept == 'fixed_intercept':
		pdb.set_trace()
		x0 = np.ones(X.shape[1])*-16
		opti=scipy.optimize.fmin_l_bfgs_b(ldsc_np_loss_and_gradient_fxn_no_intercept,x0, args=(chi_sq, ld_scores, snp_weights),approx_grad=False, factr=factr)
		print('Warnflag: ' + str(opti[2]['warnflag']))
		opt_val = opti[0]
		t_opt_val = softplus_np(opt_val)
		opt_intercept = 1.0
		opt_h2 = np.copy(t_opt_val)
		opt_h2_gene = opt_h2[1:]
	return opt_intercept, opt_h2, opt_h2_gene


def compute_genome_wide_heritability_estimates_lbfgs(ld_scores, chi_sq, snp_weights, learn_intercept='learn_intercept', factr=10000000.0):
	if learn_intercept == 'learn_intercept':
		x0 = np.ones(ld_scores.shape[1] +1)*-16
		x0[0] = .541324854612918
		ld_scores_plus_intercept = np.hstack((np.ones((ld_scores.shape[0],1)), ld_scores))
		#loss, grad = ldsc_np_loss_and_gradient_fxn(x0, chi_sq, ld_scores_plus_intercept, snp_weights)
		#loss, grad = ldsc_np_loss_and_gradient_fxn_regularized(x0, chi_sq, ld_scores_plus_intercept, snp_weights, 100000.0)
		#grad2 = scipy.optimize.approx_fprime(x0, ldsc_np_loss_fxn_regularized,1.4901161193847656e-08 ,chi_sq,ld_scores_plus_intercept,snp_weights, 100000.0)

		opti=scipy.optimize.fmin_l_bfgs_b(ldsc_np_loss_and_gradient_fxn, x0, args=(chi_sq, ld_scores_plus_intercept, snp_weights), approx_grad=False, factr=factr)
		#optir=scipy.optimize.fmin_l_bfgs_b(ldsc_np_loss_and_gradient_fxn_regularized, x0, args=(chi_sq, ld_scores_plus_intercept, snp_weights,10000000000.0), approx_grad=False, factr=factr)

		print('Warnflag: ' + str(opti[2]['warnflag']))
		opt_val = opti[0]
		t_opt_val = softplus_np(opt_val)
		opt_intercept = t_opt_val[0]
		opt_h2 = t_opt_val[1:]
		opt_h2_gene = opt_h2[1:]
	elif learn_intercept == 'fixed_intercept':
		x0 = np.ones(X.shape[1])*-16
		opti=scipy.optimize.fmin_l_bfgs_b(ldsc_np_loss_and_gradient_fxn_no_intercept,x0, args=(chi_sq, ld_scores, snp_weights),approx_grad=False, factr=factr)
		print('Warnflag: ' + str(opti[2]['warnflag']))
		opt_val = opti[0]
		t_opt_val = softplus_np(opt_val)
		opt_intercept = 1.0
		opt_h2 = np.copy(t_opt_val)
		opt_h2_gene = opt_h2[1:]
	return opt_intercept, opt_h2, opt_h2_gene

def identify_outlier_window_indices_based_on_genomic_jacknife_md(jacknifed_file, tissue_names, num_jacknife_windows, jacknife_windows):
	# First fill in h2 mat
	h2_mat = np.zeros((num_jacknife_windows, len(tissue_names)))
	mapping = {}
	for indexer, tissue_name in enumerate(tissue_names):
		mapping[tissue_name] = indexer
	f = open(jacknifed_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0] not in mapping:
			continue
		tiss_index = mapping[data[0]]
		jacknife_num = int(data[1])
		h2 = float(data[2])
		h2_mat[jacknife_num, tiss_index] = h2
	f.close()

	# Now compute mahalanobis distance on h2_mat
	V = np.cov(np.transpose(h2_mat))
	VI = np.linalg.inv(V)
	mu = np.mean(h2_mat,axis=0)
	mds = []
	for jacknife_iter in range(num_jacknife_windows):
		xx = h2_mat[jacknife_iter,:]
		md = np.sqrt(np.dot(np.dot((xx-mu),VI),(xx-mu).T))
		mds.append(md)
	mds = np.asarray(mds)
	md_pvalues = 1.0 - scipy.stats.chi2.cdf(np.square(mds),len(tissue_names) -1)

	print(sum(md_pvalues<.001))

	outlier_indices = []
	for jacknife_window in range(num_jacknife_windows):
		if md_pvalues[jacknife_window] < .001:
			outlier_indices.append(jacknife_windows[jacknife_window])
	outlier_indices = np.hstack(outlier_indices)

	return outlier_indices


def extract_top_n_windows(preprocessed_tgfm_data_dir, trait_name, window_names, num_top):
	# Keep track of maximum z score per window
	window_max_z_score_arr = []
	# Loop through windows
	for window_name in window_names:

		chi_sq_file = preprocessed_tgfm_data_dir + gene_type + '_' + window_name + '_' + trait_name + '_tgfm_ldscore_chi_squared_stats.txt'



		# Extract info from input files
		chi_sq_raw = np.loadtxt(chi_sq_file,dtype=str,delimiter='\t')

		# Get window z scores
		window_z_scores = np.sqrt(chi_sq_raw[1:,1].astype(float))

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

def compute_jacknifed_covariance_matrix(jacknifed_taus):
	jacknife_mean = np.mean(jacknifed_taus,axis=0)
	diff = jacknifed_taus - jacknife_mean
	num_jacknife_samples = jacknifed_taus.shape[0]

	jacknifed_cov = np.dot(np.transpose(diff),diff)*(num_jacknife_samples-1.0)/num_jacknife_samples

	return jacknifed_cov, jacknife_mean


def load_in_existing_multivariate_non_neg_tau_sldsc_results(file_name):
	# First get number of jacknifed windows
	global_arr = []
	mini_arr = []
	cur_window = 0
	head_count = 0
	f = open(file_name)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if cur_window != int(data[1]):
			cur_window = cur_window + 1
			global_arr.append(np.asarray(mini_arr))
			mini_arr = []
		mini_arr.append(float(data[2]))
	f.close()

	head_count = 0
	anno_names = []
	f = open(file_name)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[1] == '0':
			anno_names.append(data[0])
	f.close()

	anno_names = np.asarray(anno_names)

	global_arr.append(np.asarray(mini_arr))

	#sdev_data = np.loadtxt(anno_sdev_output_file,dtype=str, delimiter='\t')
	#sdevs = sdev_data[1:,1].astype(float)

	# Get jacknifed taus
	jacknifed_taus = np.asarray(global_arr)

	# Get jacknifed covariance matrix
	jacknifed_covariance, jacknifed_mean = compute_jacknifed_covariance_matrix(jacknifed_taus)

	jacknifed_mean_se = np.sqrt(np.diag(jacknifed_covariance))


	return jacknifed_mean, jacknifed_mean_se, jacknifed_covariance, anno_names

def remove_testis_from_tissue_names(tissue_names):
	new_tissue_names = []
	for t_index, tissue_name in enumerate(tissue_names):
		if tissue_name == 'Testis':
			testis_index = t_index
			continue
		new_tissue_names.append(tissue_name)
	return testis_index, new_tissue_names

trait_name = sys.argv[1]
ukkbb_window_summary_file = sys.argv[2]
tissue_name_file = sys.argv[3]
preprocessed_tgfm_data_dir = sys.argv[4]
learn_intercept = sys.argv[5]
output_stem = sys.argv[6]
gene_type = sys.argv[7]
top_windows = sys.argv[8]
gene_model_suffix = sys.argv[9]
remove_testis = sys.argv[10]

num_jacknife_windows = 200

np.random.seed(1)

# Get names of tissues
tissue_names = get_tissue_names(tissue_name_file)

if gene_model_suffix == '_':
	gene_model_suffix = ''

if remove_testis == 'True':
	testis_index, tissue_names = remove_testis_from_tissue_names(tissue_names)


# Get array of names of windows
window_names = get_window_names(ukkbb_window_summary_file, preprocessed_tgfm_data_dir)
if top_windows == 'True':
	window_names = extract_top_n_windows(preprocessed_tgfm_data_dir, trait_name, window_names, 500)


# Load in LDSC-style data
X,chi_sq,snp_weights = load_in_ldsc_style_data(preprocessed_tgfm_data_dir, trait_name, window_names, gene_model_suffix, window_chi_sq_lb=0.0, window_chi_sq_ub=300.0)

if remove_testis == 'True':
	X = np.delete(X, (testis_index+1), 1)

# Run LDSC-style regression (on full data)
intercept, h2, gene_h2 = compute_genome_wide_heritability_estimates_lbfgs(X, chi_sq, snp_weights, learn_intercept=learn_intercept)


# Print to output
t = open(output_stem + 'mean_estimates.txt','w')
t.write('Class_name\th2\n')
t.write('Intercept\t' + str(intercept) + '\n')
t.write('Genotype\t' + str(h2[0]) + '\n')
for tiss_index, tissue_name in enumerate(tissue_names):
	t.write(tissue_name + '\t' + str(gene_h2[tiss_index]) + '\n')
t.close()

# Jacknife ldsc-style regression

# Get indices corresponding to 200 jacknife windows
jacknife_windows = np.array_split(np.arange(X.shape[0]), num_jacknife_windows)

# Open output file
t = open(output_stem + 'jacknifed_mean_estimates.txt','w')
t.write('Class_name\tjacknife_window\th2\n')


# Loop through jacknife windows
for jacknife_window_iter in range(num_jacknife_windows):
	print(X.shape)
	# Remove points from jacknife window from data
	X_jack = np.delete(X, jacknife_windows[jacknife_window_iter], axis=0)
	chi_sq_jack = np.delete(chi_sq, jacknife_windows[jacknife_window_iter])
	snp_weights_jack = np.delete(snp_weights, jacknife_windows[jacknife_window_iter])

	print(X_jack.shape)


	# Run LDSC-style regression (on jacknifed data)
	jack_intercept, jack_h2, jack_gene_h2 = compute_genome_wide_heritability_estimates_lbfgs(X_jack, chi_sq_jack, snp_weights_jack, learn_intercept=learn_intercept)

	t.write('Intercept\t' + str(jacknife_window_iter) + '\t' + str(jack_intercept) + '\n')
	t.write('Genotype\t' + str(jacknife_window_iter) + '\t' + str(jack_h2[0]) + '\n')
	for tiss_index, tissue_name in enumerate(tissue_names):
		t.write(tissue_name + '\t' + str(jacknife_window_iter) + '\t' + str(jack_gene_h2[tiss_index]) + '\n')
	t.flush()

t.close()


# Print organized output file
jacknifed_file = output_stem + 'jacknifed_mean_estimates.txt'
jacknifed_mean, jacknifed_mean_se, jacknifed_covariance, anno_names = load_in_existing_multivariate_non_neg_tau_sldsc_results(jacknifed_file)
tau_z = jacknifed_mean/jacknifed_mean_se
# Open output handle and print to output
organized_output_file = output_stem + 'results_organized.txt'
t = open(organized_output_file,'w')
t.write('Annotation_name\ttau\ttau_se\ttau_z_score\n')
for anno_iter, anno_name in enumerate(anno_names):
	t.write(anno_name + '\t' + str(jacknifed_mean[anno_iter]) + '\t' + str(jacknifed_mean_se[anno_iter]) + '\t' + str(tau_z[anno_iter]) + '\n')
t.close()

# Print tau covariance 
jacknifed_cov_file = output_stem + 'jacknifed_covariance.txt'
np.savetxt(jacknifed_cov_file, jacknifed_covariance, fmt="%s", delimiter='\t')
















'''
# Re-run regression with jacknife windows removed
jacknifed_file = output_stem + 'jacknifed_mean_estimates.txt'
outlier_window_indices = identify_outlier_window_indices_based_on_genomic_jacknife_md(jacknifed_file, tissue_names, num_jacknife_windows, jacknife_windows)

X_non_outlier = np.delete(X, outlier_window_indices, axis=0)
chi_sq_non_outlier = np.delete(chi_sq, outlier_window_indices)
snp_weights_non_outlier = np.delete(snp_weights, outlier_window_indices)


non_outlier_intercept, non_outlier_h2, non_outlier_gene_h2 = compute_genome_wide_heritability_estimates_lbfgs(X_non_outlier, chi_sq_non_outlier, snp_weights_non_outlier, learn_intercept=learn_intercept, factr=100000.0)

# Print to output
t = open(output_stem + 'outlier_removed_mean_estimates.txt','w')
t.write('Class_name\th2\n')
t.write('Intercept\t' + str(non_outlier_intercept) + '\n')
t.write('Genotype\t' + str(non_outlier_h2[0]) + '\n')
for tiss_index, tissue_name in enumerate(tissue_names):
	t.write(tissue_name + '\t' + str(non_outlier_gene_h2[tiss_index]) + '\n')
t.close()
'''
