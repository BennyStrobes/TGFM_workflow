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

def load_in_ldsc_style_data(preprocessed_tgfm_data_dir, trait_name, window_names, window_chi_sq_lb=-1, window_chi_sq_ub=100000000000000000.0):
	# Initialize output
	X = []
	y = []
	weights = []
	num_windows = 0

	# Loop through windows
	for window_name in window_names:
		# Input files for this window
		ldscore_file = preprocessed_tgfm_data_dir + window_name + '_tgfm_ldscore_annotation_file.txt'
		#ldscore_file = preprocessed_tgfm_data_dir + window_name + '_tgfm_ldscore_not_standardized_annotation_file.txt'

		chi_sq_file = preprocessed_tgfm_data_dir + window_name + '_' + trait_name + '_tgfm_ldscore_chi_squared_stats.txt'

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

def softplus_np(x): return np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0)

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



trait_name = sys.argv[1]
ukkbb_window_summary_file = sys.argv[2]
tissue_name_file = sys.argv[3]
preprocessed_tgfm_data_dir = sys.argv[4]
output_stem = sys.argv[5]


# Get names of tissues
tissue_names = get_tissue_names(tissue_name_file)

# Get array of names of windows
window_names = get_window_names(ukkbb_window_summary_file)

# Load in LDSC-style data
X,chi_sq,snp_weights = load_in_ldsc_style_data(preprocessed_tgfm_data_dir, trait_name, window_names, window_chi_sq_lb=0.0, window_chi_sq_ub=300.0)


# Run LDSC-style regression
h2, gene_h2 = compute_genome_wide_heritability_estimates(X, chi_sq, snp_weights, learn_intercept='learn_intercept', max_epochs=20000)


t = open(output_stem + 'mean_estimates.txt','w')
t.write('Class_name\th2\n')
t.write('Genotype\t' + str(h2[0]) + '\n')
for tiss_index, tissue_name in enumerate(tissue_names):
	t.write(tissue_name + '\t' + str(gene_h2[tiss_index]) + '\n')
t.close()

