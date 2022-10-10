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



def load_in_data(dir_name, trait_name,samp_size, standardize_eqtl):
	X = []
	y = []
	weights = []

	#trait_name = 'blood_WHITE_COUNT'
	#dir_name = '/n/groups/price/ben/causal_eqtl_gwas/gtex_v8_causal_eqtl_gwas_38_tissues/pseudotissue_gtex_rss_multivariate_twas/'
	for chrom_num in [1,2,3,4,5,7,8,9,10,11,12,13,14, 15,16,17,18,20,21,22]:
		if chrom_num==6:
			continue
		if chrom_num==16:
			continue
		if chrom_num==19:
			continue
		if standardize_eqtl == 'True':
			anno_file = dir_name + trait_name + '_cis_heritable_genes_' + str(chrom_num) + '_tgfm_ldsc_standardized_annotations.txt'
		elif standardize_eqtl == 'False':
			anno_file = dir_name + trait_name + '_cis_heritable_genes_' + str(chrom_num) + '_tgfm_ldsc_annotations.txt'
		else:
			print('assumption eroror: currently not implemented')
			pdb.set_trace()
		chi_sq_file = dir_name + trait_name + '_cis_heritable_genes_' + str(chrom_num) + '_tgfm_ldsc_chi_sq.txt'
		anno_raw = np.loadtxt(anno_file,dtype=str,delimiter='\t')
		chi_sq_raw = np.loadtxt(chi_sq_file,dtype=str,delimiter='\t')
		X.append(samp_size*anno_raw[1:,2:].astype(float))
		y.append(chi_sq_raw[1:,1].astype(float))
		#weights.append(1.0/anno_raw[1:,1].astype(float))
		weights.append(anno_raw[1:,1].astype(float))

	X = np.vstack(X)
	y=np.hstack(y)
	weights = np.hstack(weights)
	#reg = LinearRegression(positive=True).fit(X, y,sample_weight=weights)
	return X,y,weights


def init_linear_mapping_ldsc_model(annotation_data_dimension):
	# Initialize Neural network model
	model = tf.keras.models.Sequential()
	model.add(tf.keras.layers.Dense(units=1, activation='linear',kernel_constraint=tf.keras.constraints.NonNeg(), input_dim=annotation_data_dimension))
	return model


def ldsc_tf_loss_fxn(chi_sq, ldsc_model,X, snp_weights, intercept_variable):
	predy = tf.linalg.matmul(X, tf.math.softplus(ldsc_model))
	pred_chi_sq = (predy) + (tf.math.softplus(intercept_variable))

	log_like = (-.5)*tf.math.log(chi_sq) - tf.math.divide(chi_sq, 2.0*pred_chi_sq) - (.5*tf.math.log(2.0*pred_chi_sq))

	weighted_middle_indices_log_like = tf.math.divide(log_like, snp_weights)

	return -tf.math.reduce_sum(weighted_middle_indices_log_like)

def softplus_np(x): return np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0)


trait_name = sys.argv[1]
dir_name = sys.argv[2]
gene_version = sys.argv[3]
samp_size = float(sys.argv[4])


learn_intercept='fixed_intercept'
learn_intercept='learn_intercept'

standardize_eqtl = 'False'
standardize_eqtl = 'True'

max_epochs=10000

X,chi_sq,snp_weights = load_in_data(dir_name, trait_name, samp_size, standardize_eqtl)

chi_sq = tf.convert_to_tensor(chi_sq.reshape(len(chi_sq),1), dtype=tf.float32)
snp_weights = snp_weights.reshape(len(snp_weights),1)
X = tf.convert_to_tensor(X, dtype=tf.float32)


# Number of annotations
annotation_data_dimension = X.shape[1]

# Initialize mapping from annotations to per snp heritability
#ldsc_model = init_linear_mapping_ldsc_model(annotation_data_dimension)
ldsc_model = tf.Variable(np.ones((39,1))*-15, trainable=True,dtype=tf.float32)
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
	# Loop through windows
	print('###################################')
	print('epoch iter ' + str(epoch_iter))
	print('###################################')



	# Use tf.gradient tape to compute gradients
	with tf.GradientTape() as tape:
		#predy = ldsc_model(X, training=True)
		loss_value = ldsc_tf_loss_fxn(chi_sq, ldsc_model, X, snp_weights, log_intercept_variable)

	# Define trainable variables
	#trainable_variables = ldsc_model.trainable_weights
	#if learn_intercept == 'learn_intercept':
		#trainable_variables.append(log_intercept_variable)
	trainable_variables = [ldsc_model]
	if learn_intercept == 'learn_intercept':
		trainable_variables.append(log_intercept_variable)
	# Compute and apply gradients
	grads = tape.gradient(loss_value, trainable_variables)
	optimizer.apply_gradients(zip(grads, trainable_variables))

	print(tf.math.softplus(ldsc_model)[:,0])

pdb.set_trace()






