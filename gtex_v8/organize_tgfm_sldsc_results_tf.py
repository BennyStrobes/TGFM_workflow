import sys
sys.path.remove('/n/app/python/3.6.0/lib/python3.6/site-packages')
import os
import pdb
import numpy as np
import pickle
import time
import scipy.sparse
import scipy.optimize
import tensorflow as tf
import tensorflow_recommenders as tfrs




def get_anno_names(example_anno_file):
	f = open(example_anno_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		anno_names = np.asarray(data[3:])
		break
	f.close()
	return anno_names




def extract_tau_and_tau_se_from_log_file(sldsc_log_file, anno_names):
	f = open(sldsc_log_file)
	for line in f:
		line = line.rstrip()
		if line.startswith('Coefficients:'):
			coef = np.asarray(line.split()[1:]).astype(float)
			if len(coef) < len(anno_names):
				#line = f.next()
				line = next(f)
				data = np.asarray(line.split()).astype(float)
				coef = np.hstack([coef, data])
		if line.startswith('Coefficient SE:'):
			coef_se = np.asarray(line.split()[2:]).astype(float)
			if len(coef_se) < len(anno_names):
				#line = f.next()
				line = next(f)
				data = np.asarray(line.split()).astype(float)
				coef_se = np.hstack([coef_se, data])
	return coef, coef_se

def compute_jacknifed_covariance_matrix(jacknifed_taus):
	jacknife_mean = np.mean(jacknifed_taus,axis=0)
	diff = jacknifed_taus - jacknife_mean
	num_jacknife_samples = jacknifed_taus.shape[0]

	jacknifed_cov = np.dot(np.transpose(diff),diff)*(num_jacknife_samples-1.0)/num_jacknife_samples

	return jacknifed_cov, jacknife_mean

def extract_non_eqtl_annotations(anno_names):
	non_eqtl_indices = []
	eqtl_indices = []
	for ii, anno_name in enumerate(anno_names):
		if anno_name.endswith('L2'):
			non_eqtl_indices.append(ii)
		else:
			eqtl_indices.append(ii)
	return np.asarray(non_eqtl_indices), np.asarray(eqtl_indices)

def load_in_annotation_sdev_data(annotation_sd_file):
	aa = np.loadtxt(annotation_sd_file,dtype=str,delimiter='\t')
	return aa[1:,0], aa[1:,1].astype(float)

def load_in_mvec(anno_stem, suffix):
	for chrom_num in range(1,23):
		filer = anno_stem + '.' + str(chrom_num) + '.l2' + suffix
		tmp = np.loadtxt(filer)
		if chrom_num == 1:
			counter = np.zeros(len(tmp))
		counter = counter + tmp
	return counter

def print_organized_summary_file(output_file_name, anno_names, tau, tau_se):
	tau_z = tau/tau_se
	t = open(output_file_name,'w')
	t.write('Annotation\ttau\ttau_se\ttau_z\n')
	for ii in range(len(anno_names)):
		t.write(anno_names[ii] + '\t' + str(tau[ii]) + '\t' + str(tau_se[ii]) + '\t' + str(tau_z[ii]) + '\n')
	t.close()

def extract_expression_mediated_h2(jacknifed_taus, eqtl_start_index, m_vec):
	jacknifed_med_h2 = jacknifed_taus*m_vec

	jacknifed_geno_h2 = np.sum(jacknifed_med_h2[:,:eqtl_start_index],axis=1)
	jacknifed_expr_h2 = np.sum(jacknifed_med_h2[:,eqtl_start_index:],axis=1)

	jacknifed_expr_med_h2 = jacknifed_expr_h2/(jacknifed_expr_h2 + jacknifed_geno_h2)

	mean_expr_med_h2 = np.mean(jacknifed_expr_med_h2)

	diff = jacknifed_expr_med_h2 - mean_expr_med_h2
	num_jacknife_samples = jacknifed_expr_med_h2.shape[0]

	jacknifed_se = np.sqrt(np.dot(np.transpose(diff),diff)*(num_jacknife_samples-1.0)/num_jacknife_samples)
	
	return mean_expr_med_h2, jacknifed_se

def print_organized_h2_mediated(output_file_name, h2_med, h2_med_se):
	t = open(output_file_name,'w')
	t.write('h2_med\th2_med_se\n')
	t.write(str(h2_med) + '\t' + str(h2_med_se) + '\n')
	t.close()

def initialize_linear_model(n_anno):
	# Initialize Neural network model
	model = tf.keras.models.Sequential()
	model.add(tf.keras.layers.Dense(units=1, activation='linear', input_dim=n_anno, dtype=tf.float64))
	return model

def multivariate_normal_loss(jacknifed_tau_mean, jacknifed_tau_precision, pred_tau, reg_weight):
	diff = pred_tau - jacknifed_tau_mean
	log_like = -.5*tf.tensordot(tf.tensordot(diff, jacknifed_tau_precision,axes=1), diff, axes=1)
	loss_value = -log_like
	loss_value = loss_value + tf.math.reduce_sum(tf.math.square(pred_tau)*reg_weight)
	return loss_value

def sparse_regression_tf(jacknifed_tau_mean, jacknifed_tau_covariance, n_tissue_anno, max_epochs=1000, learning_rate=0.001, reg_weight=1e13):
	n_anno = len(jacknifed_tau_mean)
	# Initialize mapping from annotations to per snp heritability
	#sparse_linear_model = initialize_linear_model(n_genomic_anno)
	pred_tau = tf.Variable(initial_value=np.zeros(n_anno),trainable=True, name='coef', dtype=tf.float64)

	# Initialize optimizer
	optimizer = tf.keras.optimizers.Adam(learning_rate=learning_rate)


	jacknifed_tau_mean = tf.convert_to_tensor(jacknifed_tau_mean, dtype=tf.float64)
	jacknifed_tau_covariance = tf.convert_to_tensor(jacknifed_tau_covariance, dtype=tf.float64)
	jacknifed_tau_precision = tf.linalg.inv(jacknifed_tau_covariance)
	
	# Lopp through windows
	for epoch_iter in range(max_epochs):

		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			#pred_tau = sparse_linear_model(tf.convert_to_tensor(np.ones((1,n_genomic_anno)),dtype=tf.float64), training=True)
			loss_value = multivariate_normal_loss(jacknifed_tau_mean, jacknifed_tau_precision, pred_tau, reg_weight)
			# Define trainable variables
			trainable_variables = [pred_tau]
			# Compute and apply gradients
			grads = tape.gradient(loss_value, trainable_variables)
			optimizer.apply_gradients(zip(grads, trainable_variables))
	pdb.set_trace()


sldsc_output_root = sys.argv[1]
anno_stem = sys.argv[2]


sldsc_output_root + '_tf'
# Extract relevent info from log file
sldsc_log_file = sldsc_output_root + '.log'
example_anno_file = anno_stem + '.21.l2.ldscore'
anno_names = get_anno_names(example_anno_file)
tau, tau_se = extract_tau_and_tau_se_from_log_file(sldsc_log_file, anno_names)
tau_z = tau/tau_se
# Print to output
print_organized_summary_file(sldsc_output_root + 'organized_res.txt', anno_names, tau, tau_se)


# Total counts
m_vec = load_in_mvec(anno_stem, '.M')
m_5_50_vec = load_in_mvec(anno_stem, '.M_5_50')

# Print partioned h2 to output
print_organized_summary_file(sldsc_output_root + 'organized_mediated_h2.txt', anno_names, tau*m_vec, tau_se*m_vec)
print_organized_summary_file(sldsc_output_root + 'organized_mediated_h2_5_50.txt', anno_names, tau*m_5_50_vec, tau_se*m_5_50_vec)


# Load in annotation names and annotation sdevs
annotation_sd_file = anno_stem + '_annotation_sdev.txt'
annotation_names, annotation_sdev = load_in_annotation_sdev_data(annotation_sd_file)


# Load in jacknifed data
jacknifed_taus_file = sldsc_output_root + '.part_delete'
jacknifed_taus = np.loadtxt(jacknifed_taus_file)

# Compute standardized jacknifed mean tau and covariance
standardized_jacknifed_taus = jacknifed_taus*annotation_sdev
jacknifed_tau_covariance, jacknifed_tau_mean = compute_jacknifed_covariance_matrix(standardized_jacknifed_taus)
jacknifed_tau_z = jacknifed_tau_mean/np.sqrt(np.diag(jacknifed_tau_covariance))


# Extract annotations to be fixed (ie the non-eqtl annotations)
non_eqtl_annotations, eqtl_annotations = extract_non_eqtl_annotations(anno_names)
eqtl_start_index = np.min(eqtl_annotations)

# Jacknife expression mediated h2
h2_med, h2_med_se = extract_expression_mediated_h2(jacknifed_taus, eqtl_start_index, m_vec)
h2_5_50_med, h2_5_50_med_se = extract_expression_mediated_h2(jacknifed_taus, eqtl_start_index, m_5_50_vec)

print_organized_h2_mediated(sldsc_output_root + 'h2_med.txt', h2_med, h2_med_se)
print_organized_h2_mediated(sldsc_output_root + 'h2_5_50_med.txt', h2_5_50_med, h2_5_50_med_se)


sparse_regression_tf(jacknifed_tau_mean, jacknifed_tau_covariance, 93)

