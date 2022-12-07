import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import pandas as pd
import os
import pdb
from sklearn.linear_model import LinearRegression
from sparse_sldsc_from_multivariate_summary_statistics import SPARSE_SLDSC
from sparse_sldsc_from_multivariate_summary_statistics_fixed_genotype import SPARSE_SLDSC_FIXED_TERM

'''
# Only needed if want to run default susie
import rpy2
import rpy2.robjects.numpy2ri as numpy2ri
import rpy2.robjects as ro
ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
susieR_pkg = importr('susieR')
from sparse_sldsc_from_multivariate_summary_statistics import SPARSE_SLDSC
'''

def get_tissue_names(tissue_name_file):
	aa = np.loadtxt(tissue_name_file,dtype=str,delimiter='\t')
	return aa[1:,0]

def get_window_names(ukkbb_window_summary_file, preprocessed_tgfm_data_dir, gene_type):
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



def load_in_ldsc_style_data(preprocessed_tgfm_data_dir, trait_name, window_names, window_chi_sq_lb=-1, window_chi_sq_ub=100000000000000000.0, snp_variance_weight=False):
	# Initialize output
	X = []
	y = []
	weights = []
	num_windows = 0

	# Loop through windows
	for window_name in window_names:
		# Input files for this window
		ldscore_file = preprocessed_tgfm_data_dir + gene_type + '_' + window_name + '_tgfm_ldscore_annotation_file.txt'
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
		if chi_sq_raw.shape[0] == 1:
			continue
		# Filter window based on chi-squared stats
		chi_sq_stats = chi_sq_raw[1:,1].astype(float)
		if anno_raw[1:,2:].shape[0] != len(chi_sq_stats):
			continue
		if np.max(chi_sq_stats) < window_chi_sq_lb:
			print('eroror')
			continue
		if np.max(chi_sq_stats) > window_chi_sq_ub:
			print('skip')
			print(np.max(chi_sq_stats))
			print(window_name)
			continue

		num_windows = num_windows+1
		# Add window info to array
		X.append(samp_size*anno_raw[1:,2:].astype(float))
		y.append(chi_sq_stats)
		weights.append(anno_raw[1:,1].astype(float))

	# Save to nicely formatted numpy matrices
	X = np.vstack(X)
	y=np.hstack(y)
	weights = np.hstack(weights)

	if snp_variance_weight == True:
		approx_tau_hat = (np.mean(y) - 1.0)/(samp_size*np.mean(np.sum(X/samp_size,axis=1)))
		variance_weight_term = np.square(np.sum(X/samp_size,axis=1)*samp_size*approx_tau_hat + 1.0)
		weights = weights*variance_weight_term 

	return X,y,weights






def get_anno_names(tissue_names):
	annos =[]
	annos.append('Genotype')
	for tissue_name in tissue_names:
		annos.append('Expression_' + tissue_name)
	return np.asarray(annos)


def sldsc_analysis(X, chi_sq, snp_weights, num_jacknife_windows, intercept=True, fixed_intercept=0.0):
	if intercept:
		# Run S-LDSC on full data
		reg = LinearRegression().fit(X, chi_sq, sample_weight=1.0/snp_weights)
		learned_taus = reg.coef_
		learned_intercept = reg.intercept_
	
		# Get indices corresponding to 200 jacknife windows
		jacknife_windows = np.array_split(np.arange(X.shape[0]), num_jacknife_windows)

		# Initialize array to keep track of jacknifed taus
		jacknifed_taus = []
		jacknifed_intercepts = []
		# Loop through jacknife windows
		for jacknife_window_iter in range(num_jacknife_windows):
			print(jacknife_window_iter)
			# Remove points from jacknife window from data
			X_jack = np.delete(X, jacknife_windows[jacknife_window_iter], axis=0)
			chi_sq_jack = np.delete(chi_sq, jacknife_windows[jacknife_window_iter])
			snp_weights_jack = np.delete(snp_weights, jacknife_windows[jacknife_window_iter])

			# Run S-LDSC on jacknifed data
			jack_reg = LinearRegression().fit(X_jack, chi_sq_jack, sample_weight=1.0/snp_weights_jack)
			
			# add sldsc result to array
			jacknifed_taus.append(jack_reg.coef_)
			jacknifed_intercepts.append(jack_reg.intercept_)

		#Put in compact matrix
		jacknifed_taus = np.asarray(jacknifed_taus)
		jacknifed_intercepts = np.asarray(jacknifed_intercepts)

		final_intercept = np.mean(jacknifed_intercepts)
	else:
		reg = LinearRegression(fit_intercept=False).fit(X, chi_sq-fixed_intercept, sample_weight=1.0/snp_weights)
		learned_taus = reg.coef_
		
		# Get indices corresponding to 200 jacknife windows
		jacknife_windows = np.array_split(np.arange(X.shape[0]), num_jacknife_windows)

		# Initialize array to keep track of jacknifed taus
		jacknifed_taus = []
		jacknifed_intercepts = []
		# Loop through jacknife windows
		for jacknife_window_iter in range(num_jacknife_windows):
			# Remove points from jacknife window from data
			X_jack = np.delete(X, jacknife_windows[jacknife_window_iter], axis=0)
			chi_sq_jack = np.delete(chi_sq, jacknife_windows[jacknife_window_iter])
			snp_weights_jack = np.delete(snp_weights, jacknife_windows[jacknife_window_iter])

			# Run S-LDSC on jacknifed data
			jack_reg = LinearRegression(fit_intercept=False).fit(X_jack, chi_sq_jack-fixed_intercept, sample_weight=1.0/snp_weights_jack)
			
			# add sldsc result to array
			jacknifed_taus.append(jack_reg.coef_)

		#Put in compact matrix
		jacknifed_taus = np.asarray(jacknifed_taus)
		final_intercept = fixed_intercept

	# Compute jacknifed standard errors
	jacknife_mean = np.mean(jacknifed_taus,axis=0)
	jacknife_var = np.sum(np.square(jacknifed_taus - jacknife_mean),axis=0)*(num_jacknife_windows-1.0)/num_jacknife_windows
	jacknife_se = np.sqrt(jacknife_var)

	return final_intercept, jacknife_mean, jacknife_se, jacknifed_taus

def compute_jacknifed_covariance_matrix(jacknifed_taus):
	jacknife_mean = np.mean(jacknifed_taus,axis=0)
	diff = jacknifed_taus - jacknife_mean
	num_jacknife_samples = jacknifed_taus.shape[0]

	jacknifed_cov = np.dot(np.transpose(diff),diff)*(num_jacknife_samples-1.0)/num_jacknife_samples

	return jacknifed_cov, jacknife_mean


def multivariate_sldsc_shell(X, chi_sq, snp_weights, anno_names, num_jacknife_windows, predictor_sdevs, multivariate_output_stem):
	# Run S-LDSC
	sldsc_intercept, sldsc_tau, sldsc_tau_se, jacknifed_taus = sldsc_analysis(X, chi_sq, snp_weights, num_jacknife_windows, intercept=True)

	# Print results to output
	organized_results_file = multivariate_output_stem + '_organized.txt'
	t = open(organized_results_file,'w')
	# Print header and intercept
	t.write('Annotation_name\ttau\ttau_se\ttau_z_score\n')
	t.write('Intercept\t' + str(sldsc_intercept) + '\tNA\tNA\n')
	# Print each annotation
	for anno_iter, anno_name in enumerate(anno_names):
		t.write(anno_name + '\t' + str(sldsc_tau[anno_iter]) + '\t' + str(sldsc_tau_se[anno_iter]) + '\t' + str(sldsc_tau[anno_iter]/sldsc_tau_se[anno_iter]) + '\n')
	t.close()

	# Save jacknifed taus to output
	jacknifed_taus_output_file = multivariate_output_stem + '_raw_jacknifed_taus.txt'
	np.savetxt(jacknifed_taus_output_file, jacknifed_taus, fmt='%s', delimiter='\t')

	# Get jacknifed covariance matrix
	jacknifed_covariance, jacknifed_mean = compute_jacknifed_covariance_matrix(jacknifed_taus)
	# Save jacknifed cov to output
	jacknifed_cov_output_file = multivariate_output_stem + '_raw_jacknifed_covariance.txt'
	np.savetxt(jacknifed_cov_output_file, jacknifed_covariance, fmt='%s', delimiter='\t')

	# Standardized jacknifed taus
	standardized_jacknifed_taus = jacknifed_taus*predictor_sdevs

	# Save jacknifed taus to output
	stand_jacknifed_taus_output_file = multivariate_output_stem + '_standardized_jacknifed_taus.txt'
	np.savetxt(stand_jacknifed_taus_output_file, standardized_jacknifed_taus, fmt='%s', delimiter='\t')

	# Get jacknifed covariance matrix
	stand_jacknifed_covariance, stand_jacknifed_mean = compute_jacknifed_covariance_matrix(standardized_jacknifed_taus)
	# Save jacknifed cov to output
	stand_jacknifed_cov_output_file = multivariate_output_stem + '_standardized_jacknifed_covariance.txt'
	np.savetxt(stand_jacknifed_cov_output_file, stand_jacknifed_covariance, fmt='%s', delimiter='\t')


	return

def multivariate_standardized_sldsc_shell(X, chi_sq, snp_weights, anno_names, num_jacknife_windows, multivariate_output_stem, predictor_sdevs):
	# Run S-LDSC
	sldsc_intercept, sldsc_tau, sldsc_tau_se, jacknifed_taus = sldsc_analysis(X, chi_sq, snp_weights, num_jacknife_windows, intercept=True)

	# Print results to output
	organized_results_file = multivariate_output_stem + '_organized.txt'
	t = open(organized_results_file,'w')
	# Print header and intercept
	t.write('Annotation_name\ttau\ttau_se\ttau_z_score\n')
	t.write('Intercept\t' + str(sldsc_intercept) + '\tNA\tNA\n')
	# Print each annotation
	for anno_iter, anno_name in enumerate(anno_names):
		t.write(anno_name + '\t' + str(sldsc_tau[anno_iter]) + '\t' + str(sldsc_tau_se[anno_iter]) + '\t' + str(sldsc_tau[anno_iter]/sldsc_tau_se[anno_iter]) + '\n')
	t.close()

	# Print results to output
	organized_unscaled_results_file = multivariate_output_stem + '_original_scale_organized.txt'
	t = open(organized_unscaled_results_file,'w')
	# Print header and intercept
	t.write('Annotation_name\ttau\ttau_se\ttau_z_score\n')
	t.write('Intercept\t' + str(sldsc_intercept) + '\tNA\tNA\n')
	# Print each annotation
	for anno_iter, anno_name in enumerate(anno_names):
		t.write(anno_name + '\t' + str(sldsc_tau[anno_iter]/predictor_sdevs[anno_iter]) + '\t' + str(sldsc_tau_se[anno_iter]/predictor_sdevs[anno_iter]) + '\t' + str(sldsc_tau[anno_iter]/sldsc_tau_se[anno_iter]) + '\n')
	t.close()

	# Save jacknifed taus to output
	jacknifed_taus_output_file = multivariate_output_stem + '_raw_jacknifed_taus.txt'
	np.savetxt(jacknifed_taus_output_file, jacknifed_taus, fmt='%s', delimiter='\t')

	# Get jacknifed covariance matrix
	jacknifed_covariance = compute_jacknifed_covariance_matrix(jacknifed_taus)
	# Save jacknifed cov to output
	jacknifed_cov_output_file = multivariate_output_stem + '_jacknifed_covariance.txt'
	np.savetxt(jacknifed_cov_output_file, jacknifed_covariance, fmt='%s', delimiter='\t')

	return sldsc_intercept, sldsc_tau, sldsc_tau_se


def univariate_learn_intercept_sldsc_shell(X, chi_sq, snp_weights, anno_names, num_jacknife_windows, univariate_output_stem):
	# Open output file
	organized_results_file = univariate_output_stem + '_organized.txt'
	t = open(organized_results_file, 'w')
	# Print header 
	t.write('Annotation_name\ttau\ttau_se\ttau_z_score\n')
	t.write('Intercept\t' + 'NA' + '\tNA\tNA\n')

	# Loop through annotations
	# Do sldsc analysis in each annotation independently
	for anno_iter, anno_name in enumerate(anno_names):
		print(anno_name)

		# Run S-LDSC on data from this annotation alone
		sldsc_intercept, sldsc_tau, sldsc_tau_se, jacknifed_taus = sldsc_analysis(X[:,anno_iter:(anno_iter+1)], chi_sq, snp_weights, num_jacknife_windows, intercept=True)

		# Print results to output
		t.write(anno_name + '\t' + str(sldsc_tau[0]) + '\t' + str(sldsc_tau_se[0]) + '\t' + str(sldsc_tau[0]/sldsc_tau_se[0]) + '\n')

		# Save jacknifed taus to output
		jacknifed_taus_output_file = univariate_output_stem + '_' + anno_name + '_raw_jacknifed_taus.txt'
		np.savetxt(jacknifed_taus_output_file, jacknifed_taus, fmt='%s', delimiter='\t')

	t.close()


def univariate_sldsc_shell(X, chi_sq, snp_weights, anno_names, num_jacknife_windows, univariate_output_stem, mv_sldsc_intercept):
	# Open output file
	organized_results_file = univariate_output_stem + '_organized.txt'
	t = open(organized_results_file, 'w')
	# Print header 
	t.write('Annotation_name\ttau\ttau_se\ttau_z_score\n')
	t.write('Intercept\t' + str(mv_sldsc_intercept) + '\tNA\tNA\n')

	# Loop through annotations
	# Do sldsc analysis in each annotation independently
	for anno_iter, anno_name in enumerate(anno_names):

		# Run S-LDSC on data from this annotation alone
		sldsc_intercept, sldsc_tau, sldsc_tau_se, jacknifed_taus = sldsc_analysis(X[:,anno_iter:(anno_iter+1)], chi_sq, snp_weights, num_jacknife_windows, intercept=False, fixed_intercept=mv_sldsc_intercept)

		# Print results to output
		t.write(anno_name + '\t' + str(sldsc_tau[0]) + '\t' + str(sldsc_tau_se[0]) + '\t' + str(sldsc_tau[0]/sldsc_tau_se[0]) + '\n')

		# Save jacknifed taus to output
		jacknifed_taus_output_file = univariate_output_stem + '_' + anno_name + '_raw_jacknifed_taus.txt'
		np.savetxt(jacknifed_taus_output_file, jacknifed_taus, fmt='%s', delimiter='\t')

	t.close()

def load_in_marginal_tau_sldsc_results(marginal_tau_results_file):
	tmp = np.loadtxt(marginal_tau_results_file, dtype=str,delimiter='\t')
	marginal_tau = tmp[2:,1].astype(float)
	marginal_tau_se = tmp[2:,2].astype(float)
	marginal_tau_z = tmp[2:,3].astype(float)
	return marginal_tau, marginal_tau_se, marginal_tau_z



def load_in_multivariate_tau_sldsc_results(file_name):
	aa = np.loadtxt(file_name, dtype=str, delimiter='\t')
	tau = aa[2:,1].astype(float)
	return tau



def standardize_predictors(X):
	# Compute sdev of each predictor
	predictor_sdev = np.std(X,axis=0)

	# Compute standardized predictor matrix
	X_stand = X/predictor_sdev

	return X_stand, predictor_sdev

def load_in_existing_multivariate_non_neg_tau_sldsc_results(trait_name, existing_dir, anno_sdev_output_file):
	file_name = existing_dir + 'tgfm_ldsc_style_heritability_' + trait_name + '_cis_heritable_gene_learn_intercept_top_window_False__remove_testis_True_jacknifed_mean_estimates.txt'
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
		if data[0] == 'Intercept':
			continue
		if cur_window != int(data[1]):
			cur_window = cur_window + 1
			global_arr.append(np.asarray(mini_arr))
			mini_arr = []
		mini_arr.append(float(data[2]))
	f.close()
	global_arr.append(np.asarray(mini_arr))

	sdev_data = np.loadtxt(anno_sdev_output_file,dtype=str, delimiter='\t')
	sdevs = sdev_data[1:,1].astype(float)

	# Get jacknifed taus
	jacknifed_taus = np.asarray(global_arr)
	standardized_jacknifed_taus = jacknifed_taus*sdevs

	# Get jacknifed covariance matrix
	jacknifed_covariance, jacknifed_mean = compute_jacknifed_covariance_matrix(standardized_jacknifed_taus)


	return jacknifed_mean, jacknifed_covariance

def print_annotation_names_and_sdevs(anno_names, predictor_sdevs, anno_sdev_output_file):
	t = open(anno_sdev_output_file,'w')
	t.write('annotation_name\tstandard_deviation\n')

	if len(predictor_sdevs) != len(anno_names):
		print('assumptino eroror')
		pdb.set_trace()
	for anno_iter in range(len(anno_names)):
		t.write(anno_names[anno_iter] + '\t' + str(predictor_sdevs[anno_iter]) + '\n')
	t.close()

def load_in_anno_names_and_sdev(anno_file):
	aa = np.loadtxt(anno_file,dtype=str, delimiter='\t')
	anno_names = aa[1:,0]
	anno_sdevs = aa[1:,1].astype(float)
	return anno_names, anno_sdevs

def print_sparse_sldsc_pmces(sparse_sldsc_obj, anno_names, anno_sdevs, fixed_geno, pmces_output_file):
	if fixed_geno == False:
		raw_beta_mu = sparse_sldsc_obj.beta_mu/anno_sdevs
		pmces = np.sum(raw_beta_mu*sparse_sldsc_obj.beta_phi,axis=0)
	else:
		raw_beta_mu = sparse_sldsc_obj.beta_mu/anno_sdevs[1:]
		pmces_non_geno = np.sum(raw_beta_mu*sparse_sldsc_obj.beta_phi,axis=0)
		pmces_geno = [sparse_sldsc_obj.genotype_beta_mu/anno_sdevs[0]]
		pmces = np.hstack((pmces_geno, pmces_non_geno))

	# Error check
	if len(pmces) != len(anno_names):
		print('assumption erorro')
		pdb.set_trace()

	t = open(pmces_output_file, 'w')
	t.write('annotation\tpmces_tau\n')
	for anno_iter in range(len(anno_names)):
		t.write(anno_names[anno_iter] + '\t' + str(pmces[anno_iter]) + '\n')
	t.close()

def print_sparse_sldsc_posterior_phi(sparse_sldsc_obj, anno_names, anno_sdevs, fixed_geno, posterior_phi_output_file):
	if fixed_geno == False:
		np.savetxt(posterior_phi_output_file, sparse_sldsc_obj.beta_phi, fmt="%s", delimiter='\t')
	else:
		orig_L, orig_num_anno = sparse_sldsc_obj.beta_phi.shape
		tmp_phi = np.vstack((np.zeros((1, orig_num_anno)), sparse_sldsc_obj.beta_phi))
		tmp_vec = np.zeros(((orig_L + 1), 1))
		tmp_vec[0,0] = 1
		final_phi = np.hstack((tmp_vec, tmp_phi))
		np.savetxt(posterior_phi_output_file, final_phi, fmt="%s", delimiter='\t')

def print_sparse_sldsc_posterior_beta_mu(sparse_sldsc_obj, anno_names, anno_sdevs, fixed_geno, posterior_beta_mu_output_file):
	if fixed_geno == False:
		raw_beta_mu = sparse_sldsc_obj.beta_mu/anno_sdevs
		np.savetxt(posterior_beta_mu_output_file, raw_beta_mu, fmt="%s", delimiter='\t')
	else:
		raw_beta_mu = sparse_sldsc_obj.beta_mu/anno_sdevs[1:]
		orig_L, orig_num_anno = sparse_sldsc_obj.beta_mu.shape
		tmp_mu = np.vstack((np.zeros((1, orig_num_anno)), raw_beta_mu))
		tmp_vec = np.zeros(((orig_L + 1), 1))
		tmp_vec[0,0] = sparse_sldsc_obj.genotype_beta_mu/anno_sdevs[0]
		final_mu = np.hstack((tmp_vec, tmp_mu))
		np.savetxt(posterior_beta_mu_output_file, final_mu, fmt="%s", delimiter='\t')

def print_sparse_sldsc_posterior_beta_var(sparse_sldsc_obj, anno_names, anno_sdevs, fixed_geno, posterior_beta_var_output_file):
	if fixed_geno == False:
		raw_beta_var = np.square(np.sqrt(sparse_sldsc_obj.beta_var)/anno_sdevs)
		np.savetxt(posterior_beta_var_output_file, raw_beta_var, fmt="%s", delimiter='\t')
	else:
		raw_beta_var = np.square(np.sqrt(sparse_sldsc_obj.beta_var)/anno_sdevs[1:])
		orig_L, orig_num_anno = sparse_sldsc_obj.beta_var.shape
		tmp_var = np.vstack((np.zeros((1, orig_num_anno)), raw_beta_var))
		tmp_vec = np.zeros(((orig_L + 1), 1))
		tmp_vec[0,0] = np.square(np.sqrt(sparse_sldsc_obj.genotype_beta_var)/anno_sdevs[0])
		final_var = np.hstack((tmp_vec, tmp_var))
		np.savetxt(posterior_beta_var_output_file, final_var, fmt="%s", delimiter='\t')

def print_sparse_sldsc_posterior_component_variance(sparse_sldsc_obj, anno_names, anno_sdevs, fixed_geno, posterior_component_var_output_file):
	if fixed_geno == False:
		np.savetxt(posterior_component_var_output_file, sparse_sldsc_obj.component_variances, fmt="%s", delimiter='\t')
	else:
		component_variances = np.hstack(([sparse_sldsc_obj.genotype_component_variance], sparse_sldsc_obj.component_variances))
		np.savetxt(posterior_component_var_output_file, component_variances, fmt="%s", delimiter='\t')

def save_sparse_sldsc_results(sparse_sldsc_obj, anno_file, sparse_sldsc_output_stem, fixed_geno=False):
	# Load in anno names and anno_sdev
	anno_names, anno_sdevs = load_in_anno_names_and_sdev(anno_file)

	# Print PMCES
	pmces_output_file = sparse_sldsc_output_stem + 'pmces.txt'
	print_sparse_sldsc_pmces(sparse_sldsc_obj, anno_names, anno_sdevs, fixed_geno, pmces_output_file)

	# Print phi
	posterior_phi_output_file = sparse_sldsc_output_stem + 'posterior_phi.txt'
	print_sparse_sldsc_posterior_phi(sparse_sldsc_obj, anno_names, anno_sdevs, fixed_geno, posterior_phi_output_file)

	# Print beta mu
	posterior_beta_mu_output_file = sparse_sldsc_output_stem + 'posterior_beta_mu.txt'
	print_sparse_sldsc_posterior_beta_mu(sparse_sldsc_obj, anno_names, anno_sdevs, fixed_geno, posterior_beta_mu_output_file)

	# Print beta var
	posterior_beta_var_output_file = sparse_sldsc_output_stem + 'posterior_beta_var.txt'
	print_sparse_sldsc_posterior_beta_var(sparse_sldsc_obj, anno_names, anno_sdevs, fixed_geno, posterior_beta_var_output_file)

	# Print standardized component variances
	posterior_component_var_output_file = sparse_sldsc_output_stem + 'posterior_component_variance.txt'
	print_sparse_sldsc_posterior_component_variance(sparse_sldsc_obj, anno_names, anno_sdevs, fixed_geno, posterior_component_var_output_file)

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
tgfm_heritability_results_dir = sys.argv[5]
learn_intercept = sys.argv[6]
output_stem = sys.argv[7]
gene_type = sys.argv[8]


remove_testis = 'True'

num_jacknife_windows = 200

np.random.seed(1)

# Get names of tissues
tissue_names = get_tissue_names(tissue_name_file)

if remove_testis == 'True':
	testis_index, tissue_names = remove_testis_from_tissue_names(tissue_names)
	output_stem = output_stem + 'remove_testis_'

# Get array of names of windows
window_names = get_window_names(ukkbb_window_summary_file, preprocessed_tgfm_data_dir, gene_type)


# Get annotation names
anno_names = get_anno_names(tissue_names)

'''
# Load in LDSC-style data
X,chi_sq,snp_weights = load_in_ldsc_style_data(preprocessed_tgfm_data_dir, trait_name, window_names, window_chi_sq_lb=0.0, window_chi_sq_ub=320.0, snp_variance_weight=True)

if remove_testis == 'True':
	X = np.delete(X, (testis_index+1), 1)

# Multivariate standardized updates
# Standardize predictors
X_stand, predictor_sdevs = standardize_predictors(X)
# Print to output
'''
anno_sdev_output_file = output_stem + 'annotation_sdev.txt'
#print_annotation_names_and_sdevs(anno_names, predictor_sdevs, anno_sdev_output_file)


# Multivariate Updates
multivariate_output_stem = output_stem + 'multivariate_results'
#multivariate_sldsc_shell(X, chi_sq, snp_weights, anno_names, num_jacknife_windows, predictor_sdevs, multivariate_output_stem)

'''
# Load in standardized jacknifed taus
standardarized_jacknifed_taus_file = multivariate_output_stem + '_standardized_jacknifed_taus.txt'
standardized_jacknifed_taus = np.loadtxt(standardarized_jacknifed_taus_file)

# Compute standardized jacknifed mean tau and covariance
standardized_jacknifed_covariance, standardized_jacknifed_mean = compute_jacknifed_covariance_matrix(standardized_jacknifed_taus)

# Run sparse regression (note this is working in standardized space)
sparse_sldsc_obj = SPARSE_SLDSC(max_iter=1000)
sparse_sldsc_obj.fit(tau=standardized_jacknifed_mean, tau_cov=standardized_jacknifed_covariance)

# Save sparse sldsc results to output
sparse_sldsc_output_stem = multivariate_output_stem + '_sparse_sldsc_results_'
save_sparse_sldsc_results(sparse_sldsc_obj, anno_sdev_output_file, sparse_sldsc_output_stem)

# Run sparse regression with fixed genotype (note this is working in standardized space)
sparse_sldsc_obj_fixed_geno = SPARSE_SLDSC_FIXED_TERM(max_iter=1000)
sparse_sldsc_obj_fixed_geno.fit(tau=standardized_jacknifed_mean, tau_cov=standardized_jacknifed_covariance)

# Save sparse sldsc with fixed genotype results to output
sparse_sldsc_fixed_geno_output_stem = multivariate_output_stem + '_sparse_sldsc_fixed_genotype_results_'
save_sparse_sldsc_results(sparse_sldsc_obj_fixed_geno, anno_sdev_output_file, sparse_sldsc_fixed_geno_output_stem, fixed_geno=True)
'''


############################################################ 
# Explorations into non-negative analysis
###########################################################
multivariate_non_neg_tau, multivariate_non_neg_tau_cov = load_in_existing_multivariate_non_neg_tau_sldsc_results(trait_name, tgfm_heritability_results_dir, anno_sdev_output_file)
sparse_sldsc_obj_fixed_geno = SPARSE_SLDSC_FIXED_TERM(max_iter=1000)
sparse_sldsc_obj_fixed_geno.fit(tau=multivariate_non_neg_tau, tau_cov=multivariate_non_neg_tau_cov)


# Save sparse sldsc with fixed genotype results to output
sparse_sldsc_fixed_geno_output_stem = multivariate_output_stem + '_sparse_loglss_sldsc_fixed_genotype_results_'
save_sparse_sldsc_results(sparse_sldsc_obj_fixed_geno, anno_sdev_output_file, sparse_sldsc_fixed_geno_output_stem, fixed_geno=True)















###############
# OLD
###############



# Run regression
#multivariate_standardized_output_stem = output_stem + 'multivariate_standardized_results'
#mv_stand_sldsc_intercept, mv_stand_sldsc_tau, mv_stand_sldsc_tau_se = multivariate_standardized_sldsc_shell(X_stand, chi_sq, snp_weights, anno_names, num_jacknife_windows, multivariate_standardized_output_stem, predictor_sdevs)



############################################################ 
# Explorations into non-negative analysis
###########################################################
#multivariate_non_neg_tau, multivariate_non_neg_tau_cov = load_in_existing_multivariate_non_neg_tau_sldsc_results(trait_name, '/n/scratch3/users/b/bes710/causal_eqtl_gwas/gtex/tgfm_heritability_results/')

#sparse_sldsc_obj = SPARSE_SLDSC(max_iter=1000)
#sparse_sldsc_obj.fit(tau=multivariate_non_neg_tau, tau_cov=multivariate_non_neg_tau_cov)

############################################################ 
# Explorations into standardized multivariate summary statistic analysis
###########################################################
# Multivariate taus
#multivariate_standardized_tau = load_in_multivariate_tau_sldsc_results(multivariate_standardized_output_stem + '_organized.txt')
# Tau covariance matrix
#multivariate_standardized_tau_covariance = np.loadtxt(multivariate_standardized_output_stem + '_jacknifed_covariance.txt')

#sparse_sldsc_obj = SPARSE_SLDSC_FIXED_TERM(max_iter=1000)
#sparse_sldsc_obj.fit(tau=multivariate_standardized_tau, tau_cov=multivariate_standardized_tau_covariance)
#sparse_sldsc_obj.fit(tau=jacknifed_mean, tau_cov=jacknifed_covariance)
#pdb.set_trace()


# Run Sparse SLDSC
#sparse_sldsc_obj = SPARSE_SLDSC(max_iter=1000)
#sparse_sldsc_obj.fit(tau=multivariate_standardized_tau, tau_cov=multivariate_standardized_tau_covariance)

#sparse_sldsc_obj = SPARSE_SLDSC()
#sparse_sldsc_obj2.fit(tau=multivariate_standardized_tau[1:], tau_cov=multivariate_standardized_tau_covariance[1:,:][:,1:])








############################################################ 
# Explorations into multivariate summary statistic analysis
###########################################################
# Load in some data just generated by the above functions

# multivariate taus
#multivariate_tau = load_in_multivariate_tau_sldsc_results(multivariate_output_stem + '_organized.txt')
# tau covariance matrix
#multivariate_tau_covariance = np.loadtxt(multivariate_output_stem + '_jacknifed_covariance.txt')

# Run Sparse SLDSC
#sparse_sldsc_obj = SPARSE_SLDSC()
#sparse_sldsc_obj.fit(tau=multivariate_tau, tau_cov=multivariate_tau_covariance)

# Save sparse sldsc results
#save_sparse_sldsc_results(sparse_sldsc_obj, multivariate_output_stem + '_sparse_sldsc_')





















############################################################ 
#Explorations into univariate summary statistic analysis
# Does not seem to be currently working. Presumably because each effect explains a relatively large variance?
###########################################################
'''

# Compute LD-scores LD (LDSLD)
#LDSLD = np.corrcoef(np.transpose(X))
ldscores_ld_output_file = output_stem + 'ld_scores_ld.txt'
np.savetxt(ldscores_ld_output_file, LDSLD, fmt="%s", delimiter='\t')

# Univariate Updates (each annotation independently)
univariate_output_stem = output_stem + 'univariate_results'
#univariate_sldsc_shell(X, chi_sq, snp_weights, anno_names, num_jacknife_windows, univariate_output_stem, mv_sldsc_intercept)

# Univariate Updates (each annotation independently)
univariate_learn_intercept_output_stem = output_stem + 'univariate_learn_intercept_results'
#univariate_learn_intercept_sldsc_shell(X, chi_sq, snp_weights, anno_names, num_jacknife_windows, univariate_learn_intercept_output_stem)



# Load in some data just generated by the above functions
# marginal taus
marginal_tau, marginal_tau_se, marginal_tau_z = load_in_marginal_tau_sldsc_results(univariate_learn_intercept_output_stem + '_organized.txt')
# LDSLD
LDSLD = np.loadtxt(ldscores_ld_output_file)


m=len(marginal_tau)

#susie_variant_obj = susieR_pkg.susie_rss(bhat=marginal_tau.reshape((m,1)), shat=marginal_tau_se.reshape((m,1)), R=LDSLD, n=200000.0)
susie_variant_obj = susieR_pkg.susie_rss(z=marginal_tau_z.reshape((m,1)), R=LDSLD, n=200000.0, estimate_residual_variance=True, var_y=np.var(chi_sq))
'''

