import numpy as np 
import os
import sys
import pdb
import math
import scipy.stats
import scipy.optimize
import scipy.special


def log_sum(lnbfs, weights=None):
	# Taken from sumstats package
	"""Sum of a sequence of bayes factors in logarithmic space
	# Basically this np.log(np.sum(np.exp(lnbfs)))
	Parameters
	----------
	lnbfs
		sequence of (natural) log bayes factors
	weights
		sequence of weights for `lnbfs`
	
	Returns
	-------
	float
		the logarithm of the sum of the bayes factors
	"""
	
	lnbfs = tuple(lnbfs)
	if not weights:
		weights = (1,) * len(lnbfs)
	max_lnbf = max(lnbfs)
	try:
		return (
			max_lnbf + math.log(
				math.fsum(
				   math.exp(lnbf - max_lnbf) * weight
				   for lnbf, weight in zip(lnbfs, weights)
				)
			)
		)
	except ValueError:
		if len(lnbfs) == 2:
			return min(lnbfs)
		else:
			raise RuntimeError(
				'The sum of absolute bayes factors may have been rounded to '
				'zero due to a pathalogically high maximum'
			)


def log_diff(x,y):
	my_max = np.max((x,y))
	res = my_max + np.log(np.exp(x-my_max) - np.exp(y-my_max))
	return res


def get_coloc_object(study_name, snp_names, sample_size, beta_vec, var_beta_vec):
	dicti = {}
	dicti['study_name'] = study_name
	dicti['snp_names'] = snp_names
	dicti['sample_size'] = sample_size
	dicti['beta'] = beta_vec
	dicti['var_beta'] = var_beta_vec
	return dicti

def get_causal_coloc_prob_from_pph4_vec(pph4s):
	num_tissues = len(pph4s)
	not_pph4s = 1.0 - pph4s
	log_pph4s = np.log(pph4s)
	log_not_pph4s = np.log(not_pph4s)
	causal_prob_un_normalized_arr = []
	for tissue_num in range(num_tissues):
		log_prob = log_pph4s[tissue_num] + np.sum(np.delete(log_not_pph4s, tissue_num))
		causal_prob_un_normalized_arr.append(log_prob)
	causal_prob_un_normalized_arr = np.asarray(causal_prob_un_normalized_arr)

	causal_prob = np.exp(causal_prob_un_normalized_arr - log_sum(causal_prob_un_normalized_arr))
	return causal_prob

def get_causal_coloc_prob_from_pph4_vec_with_prior(pph4s, prior):
	num_tissues = len(pph4s)
	not_pph4s = 1.0 - pph4s
	log_pph4s = np.log(pph4s)
	log_not_pph4s = np.log(not_pph4s)
	causal_prob_un_normalized_arr = []
	for tissue_num in range(num_tissues):
		log_prob = log_pph4s[tissue_num] + np.sum(np.delete(log_not_pph4s, tissue_num)) + np.log(prior[tissue_num])
		causal_prob_un_normalized_arr.append(log_prob)
	causal_prob_un_normalized_arr = np.asarray(causal_prob_un_normalized_arr)

	causal_prob = np.exp(causal_prob_un_normalized_arr - log_sum(causal_prob_un_normalized_arr))
	return causal_prob

def get_causal_coloc_prob_from_pph4_vec_v2(pph4s):
	num_tissues = len(pph4s)
	not_pph4s = 1.0 - pph4s
	log_pph4s = np.log(pph4s)
	log_not_pph4s = np.log(not_pph4s)
	causal_prob_un_normalized_arr = []
	for tissue_num in range(num_tissues):
		log_prob = np.sum(np.delete(log_not_pph4s, tissue_num))
		causal_prob_un_normalized_arr.append(log_prob)
	causal_prob_un_normalized_arr = np.asarray(causal_prob_un_normalized_arr)

	causal_prob = np.exp(causal_prob_un_normalized_arr - log_sum(causal_prob_un_normalized_arr))
	return causal_prob

def get_causal_coloc_prob_from_pph4_vec_v2_with_prior(pph4s, prior):
	num_tissues = len(pph4s)
	not_pph4s = 1.0 - pph4s
	log_pph4s = np.log(pph4s)
	log_not_pph4s = np.log(not_pph4s)
	causal_prob_un_normalized_arr = []
	for tissue_num in range(num_tissues):
		log_prob = np.sum(np.delete(log_not_pph4s, tissue_num)) + np.log(prior[tissue_num])
		causal_prob_un_normalized_arr.append(log_prob)
	causal_prob_un_normalized_arr = np.asarray(causal_prob_un_normalized_arr)

	causal_prob = np.exp(causal_prob_un_normalized_arr - log_sum(causal_prob_un_normalized_arr))
	return causal_prob


def meta_analysis(effects, se, method='random', weights=None):
	# From Omer Weissbrod
	assert method in ['fixed', 'random']
	d = effects
	variances = se**2

	#compute random-effects variance tau2
	vwts = 1.0 / variances
	fixedsumm = vwts.dot(d) / vwts.sum()
	Q = np.sum(((d - fixedsumm)**2) / variances)
	df = len(d)-1
	tau2 = np.maximum(0, (Q-df) / (vwts.sum() - vwts.dot(vwts) / vwts.sum()))

	#defing weights
	if weights is None:
		if method == 'fixed':
			wt = 1.0 / variances
		else:
			wt = 1.0 / (variances + tau2)
	else:
		wt = weights

	#compute summtest
	summ = wt.dot(d) / wt.sum()
	if method == 'fixed':
		varsum = np.sum(wt*wt*variances) / (np.sum(wt)**2)
	else:
		varsum = np.sum(wt*wt*(variances+tau2)) / (np.sum(wt)**2)
	###summtest = summ / np.sqrt(varsum)

	summary=summ
	se_summary=np.sqrt(varsum)

	return summary, se_summary

def get_approx_log_bf_estimates(coloc_object, sd_prior=0.15):
	z = coloc_object['beta']/np.sqrt(coloc_object['var_beta'])
	r = np.square(sd_prior)/(np.square(sd_prior) + coloc_object['var_beta'])
	lABF = 0.5*(np.log(1.0 - r) + (r*np.square(z)))
	return lABF

def get_snp_pph4_from_log_bf(lbf_1, lbf_2):
	internal_sum_lbf = lbf_1 + lbf_2
	denom_log_abf = log_sum(internal_sum_lbf)
	snp_pph4 = np.exp(internal_sum_lbf - denom_log_abf)
	return snp_pph4

def get_log_bayes_sums(l1, l2):
	lsum = l1 + l2

	lb_sum_h1 = log_sum(l1)
	lb_sum_h2 = log_sum(l2)
	lb_sum_h3 =  log_diff(log_sum(l1) + log_sum(l2), log_sum(lsum))
	lb_sum_h4 = log_sum(lsum)
	return lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4

def run_coloc_with_precomputed_log_bayes_sums(lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4, p0=1.0, p1=1e-4, p2=1e-4, p12=1e-5):
	lb_sum_h0 = 0.0

	joint_sum_h0 = lb_sum_h0 + 0.0
	joint_sum_h1 = lb_sum_h1 + np.log(p1) - np.log(p0)
	joint_sum_h2 = lb_sum_h2 + np.log(p2) - np.log(p0)
	joint_sum_h3 = lb_sum_h3 + np.log(p1) + np.log(p2) - (2.0*np.log(p0))
	joint_sum_h4 = lb_sum_h4 + np.log(p12) - np.log(p0)

	all_abf = np.asarray([joint_sum_h0, joint_sum_h1, joint_sum_h2, joint_sum_h3, joint_sum_h4])
	denom_log_abf = log_sum(all_abf)
	pps = np.exp(all_abf - denom_log_abf)
	return pps

def joint_neg_log_pph4(w_init, trait_approx_lbf, eqtl_approx_lbfs_arr):
	weighted_eqtl_approx_lbf = np.dot(w_init, eqtl_approx_lbfs_arr)
	# Get log bayes sums (basically summary stats for the gene
	lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4 = get_log_bayes_sums(weighted_eqtl_approx_lbf, trait_approx_lbf)

	# Keep track of Coloc probabilities
	pph_vec = run_coloc_with_precomputed_log_bayes_sums(lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4)

	return -np.log(pph_vec[-1])

def joint_neg_log_pph4_unconstrained_simplex(y_init, trait_approx_lbf, eqtl_approx_lbfs_arr):
	w_init = reverse_stickbreaking_simplex_transform(y_init)
	neg_log_like = joint_neg_log_pph4(w_init, trait_approx_lbf, eqtl_approx_lbfs_arr)
	return neg_log_like

def joint_pph4(w_init, trait_approx_lbf, eqtl_approx_lbfs_arr):
	weighted_eqtl_approx_lbf = np.dot(w_init, eqtl_approx_lbfs_arr)
	# Get log bayes sums (basically summary stats for the gene
	lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4 = get_log_bayes_sums(weighted_eqtl_approx_lbf, trait_approx_lbf)

	# Keep track of Coloc probabilities
	pph_vec = run_coloc_with_precomputed_log_bayes_sums(lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4)

	return pph_vec[-1]

def joint_pph(w_init, trait_approx_lbf, eqtl_approx_lbfs_arr):
	weighted_eqtl_approx_lbf = np.dot(w_init, eqtl_approx_lbfs_arr)
	# Get log bayes sums (basically summary stats for the gene
	lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4 = get_log_bayes_sums(weighted_eqtl_approx_lbf, trait_approx_lbf)

	# Keep track of Coloc probabilities
	pph_vec = run_coloc_with_precomputed_log_bayes_sums(lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4)

	return pph_vec

def joint_neg_log_pph4_unconstrained_simplex_meta_analysis(y_init, trait_beta, trait_var_beta, eqtl_beta_arr, eqtl_var_beta_arr):
	w_init = reverse_stickbreaking_simplex_transform(y_init)
	num_snps = eqtl_beta_arr.shape[1]
	meta_analyzed_betas = []
	meta_analyzed_std_errs = []
	for snp_num in range(num_snps):
		ma_beta, ma_se = meta_analysis(eqtl_beta_arr[:,snp_num], np.sqrt(eqtl_var_beta_arr[:,snp_num]), method='fixed', weights=w_init)
		meta_analyzed_betas.append(ma_beta)
		meta_analyzed_std_errs.append(ma_se)
	meta_analyzed_betas = np.asarray(meta_analyzed_betas)
	meta_analyzed_varz = np.square(np.asarray(meta_analyzed_std_errs))
	trait_approx_lbf = get_approx_log_bf_estimates({'beta':trait_beta, 'var_beta':trait_var_beta})
	ma_eqtl_approx_lbf = get_approx_log_bf_estimates({'beta':meta_analyzed_betas, 'var_beta':meta_analyzed_varz})



	lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4 = get_log_bayes_sums(ma_eqtl_approx_lbf, trait_approx_lbf)

	# Keep track of Coloc probabilities
	pph_vec = run_coloc_with_precomputed_log_bayes_sums(lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4)

	return -np.log(pph_vec[-1])

def joint_pph4_unconstrained_simplex_meta_analysis(y_init, trait_beta, trait_var_beta, eqtl_beta_arr, eqtl_var_beta_arr):
	w_init = reverse_stickbreaking_simplex_transform(y_init)
	num_snps = eqtl_beta_arr.shape[1]
	meta_analyzed_betas = []
	meta_analyzed_std_errs = []
	for snp_num in range(num_snps):
		ma_beta, ma_se = meta_analysis(eqtl_beta_arr[:,snp_num], np.sqrt(eqtl_var_beta_arr[:,snp_num]), method='fixed', weights=w_init)
		meta_analyzed_betas.append(ma_beta)
		meta_analyzed_std_errs.append(ma_se)
	meta_analyzed_betas = np.asarray(meta_analyzed_betas)
	meta_analyzed_varz = np.square(np.asarray(meta_analyzed_std_errs))
	trait_approx_lbf = get_approx_log_bf_estimates({'beta':trait_beta, 'var_beta':trait_var_beta})
	ma_eqtl_approx_lbf = get_approx_log_bf_estimates({'beta':meta_analyzed_betas, 'var_beta':meta_analyzed_varz})



	lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4 = get_log_bayes_sums(ma_eqtl_approx_lbf, trait_approx_lbf)

	# Keep track of Coloc probabilities
	pph_vec = run_coloc_with_precomputed_log_bayes_sums(lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4)

	return pph_vec




# Go from K dimensional x (simplex) space to K-1 dimensional unconstrained space
def forward_stickbreaking_simplex_transform(x):
	K = len(x)
	y = []
	for k in range(K-1):
		z_k = x[k]/(1.0 - np.sum(x[:k]))
		y_k = scipy.special.logit(z_k) - np.log(1.0/(K - k - 1.0))
		y.append(y_k)
	return np.asarray(y)

# Go from K-1 dimensional unconstrained space to K dimensional x (simplex) space
def reverse_stickbreaking_simplex_transform(y):
	x = []
	K = len(y) + 1
	for k in range(K-1):
		z_k = scipy.special.expit(y[k] + np.log(1.0/(K - k - 1.0)))
		x_k = (1.0 - sum(x))*z_k
		x.append(x_k)
	x_K = 1.0 - sum(x)
	x.append(x_K)
	return np.asarray(x)

def competitive_coloc_study_posterior_probability_estimation(trait_approx_lbf, eqtl_approx_lbfs_arr, observed_pph4s):
	num_studies = eqtl_approx_lbfs_arr.shape[0]
	#w_init = observed_pph4s/np.sum(observed_pph4s)
	w_init = np.zeros(len(observed_pph4s)) + (.01/len(observed_pph4s))
	w_init[np.argmax(observed_pph4s)] = .99
	print(np.sum(w_init))

	y_init = forward_stickbreaking_simplex_transform(w_init)
	#w_init2 = reverse_stickbreaking_simplex_transform(y_init)
	#pdb.set_trace()
	#log_pph4 = joint_neg_log_pph4(w_init, trait_approx_lbf, eqtl_approx_lbfs_arr)
	val = scipy.optimize.fmin_l_bfgs_b(joint_neg_log_pph4_unconstrained_simplex, y_init, args=(trait_approx_lbf,eqtl_approx_lbfs_arr), approx_grad=True)

	w_opt = reverse_stickbreaking_simplex_transform(val[0])
	competitive_pps = joint_pph(w_opt, trait_approx_lbf, eqtl_approx_lbfs_arr)
	return competitive_pps, w_opt, val[2]['warnflag']


def competitive_coloc_study_posterior_probability_estimation_via_meta_analysis(trait_beta, trait_var_beta, eqtl_beta_arr, eqtl_var_beta_arr, observed_pph4s):
	num_studies = eqtl_beta_arr.shape[0]
	w_init = np.zeros(len(observed_pph4s)) + (.01/len(observed_pph4s))
	w_init[np.argmax(observed_pph4s)] = .99

	#w_init = np.zeros(len(observed_pph4s)) 
	#w_init[0] = 1.0

	y_init = forward_stickbreaking_simplex_transform(w_init)

	neg_log_pph4 = joint_neg_log_pph4_unconstrained_simplex_meta_analysis(y_init, trait_beta, trait_var_beta, eqtl_beta_arr, eqtl_var_beta_arr)

	val = scipy.optimize.fmin_l_bfgs_b(joint_neg_log_pph4_unconstrained_simplex_meta_analysis, y_init, args=( trait_beta, trait_var_beta, eqtl_beta_arr, eqtl_var_beta_arr), approx_grad=True)


	w_opt = reverse_stickbreaking_simplex_transform(val[0])

	joint_pph4_vec = joint_pph4_unconstrained_simplex_meta_analysis(val[0], trait_beta, trait_var_beta, eqtl_beta_arr, eqtl_var_beta_arr)

	return joint_pph4_vec, w_opt, val[2]['warnflag']








