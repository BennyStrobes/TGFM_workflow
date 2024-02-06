import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import pickle
import tgfm
import rpy2
import rpy2.robjects.numpy2ri as numpy2ri
import rpy2.robjects as ro
import scipy.stats
ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
susieR_pkg = importr('susieR')


def extract_tissue_names(gtex_pseudotissue_file, remove_testis=False):
	f = open(gtex_pseudotissue_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0] == 'Testis' and remove_testis:
			continue
		arr.append(data[0])
	f.close()
	return np.asarray(arr)


def standardize_eqtl_pmces_old(eqtl_pmces, variant_ld):

	#gene_variances = np.diag(np.dot(np.dot(eqtl_pmces, variant_ld), np.transpose(eqtl_pmces)))
	gene_variances = np.sum(np.dot(eqtl_pmces, variant_ld)*eqtl_pmces, axis=1)
	pdb.set_trace()

	# compute gene variance
	n_genes = eqtl_pmces.shape[0]

	# standardize
	for gene_iter in range(n_genes):
		eqtl_pmces[gene_iter,:] = eqtl_pmces[gene_iter,:]/np.sqrt(gene_variances[gene_iter])

	return eqtl_pmces

def standardize_eqtl_pmces(eqtl_pmces, gene_variances, reference_ld):
	# compute gene variance
	n_genes = eqtl_pmces.shape[0]

	# standardize
	for gene_iter in range(n_genes):
		#np.dot(np.dot(eqtl_pmces[gene_iter,:], reference_ld), (eqtl_pmces[gene_iter,:]))
		eqtl_pmces[gene_iter,:] = eqtl_pmces[gene_iter,:]/np.sqrt(gene_variances[gene_iter])

	return eqtl_pmces


def extract_full_gene_variant_ld(standardized_eqtl_effects, variant_ld):
	expression_covariance = np.dot(np.dot(standardized_eqtl_effects, variant_ld), np.transpose(standardized_eqtl_effects))
	np.fill_diagonal(expression_covariance, 1.0)
	dd = np.diag(1.0/np.sqrt(np.diag(expression_covariance)))
	ge_ld = np.dot(np.dot(dd, expression_covariance),dd)
	gene_variant_ld = np.dot(standardized_eqtl_effects,variant_ld) # Ngenes X n_variants
	top = np.hstack((ge_ld, gene_variant_ld))
	bottom = np.hstack((np.transpose(gene_variant_ld), variant_ld))
	full_ld = np.vstack((top,bottom))
	return full_ld

def get_credible_set_genes(phi_comp, cs_thresh):
	ordered_genes = np.argsort(-phi_comp)
	cs_genes = []
	cs_counter = 0.0
	for gene_index in ordered_genes:
		if cs_counter < cs_thresh:
			cs_genes.append(gene_index)
		cs_counter = cs_counter + phi_comp[gene_index]
	cs_genes = np.asarray(cs_genes)
	return cs_genes

def load_in_log_priors(log_prior_prob_file, variant_names, gene_names):
	prior_prob_raw = np.loadtxt(log_prior_prob_file, dtype=str, delimiter='\t')[1:,:]
	n_var = len(variant_names)
	variant_log_prob = prior_prob_raw[:n_var,1].astype(float)
	gene_log_prob = prior_prob_raw[n_var:,1].astype(float)
	temp_var_names = prior_prob_raw[:n_var,0]
	temp_gene_names = prior_prob_raw[n_var:,0]

	# error check
	if np.array_equal(temp_var_names, variant_names) == False:
		print('assumption eroror')
		pdb.set_trace()
	if len(gene_names) > 0:
		if np.array_equal(temp_gene_names, gene_names) == False:
			print('assumption eroror')
			pdb.set_trace()
	
	return variant_log_prob, gene_log_prob


def extract_valid_joint_susie_components_from_full_ld(alpha_phi, beta_phi, full_ld, ld_thresh, subset_n = 100):
	num_components = alpha_phi.shape[0]
	valid_components = []

	num_genes = alpha_phi.shape[1]
	num_variants = beta_phi.shape[1]

	for component_num in range(num_components):
		cs_predictors = get_credible_set_genes(np.hstack((alpha_phi[component_num,:], beta_phi[component_num,:])), .95)

		if subset_n > len(cs_predictors):
			# absolute ld among genes and variants in credible set
			if np.min(np.abs(full_ld[cs_predictors,:][:, cs_predictors])) > ld_thresh:
				valid_components.append(component_num)
		else:
			# First run subsetted analysis
			cs_predictors_subset = np.random.choice(cs_predictors, size=subset_n, replace=False, p=None)
			if np.min(np.abs(full_ld[cs_predictors_subset,:][:, cs_predictors_subset])) > ld_thresh:
				if np.min(np.abs(full_ld[cs_predictors,:][:, cs_predictors])) > ld_thresh:
					valid_components.append(component_num)

	return valid_components

def extract_middle_genetic_elements(ordered_genes, middle_gene_indices, ordered_variants, middle_variant_indices):
	dicti = {}
	for gene_name in ordered_genes[middle_gene_indices]:
		dicti[gene_name] = 1
	for variant_name in ordered_variants[middle_variant_indices.astype(int)]:
		dicti[variant_name] = 1
	return dicti

def get_tissues_from_full_gene_names(gene_names):
	tissues = []
	for gene_name in gene_names:
		gene_info = gene_name.split('_')
		tissue = '_'.join(gene_info[1:])
		tissues.append(tissue)
	return np.asarray(tissues)

def get_probability_coming_from_each_tissue(ordered_tissue_names, tissue_to_position_mapping, window_gene_names, window_gene_probs):
	tiss_probs = np.zeros(len(ordered_tissue_names))
	window_tissue_names = get_tissues_from_full_gene_names(window_gene_names)
	for gene_index, window_tissue_name in enumerate(window_tissue_names):
		tissue_position = tissue_to_position_mapping[window_tissue_name]
		tiss_probs[tissue_position] = tiss_probs[tissue_position] + window_gene_probs[gene_index]
	return tiss_probs

def run_susie_debug(tgfm_obj, gene_variant_full_ld,tgfm_data):
	z_vec = np.hstack((tgfm_obj.nominal_twas_z,tgfm_data['gwas_beta']/tgfm_data['gwas_beta_se']))
	susie_variant_obj_orig = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20)
	susie_alpha = susie_variant_obj_orig.rx2('alpha')
	pdb.set_trace()

def elbo_calc(z_vec, LD, samp_size, alpha, mu, mu2, KL_terms):
	bb = alpha*mu
	b_bar = np.sum(bb,axis=0)
	postb2 = alpha*mu2
	elbo_term1 = samp_size -1
	elbo_term2 = -2.0*np.sum(np.sqrt(samp_size-1)*b_bar*z_vec)
	elbo_term3 = np.sum(b_bar*np.dot((samp_size-1.0)*LD, b_bar))
	elbo_term4 = - np.sum(np.dot(bb, (samp_size-1.0)*LD)*bb)
	elbo_term5 = np.sum(np.dot(np.diag(LD*(samp_size-1)), np.transpose(postb2)))
	elbo_term7 = (-samp_size/2.0)*np.log(2.0*np.pi)

	elbo = elbo_term7 - .5*(elbo_term1 + elbo_term2 + elbo_term3 + elbo_term4 + elbo_term5) - np.sum(KL_terms)
	return elbo

def compute_pips(alpha_mat):
	LL = alpha_mat.shape[0]
	n_elements = alpha_mat.shape[1]
	anti_pips = np.ones(n_elements)

	for component_iter in range(LL):
		anti_pips = anti_pips*(1.0 - alpha_mat[component_iter,:])
	pips = 1.0 - anti_pips
	return pips

def update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_obj):
	# Alphas
	susie_alpha = susie_obj.rx2('alpha')
	tgfm_obj.alpha_phi = susie_alpha[:,:(tgfm_obj.G)]
	tgfm_obj.beta_phi = susie_alpha[:,(tgfm_obj.G):]

	# Mus
	susie_mu = susie_obj.rx2('mu')
	tgfm_obj.alpha_mu = susie_mu[:,:(tgfm_obj.G)]
	tgfm_obj.beta_mu = susie_mu[:,(tgfm_obj.G):]

	# susie_mu_var
	susie_mu_var = susie_obj.rx2('mu2') - np.square(susie_mu)
	tgfm_obj.alpha_var = susie_mu_var[:,:(tgfm_obj.G)]
	tgfm_obj.beta_var = susie_mu_var[:,(tgfm_obj.G):]

	# Compute pips
	pips = compute_pips(susie_obj.rx2('alpha'))
	tgfm_obj.alpha_pip = pips[:(tgfm_obj.G)]
	tgfm_obj.beta_pip = pips[(tgfm_obj.G):]

	lbf = susie_obj.rx2('lbf_variable')
	tgfm_obj.alpha_lbf = lbf[:, :(tgfm_obj.G)]
	tgfm_obj.beta_lbf = lbf[:, (tgfm_obj.G):]

	return tgfm_obj

def tgfm_inference_shell(tgfm_data, gene_log_prior, var_log_prior, gene_variant_full_ld, init_method, est_resid_var_bool):
	# Hacky: Initialize old TGFM object using only one iter of optimization
	tgfm_obj = tgfm.TGFM(L=15, estimate_prior_variance=False, gene_init_log_pi=gene_log_prior, variant_init_log_pi=var_log_prior, convergence_thresh=1e-5, max_iter=1)
	tgfm_obj.fit(twas_data_obj=tgfm_data)
	# More hack: need to redo twas z
	variant_z = tgfm_data['gwas_beta']/tgfm_data['gwas_beta_se']
	new_gene_z = np.dot(tgfm_data['gene_eqtl_pmces'], variant_z)

	#tgfm_obj_temp = tgfm.TGFM(L=20, estimate_prior_variance=True, gene_init_log_pi=gene_log_prior, variant_init_log_pi=var_log_prior, convergence_thresh=1e-8, max_iter=10)
	#tgfm_obj_temp.fit(twas_data_obj=tgfm_data)

	tgfm_obj.nominal_twas_z = new_gene_z

	# Create vector of concatenated z-scores
	z_vec = np.hstack((new_gene_z,variant_z))
	#tgfm_obj_temp_elbo = elbo_calc(z_vec, gene_variant_full_ld, tgfm_data['gwas_sample_size'], np.hstack((tgfm_obj_temp.alpha_phi, tgfm_obj_temp.beta_phi)), np.hstack((tgfm_obj_temp.alpha_mu, tgfm_obj_temp.beta_mu)), np.square(np.hstack((tgfm_obj_temp.alpha_mu, tgfm_obj_temp.beta_mu))) + np.hstack((tgfm_obj_temp.alpha_var, tgfm_obj_temp.beta_var)), tgfm_obj_temp.KL_terms)


	# Create concatenated vector of prior probs
	prior_probs = np.hstack((np.exp(gene_log_prior), np.exp(var_log_prior)))



	if init_method == 'best':
		rpy2.robjects.r['options'](warn=1)

		# Run susie with null initialization
		susie_null_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=15, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool)
		tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_null_init)

		# Only run alternative variant-only inititialization if we have high pip genes
		if np.max(tgfm_obj.alpha_pip) > .2:
			# Run susie with only variants
			p_var_only = np.ones(len(z_vec))
			p_var_only[:len(tgfm_obj.nominal_twas_z)] = 0.0
			p_var_only = p_var_only/np.sum(p_var_only)
			susie_variant_only = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=15, prior_weights=p_var_only.reshape((len(p_var_only),1)), estimate_residual_variance=est_resid_var_bool)

			# Run susie with variant initialization
			init_obj = {'alpha':susie_variant_only.rx2('alpha'), 'mu':susie_variant_only.rx2('mu'),'mu2':susie_variant_only.rx2('mu2')}
			init_obj2 = ro.ListVector(init_obj)
			init_obj2.rclass = rpy2.robjects.StrVector(("list", "susie"))
			susie_variant_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=15, s_init=init_obj2, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool)

			# Select model with largest elbo
			susie_null_init_elbo = susie_null_init.rx2('elbo')[-1]
			susie_variant_init_elbo = susie_variant_init.rx2('elbo')[-1]

			if susie_variant_init_elbo > susie_null_init_elbo:
				# Variant init wins
				tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_variant_init)
			else:
				# Null init wins
				tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_null_init)
	elif init_method == 'refine_best':
		rpy2.robjects.r['options'](warn=1)
		# Run susie with only variants
		p_var_only = np.ones(len(z_vec))
		p_var_only[:len(tgfm_obj.nominal_twas_z)] = 0.0
		p_var_only = p_var_only/np.sum(p_var_only)
		susie_variant_only = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=p_var_only.reshape((len(p_var_only),1)), estimate_residual_variance=est_resid_var_bool)
	

		'''
		susie_alt_obj = susie_alt.SUSIE_ALT(L=20, estimate_prior_variance=False, estimate_prior_prob=True, alpha_0=1e-2,convergence_thresh=1e-5, max_iter=100, ard_a_prior=0.0, ard_b_prior=0.0)
		se = np.ones(len(z_vec))/np.sqrt(100000)
		beta = z_vec*se
		susie_alt_obj.fit(beta=beta, beta_se=se, LD=gene_variant_full_ld)
		'''
		
		# Run susie with null initialization
		susie_null_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool,tol=1e-6)

		# Run susie with variant initialization
		init_obj = {'alpha':susie_variant_only.rx2('alpha'), 'mu':susie_variant_only.rx2('mu'),'mu2':susie_variant_only.rx2('mu2')}
		init_obj2 = ro.ListVector(init_obj)
		init_obj2.rclass = rpy2.robjects.StrVector(("list", "susie"))
		susie_variant_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, s_init=init_obj2, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool)
		
		#res=susie_inf.susie(z_vec, 1.0, n=100000,L=20, LD=gene_variant_full_ld)
		susie_variant_refine_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, s_init=init_obj2, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool, refine=True)

		# Select model with largest elbo
		susie_null_init_elbo = susie_null_init.rx2('elbo')[-1]
		susie_variant_init_elbo = susie_variant_init.rx2('elbo')[-1]
		susie_variant_refine_init_elbo = susie_variant_refine_init.rx2('elbo')[-1]


		'''
		# Run susie with only genes
		p_gene_only = np.zeros(len(z_vec))
		p_gene_only[:len(tgfm_obj.nominal_twas_z)] = 1.0
		p_gene_only = p_gene_only/np.sum(p_gene_only)
		susie_gene_only = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=p_gene_only.reshape((len(p_gene_only),1)), estimate_residual_variance=est_resid_var_bool)

		init_obj = {'alpha':susie_gene_only.rx2('alpha'), 'mu':susie_gene_only.rx2('mu'),'mu2':susie_gene_only.rx2('mu2')}
		init_obj2 = ro.ListVector(init_obj)
		init_obj2.rclass = rpy2.robjects.StrVector(("list", "susie"))
		susie_variant_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, s_init=init_obj2, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool)
		'''
		#res=susie_inf.susie(z_vec, 1.0, n=100000,L=20, LD=gene_variant_full_ld)

		if susie_variant_init_elbo > susie_null_init_elbo and susie_variant_init_elbo > susie_variant_refine_init_elbo:
			# Variant init wins
			tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_variant_init)
		elif susie_null_init_elbo > susie_variant_init_elbo and susie_null_init_elbo > susie_variant_refine_init_elbo:
			# Null init wins
			tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_null_init)
		else:
			# Refine init wins
			tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_variant_refine_init)
	elif init_method == 'null':
		# Run susie with null initialization
		susie_null_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool)
		tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_null_init)
	elif init_method == 'variant_only':
		# Run susie with only variants
		p_var_only = np.ones(len(z_vec))
		p_var_only[:len(tgfm_obj.nominal_twas_z)] = 0.0
		p_var_only = p_var_only/np.sum(p_var_only)
		susie_variant_only = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=p_var_only.reshape((len(p_var_only),1)), estimate_residual_variance=est_resid_var_bool)
	
		# Run susie with variant initialization
		init_obj = {'alpha':susie_variant_only.rx2('alpha'), 'mu':susie_variant_only.rx2('mu'),'mu2':susie_variant_only.rx2('mu2')}
		init_obj2 = ro.ListVector(init_obj)
		init_obj2.rclass = rpy2.robjects.StrVector(("list", "susie"))
		susie_variant_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, s_init=init_obj2, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool)
		tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_variant_init)
	else:
		print('assumption errror: susie initialization method ' + init_method + ' not recognized')
		pdb.set_trace()
	return tgfm_obj

def correct_ld_mat_for_af_standardization(ld_mat):
	n_snps = ld_mat.shape[0]
	correction = 1.0/np.diag(ld_mat)
	for snp_iter in range(n_snps):
		ld_mat[:,snp_iter] = ld_mat[:,snp_iter]*np.sqrt(correction[snp_iter])
		ld_mat[snp_iter,:] = ld_mat[snp_iter,:]*np.sqrt(correction[snp_iter])
	return ld_mat


# Convert gwas summary statistics to *STANDARDIZED* effect sizes
# Following SuSiE code found in these two places:
########1. https://github.com/stephenslab/susieR/blob/master/R/susie_rss.R  (LINES 277-279)
########2. https://github.com/stephenslab/susieR/blob/master/R/susie_ss.R (LINES 148-156 AND 203-205)
def convert_to_standardized_summary_statistics(gwas_beta_raw, gwas_beta_se_raw, gwas_sample_size, R, sigma2=1.0):
	gwas_z_raw = gwas_beta_raw/gwas_beta_se_raw

	XtX = (gwas_sample_size-1)*R
	Xty = np.sqrt(gwas_sample_size-1)*gwas_z_raw
	var_y = 1

	dXtX = np.diag(XtX)
	csd = np.sqrt(dXtX/(gwas_sample_size-1))
	csd[csd == 0] = 1

	XtX = (np.transpose((1/csd) * XtX) / csd)
	Xty = Xty / csd

	dXtX2 = np.diag(XtX)

	beta_scaled = (1/dXtX2)*Xty
	beta_se_scaled = np.sqrt(sigma2/dXtX2)

	return beta_scaled, beta_se_scaled


def extract_log_prior_probabilities_from_summary_file(log_prior_file, existing_var_names, existing_gene_names):
	var_names = []
	gene_names = []
	var_probs = []
	gene_probs = []

	f = open(log_prior_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		element_name = data[0]
		log_prob = float(data[1])
		if element_name.startswith('ENSG'):
			gene_probs.append(log_prob)
			gene_names.append(element_name)
		else:
			var_probs.append(log_prob)
			var_names.append(element_name)
	f.close()


	# Quick error checking
	var_names = np.asarray(var_names)
	gene_names = np.asarray(gene_names)
	if np.array_equal(var_names, existing_var_names) == False:
		print('assumption eorroror')
		pdb.set_trace()
	if np.array_equal(gene_names, existing_gene_names) == False:
		print('assumption eorroror')
		pdb.set_trace()


	return np.asarray(var_probs), np.asarray(gene_probs)


def extract_log_prior_probabilities(variant_names, gene_names, ln_pi_method_name):
	n_snps = len(variant_names)
	n_genes = len(gene_names)
	n_ele = n_snps + n_genes

	ele_prob = 1.0/n_ele

	snp_prior = np.ones(n_snps)*ele_prob
	gene_prior = np.ones(n_genes)*ele_prob

	return np.log(snp_prior), np.log(gene_prior)

def filter_genes_in_susie_object(susie_obj, valid_gt_indices):
	new_susie_obj = []
	for gt_index in valid_gt_indices:
		new_susie_obj.append(susie_obj[gt_index])
	return new_susie_obj

def fill_in_causal_effect_size_matrix(null_mat, bs_eqtls_pmces_sparse):
	null_mat[bs_eqtls_pmces_sparse[:,0].astype(int), bs_eqtls_pmces_sparse[:,1].astype(int)] = bs_eqtls_pmces_sparse[:,2]
	return null_mat


def filter_tgfm_data_structure_to_remove_tissue0_gene_tissue_pairs(tgfm_data):
	################
	# First get valid gt indices by removing gt-pairs corresponding to tissue 0
	valid_gt_indices = []
	gt_old_to_new = []
	counter = 0
	for ii, gt_pair in enumerate(tgfm_data['genes']):
		if gt_pair.split('_')[1] == 'tissue0':
			gt_old_to_new.append(-2.0)
			continue
		gt_old_to_new.append(float(counter))
		valid_gt_indices.append(ii)
		counter = counter + 1
	if len(valid_gt_indices) == 0:
		return tgfm_data, True
	valid_gt_indices = np.asarray(valid_gt_indices)
	gt_old_to_new = np.asarray(gt_old_to_new)

	################
	# Filter 'genes' key
	tgfm_data['genes'] = tgfm_data['genes'][valid_gt_indices]

	#################
	# Filter 'sparse_sampled_gene_eqtl_pmces' key
	n_bs = len(tgfm_data['sparse_sampled_gene_eqtl_pmces'])
	new_sparse_sampled_gene_eqtl_pmces = []
	for bs_iter in range(n_bs):
		bs_eqtl_effects_sparse_old = tgfm_data['sparse_sampled_gene_eqtl_pmces'][bs_iter]
		bs_eqtl_effects_sparse_new = []
		for row_iter in range(bs_eqtl_effects_sparse_old.shape[0]):
			if gt_old_to_new[int(bs_eqtl_effects_sparse_old[row_iter,0])] == -2.0:
				continue
			new = np.zeros(3)
			new[0] = gt_old_to_new[int(bs_eqtl_effects_sparse_old[row_iter,0])]
			new[1] = bs_eqtl_effects_sparse_old[row_iter,1]
			new[2] = bs_eqtl_effects_sparse_old[row_iter,2]
			bs_eqtl_effects_sparse_new.append(new)
		bs_eqtl_effects_sparse_new = np.asarray(bs_eqtl_effects_sparse_new)
		new_sparse_sampled_gene_eqtl_pmces.append(bs_eqtl_effects_sparse_new)
	## error checking
	#old_mat = fill_in_causal_effect_size_matrix(np.zeros((len(gt_old_to_new), len(tgfm_data['variants']))), tgfm_data['sparse_sampled_gene_eqtl_pmces'][10])
	#new_mat = fill_in_causal_effect_size_matrix(np.zeros((len(valid_gt_indices), len(tgfm_data['variants']))), new_sparse_sampled_gene_eqtl_pmces[10])
	#print(np.array_equal(old_mat[valid_gt_indices,:], new_mat))
	#pdb.set_trace()
	tgfm_data['sparse_sampled_gene_eqtl_pmces'] = new_sparse_sampled_gene_eqtl_pmces
	
	#################
	# Filter 'middle_gene_indices' key
	new_middle_genes_indices = []
	for old_index in tgfm_data['middle_gene_indices']:
		new_index = gt_old_to_new[old_index]
		if new_index == -2.0:
			continue
		new_middle_genes_indices.append(int(new_index))
	new_middle_genes_indices = np.asarray(new_middle_genes_indices)
	tgfm_data['middle_gene_indices'] = new_middle_genes_indices.astype(int)
	
	#################
	# Filter 'gene_eqtl_pmces' key
	tgfm_data['gene_eqtl_pmces'] = tgfm_data['gene_eqtl_pmces'][valid_gt_indices,:]

	#################
	# Filter 'gene_variances' key
	tgfm_data['gene_variances'] = tgfm_data['gene_variances'][valid_gt_indices]

	#################
	# Filter 'tss' key
	tgfm_data['tss'] = tgfm_data['tss'][valid_gt_indices]

	#################
	# Filter 'gene_susie_mu' key
	new_gene_susie_mu = filter_genes_in_susie_object(tgfm_data['gene_susie_mu'], valid_gt_indices)
	tgfm_data['gene_susie_mu'] = new_gene_susie_mu

	#################
	# Filter 'gene_susie_alpha' key
	new_gene_susie_alpha = filter_genes_in_susie_object(tgfm_data['gene_susie_alpha'], valid_gt_indices)
	tgfm_data['gene_susie_alpha'] = new_gene_susie_alpha

	#################
	# Filter 'gene_susie_mu_var' key
	new_gene_susie_mu_var = filter_genes_in_susie_object(tgfm_data['gene_susie_mu_var'], valid_gt_indices)
	tgfm_data['gene_susie_mu_var'] = new_gene_susie_mu_var

	#################
	# Filter 'gene_susie_indices' key
	new_gene_susie_indices = filter_genes_in_susie_object(tgfm_data['gene_susie_indices'], valid_gt_indices)
	tgfm_data['gene_susie_indices'] = new_gene_susie_indices

	return tgfm_data, False


######################
# Command line args
######################
tgfm_input_summary_file = sys.argv[1]
tgfm_output_stem = sys.argv[2]
init_method = sys.argv[3]
est_resid_var_str = sys.argv[4]
ln_pi_method_name = sys.argv[5]
tgfm_tissues = sys.argv[6]





window_pvalue_thresh = 1e-5



if est_resid_var_str == 'False':
	est_resid_var_bool = False
elif est_resid_var_str == 'True':
	est_resid_var_bool = True 
else:
	print('assumption eroror')
	pdb.set_trace()



# Open PIP file handle
pip_output_file = tgfm_output_stem + '_tgfm_pip_summary.txt'
t_pip = open(pip_output_file,'w')
t_pip.write('window_name\tinclusion_elements\tinclusion_probabilities\n')




# Now loop through windows
# In each window run TGFM independently
# Loop through trait components
tgfm_input_data = np.loadtxt(tgfm_input_summary_file,dtype=str,delimiter='\t')
tgfm_input_data = tgfm_input_data[1:,:]



# Get n_windows on this run
n_windows = tgfm_input_data.shape[0]

for window_iter in range(n_windows):
	data = tgfm_input_data[window_iter, :]

	##############################
	# Extract relevent fields
	###############################
	window_name = data[0]
	print(window_name)

	ld_file = data[1]
	tgfm_input_pkl = data[2]

	##############################
	# Load in Data
	###############################
	# Load in tgfm input data
	g = open(tgfm_input_pkl, "rb")
	tgfm_data = pickle.load(g)
	g.close()

	# Skip windows with no genes
	if len(tgfm_data['genes']) == 0:
		print('skipped because of no genes')
		t_pip.write(window_name + '\tNA\tNA\n')
		continue


	if tgfm_tissues == 'no_t0':
		tgfm_data, filter_error_bool = filter_tgfm_data_structure_to_remove_tissue0_gene_tissue_pairs(tgfm_data)
		if filter_error_bool:
			print('skipped because of no genes')
			t_pip.write(window_name + '\tNA\tNA\n')
			continue

	# Extract gwas p
	gwas_z = tgfm_data['gwas_beta']/tgfm_data['gwas_beta_se']
	gwas_p = scipy.stats.norm.sf(abs(gwas_z))*2.0

	# Ignore windows with no pvalues less than some threshold
	if np.min(gwas_p) > window_pvalue_thresh:
		print('skipped because of window pvalue threshold')
		t_pip.write(window_name + '\tNA\tNA\n')
		continue

	# Load in LD
	ld_mat = np.load(ld_file)
	# Add ld to tgfm_data obj
	tgfm_data['reference_ld'] = ld_mat


	# Extract log prior probabilities from summary file
	#log_prior_file = ln_pi_file_stem + '_' + ln_pi_method_name + '.txt'
	#var_log_prior, gene_log_prior = extract_log_prior_probabilities_from_summary_file(log_prior_file, tgfm_data['variants'], tgfm_data['genes'])
	var_log_prior, gene_log_prior = extract_log_prior_probabilities(tgfm_data['variants'], tgfm_data['genes'], ln_pi_method_name)


	# Extract full ld between genes, variants, and gene-variants
	gene_variant_full_ld = extract_full_gene_variant_ld(tgfm_data['gene_eqtl_pmces'], tgfm_data['reference_ld'])


	##############################
	# Run TGFM
	###############################
	tgfm_obj = tgfm_inference_shell(tgfm_data, gene_log_prior, var_log_prior, gene_variant_full_ld, init_method, est_resid_var_bool)

	##############################
	# Organize TGFM data and print to results
	###############################
	# Extract components that pass purity filter
	valid_tgfm_components = extract_valid_joint_susie_components_from_full_ld(tgfm_obj.alpha_phi, tgfm_obj.beta_phi, gene_variant_full_ld, .5)

	# Extract names of genetic elements
	genetic_element_names = np.hstack((tgfm_data['genes'], tgfm_data['variants']))
	# Extract dictionary list of genetic elements in the middel of this window
	middle_genetic_elements = extract_middle_genetic_elements(tgfm_data['genes'], tgfm_data['middle_gene_indices'], tgfm_data['variants'], tgfm_data['middle_variant_indices'])


	# First print PIP results and then print component results
	# Extract genetic element pips
	genetic_element_pips = np.hstack((tgfm_obj.alpha_pip, tgfm_obj.beta_pip))

	# Extract genetic elements and pips only corresponding to middle genetic elements
	middle_pips = []
	middle_names = []
	for genetic_element_iter, genetic_element_name in enumerate(genetic_element_names):
		if genetic_element_name not in middle_genetic_elements:
			continue
		if genetic_element_pips[genetic_element_iter] < .01:
			continue
		middle_pips.append(genetic_element_pips[genetic_element_iter])
		middle_names.append(genetic_element_name)
	middle_pips = np.asarray(middle_pips)
	middle_names = np.asarray(middle_names)

	# Sort middle pips
	indices = np.argsort(-middle_pips)
	ordered_middle_pips = middle_pips[indices]
	ordered_middle_names = middle_names[indices]

	# Write to credible set output
	t_pip.write(window_name + '\t')
	t_pip.write(';'.join(ordered_middle_names) + '\t')
	t_pip.write(';'.join(ordered_middle_pips.astype(str)) + '\n')

	t_pip.flush()

	# Save all TGFM results to pkl
	tgfm_results = {}
	tgfm_results['variants'] = tgfm_data['variants']
	tgfm_results['genes'] = tgfm_data['genes']
	tgfm_results['alpha_phi'] = tgfm_obj.alpha_phi
	tgfm_results['beta_phi'] = tgfm_obj.beta_phi
	tgfm_results['alpha_lbf'] = tgfm_obj.alpha_lbf
	tgfm_results['beta_lbf'] = tgfm_obj.beta_lbf
	tgfm_results['alpha_mu'] = tgfm_obj.alpha_mu
	tgfm_results['beta_mu'] = tgfm_obj.beta_mu
	tgfm_results['alpha_var'] = tgfm_obj.alpha_var
	tgfm_results['beta_var'] = tgfm_obj.beta_var
	tgfm_results['alpha_pip'] = tgfm_obj.alpha_pip
	tgfm_results['beta_pip'] = tgfm_obj.beta_pip
	tgfm_results['component_variances'] = tgfm_obj.component_variances
	tgfm_results['valid_components'] = valid_tgfm_components
	tgfm_results['nominal_twas_z'] = tgfm_obj.nominal_twas_z
	tgfm_results['middle_variant_indices'] = tgfm_data['middle_variant_indices']
	tgfm_results['middle_gene_indices'] = tgfm_data['middle_gene_indices']

	# Write pickle file
	window_tgfm_output_file = tgfm_output_stem + '_' + window_name + '_results.pkl'
	g = open(window_tgfm_output_file, "wb")
	pickle.dump(tgfm_results, g)
	g.close()

# Close file handles
t_pip.close()
