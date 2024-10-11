import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import pickle
import tgfm
import bootstrapped_tgfm
import rpy2
import rpy2.robjects.numpy2ri as numpy2ri
import rpy2.robjects as ro
import scipy.stats
ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
susieR_pkg = importr('susieR')


def extract_tissue_names(gtex_pseudotissue_file,ignore_tissues, remove_testis=False):
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
		if data[0] == ignore_tissues:
			continue
		arr.append(data[0])
	f.close()
	return np.asarray(arr)

def standardize_eqtl_pmces_old(eqtl_pmces, variant_ld):

	gene_variances = np.diag(np.dot(np.dot(eqtl_pmces, variant_ld), np.transpose(eqtl_pmces)))
	# compute gene variance
	n_genes = eqtl_pmces.shape[0]

	# standardize
	for gene_iter in range(n_genes):
		eqtl_pmces[gene_iter,:] = eqtl_pmces[gene_iter,:]/np.sqrt(gene_variances[gene_iter])

	return eqtl_pmces

def standardize_eqtl_pmces(eqtl_pmces, gene_variances, reference_ld):

	pmces_gene_variances = np.diag(np.dot(np.dot(eqtl_pmces, reference_ld), np.transpose(eqtl_pmces)))

	# compute gene variance
	n_genes = eqtl_pmces.shape[0]
	predictor_variances = []

	# standardize
	for gene_iter in range(n_genes):		
		expected_total_genetic_var = gene_variances[gene_iter]
		genetic_gene_pmces_var = pmces_gene_variances[gene_iter]
		genetic_gene_noise_var = expected_total_genetic_var - genetic_gene_pmces_var

		unbiased_genetic_gene_var_est = genetic_gene_pmces_var - genetic_gene_noise_var
		

		att_lambda = unbiased_genetic_gene_var_est/genetic_gene_pmces_var

		standardized_predictor_variance = expected_total_genetic_var/genetic_gene_pmces_var
		standardized_attenuated_predictor_variance = standardized_predictor_variance*np.square(att_lambda)

		eqtl_pmces[gene_iter,:] = eqtl_pmces[gene_iter,:]/np.sqrt(gene_variances[gene_iter])
		predictor_variances.append(standardized_attenuated_predictor_variance)



	return eqtl_pmces, np.asarray(predictor_variances)


def extract_full_gene_variant_ld(standardized_eqtl_effects, variant_ld):
	gene_variant_ld = np.dot(standardized_eqtl_effects,variant_ld) # Ngenes X n_variants
	expression_covariance = np.dot(gene_variant_ld, np.transpose(standardized_eqtl_effects))
	np.fill_diagonal(expression_covariance, 1.0)
	dd = np.diag(1.0/np.sqrt(np.diag(expression_covariance)))
	ge_ld = np.dot(np.dot(dd, expression_covariance),dd)
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


def extract_valid_joint_susie_components_from_full_ld(alpha_phi, beta_phi, full_ld, ld_thresh):
	num_components = alpha_phi.shape[0]
	valid_components = []

	num_genes = alpha_phi.shape[1]
	num_variants = beta_phi.shape[1]

	for component_num in range(num_components):
		cs_predictors = get_credible_set_genes(np.hstack((alpha_phi[component_num,:], beta_phi[component_num,:])), .95)
		cs_genes = cs_predictors[cs_predictors < num_genes]
		cs_variants = cs_predictors[cs_predictors >= num_genes] - num_genes

		# absolute ld among genes and variants in credible set
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

def temp_elbo_calc(z_vec, LD, samp_size, alpha, mu, mu2, KL_terms):
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

	return tgfm_obj


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

def fill_in_causal_effect_size_matrix(null_mat, bs_eqtls_pmces_sparse):
	null_mat[bs_eqtls_pmces_sparse[:,0].astype(int), bs_eqtls_pmces_sparse[:,1].astype(int)] = bs_eqtls_pmces_sparse[:,2]
	return null_mat

def extract_susie_obj_from_this_bootstrap(alpha_phis, beta_phis, bs_iter, num_genes, num_snps, num_components):
	final = np.zeros((num_components, num_genes+num_snps))
	for l_iter in range(num_components):
		final[l_iter,:] = np.hstack((alpha_phis[l_iter][bs_iter,:], beta_phis[l_iter][bs_iter,:]))
	return final

def extract_susie_kl_terms_from_this_bootstrap(KL_terms, bs_iter, num_components):
	final = []
	for l_iter in range(num_components):
		final.append(KL_terms[l_iter][bs_iter])
	return np.asarray(final)


def compute_elbo_for_bootstrapped_tgfm_obj_shell(tgfm_obj, variant_z_vec, variant_ld_mat, gwas_sample_size):
	elbos = []
	for bs_iter in range(tgfm_obj.n_bs):
		# Extract causal eqtl effects for this bootstrap
		bs_eqtls_pmces = np.zeros((tgfm_obj.G, tgfm_obj.K))
		bs_eqtls_pmces_sparse = np.load(tgfm_obj.pmces_files[bs_iter])
		bs_eqtls_pmces = fill_in_causal_effect_size_matrix(bs_eqtls_pmces, bs_eqtls_pmces_sparse)

		gene_variant_z_vec = np.hstack((tgfm_obj.nominal_twas_z[bs_iter], tgfm_obj.gwas_variant_z))

		gene_variant_full_ld = extract_full_gene_variant_ld(bs_eqtls_pmces, tgfm_data['reference_ld'])

		bs_phi = extract_susie_obj_from_this_bootstrap(tgfm_obj.alpha_phis, tgfm_obj.beta_phis, bs_iter, tgfm_obj.G, tgfm_obj.K, tgfm_obj.L)
		bs_mu = extract_susie_obj_from_this_bootstrap(tgfm_obj.alpha_mus, tgfm_obj.beta_mus, bs_iter, tgfm_obj.G, tgfm_obj.K, tgfm_obj.L)
		bs_var = extract_susie_obj_from_this_bootstrap(tgfm_obj.alpha_vars, tgfm_obj.beta_vars, bs_iter, tgfm_obj.G, tgfm_obj.K, tgfm_obj.L)
		bs_mu2 = np.square(bs_mu) + bs_var
		bs_kl = extract_susie_kl_terms_from_this_bootstrap(tgfm_obj.KL_terms, bs_iter, tgfm_obj.L)

		elbo = elbo_calc(gene_variant_z_vec, gene_variant_full_ld, gwas_sample_size, bs_phi, bs_mu, bs_mu2, bs_kl)
		elbos.append(elbo)
	return np.asarray(elbos)

def extract_valid_tgfm_sampler_components(tgfm_data, tgfm_obj, subset_n = 100, ld_thresh=0.5):
	bs_eqtls_pmces = np.zeros((tgfm_obj.G, tgfm_obj.K))
	valid_components = []
	component_marginal_variant_preds = []
	for bs_iter in range(tgfm_obj.n_bs):
		# Initialize component array for this bootstrap
		bs_valid_components = []
		bs_component_marginal_variant_preds = []

		# Extract eqtl causal effects for this bootstrap
		bs_eqtls_pmces = bs_eqtls_pmces*0.0
		bs_eqtls_pmces_sparse = tgfm_obj.sparse_sampled_gene_eqtl_pmces[bs_iter]
		bs_eqtls_pmces = fill_in_causal_effect_size_matrix(bs_eqtls_pmces, bs_eqtls_pmces_sparse)
		gene_variant_z_vec = np.hstack((tgfm_obj.nominal_twas_z[bs_iter], tgfm_obj.gwas_variant_z))
		gene_variant_full_ld = extract_full_gene_variant_ld(bs_eqtls_pmces, tgfm_data['reference_ld'])


		for l_iter in range(tgfm_obj.L):
			cs_predictors = get_credible_set_genes(np.hstack((tgfm_obj.alpha_phis[l_iter][bs_iter, :], tgfm_obj.beta_phis[l_iter][bs_iter, :])), .95)

			if subset_n > len(cs_predictors):
				# absolute ld among genes and variants in credible set
				if np.min(np.abs(gene_variant_full_ld[cs_predictors,:][:, cs_predictors])) > ld_thresh:
					bs_valid_components.append(l_iter)

					causal_effect_vec = np.hstack((tgfm_obj.alpha_phis[l_iter][bs_iter, :]*tgfm_obj.alpha_mus[l_iter][bs_iter, :], tgfm_obj.beta_phis[l_iter][bs_iter, :]*tgfm_obj.beta_mus[l_iter][bs_iter, :]))
					pred_marginal = np.dot(gene_variant_full_ld, causal_effect_vec)
					var_pred_marginal = pred_marginal[(tgfm_obj.G):]
					bs_component_marginal_variant_preds.append(var_pred_marginal)
			else:
				# First run subsetted analysis
				cs_predictors_subset = np.random.choice(cs_predictors, size=subset_n, replace=False, p=None)
				if np.min(np.abs(gene_variant_full_ld[cs_predictors_subset,:][:, cs_predictors_subset])) > ld_thresh:
					if np.min(np.abs(gene_variant_full_ld[cs_predictors,:][:, cs_predictors])) > ld_thresh:
						bs_valid_components.append(l_iter)

						causal_effect_vec = np.hstack((tgfm_obj.alpha_phis[l_iter][bs_iter, :]*tgfm_obj.alpha_mus[l_iter][bs_iter, :], tgfm_obj.beta_phis[l_iter][bs_iter, :]*tgfm_obj.beta_mus[l_iter][bs_iter, :]))
						pred_marginal = np.dot(gene_variant_full_ld, causal_effect_vec)
						var_pred_marginal = pred_marginal[(tgfm_obj.G):]
						bs_component_marginal_variant_preds.append(var_pred_marginal)
		valid_components.append(np.asarray(bs_valid_components))
		component_marginal_variant_preds.append(bs_component_marginal_variant_preds)
	return valid_components, component_marginal_variant_preds


def merge_two_bootstrapped_tgfms_based_on_elbo(tgfm_obj, tgfm_obj2, variant_z_vec, variant_ld_mat, gwas_sample_size):
	for bs_iter in range(tgfm_obj.n_bs):
		# Extract causal eqtl effects for this bootstrap
		bs_eqtls_pmces = np.zeros((tgfm_obj.G, tgfm_obj.K))
		#bs_eqtls_pmces_sparse = np.load(tgfm_obj.pmces_files[bs_iter])
		bs_eqtls_pmces_sparse = tgfm_obj.sparse_sampled_gene_eqtl_pmces[bs_iter]
		bs_eqtls_pmces = fill_in_causal_effect_size_matrix(bs_eqtls_pmces, bs_eqtls_pmces_sparse)
		gene_variant_z_vec = np.hstack((tgfm_obj.nominal_twas_z[bs_iter], tgfm_obj.gwas_variant_z))
		gene_variant_full_ld = extract_full_gene_variant_ld(bs_eqtls_pmces, tgfm_data['reference_ld'])

		# Compute ELBO for object 1
		bs_phi = extract_susie_obj_from_this_bootstrap(tgfm_obj.alpha_phis, tgfm_obj.beta_phis, bs_iter, tgfm_obj.G, tgfm_obj.K, tgfm_obj.L)
		bs_mu = extract_susie_obj_from_this_bootstrap(tgfm_obj.alpha_mus, tgfm_obj.beta_mus, bs_iter, tgfm_obj.G, tgfm_obj.K, tgfm_obj.L)
		bs_var = extract_susie_obj_from_this_bootstrap(tgfm_obj.alpha_vars, tgfm_obj.beta_vars, bs_iter, tgfm_obj.G, tgfm_obj.K, tgfm_obj.L)
		bs_mu2 = np.square(bs_mu) + bs_var
		bs_kl = extract_susie_kl_terms_from_this_bootstrap(tgfm_obj.KL_terms, bs_iter, tgfm_obj.L)
		elbo = elbo_calc(gene_variant_z_vec, gene_variant_full_ld, gwas_sample_size, bs_phi, bs_mu, bs_mu2, bs_kl)

		# Compute ELBO for object 1
		bs_phi_2 = extract_susie_obj_from_this_bootstrap(tgfm_obj2.alpha_phis, tgfm_obj2.beta_phis, bs_iter, tgfm_obj2.G, tgfm_obj2.K, tgfm_obj2.L)
		bs_mu_2 = extract_susie_obj_from_this_bootstrap(tgfm_obj2.alpha_mus, tgfm_obj2.beta_mus, bs_iter, tgfm_obj2.G, tgfm_obj2.K, tgfm_obj2.L)
		bs_var_2 = extract_susie_obj_from_this_bootstrap(tgfm_obj2.alpha_vars, tgfm_obj2.beta_vars, bs_iter, tgfm_obj2.G, tgfm_obj2.K, tgfm_obj2.L)
		bs_mu2_2 = np.square(bs_mu_2) + bs_var_2
		bs_kl_2 = extract_susie_kl_terms_from_this_bootstrap(tgfm_obj2.KL_terms, bs_iter, tgfm_obj2.L)
		elbo_2 = elbo_calc(gene_variant_z_vec, gene_variant_full_ld, gwas_sample_size, bs_phi_2, bs_mu_2, bs_mu2_2, bs_kl_2)

		# Only need to update tgfm_obj if this is the case
		if elbo_2 > elbo:
			for l_iter in range(tgfm_obj.L):
				# Update Phi
				tgfm_obj.alpha_phis[l_iter][bs_iter, :] = np.copy(tgfm_obj2.alpha_phis[l_iter][bs_iter, :])
				tgfm_obj.beta_phis[l_iter][bs_iter, :] = np.copy(tgfm_obj2.beta_phis[l_iter][bs_iter, :])
				# Update mu
				tgfm_obj.alpha_mus[l_iter][bs_iter, :] = np.copy(tgfm_obj2.alpha_mus[l_iter][bs_iter, :])
				tgfm_obj.beta_mus[l_iter][bs_iter, :] = np.copy(tgfm_obj2.beta_mus[l_iter][bs_iter, :])
				# Update mu_var
				tgfm_obj.alpha_vars[l_iter][bs_iter, :] = np.copy(tgfm_obj2.alpha_vars[l_iter][bs_iter, :])
				tgfm_obj.beta_vars[l_iter][bs_iter, :] = np.copy(tgfm_obj2.beta_vars[l_iter][bs_iter, :])
				# Update LBF
				tgfm_obj.alpha_lbfs[l_iter][bs_iter, :] = np.copy(tgfm_obj2.alpha_lbfs[l_iter][bs_iter, :])
				tgfm_obj.beta_lbfs[l_iter][bs_iter, :] = np.copy(tgfm_obj2.beta_lbfs[l_iter][bs_iter, :])
	# Recompute PIPs
	tgfm_obj.compute_pips()

	return tgfm_obj




def tgfm_inference_shell(tgfm_data, gene_log_prior, var_log_prior, init_method, bootstrap_prior):
	if init_method == 'null':
		tgfm_obj = bootstrapped_tgfm.TGFM(L=7, estimate_prior_variance=True, gene_init_log_pi=gene_log_prior, variant_init_log_pi=var_log_prior, convergence_thresh=1e-5, max_iter=5)
		tgfm_obj.fit(twas_data_obj=tgfm_data)
	elif init_method == 'best':
		# Run susie with only variants
		z_vec = tgfm_data['gwas_beta']/tgfm_data['gwas_beta_se']

		# Now run tgfm sampler with null init
		tgfm_obj = bootstrapped_tgfm.TGFM(L=10, estimate_prior_variance=True, gene_init_log_pi=gene_log_prior, variant_init_log_pi=var_log_prior, convergence_thresh=1e-5, max_iter=5, bootstrap_prior=bootstrap_prior)
		tgfm_obj.fit(twas_data_obj=tgfm_data)
		#elbo_null_init = compute_elbo_for_bootstrapped_tgfm_obj_shell(tgfm_obj, z_vec, ld_mat, tgfm_data['gwas_sample_size'])

		if np.max(np.max(tgfm_obj.expected_alpha_pips)) > .2:
			print('extra')
			# Create initialization alpha, mu, and mu_var
			susie_variant_only = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=tgfm_data['reference_ld'], n=tgfm_data['gwas_sample_size'], L=10, estimate_residual_variance=False)

			num_snps = len(z_vec)
			num_genes = len(tgfm_data['genes'])
			alpha_init = np.zeros((10, num_genes + num_snps))
			mu_init = np.zeros((10, num_genes + num_snps))
			mu_var_init = np.zeros((10, num_genes + num_snps))
			alpha_init[:, num_genes:] = susie_variant_only.rx2('alpha')
			mu_init[:, num_genes:] = susie_variant_only.rx2('mu')
			variant_only_mu_var = susie_variant_only.rx2('mu2') - np.square(susie_variant_only.rx2('mu'))
			mu_var_init[:, num_genes:] = variant_only_mu_var
			del susie_variant_only

			# Run tgfm sampler with variant init
			tgfm_obj_variant_init = bootstrapped_tgfm.TGFM(L=10, estimate_prior_variance=True, gene_init_log_pi=gene_log_prior, variant_init_log_pi=var_log_prior, convergence_thresh=1e-5, max_iter=5, bootstrap_prior=bootstrap_prior)
			tgfm_obj_variant_init.fit(twas_data_obj=tgfm_data, phi_init=alpha_init, mu_init=mu_init, mu_var_init=mu_var_init)

			tgfm_obj = merge_two_bootstrapped_tgfms_based_on_elbo(tgfm_obj, tgfm_obj_variant_init, z_vec, tgfm_data['reference_ld'], tgfm_data['gwas_sample_size'])
		
	return tgfm_obj


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

def extract_log_prior_probabilities_for_tglr_bs_nn_pmces(prior_file, variant_names, gene_names):
	element_name_to_h2 = {}
	head_count = 0
	f = open(prior_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count +1
			continue
		ele_name = data[0]
		ele_h2 = float(data[1])
		if ele_h2 == 0:
			ele_h2 = 1e-100
		element_name_to_h2[ele_name] = ele_h2
	f.close()

	snp_h2 = np.zeros(len(variant_names)) + element_name_to_h2['variant']
	gene_h2 = []
	for gene_name in gene_names:
		tissue_name = '_'.join(gene_name.split('_')[1:])
		gene_h2.append(element_name_to_h2[tissue_name])
	gene_h2 = np.asarray(gene_h2)

	# normalize to get probs
	normalizer = np.sum(snp_h2) + np.sum(gene_h2)

	variant_probs = snp_h2/normalizer
	gene_probs = gene_h2/normalizer


	return np.log(variant_probs), np.log(gene_probs)

def extract_log_prior_probabilities_for_tglr_bs_nn_sampler(prior_file, variants, genes):
	element_name_to_h2 = {}
	head_count = 0
	f = open(prior_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count +1
			continue
		ele_name = data[0]
		bs_probs = np.asarray(data[2].split(';')).astype(float)
		new_bs_probs = []
		for bs_prob in bs_probs:
			if bs_prob == 0:
				new_bs_probs.append(1e-100)
			else:
				new_bs_probs.append(bs_prob)
		new_bs_probs = np.asarray(new_bs_probs)
		element_name_to_h2[ele_name] = new_bs_probs
	f.close()

	n_var = len(variants)
	n_genes = len(genes)
	n_bs = len(element_name_to_h2['variant'])
	var_probs = np.zeros((n_var, n_bs))
	gene_probs = np.zeros((n_genes, n_bs))

	for var_iter in range(n_var):
		var_probs[var_iter, :] = element_name_to_h2['variant']
	for gene_iter, gene_name in enumerate(genes):
		tissue_name = '_'.join(gene_name.split('_')[1:])
		gene_probs[gene_iter, :] = element_name_to_h2[tissue_name]
	
	# Normalize rows
	normalizers = np.sum(gene_probs,axis=0) + np.sum(var_probs,axis=0)
	norm_var_probs = var_probs/normalizers
	norm_gene_probs = gene_probs/normalizers

	return np.log(norm_var_probs), np.log(norm_gene_probs)


def extract_log_prior_probabilities_for_iterative_prior_sampler_sampler(prior_file, variants, genes):
	element_name_to_h2 = {}
	head_count = 0
	f = open(prior_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count +1
			continue
		ele_name = data[0]
		bs_probs = np.asarray(data[3].split(';')).astype(float)
		new_bs_probs = []
		for bs_prob in bs_probs:
			if bs_prob == 0:
				new_bs_probs.append(1e-100)
			else:
				new_bs_probs.append(bs_prob)
		new_bs_probs = np.asarray(new_bs_probs)
		element_name_to_h2[ele_name] = new_bs_probs
	f.close()

	n_var = len(variants)
	n_genes = len(genes)
	n_bs = len(element_name_to_h2['variant'])
	var_probs = np.zeros((n_var, n_bs))
	gene_probs = np.zeros((n_genes, n_bs))

	for var_iter in range(n_var):
		var_probs[var_iter, :] = element_name_to_h2['variant']
	for gene_iter, gene_name in enumerate(genes):
		tissue_name = '_'.join(gene_name.split('_')[1:])
		gene_probs[gene_iter, :] = element_name_to_h2[tissue_name]
	
	# Normalize rows
	normalizers = np.sum(gene_probs,axis=0) + np.sum(var_probs,axis=0)
	norm_var_probs = var_probs/normalizers
	norm_gene_probs = gene_probs/normalizers

	return np.log(norm_var_probs), np.log(norm_gene_probs)

def extract_log_prior_probabilities_for_iterative_prior_pmces(prior_file, variants, genes):
	element_name_to_h2 = {}
	head_count = 0
	f = open(prior_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count +1
			continue
		ele_name = data[0]
		bs_prob = float(data[1])
		if bs_prob == 0:
			bs_prob = 1e-100
		element_name_to_h2[ele_name] = bs_prob
	f.close()

	n_var = len(variants)
	n_genes = len(genes)
	#n_bs = len(element_name_to_h2['variant'])
	var_probs = np.zeros(n_var)
	gene_probs = np.zeros(n_genes)

	for var_iter in range(n_var):
		var_probs[var_iter] = element_name_to_h2['variant']
	for gene_iter, gene_name in enumerate(genes):
		tissue_name = '_'.join(gene_name.split('_')[1:])
		gene_probs[gene_iter] = element_name_to_h2[tissue_name]
	
	# Normalize rows
	normalizers = np.sum(gene_probs,axis=0) + np.sum(var_probs,axis=0)
	norm_var_probs = var_probs/normalizers
	norm_gene_probs = gene_probs/normalizers

	return np.log(norm_var_probs), np.log(norm_gene_probs)


def component_in_middle_of_window(alpha_phi_vec, beta_phi_vec, middle_gene_arr, middle_variant_arr):
	middle_gene_indices_dicti = {}
	for middle_gene_index in middle_gene_arr:
		middle_gene_indices_dicti[middle_gene_index] = 1
	middle_variant_indices_dicti = {}
	for middle_variant_index in middle_variant_arr:
		middle_variant_indices_dicti[middle_variant_index] = 1

	booler = False
	# Gene wins
	if np.max(alpha_phi_vec) > np.max(beta_phi_vec):
		best_index = np.argmax(alpha_phi_vec)
		if best_index in middle_gene_indices_dicti:
			booler = True
	else:
		best_index = np.argmax(beta_phi_vec)
		if best_index in middle_variant_indices_dicti:
			booler = True
	return booler



def filter_genes_in_susie_object(susie_obj, valid_gt_indices):
	new_susie_obj = []
	for gt_index in valid_gt_indices:
		new_susie_obj.append(susie_obj[gt_index])
	return new_susie_obj



def filter_tgfm_data_structure_to_remove_specified_tissue_gene_tissue_pairs(tgfm_data, ignore_tissues):
	################
	# First get valid gt indices by removing gt-pairs corresponding to tissue 0
	valid_gt_indices = []
	gt_old_to_new = []
	counter = 0
	for ii, gt_pair in enumerate(tgfm_data['genes']):
		if '_'.join(gt_pair.split('_')[1:]) == ignore_tissues:
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
	#pdb.set_trace()
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
	tgfm_data['full_gene_variances'] = tgfm_data['full_gene_variances'][valid_gt_indices]

	#################
	# Filter 'tss' key
	tgfm_data['tss'] = tgfm_data['tss'][valid_gt_indices]


	#################
	# Filter 'gene_susie_indices' key
	new_gene_susie_indices = filter_genes_in_susie_object(tgfm_data['valid_susie_components'], valid_gt_indices)
	tgfm_data['valid_susie_components'] = new_gene_susie_indices

	return tgfm_data, False

def get_agg_non_med_probs(beta_phis):
	n_comp = len(beta_phis)
	n_bs = beta_phis[0].shape[0]

	agg_nm_probs = np.zeros((n_comp, n_bs))

	for comp_iter in range(n_comp):
		agg_nm_probs[comp_iter, :] = np.sum(beta_phis[comp_iter],axis=1)
	return agg_nm_probs





######################
# Command line args
######################
trait_name = sys.argv[1]
tgfm_input_summary_file = sys.argv[2]
tgfm_output_stem = sys.argv[3]
job_number = int(sys.argv[4])
num_jobs = int(sys.argv[5])
init_method = sys.argv[6]
est_resid_var_str = sys.argv[7]
ln_pi_method_name = sys.argv[8]
gtex_pseudotissue_file = sys.argv[9]
prior_dir= sys.argv[10]
ignore_tissues = sys.argv[11]


window_pvalue_thresh = 1e-5



if est_resid_var_str == 'False':
	est_resid_var_bool = False
elif est_resid_var_str == 'True':
	est_resid_var_bool = True 
else:
	print('assumption eroror')
	pdb.set_trace()



# Extract ordered tissue information
# Extract ordered tissue information
ordered_tissue_names = extract_tissue_names(gtex_pseudotissue_file, ignore_tissues, remove_testis=True)
tissue_to_position_mapping = {}
for i, val in enumerate(ordered_tissue_names):
	tissue_to_position_mapping[val] = i



# Append job number and num jobs to TGFM ouptut stem
new_tgfm_output_stem = tgfm_output_stem + '_' + str(job_number) + '_' + str(num_jobs)
print(new_tgfm_output_stem)
# Open PIP file handle
pip_output_file = new_tgfm_output_stem + '_tgfm_pip_summary.txt'
t_pip = open(pip_output_file,'w')
t_pip.write('window_name\tinclusion_elements\tinclusion_probabilities\n')




# Now loop through windows
# In each window run TGFM independently
# Loop through trait components
tgfm_input_data = np.loadtxt(tgfm_input_summary_file,dtype=str,delimiter='\t')
tgfm_input_data = tgfm_input_data[1:,:]

# Subset to just windows in this parallel run
tgfm_input_data_parr = np.array_split(tgfm_input_data, num_jobs)[job_number]

# Get n_windows on this run
n_windows = tgfm_input_data_parr.shape[0]



# Now loop through windows
# In each window run TGFM independently
# Loop through trait components
tgfm_input_data = np.loadtxt(tgfm_input_summary_file,dtype=str,delimiter='\t')
tgfm_input_data = tgfm_input_data[1:,:]

# Subset to just windows in this parallel run
tgfm_input_data_parr = np.array_split(tgfm_input_data, num_jobs)[job_number]

# Get n_windows on this run
n_windows = tgfm_input_data_parr.shape[0]

for window_iter in range(n_windows):
	data = tgfm_input_data_parr[window_iter, :]

	##############################
	# Extract relevent fields
	###############################
	window_name = data[0]
	print(window_name)

	ld_file = data[1]
	tgfm_input_pkl = data[2]
	tgfm_trait_input_pkl = data[3]

	##############################
	# Load in Data
	###############################
	# Load in tgfm trait input data
	g = open(tgfm_trait_input_pkl, "rb")
	tgfm_trait_data = pickle.load(g)
	g.close()

	# Extract index in tgfm_trait_data corresponding to current gwas study
	trait_index = np.where(tgfm_trait_data['gwas_study_names'] == trait_name)[0][0]

	# Extract gwas betas and standard errrors
	gwas_beta = tgfm_trait_data['gwas_beta'][trait_index,:]
	gwas_beta_se = tgfm_trait_data['gwas_beta_se'][trait_index,:]
	gwas_sample_size = int(tgfm_trait_data['gwas_sample_size'][trait_index])
	gwas_z = gwas_beta/gwas_beta_se
	gwas_p = scipy.stats.norm.sf(abs(gwas_z))*2.0

	# Ignore windows with no pvalues less than some threshold
	if np.min(gwas_p) > window_pvalue_thresh:
		print('skipped because of window pvalue threshold')
		continue
	# Hacky fix
	if trait_name == 'blood_RED_COUNT' and window_name == '10:44014743:47014743':
		continue
	if trait_name == 'blood_RED_COUNT' and window_name == '10:45014743:48014743':
		continue
	if trait_name == 'blood_RED_COUNT' and window_name == '10:46014743:49014743':
		continue


	# Load in tgfm input data
	g = open(tgfm_input_pkl, "rb")
	tgfm_data = pickle.load(g)
	g.close()

	if ignore_tissues != 'None':
		tgfm_data, filter_error_bool = filter_tgfm_data_structure_to_remove_specified_tissue_gene_tissue_pairs(tgfm_data, ignore_tissues)
		if filter_error_bool:
			print('skipped because of no genes')
			continue

	# Load in LD
	tgfm_data['reference_ld'] = np.load(ld_file)
	# Add ld to tgfm_data obj
	#tgfm_data['reference_ld'] = ld_mat

	# Skip windows with no genes
	if len(tgfm_data['genes']) == 0:
		continue

	bootstrap_prior = False
	if ln_pi_method_name == 'variant_gene':
		var_log_prior = tgfm_trait_data['variant_gene_ln_prior_variant'][trait_index,:]
		gene_log_prior = tgfm_trait_data['variant_gene_ln_prior_gene'][trait_index,:]
	elif ln_pi_method_name == 'sparse_variant_gene_tissue':
		var_log_prior = tgfm_trait_data['sparse_variant_gene_tissue_ln_prior_variant'][trait_index,:]
		gene_log_prior = tgfm_trait_data['sparse_variant_gene_tissue_ln_prior_gene'][trait_index,:]
	elif ln_pi_method_name == 'iterative_variant_gene_tissue':
		log_prior_file = tgfm_output_stem.split('_iterative_variant')[0] + '_variant_gene_iterative_emperical_distribution_prior_' + window_name + '.txt'
		var_log_prior, gene_log_prior = extract_log_prior_probabilities_from_summary_file(log_prior_file, tgfm_data['variants'], tgfm_data['genes'])
	elif ln_pi_method_name == 'uniform':
		var_log_prior = tgfm_trait_data['uniform_ln_prior_variant'][trait_index,:]
		gene_log_prior = tgfm_trait_data['uniform_ln_prior_gene'][trait_index,:]
	elif ln_pi_method_name == 'tglr_bootstrapped_nonnegative_pmces':
		tmp_sldsc_dir = '/n/scratch3/users/b/bes710/causal_eqtl_gwas/gtex/tgfm_sldsc_results/'
		prior_file = tmp_sldsc_dir + trait_name + '_baseline_no_qtl_component_gene_no_testis_pmces_gene_adj_ld_scores_nonnegative_eqtl_bootstrapped_sldsc_per_element_h2.txt'
		var_log_prior, gene_log_prior = extract_log_prior_probabilities_for_tglr_bs_nn_pmces(prior_file, tgfm_data['variants'], tgfm_data['genes'])
	elif ln_pi_method_name == 'tglr_bootstrapped_nonnegative_sampler':
		tmp_sldsc_dir = '/n/scratch3/users/b/bes710/causal_eqtl_gwas/gtex/tgfm_sldsc_results/'
		prior_file = tmp_sldsc_dir + trait_name + '_baseline_no_qtl_component_gene_no_testis_pmces_gene_adj_ld_scores_nonnegative_eqtl_bootstrapped_sldsc_per_element_h2.txt'
		var_log_prior, gene_log_prior = extract_log_prior_probabilities_for_tglr_bs_nn_sampler(prior_file, tgfm_data['variants'], tgfm_data['genes'])
		bootstrap_prior = True
	elif ln_pi_method_name == 'uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler':
		log_prior_file = tgfm_output_stem.split('susie')[0] + 'susie_pmces_uniform_iterative_variant_gene_prior_v2_pip_level_bootstrapped.txt'
		log_prior_file = prior_dir + log_prior_file.split('/')[-1]
		#log_prior_file = prior_dir + 'tgfm_results_' + trait_name + '_'
		var_log_prior, gene_log_prior = extract_log_prior_probabilities_for_iterative_prior_sampler_sampler(log_prior_file, tgfm_data['variants'], tgfm_data['genes'])
		bootstrap_prior = True
	elif ln_pi_method_name == 'uniform_pmces_iterative_variant_gene_tissue_pip_level_pmces':
		log_prior_file = tgfm_output_stem.split('_component_gene')[0] + '_component_gene_susie_pmces_uniform_iterative_variant_gene_prior_v2_pip_level_bootstrapped.txt'
		var_log_prior, gene_log_prior = extract_log_prior_probabilities_for_iterative_prior_pmces(log_prior_file, tgfm_data['variants'], tgfm_data['genes'])
	else:
		print('assumption erororo: ' + str(ln_pi_method_name) + ' is not yet implemented')
		pdb.set_trace()


	# Standardize gwas summary statistics
	gwas_beta_scaled, gwas_beta_se_scaled = convert_to_standardized_summary_statistics(gwas_beta, gwas_beta_se, float(gwas_sample_size), tgfm_data['reference_ld'])

	# Extract full ld between genes, variants, and gene-variants
	#gene_variant_full_ld = extract_full_gene_variant_ld(tgfm_data['gene_eqtl_pmces'], tgfm_data['reference_ld'])

	# Add extra info to tgfm_data
	tgfm_data['gwas_beta'] = gwas_beta_scaled
	tgfm_data['gwas_beta_se'] = gwas_beta_se_scaled
	tgfm_data['gwas_sample_size'] = gwas_sample_size

	##############################
	# Run TGFM
	###############################
	tgfm_obj = tgfm_inference_shell(tgfm_data, gene_log_prior, var_log_prior, init_method, bootstrap_prior)

	print('organizing')

	##############################
	# Organize TGFM data and print to results
	###############################
	# Extract names of genetic elements
	genetic_element_names = np.hstack((tgfm_data['genes'], tgfm_data['variants']))
	# Extract dictionary list of genetic elements in the middel of this window
	middle_genetic_elements = extract_middle_genetic_elements(tgfm_data['genes'], tgfm_data['middle_gene_indices'], tgfm_data['variants'], tgfm_data['middle_variant_indices'])
	# Extract genetic element pips
	genetic_element_pips = np.hstack((tgfm_obj.expected_alpha_pips, tgfm_obj.expected_beta_pips))

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

	# Note could compute expected mediated probability

	# Write to credible set output
	t_pip.write(window_name + '\t')
	t_pip.write(';'.join(ordered_middle_names) + '\t')
	t_pip.write(';'.join(ordered_middle_pips.astype(str)) + '\n')

	# Save all TGFM results to pkl
	'''
	tgfm_results = {}
	tgfm_results['variants'] = tgfm_data['variants']
	tgfm_results['genes'] = tgfm_data['genes']
	tgfm_results['alpha_phis'] = tgfm_obj.alpha_phis
	tgfm_results['beta_phis'] = tgfm_obj.beta_phis
	tgfm_results['alpha_lbfs'] = tgfm_obj.alpha_lbfs
	tgfm_results['beta_lbfs'] = tgfm_obj.beta_lbfs
	tgfm_results['alpha_pips'] = tgfm_obj.alpha_pips
	tgfm_results['beta_pips'] = tgfm_obj.beta_pips
	tgfm_results['expected_alpha_pips'] = tgfm_obj.expected_alpha_pips
	tgfm_results['expected_beta_pips'] = tgfm_obj.expected_beta_pips
	tgfm_results['valid_components'] = valid_tgfm_sampler_components
	tgfm_results['nominal_twas_z'] = tgfm_obj.nominal_twas_z
	'''
	agg_nm_probs = get_agg_non_med_probs(tgfm_obj.beta_phis)

	tgfm_results = {}
	tgfm_results['variants'] = tgfm_data['variants']
	tgfm_results['genes'] = tgfm_data['genes']
	tgfm_results['alpha_phis'] = tgfm_obj.alpha_phis
	tgfm_results['alpha_mus'] = tgfm_obj.alpha_mus
	tgfm_results['alpha_vars'] = tgfm_obj.alpha_vars
	tgfm_results['beta_phis'] = tgfm_obj.beta_phis
	#tgfm_results['alpha_lbfs'] = tgfm_obj.alpha_lbfs
	#tgfm_results['beta_lbfs'] = tgfm_obj.beta_lbfs
	tgfm_results['alpha_pips'] = tgfm_obj.alpha_pips
	#tgfm_results['beta_pips'] = tgfm_obj.beta_pips
	tgfm_results['expected_alpha_pips'] = tgfm_obj.expected_alpha_pips
	tgfm_results['expected_beta_pips'] = tgfm_obj.expected_beta_pips
	#tgfm_results['valid_components'] = valid_tgfm_sampler_components
	#tgfm_results['valid_middle_components'] = valid_middle_tgfm_components
	tgfm_results['nominal_twas_z'] = tgfm_obj.nominal_twas_z
	tgfm_results['middle_variant_indices'] = tgfm_data['middle_variant_indices']
	tgfm_results['middle_gene_indices'] = tgfm_data['middle_gene_indices']
	tgfm_results['aggregated_nm_probs'] = agg_nm_probs
	

	# Write pickle file
	window_tgfm_output_file = tgfm_output_stem + '_' + window_name + '_results.pkl'
	g = open(window_tgfm_output_file, "wb")
	pickle.dump(tgfm_results, g)
	g.close()



# Close file handles
t_pip.close()

