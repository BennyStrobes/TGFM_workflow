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
import susie_inf
ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
susieR_pkg = importr('susieR')


def standardize_eqtl_pmces_old(eqtl_pmces, variant_ld):

	#gene_variances = np.diag(np.dot(np.dot(eqtl_pmces, variant_ld), np.transpose(eqtl_pmces)))
	gene_variances = np.sum(np.dot(eqtl_pmces, variant_ld)*eqtl_pmces, axis=1)

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

	tgfm_obj.elbo = susie_obj.rx2('elbo')[-1]

	return tgfm_obj


def update_tgfm_obj_with_susie_inf_res_obj(tgfm_obj, inf_model):
	# Alphas
	susie_alpha = np.transpose(inf_model['PIP'])
	tgfm_obj.alpha_phi = susie_alpha[:,:(tgfm_obj.G)]
	tgfm_obj.beta_phi = susie_alpha[:,(tgfm_obj.G):]

	# Mus
	susie_mu = np.transpose(inf_model['mu'])
	tgfm_obj.alpha_mu = susie_mu[:,:(tgfm_obj.G)]
	tgfm_obj.beta_mu = susie_mu[:,(tgfm_obj.G):]

	return tgfm_obj

def calculate_emperical_null(susie_variant_obj, tgfm_data, gene_variant_full_ld, null_perm_output, variant_z, new_gene_z):
	n_genes = len(tgfm_data['genes'])
	gene_names = tgfm_data['genes']
	n_components = susie_variant_obj.rx2('alpha').shape[0]
	t = open(null_perm_output,'w')
	t.write('component_number\tgene_name\tprob\temperical_p\n')

	for component_iter in range(n_components):
		# Ignore components where no highly probable genes
		if np.max(susie_variant_obj.rx2('alpha')[component_iter,:n_genes]) < .5:
			continue
		gene_index = np.argmax(susie_variant_obj.rx2('alpha')[component_iter,:n_genes])
		gene_prob = susie_variant_obj.rx2('alpha')[component_iter, gene_index]
		gene_name = gene_names[gene_index]

		n_perm = 50000
		n_var = tgfm_data['gene_eqtl_pmces'].shape[1]
		n_genes = tgfm_data['gene_eqtl_pmces'].shape[0]
		n_eqtl_per_gene = 5
		# Simulate unstandardized causal eqtl effects
		perm_causal_effects_unstandardized = np.zeros((n_perm, n_var))
		all_indices = np.arange(n_var)
		for perm_iter in range(n_perm):
			causal_indices = np.random.choice(all_indices, size=n_eqtl_per_gene, replace=False)
			perm_causal_effects_unstandardized[perm_iter, causal_indices] = np.random.normal(loc=0.0, scale=1.0,size=n_eqtl_per_gene)
		# simulate standardized causal eqtl effects
		variant_ld = gene_variant_full_ld[n_genes:,n_genes:]
		perm_causal_effects_standardized = standardize_eqtl_pmces_old(perm_causal_effects_unstandardized, variant_ld)

		perm_z = np.dot(perm_causal_effects_standardized, variant_z)

		variant_perm_ld = np.transpose(np.dot(perm_causal_effects_standardized, variant_ld))
		gene_perm_ld = np.transpose(np.dot(np.dot(perm_causal_effects_standardized, variant_ld), np.transpose(tgfm_data['gene_eqtl_pmces'])))

		element_perm_ld = np.vstack((gene_perm_ld, variant_perm_ld))

		standard_error = 1.0/np.sqrt(tgfm_data['gwas_sample_size'])

		full_z = np.hstack((new_gene_z, variant_z))
		resid_z = np.copy(full_z)

		resid_perm_z = np.copy(perm_z)


		for component_iter_prime in range(n_components):
			if component_iter_prime == component_iter:
				continue
			component_pred = susie_variant_obj.rx2('alpha')[component_iter_prime,:]*susie_variant_obj.rx2('mu')[component_iter_prime,:]
			resid_z = resid_z - (np.dot(gene_variant_full_ld, component_pred)/standard_error)

			resid_perm_z = resid_perm_z - (np.dot(component_pred, element_perm_ld)/standard_error)


		obs_z = np.abs(resid_z[gene_index])
		pvalue = np.sum(np.abs(resid_perm_z) > obs_z)/len(resid_perm_z)

		t.write(str(component_iter) + '\t' + gene_name + '\t' + str(gene_prob) + '\t' + str(pvalue) + '\n')
	t.close()
	return

def extract_causal_variants_and_genes_in_window(trait_dir, window_stem, window_variants, window_genes):
	causal_variants = {}
	causal_genes = {}
	cv_file = trait_dir + window_stem + '_non_mediated_variant_causal_effect_sizes.txt'
	f = open(cv_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		var_id = data[0]
		effect_size = data[1]
		if effect_size == '0.0':
			continue
		causal_variants[var_id] = float(effect_size)
	f.close()
	cg_file = trait_dir + window_stem + '_expression_mediated_gene_causal_effect_sizes.txt'
	f = open(cg_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		gene_name = data[0]
		effect_sizes = np.asarray(data[1:]).astype(float)
		n_tiss = len(effect_sizes)
		for tiss_iter in range(n_tiss):
			if effect_sizes[tiss_iter] == 0.0:
				continue
			full_gene_name = gene_name + '_tissue' + str(tiss_iter)
			causal_genes[full_gene_name] = effect_sizes[tiss_iter]
	f.close()

	window_causal_gene_names = []
	window_causal_gene_indices = []
	window_causal_variant_names = []
	window_causal_variant_indices = []
	counter = 0
	for window_gene in window_genes:
		if window_gene in causal_genes:
			window_causal_gene_names.append(window_gene)
			window_causal_gene_indices.append(counter)
		counter = counter + 1
	for window_variant in window_variants:
		if window_variant in causal_variants:
			window_causal_variant_names.append(window_variant)
			window_causal_variant_indices.append(counter)
		counter = counter + 1	

	return np.asarray(window_causal_gene_names), np.asarray(window_causal_gene_indices), np.asarray(window_causal_variant_names), np.asarray(window_causal_variant_indices)	



def tgfm_inference_shell(tgfm_data, gene_log_prior, var_log_prior, gene_variant_full_ld, init_method, est_resid_var_bool, inference_output_stem):
	# Hacky: Initialize old TGFM object using only one iter of optimization
	tgfm_obj = tgfm.TGFM(L=20, estimate_prior_variance=False, gene_init_log_pi=gene_log_prior, variant_init_log_pi=var_log_prior, convergence_thresh=1e-5, max_iter=1)
	tgfm_obj.fit(twas_data_obj=tgfm_data)
	tgfm_obj2 = tgfm.TGFM(L=20, estimate_prior_variance=False, gene_init_log_pi=gene_log_prior, variant_init_log_pi=var_log_prior, convergence_thresh=1e-5, max_iter=1)
	tgfm_obj2.fit(twas_data_obj=tgfm_data)
	# More hack: need to redo twas z
	variant_z = tgfm_data['gwas_beta']/tgfm_data['gwas_beta_se']
	new_gene_z = np.dot(tgfm_data['gene_eqtl_pmces'], variant_z)

	#tgfm_obj_temp = tgfm.TGFM(L=20, estimate_prior_variance=True, gene_init_log_pi=gene_log_prior, variant_init_log_pi=var_log_prior, convergence_thresh=1e-5, max_iter=10)
	#tgfm_obj_temp.fit(twas_data_obj=tgfm_data)

	tgfm_obj.nominal_twas_z = new_gene_z

	# Create vector of concatenated z-scores
	z_vec = np.hstack((new_gene_z,variant_z))

	# Create concatenated vector of prior probs
	prior_probs = np.hstack((np.exp(gene_log_prior), np.exp(var_log_prior)))

	if init_method == 'best':
		rpy2.robjects.r['options'](warn=1)
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
		#susie_variant_init2 = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, s_init=init_obj2, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool, refine=True)

		# Run susie with null initialization
		susie_null_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool)
		#susie_null_init2 = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool, refine=True)

		# Select model with largest elbo
		susie_null_init_elbo = susie_null_init.rx2('elbo')[-1]
		susie_variant_init_elbo = susie_variant_init.rx2('elbo')[-1]


		if susie_variant_init_elbo > susie_null_init_elbo:
			# Variant init wins
			tgfm_obj_final = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_variant_init)
			tgfm_obj_variant = update_tgfm_obj_with_susie_res_obj(tgfm_obj2, susie_variant_only)
		else:
			# Null init wins
			tgfm_obj_final = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_null_init)
			tgfm_obj_variant = update_tgfm_obj_with_susie_res_obj(tgfm_obj2, susie_variant_only)

	elif init_method == 'debug':
		# Extract causal variants in window
		window_stem = inference_output_stem.split('/')[-1].split('_eqtl_ss_inf')[0]
		trait_dir = '/n/scratch3/users/b/bes710/causal_eqtl_gwas/tgfm_inf_simulation/simulated_trait/'
		causal_genes_in_window, causal_gene_indices_in_window, causal_variants_in_window, causal_variant_indices_in_window = extract_causal_variants_and_genes_in_window(trait_dir, window_stem, tgfm_data['variants'], tgfm_data['genes'])
		print('hello')
		# Run 'best' version of tgfm
		rpy2.robjects.r['options'](warn=1)
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
		#susie_variant_init2 = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, s_init=init_obj2, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool, refine=True)

		# Run susie with null initialization
		susie_null_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool)
		#susie_null_init2 = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool, refine=True)

		# Select model with largest elbo
		susie_null_init_elbo = susie_null_init.rx2('elbo')[-1]
		susie_variant_init_elbo = susie_variant_init.rx2('elbo')[-1]


		pdb.set_trace()



	elif init_method == 'null_perm':
		null_perm_output = inference_output_stem + 'null_perm_window_res.txt'
		rpy2.robjects.r['options'](warn=1)
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
		#susie_variant_init2 = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, s_init=init_obj2, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool, refine=True)

		# Run susie with null initialization
		susie_null_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool)
		#susie_null_init2 = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool, refine=True)

		# Select model with largest elbo
		susie_null_init_elbo = susie_null_init.rx2('elbo')[-1]
		susie_variant_init_elbo = susie_variant_init.rx2('elbo')[-1]

		if susie_variant_init_elbo > susie_null_init_elbo:
			# Variant init wins
			tgfm_obj_final = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_variant_init)
			tgfm_obj_variant = update_tgfm_obj_with_susie_res_obj(tgfm_obj2, susie_variant_only)
			calculate_emperical_null(susie_variant_init, tgfm_data, gene_variant_full_ld, null_perm_output, variant_z, new_gene_z)
		else:
			# Null init wins
			tgfm_obj_final = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_null_init)
			tgfm_obj_variant = update_tgfm_obj_with_susie_res_obj(tgfm_obj2, susie_variant_only)
			calculate_emperical_null(susie_null_init, tgfm_data, gene_variant_full_ld, null_perm_output, variant_z, new_gene_z)

	elif init_method == 'susie_inf':
		rpy2.robjects.r['options'](warn=1)
		# Run susie with only variants
		p_var_only = np.ones(len(z_vec))
		p_var_only[:len(tgfm_obj.nominal_twas_z)] = 0.0
		p_var_only = p_var_only/np.sum(p_var_only)
		susie_variant_only = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=p_var_only.reshape((len(p_var_only),1)), estimate_residual_variance=est_resid_var_bool)
	
		inf_model = susie_inf.susie(z_vec, 1.0, n=tgfm_data['gwas_sample_size'], L=20, LD=gene_variant_full_ld)

		tgfm_obj_variant = update_tgfm_obj_with_susie_res_obj(tgfm_obj2, susie_variant_only)
		tgfm_obj_final = update_tgfm_obj_with_susie_inf_res_obj(tgfm_obj, inf_model)

	elif init_method == 'refine_best':
		rpy2.robjects.r['options'](warn=1)
		# Run susie with only variants
		p_var_only = np.ones(len(z_vec))
		p_var_only[:len(tgfm_obj.nominal_twas_z)] = 0.0
		p_var_only = p_var_only/np.sum(p_var_only)
		susie_variant_only = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=p_var_only.reshape((len(p_var_only),1)), estimate_residual_variance=est_resid_var_bool)
	

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
	elif init_method == 'finemap':
		# Prepare input data for FINEMAP
		tmp_output_filestem = inference_output_stem + '_finemap_tmper_'
		tmp_ld_file = tmp_output_filestem + '.ld'
		np.savetxt(tmp_ld_file, gene_variant_full_ld, delimiter=" ")
		tmp_z_file = tmp_output_filestem + '.z'
		t = open(tmp_z_file,'w')
		t.write('rsid chromosome position allele1 allele2 maf beta se\n')
		se = 1.0/np.sqrt(100000.0)
		for var_iter in range(len(z_vec)):
			zz = z_vec[var_iter]
			beta = zz*se
			t.write('rs' + str(var_iter+1) + ' 1 1 A T .3 ' + str(beta) + ' ' + str(se) + '\n')
		t.close()

		tmp_k_file = tmp_output_filestem + '.k'
		probs = np.ones(10)/10.0
		t = open(tmp_k_file,'w')
		t.write(' '.join(probs.astype(str)) + '\n')
		t.close()

		'''
		tmp_k_file = tmp_output_filestem + '.k'
		probs = np.zeros(10) + 1e-5
		probs[1] = probs[1] + 0.155
		probs[2] = probs[2] +0.397
		probs[3] = probs[3] + 0.445
		probs[4] = probs[4] + 0.00308
		probs = probs/np.sum(probs)
		t = open(tmp_k_file,'w')
		t.write(' '.join(probs.astype(str)) + '\n')
		t.close()

		'''
		window_stem = inference_output_stem.split('/')[-1].split('_eqtl_ss_inf')[0]
		trait_dir = '/n/scratch3/users/b/bes710/causal_eqtl_gwas/tgfm_inf_simulation/simulated_trait/'
		causal_genes_in_window, causal_gene_indices_in_window, causal_variants_in_window, causal_variant_indices_in_window = extract_causal_variants_and_genes_in_window(trait_dir, window_stem, tgfm_data['variants'], tgfm_data['genes'])

		tmp_master_file = tmp_output_filestem + '.master'
		t = open(tmp_master_file,'w')
		t.write('z;ld;k;snp;config;cred;log;n_samples\n')
		t.write(tmp_z_file + ';' + tmp_ld_file + ';' + tmp_k_file + ';' + tmp_output_filestem + '.snp;' + tmp_output_filestem + '.config;' + tmp_output_filestem + '.cred;' + tmp_output_filestem + '.log;' + '100000\n')
		t.close()

		finemap_executable='/n/groups/price/ben/tools/finemap/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64'

		pdb.set_trace()
		#os.system(finemap_executable + ' --in-files ' + tmp_master_file + ' --sss --n-causal-snps 10')
		os.system(finemap_executable + ' --in-files ' + tmp_master_file + ' --prior-k --sss --n-conv-sss 1000')
		pdb.set_trace()

	else:
		print('assumption errror: susie initialization method ' + init_method + ' not recognized')
		pdb.set_trace()
	return tgfm_obj_final, tgfm_obj_variant

######################
# Command line args
######################
tgfm_input_file = sys.argv[1]
tgfm_output_stem = sys.argv[2]
ln_pi_method_name = sys.argv[3]
init_method = sys.argv[4]
est_resid_var_str = sys.argv[5]

if est_resid_var_str == 'False':
	est_resid_var_bool = False
elif est_resid_var_str == 'True':
	est_resid_var_bool = True 
else:
	print('assumption eroror')
	pdb.set_trace()



# Extract ordered tissue information
tissue_to_position_mapping = {}
ordered_tissue_names = []
for i in range(10):
	tissue_to_position_mapping['tissue' + str(i)] = i
	ordered_tissue_names.append('tissue' + str(i))
ordered_tissue_names = np.asarray(ordered_tissue_names)


# Open cs output file handle
component_cs_output_file = tgfm_output_stem + '_tgfm_component_cs_summary.txt'
t_cs = open(component_cs_output_file,'w')
t_cs.write('window_name\tcomponent_index\tgene_mediated_probability\tinclusion_elements\tinclusion_probabilities\n')
# Open tissue prob output file handle
component_tissue_output_file = tgfm_output_stem + '_tgfm_component_tissue_prob_summary.txt'
t_tiss = open(component_tissue_output_file,'w')
t_tiss.write('window_name\tcomponent_index\tgene_mediated_probability')
for tissue_name in ordered_tissue_names:
	t_tiss.write('\t' + tissue_name)
t_tiss.write('\n')
# Open gene prob output file handle
component_gene_output_file = tgfm_output_stem + '_tgfm_component_gene_prob_summary.txt'
t_gene = open(component_gene_output_file,'w')
t_gene.write('window_name\tcomponent_index\tgene_names\tn_components_per_gene\tgene_mediated_probability\n')



# Now loop through windows
# In each window run TGFM independently
f = open(tgfm_input_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue

	##############################
	# Extract relevent fields
	###############################
	window_name = data[0]
	if window_name != 'chr1_149734349_152734349':
		continue

	ld_file = data[1]
	tgfm_input_pkl = data[2]
	log_prior_prob_file = data[3] + '_' + ln_pi_method_name + '.txt'

	##############################
	# Load in Data
	###############################
	# Load in LD
	ld_mat = np.load(ld_file)
	# Load in tgfm input data
	g = open(tgfm_input_pkl, "rb")
	tgfm_data = pickle.load(g)
	g.close()
	# Add ld to tgfm_data obj
	tgfm_data['reference_ld'] = ld_mat

	# Skip windows with no genes
	if len(tgfm_data['genes']) == 0:
		continue

	# Load in log_priors
	var_log_prior, gene_log_prior = load_in_log_priors(log_prior_prob_file, tgfm_data['variants'], tgfm_data['genes'])

	# Standardize eqtl PMCES
	tgfm_data['gene_eqtl_pmces'] = standardize_eqtl_pmces(tgfm_data['gene_eqtl_pmces'], tgfm_data['gene_variances'], tgfm_data['reference_ld'])
	#tgfm_data['gene_eqtl_pmces'] = standardize_eqtl_pmces_old(tgfm_data['gene_eqtl_pmces'], tgfm_data['reference_ld'])

	# Extract full ld between genes, variants, and gene-variants
	gene_variant_full_ld = extract_full_gene_variant_ld(tgfm_data['gene_eqtl_pmces'], tgfm_data['reference_ld'])

	##############################
	# Run TGFM
	###############################
	inference_output_stem = tgfm_output_stem + '_' + window_name + '_'
	tgfm_obj, tgfm_obj_variant_only = tgfm_inference_shell(tgfm_data, gene_log_prior, var_log_prior, gene_variant_full_ld, init_method, est_resid_var_bool, inference_output_stem)

	##############################
	# Organize TGFM data and print to results
	###############################
	# Extract components that pass purity filter
	valid_tgfm_components = extract_valid_joint_susie_components_from_full_ld(tgfm_obj.alpha_phi, tgfm_obj.beta_phi, gene_variant_full_ld, .5)
	# Extract names of genetic elements
	genetic_element_names = np.hstack((tgfm_data['genes'], tgfm_data['variants']))
	# Extract dictionary list of genetic elements in the middel of this window
	middle_genetic_elements = extract_middle_genetic_elements(tgfm_data['genes'], tgfm_data['middle_gene_indices'], tgfm_data['variants'], tgfm_data['middle_variant_indices'])

	#run_susie_debug(tgfm_obj, gene_variant_full_ld,tgfm_data)

	# loop through TGFM components for this window
	for tgfm_component in valid_tgfm_components:
		# Get probability component is mediated by gene expression in any cis tissue, gene
		mediated_probability = np.sum(tgfm_obj.alpha_phi[tgfm_component,:])
		# Get probability coming from each tissue
		tissue_mediated_probabilities = get_probability_coming_from_each_tissue(ordered_tissue_names, tissue_to_position_mapping, tgfm_data['genes'], tgfm_obj.alpha_phi[tgfm_component,:])
		# Get probability of each element (concatenated across genes and variants)
		element_probabilities = np.hstack((tgfm_obj.alpha_phi[tgfm_component,:], tgfm_obj.beta_phi[tgfm_component,:]))
		# Get indices of each element included in 95% cs
		cs_element_indices = get_credible_set_genes(element_probabilities, .95)
		# Get probabilities of each element in 95% cs
		cs_element_prob = element_probabilities[cs_element_indices]
		# Get cs genetic element names
		cs_element_names = genetic_element_names[cs_element_indices]
		# Get top element
		top_element_name = cs_element_names[0]
		# Ignore components for this window not in middle
		if top_element_name not in middle_genetic_elements:
			continue
		# Write to credible set output
		t_cs.write(window_name + '\t' + str(tgfm_component) + '\t' + str(mediated_probability) + '\t')
		t_cs.write(';'.join(cs_element_names) + '\t')
		t_cs.write(';'.join(cs_element_prob.astype(str)) + '\n')
		# Write to tissue output
		t_tiss.write(window_name + '\t' + str(tgfm_component) + '\t' + str(mediated_probability) + '\t')
		t_tiss.write('\t'.join(tissue_mediated_probabilities.astype(str)) + '\n')
		# Write to gene output
		t_gene.write(window_name + '\t' + str(tgfm_component) + '\t' + ';'.join(tgfm_data['genes']) + '\t' + 'NaN' + '\t' + ';'.join(tgfm_obj.alpha_phi[tgfm_component,:].astype(str)) + '\n')
		t_cs.flush()
		t_tiss.flush()
		t_gene.flush()
	# Save all TGFM results to pkl
	tgfm_results = {}
	tgfm_results['variants'] = tgfm_data['variants']
	tgfm_results['genes'] = tgfm_data['genes']
	tgfm_results['alpha_phi'] = tgfm_obj.alpha_phi
	tgfm_results['beta_phi'] = tgfm_obj.beta_phi
	tgfm_results['alpha_mu'] = tgfm_obj.alpha_mu
	tgfm_results['beta_mu'] = tgfm_obj.beta_mu
	tgfm_results['alpha_var'] = tgfm_obj.alpha_var
	tgfm_results['beta_var'] = tgfm_obj.beta_var
	tgfm_results['component_variances'] = tgfm_obj.component_variances
	tgfm_results['elbo'] = tgfm_obj.elbo

	# Write pickle file
	window_tgfm_output_file = tgfm_output_stem + '_' + window_name + '_results.pkl'
	g = open(window_tgfm_output_file, "wb")
	pickle.dump(tgfm_results, g)
	g.close()

	# Save all TGFM results to pkl
	tgfm_results2 = {}
	tgfm_results2['variants'] = tgfm_data['variants']
	tgfm_results2['genes'] = tgfm_data['genes']
	tgfm_results2['alpha_phi'] = tgfm_obj_variant_only.alpha_phi
	tgfm_results2['beta_phi'] = tgfm_obj_variant_only.beta_phi
	tgfm_results2['alpha_mu'] = tgfm_obj_variant_only.alpha_mu
	tgfm_results2['beta_mu'] = tgfm_obj_variant_only.beta_mu
	tgfm_results2['alpha_var'] = tgfm_obj_variant_only.alpha_var
	tgfm_results2['beta_var'] = tgfm_obj_variant_only.beta_var
	tgfm_results2['component_variances'] = tgfm_obj_variant_only.component_variances
	tgfm_results2['elbo'] = tgfm_obj_variant_only.elbo

	# Write pickle file
	window_tgfm_output_file = tgfm_output_stem + '_' + window_name + '_variant_only_results.pkl'
	g = open(window_tgfm_output_file, "wb")
	pickle.dump(tgfm_results2, g)
	g.close()

f.close()

# Close file handles
t_cs.close()
t_gene.close()
t_tiss.close()
