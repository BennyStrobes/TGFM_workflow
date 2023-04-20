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
import susie_alt
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

	return tgfm_obj

def tgfm_inference_shell(tgfm_data, gene_log_prior, var_log_prior, gene_variant_full_ld, init_method, est_resid_var_bool, tmp_output_filestem):
	# Hacky: Initialize old TGFM object using only one iter of optimization
	tgfm_obj = tgfm.TGFM(L=20, estimate_prior_variance=False, gene_init_log_pi=gene_log_prior, variant_init_log_pi=var_log_prior, convergence_thresh=1e-5, max_iter=1)
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



	#susie_null_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool)

	#elbo = temp_elbo_calc(z_vec, gene_variant_full_ld, tgfm_data['gwas_sample_size'], susie_null_init.rx2('alpha'), susie_null_init.rx2('mu'), susie_null_init.rx2('mu2'), susie_null_init.rx2('KL'))

	if init_method == 'best':
		rpy2.robjects.r['options'](warn=1)
		# Run susie with only variants
		p_var_only = np.ones(len(z_vec))
		p_var_only[:len(tgfm_obj.nominal_twas_z)] = 0.0
		p_var_only = p_var_only/np.sum(p_var_only)
		susie_variant_only = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=p_var_only.reshape((len(p_var_only),1)), estimate_residual_variance=est_resid_var_bool)
	
		#elbo = temp_elbo_calc(z_vec, LD, tgfm_data['gwas_sample_size'], susie_variant_only.rx2('alpha'), susie_variant_only.rx2('mu'), susie_variant_only.rx2('mu2'), susie_variant_only.rx2('KL'))


		# Run susie with variant initialization
		init_obj = {'alpha':susie_variant_only.rx2('alpha'), 'mu':susie_variant_only.rx2('mu'),'mu2':susie_variant_only.rx2('mu2')}
		init_obj2 = ro.ListVector(init_obj)
		init_obj2.rclass = rpy2.robjects.StrVector(("list", "susie"))
		susie_variant_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, s_init=init_obj2, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool)
		#elbo1 = elbo_calc(z_vec, gene_variant_full_ld, tgfm_data['gwas_sample_size'], susie_variant_init.rx2('alpha'), susie_variant_init.rx2('mu'), susie_variant_init.rx2('mu2'), susie_variant_init.rx2('KL'))

		#susie_variant_init2 = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, s_init=init_obj2, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool, refine=True)

		# Run susie with null initialization
		susie_null_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool)
		#susie_null_init2 = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=est_resid_var_bool, refine=True)
		#elbo2 = elbo_calc(z_vec, gene_variant_full_ld, tgfm_data['gwas_sample_size'], susie_null_init.rx2('alpha'), susie_null_init.rx2('mu'), susie_null_init.rx2('mu2'), susie_null_init.rx2('KL'))

		# Select model with largest elbo
		susie_null_init_elbo = susie_null_init.rx2('elbo')[-1]
		susie_variant_init_elbo = susie_variant_init.rx2('elbo')[-1]

		if susie_variant_init_elbo > susie_null_init_elbo:
			# Variant init wins
			tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_variant_init)
		else:
			# Null init wins
			tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_null_init)
	elif init_method == 'null_perm':
		n_perm = 100000
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
			tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_variant_init)
		else:
			# Null init wins
			tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_null_init)


		pdb.set_trace()
	elif init_method == 'finemap':
		# Prepare input data for FINEMAP
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


		tmp_master_file = tmp_output_filestem + '.master'
		t = open(tmp_master_file,'w')
		t.write('z;ld;k;snp;config;cred;log;n_samples\n')
		t.write(tmp_z_file + ';' + tmp_ld_file + ';' + tmp_k_file + ';' + tmp_output_filestem + '.snp;' + tmp_output_filestem + '.config;' + tmp_output_filestem + '.cred;' + tmp_output_filestem + '.log;' + '100000\n')
		t.close()

		finemap_executable='/n/groups/price/ben/tools/finemap/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64'

		pdb.set_trace()
		#os.system(finemap_executable + ' --in-files ' + tmp_master_file + ' --sss --n-causal-snps 10')
		os.system(finemap_executable + ' --in-files ' + tmp_master_file + ' --prior-k --sss --prob-conv-sss-tol 0.000000000000000001')
		pdb.set_trace()

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
'''
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
'''

# Open PIP file handle
pip_output_file = tgfm_output_stem + '_tgfm_pip_summary.txt'
t_pip = open(pip_output_file,'w')
t_pip.write('window_name\tinclusion_elements\tinclusion_probabilities\n')




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
	print(window_name)

	ld_file = data[1]
	tgfm_input_pkl = data[2]
	log_prior_prob_file = data[3] + '_' + ln_pi_method_name + '.txt'

	##############################
	# Load in Data
	###############################
	'''
	# Load in LD
	ld_mat = np.load(ld_file)
	'''
	# Load in tgfm input data
	g = open(tgfm_input_pkl, "rb")
	tgfm_data = pickle.load(g)
	g.close()
	'''
	# Add ld to tgfm_data obj
	tgfm_data['reference_ld'] = ld_mat
	'''

	# Skip windows with no genes
	if len(tgfm_data['genes']) == 0:
		continue
	'''

	# Load in log_priors
	var_log_prior, gene_log_prior = load_in_log_priors(log_prior_prob_file, tgfm_data['variants'], tgfm_data['genes'])
	'''

	# Standardize eqtl PMCES
	#tgfm_data['gene_eqtl_pmces'] = standardize_eqtl_pmces(tgfm_data['gene_eqtl_pmces'], tgfm_data['gene_variances'], tgfm_data['reference_ld'])
	'''
	tgfm_data['gene_eqtl_pmces'] = standardize_eqtl_pmces_old(tgfm_data['gene_eqtl_pmces'], tgfm_data['reference_ld'])
	'''
	# Extract full ld between genes, variants, and gene-variants
	#gene_variant_full_ld = extract_full_gene_variant_ld(tgfm_data['gene_eqtl_pmces'], tgfm_data['reference_ld'])

	##############################
	# Run TGFM
	###############################
	#tmp_output_file = tgfm_output_stem + '_' + window_name + '_temp_output_'
	#tgfm_obj = tgfm_inference_shell(tgfm_data, gene_log_prior, var_log_prior, gene_variant_full_ld, init_method, est_resid_var_bool, tmp_output_file)

	##############################
	# Organize TGFM data and print to results
	###############################
	# Extract names of genetic elements
	genetic_element_names = np.hstack((tgfm_data['genes'], tgfm_data['variants']))
	# Extract dictionary list of genetic elements in the middel of this window
	middle_genetic_elements = extract_middle_genetic_elements(tgfm_data['genes'], tgfm_data['middle_gene_indices'], tgfm_data['variants'], tgfm_data['middle_variant_indices'])

	# Write pickle file
	window_tgfm_output_file = tgfm_output_stem + '_' + window_name + '_results.pkl'
	# Load in tgfm input data
	g = open(window_tgfm_output_file, "rb")
	tgfm_res = pickle.load(g)
	g.close()


	# First print PIP results and then print component results
	# Extract genetic element pips
	genetic_element_pips = compute_pips(np.hstack((tgfm_res['alpha_phi'], tgfm_res['beta_phi'])))

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


f.close()

# Close file handles

t_pip.close()
