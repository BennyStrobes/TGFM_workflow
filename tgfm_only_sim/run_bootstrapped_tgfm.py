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
ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
susieR_pkg = importr('susieR')


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

def merge_two_bootstrapped_tgfms_based_on_elbo(tgfm_obj, tgfm_obj2, variant_z_vec, variant_ld_mat, gwas_sample_size):
	for bs_iter in range(tgfm_obj.n_bs):
		# Extract causal eqtl effects for this bootstrap
		bs_eqtls_pmces = np.zeros((tgfm_obj.G, tgfm_obj.K))
		bs_eqtls_pmces_sparse = np.load(tgfm_obj.pmces_files[bs_iter])
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
	# Recompute PIPs
	tgfm_obj.compute_pips()

	return tgfm_obj




def tgfm_inference_shell(tgfm_data, gene_log_prior, var_log_prior, ld_mat, init_method):
	if init_method == 'null':
		tgfm_obj = bootstrapped_tgfm.TGFM(L=20, estimate_prior_variance=True, gene_init_log_pi=gene_log_prior, variant_init_log_pi=var_log_prior, convergence_thresh=1e-5, max_iter=10)
		tgfm_obj.fit(twas_data_obj=tgfm_data)
	elif init_method == 'best':
		# Run susie with only variants
		z_vec = tgfm_data['gwas_beta']/tgfm_data['gwas_beta_se']

		# Now run tgfm sampler with null init
		tgfm_obj = bootstrapped_tgfm.TGFM(L=20, estimate_prior_variance=True, gene_init_log_pi=gene_log_prior, variant_init_log_pi=var_log_prior, convergence_thresh=1e-5, max_iter=10)
		tgfm_obj.fit(twas_data_obj=tgfm_data)
		#elbo_null_init = compute_elbo_for_bootstrapped_tgfm_obj_shell(tgfm_obj, z_vec, ld_mat, tgfm_data['gwas_sample_size'])

		if np.max(np.max(tgfm_obj.expected_alpha_pips)) > .5:
			# Create initialization alpha, mu, and mu_var
			susie_variant_only = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=ld_mat, n=tgfm_data['gwas_sample_size'], L=20, estimate_residual_variance=False)

			num_snps = len(z_vec)
			num_genes = len(tgfm_data['genes'])
			alpha_init = np.zeros((20, num_genes + num_snps))
			mu_init = np.zeros((20, num_genes + num_snps))
			mu_var_init = np.zeros((20, num_genes + num_snps))
			alpha_init[:, num_genes:] = susie_variant_only.rx2('alpha')
			mu_init[:, num_genes:] = susie_variant_only.rx2('mu')
			variant_only_mu_var = susie_variant_only.rx2('mu2') - np.square(susie_variant_only.rx2('mu'))
			mu_var_init[:, num_genes:] = variant_only_mu_var

			# Run tgfm sampler with variant init
			tgfm_obj_variant_init = bootstrapped_tgfm.TGFM(L=20, estimate_prior_variance=True, gene_init_log_pi=gene_log_prior, variant_init_log_pi=var_log_prior, convergence_thresh=1e-5, max_iter=10)
			tgfm_obj_variant_init.fit(twas_data_obj=tgfm_data, phi_init=alpha_init, mu_init=mu_init, mu_var_init=mu_var_init)
			#elbo_variant_init = compute_elbo_for_bootstrapped_tgfm_obj_shell(tgfm_obj_variant_init, z_vec, ld_mat, tgfm_data['gwas_sample_size'])

			tgfm_obj = merge_two_bootstrapped_tgfms_based_on_elbo(tgfm_obj, tgfm_obj_variant_init, z_vec, ld_mat, tgfm_data['gwas_sample_size'])
		
	return tgfm_obj



######################
# Command line args
######################
tgfm_input_file = sys.argv[1]
tgfm_output_stem = sys.argv[2]
ln_pi_method_name = sys.argv[3]
init_method = sys.argv[4]
est_resid_var = sys.argv[5]


# Extract ordered tissue information
tissue_to_position_mapping = {}
ordered_tissue_names = []
for i in range(10):
	tissue_to_position_mapping['tissue' + str(i)] = i
	ordered_tissue_names.append('tissue' + str(i))
ordered_tissue_names = np.asarray(ordered_tissue_names)


# Open cs output file handle
component_cs_output_file = tgfm_output_stem + '_tgfm_component_bs_cs_summary.txt'
t_cs = open(component_cs_output_file,'w')
t_cs.write('window_name\tinclusion_elements\tinclusion_probabilities\n')



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


	# Extract full ld between genes, variants, and gene-variants
	#gene_variant_full_ld = extract_full_gene_variant_ld(tgfm_data['gene_eqtl_pmces'], tgfm_data['reference_ld'])

	##############################
	# Run TGFM
	###############################
	tgfm_obj = tgfm_inference_shell(tgfm_data, gene_log_prior, var_log_prior, ld_mat, init_method)

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
	t_cs.write(window_name + '\t')
	t_cs.write(';'.join(ordered_middle_names) + '\t')
	t_cs.write(';'.join(ordered_middle_pips.astype(str)) + '\n')

	# Save all TGFM results to pkl
	tgfm_results = {}
	tgfm_results['variants'] = tgfm_data['variants']
	tgfm_results['genes'] = tgfm_data['genes']
	tgfm_results['alpha_phis'] = tgfm_obj.alpha_phis
	tgfm_results['beta_phis'] = tgfm_obj.beta_phis
	tgfm_results['alpha_pips'] = tgfm_obj.alpha_pips
	tgfm_results['beta_pips'] = tgfm_obj.beta_pips
	tgfm_results['expected_alpha_pips'] = tgfm_obj.expected_alpha_pips
	tgfm_results['expected_beta_pips'] = tgfm_obj.expected_beta_pips


	# Write pickle file
	window_tgfm_output_file = tgfm_output_stem + '_' + window_name + '_results.pkl'
	g = open(window_tgfm_output_file, "wb")
	pickle.dump(tgfm_results, g)
	g.close()

f.close()

# Close file handles
t_cs.close()

