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

	return tgfm_obj


def tgfm_inference_shell(tgfm_data, gene_log_prior, var_log_prior, ld_mat, init_method):
	tgfm_obj = bootstrapped_tgfm.TGFM(L=20, estimate_prior_variance=True, gene_init_log_pi=gene_log_prior, variant_init_log_pi=var_log_prior, convergence_thresh=1e-5, max_iter=5)
	tgfm_obj.fit(twas_data_obj=tgfm_data)
	pdb.set_trace()


def tgfm_inference_shell_old(tgfm_data, gene_log_prior, var_log_prior, gene_variant_full_ld, init_method, attenuated_gene_variances):
	# Hacky: Initialize old TGFM object using only one iter of optimization
	tgfm_obj = tgfm.TGFM(L=20, estimate_prior_variance=False, gene_init_log_pi=gene_log_prior, variant_init_log_pi=var_log_prior, convergence_thresh=1e-5, max_iter=1)
	tgfm_obj.fit(twas_data_obj=tgfm_data)
	# More hack: need to redo twas z
	variant_z = tgfm_data['gwas_beta']/tgfm_data['gwas_beta_se']
	new_gene_z = np.dot(tgfm_data['gene_eqtl_pmces'], variant_z)


	#tgfm_obj_temp = tgfm.TGFM(L=20, estimate_prior_variance=True, gene_init_log_pi=gene_log_prior, variant_init_log_pi=var_log_prior, convergence_thresh=1e-5, max_iter=10, single_variance_component=False)
	#tgfm_obj_temp.fit(twas_data_obj=tgfm_data)

	tgfm_obj.nominal_twas_z = new_gene_z

	# Create vector of concatenated z-scores
	z_vec = np.hstack((new_gene_z,variant_z))

	# Create concatenated vector of prior probs
	prior_probs = np.hstack((np.exp(gene_log_prior), np.exp(var_log_prior)))

	pdb.set_trace()
	temp_se = np.ones(len(z_vec))/np.sqrt(100000.0)
	temp_se[:len(new_gene_z)] = temp_se[:len(new_gene_z)]*(1.0/np.sqrt(attenuated_gene_variances))
	susie_null_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), shat=temp_se.reshape((len(temp_se), 1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=False)



	if init_method == 'best':
		# Run susie with only variants
		p_var_only = np.ones(len(z_vec))
		p_var_only[:len(tgfm_obj.nominal_twas_z)] = 0.0
		p_var_only = p_var_only/np.sum(p_var_only)
		susie_variant_only = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=p_var_only.reshape((len(p_var_only),1)), estimate_residual_variance=False)
	
		# Run susie with variant initialization
		init_obj = {'alpha':susie_variant_only.rx2('alpha'), 'mu':susie_variant_only.rx2('mu'),'mu2':susie_variant_only.rx2('mu2')}
		init_obj2 = ro.ListVector(init_obj)
		init_obj2.rclass = rpy2.robjects.StrVector(("list", "susie"))
		susie_variant_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, s_init=init_obj2, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=False)

		# Run susie with null initialization
		susie_null_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=False)

		# Select model with largest elbo
		susie_null_init_elbo = susie_null_init.rx2('elbo')[-1]
		susie_variant_init_elbo = susie_variant_init.rx2('elbo')[-1]

		if susie_variant_init_elbo > susie_null_init_elbo:
			# Variant init wins
			tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_variant_init)
		else:
			# Null init wins
			tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_null_init)
	elif init_method == 'null':
		# Run susie with null initialization
		susie_null_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=False)
		tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_null_init)
	elif init_method == 'variant_only':
		# Run susie with only variants
		p_var_only = np.ones(len(z_vec))
		p_var_only[:len(tgfm_obj.nominal_twas_z)] = 0.0
		p_var_only = p_var_only/np.sum(p_var_only)
		susie_variant_only = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, prior_weights=p_var_only.reshape((len(p_var_only),1)), estimate_residual_variance=False)
	
		# Run susie with variant initialization
		init_obj = {'alpha':susie_variant_only.rx2('alpha'), 'mu':susie_variant_only.rx2('mu'),'mu2':susie_variant_only.rx2('mu2')}
		init_obj2 = ro.ListVector(init_obj)
		init_obj2.rclass = rpy2.robjects.StrVector(("list", "susie"))
		susie_variant_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=20, s_init=init_obj2, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=False)
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

	ld_file = data[1]
	tgfm_input_pkl = data[2]
	log_prior_prob_file = data[3] + '_' + ln_pi_method_name + '.txt'

	if window_name != 'chr1_109734349_112734349':
		continue


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

	# Write pickle file
	window_tgfm_output_file = tgfm_output_stem + '_' + window_name + '_results.pkl'
	g = open(window_tgfm_output_file, "wb")
	pickle.dump(tgfm_results, g)
	g.close()

f.close()

# Close file handles
t_cs.close()
t_gene.close()
t_tiss.close()
