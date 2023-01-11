import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import pickle
import tgfm



def load_in_heritability_model(h2_model_file, num_variant_anno):
	tmp = np.loadtxt(h2_model_file, dtype=str, delimiter='\t')
	variant_anno_coef = tmp[1:(num_variant_anno+1),1].astype(float)
	gene_coef = tmp[(num_variant_anno+1):,1,].astype(float)
	return variant_anno_coef, gene_coef

def extract_tissue_names(gtex_pseudotissue_file):
	f = open(gtex_pseudotissue_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0])
	f.close()
	return np.asarray(arr)

def extract_expected_non_mediated_variant_variance(variant_non_mediated_heritability_model, annotation_mat):
	return np.dot(annotation_mat, variant_non_mediated_heritability_model)


def get_tissues_from_full_gene_names(gene_names):
	tissues = []
	for gene_name in gene_names:
		gene_info = gene_name.split('_')
		tissue = '_'.join(gene_info[1:])
		tissues.append(tissue)
	return np.asarray(tissues)

def extract_expected_gene_mediated_variance(full_gene_names, tissue_to_position_mapping, gene_mediated_heritability_model):
	tissue_names = get_tissues_from_full_gene_names(full_gene_names)
	variances = []
	for tissue_name in tissue_names:
		per_gene_variance = gene_mediated_heritability_model[tissue_to_position_mapping[tissue_name]]
		variances.append(per_gene_variance)

	return np.asarray(variances)

def extract_per_gene_eqtl_pmces_from_susie_eqtl_data(susie_mus, susie_alphas, n_genes):
	eqtl_pmces_arr = []	
	for g_index in range(n_genes):
		eqtl_pmces = np.sum((susie_mus[g_index])*(susie_alphas[g_index]),axis=0)
		eqtl_pmces_arr.append(eqtl_pmces)

	return eqtl_pmces_arr

def compute_gene_variance(susie_mu, susie_mu_sd, susie_alpha, ld):
	gene_var = 0.0

	# Component level eqtl effect sizes for this gene		
	gene_component_effect_sizes = (susie_mu)*susie_alpha

	# eQTL effect sizes for this gene
	gene_eqtl_effect_sizes = np.sum(gene_component_effect_sizes,axis=0)


	num_susie_components = susie_mu.shape[0]
	for k_index in range(num_susie_components):
		gene_var = gene_var + np.sum((np.square(susie_mu[k_index,:]) + np.square(susie_mu_sd[k_index,:]))*np.diag(ld)*susie_alpha[k_index,:])
		eqtl_component_pmces = (susie_mu[k_index,:])*(susie_alpha[k_index,:])
		gene_var = gene_var - np.dot(np.dot(eqtl_component_pmces,ld), eqtl_component_pmces)
	gene_var = gene_var + np.dot(np.dot(gene_eqtl_effect_sizes,ld), gene_eqtl_effect_sizes)
				
	return gene_var



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


def extract_valid_joint_susie_components(alpha_phi, beta_phi, alpha_ld, beta_ld, ld_thresh):
	num_components = alpha_phi.shape[0]
	valid_components = []

	num_genes = alpha_phi.shape[1]
	num_variants = beta_phi.shape[1]

	for component_num in range(num_components):
		cs_predictors = get_credible_set_genes(np.hstack((alpha_phi[component_num,:], beta_phi[component_num,:])), .95)
		cs_genes = cs_predictors[cs_predictors < num_genes]
		cs_variants = cs_predictors[cs_predictors >= num_genes] - num_genes

		# absolute ld among genes in credible set
		cs_gene_ld = np.abs(alpha_ld[cs_genes,:][:,cs_genes])
		cs_variant_ld = np.abs(beta_ld[cs_variants,:][:,cs_variants])

		if len(cs_genes) >= 1 and np.min(cs_gene_ld) > ld_thresh:
			valid_components.append(component_num)
		elif len(cs_variants) <= 1 or np.min(cs_variant_ld) > ld_thresh:
			valid_components.append(component_num)

	return valid_components


def get_probability_coming_from_each_tissue(ordered_tissue_names, tissue_to_position_mapping, window_gene_names, window_gene_probs):
	tiss_probs = np.zeros(len(ordered_tissue_names))
	window_tissue_names = get_tissues_from_full_gene_names(window_gene_names)
	for gene_index, window_tissue_name in enumerate(window_tissue_names):
		tissue_position = tissue_to_position_mapping[window_tissue_name]
		tiss_probs[tissue_position] = tiss_probs[tissue_position] + window_gene_probs[gene_index]
	return tiss_probs

def extract_middle_genetic_elements(ordered_genes, middle_gene_indices, ordered_variants, middle_variant_indices):
	dicti = {}
	for gene_name in ordered_genes[middle_gene_indices]:
		dicti[gene_name] = 1
	for variant_name in ordered_variants[middle_variant_indices]:
		dicti[variant_name] = 1
	return dicti

def get_n_eqtl_components_per_gene(valid_susie_eqtl_components):
	n_eqtl_components_per_gene = []
	for valid_susie_eqtl_component in valid_susie_eqtl_components:
		n_eqtl_components_per_gene.append(len(valid_susie_eqtl_component))
	return np.asarray(n_eqtl_components_per_gene)

def get_ukbb_windows(window_input_data):
	window_names = []
	n_windows = window_input_data.shape[0]
	for window_iter in range(n_windows):
		# Extract data for window corresponding to window iter
		data = window_input_data[window_iter,:]

		# Name of window
		window_name = data[0] + ':' + data[1] + ':' + data[2]
		window_names.append(window_name)
	return np.asarray(window_names)


#########################
# Command line arguments
#########################
trait_name = sys.argv[1]
ukkbb_window_summary_file = sys.argv[2]
gtex_pseudotissue_file = sys.argv[3]
preprocessed_tgfm_data_dir = sys.argv[4]
tgfm_sldsc_results_dir = sys.argv[5]
samp_size = sys.argv[6]
gene_type = sys.argv[7]
tgfm_results_dir = sys.argv[8]
job_number = int(sys.argv[9])
num_jobs = int(sys.argv[10])
if gene_type == 'cis_heritable_genes':
	gene_type = 'cis_heritable_gene'





# Load in h2 model for this trait
num_variant_anno = 93
h2_model_file = tgfm_sldsc_results_dir + trait_name + '_baselineLD_no_qtl_pmces_gene_adj_ld_scores_organized_sparse_ard_no_geno_regularization_res.txt'
variant_non_mediated_heritability_model, gene_mediated_heritability_model = load_in_heritability_model(h2_model_file, num_variant_anno)

# Extract ordered tissue information
ordered_tissue_names = extract_tissue_names(gtex_pseudotissue_file)
tissue_to_position_mapping = {}
for i, val in enumerate(ordered_tissue_names):
	tissue_to_position_mapping[val] = i


# Open cs output file handle
component_cs_output_file = tgfm_results_dir + trait_name + '_tgfm_component_cs_summary_' + str(job_number) + '_' + str(num_jobs) + '.txt'
t_cs = open(component_cs_output_file,'w')
t_cs.write('window_name\tcomponent_index\tgene_mediated_probability\tinclusion_elements\tinclusion_probabilities\n')
# Open tissue prob output file handle
component_tissue_output_file = tgfm_results_dir + trait_name + '_tgfm_component_tissue_prob_summary_' + str(job_number) + '_' + str(num_jobs) + '.txt'
t_tiss = open(component_tissue_output_file,'w')
t_tiss.write('window_name\tcomponent_index\tgene_mediated_probability')
for tissue_name in ordered_tissue_names:
	t_tiss.write('\t' + tissue_name)
t_tiss.write('\n')
# Open gene prob output file handle
component_gene_output_file = tgfm_results_dir + trait_name + '_tgfm_component_gene_prob_summary_' + str(job_number) + '_' + str(num_jobs) + '.txt'
t_gene = open(component_gene_output_file,'w')
t_gene.write('window_name\tcomponent_index\tgene_names\tn_components_per_gene\tgene_mediated_probability\n')



# Now loop through windows
# In each window run TGFM independently
# Loop through trait components
window_input_data = np.loadtxt(ukkbb_window_summary_file,dtype=str,delimiter='\t')
window_input_data = window_input_data[1:,:]
n_windows = window_input_data.shape[0]
ukbb_windows = get_ukbb_windows(window_input_data)

# Subset to just windows in this parallel run
ukbb_windows_parr = np.array_split(ukbb_windows, num_jobs)[job_number]

for window_name in ukbb_windows_parr:
	print(window_name)

	# Load in data for this window
	# Trait agnostic data
	trait_agnostic_data_file = preprocessed_tgfm_data_dir + gene_type + '_' + window_name + '_tgfm_trait_agnostic_data_obj.pkl'
	f = open(trait_agnostic_data_file, "rb")
	trait_agnostic_data = pickle.load(f)
	f.close()
	# Trait shared data
	trait_shared_data_file = preprocessed_tgfm_data_dir + gene_type + '_' + window_name + '_rss_likelihood_shared_standardized_data.pkl'
	f = open(trait_shared_data_file, "rb")
	trait_shared_data = pickle.load(f)
	f.close()

	# Trait shared data
	trait_specific_data_file = preprocessed_tgfm_data_dir + gene_type + '_' + window_name + '_rss_likelihood_' + trait_name + '_standardized_data.pkl'
	f = open(trait_specific_data_file, "rb")
	trait_specific_data = pickle.load(f)
	f.close()

	# Extract expected genetic element variances
	non_med_variant_h2 = extract_expected_non_mediated_variant_variance(variant_non_mediated_heritability_model, trait_agnostic_data['annotation'])
	gene_h2 = extract_expected_gene_mediated_variance(trait_shared_data['genes'], tissue_to_position_mapping, gene_mediated_heritability_model)
	# Make negative entries zero
	non_med_variant_h2[non_med_variant_h2 < 0.0] = 1e-30
	gene_h2[gene_h2 < 0.0] = 1e-30
	#gene_h2=gene_h2*50.0
	#pdb.set_trace()

	# Convert gene variances into categorical random variables
	gene_pi_init = gene_h2/(np.sum(gene_h2) + np.sum(non_med_variant_h2))
	non_med_variant_pi_init = non_med_variant_h2/(np.sum(gene_h2) + np.sum(non_med_variant_h2))

	# Get per gene eqtl pmces
	per_gene_eqtl_pmces = trait_shared_data['gene_eqtl_pmces']

	# Re-organize TGFM data
	tgfm_data = {}
	tgfm_data['genes'] = trait_shared_data['genes']
	tgfm_data['variants'] = trait_shared_data['variants']
	tgfm_data['gwas_beta'] = trait_specific_data['standardized_gwas_beta']
	tgfm_data['gwas_beta_se'] = trait_specific_data['standardized_gwas_beta_se']
	tgfm_data['s_diag'] = trait_specific_data['s_diag']
	tgfm_data['s_inv_2_diag'] = trait_specific_data['s_inv_2_diag']
	tgfm_data['D_diag'] = trait_specific_data['D_diag']
	tgfm_data['gwas_sample_size'] = trait_specific_data['gwas_sample_size']
	tgfm_data['reference_ld'] = trait_shared_data['reference_ld']
	tgfm_data['gene_eqtl_pmces'] = per_gene_eqtl_pmces

	#tgfm_data['precomputed_a_terms'] = trait_specific_data['precomputed_a_terms']
	#tgfm_data['precomputed_gene_gene_terms'] = trait_specific_data['precomputed_gene_gene_terms']

	#tgfm_data['susie_mu'] = trait_agnostic_data['susie_mu']
	#tgfm_data['susie_alpha'] = trait_agnostic_data['susie_alpha']
	#tgfm_data['susie_mu_sd'] = trait_agnostic_data['susie_mu_sd']
	'''
	for g_index in range(len(tgfm_data['susie_mu'])):
		gene_variance = compute_gene_variance(tgfm_data['susie_mu'][g_index], tgfm_data['susie_mu_sd'][g_index], tgfm_data['susie_alpha'][g_index], tgfm_data['reference_ld'])
		tgfm_data['susie_mu'][g_index] = tgfm_data['susie_mu'][g_index]/np.sqrt(gene_variance)
		tgfm_data['susie_mu_sd'][g_index] = tgfm_data['susie_mu_sd'][g_index]/np.sqrt(gene_variance)
	'''
	#per_gene_eqtl_pmces = extract_per_gene_eqtl_pmces_from_susie_eqtl_data(tgfm_data['susie_mu'], tgfm_data['susie_alpha'], len(trait_agnostic_data['genes']))
	if len(tgfm_data['genes']) == 0:
		continue

	# Standardized eqtl pmces
	tmp = np.asarray(tgfm_data['gene_eqtl_pmces'])
	gene_variances = np.diag(np.dot(np.dot(tmp, tgfm_data['reference_ld']), np.transpose(tmp)))
	for g_index in range(len(tgfm_data['genes'])):
		tgfm_data['gene_eqtl_pmces'][g_index] = tgfm_data['gene_eqtl_pmces'][g_index]/np.sqrt(gene_variances[g_index])


	# RUN TGFM
	#tgfm_obj = tgfm.TGFM(estimate_prior_variance=True, init_pi=gene_pi_init, variant_init_pi=non_med_variant_pi_init, convergence_thresh=1e-7, max_iter=500)
	tgfm_obj = tgfm.TGFM(estimate_prior_variance=True, init_pi=gene_pi_init, variant_init_pi=non_med_variant_pi_init, convergence_thresh=1e-5, max_iter=500)
	tgfm_obj.fit(twas_data_obj=tgfm_data)

	# Re-organize TGFM results to summary format
	# Extract components that pass purity filter
	valid_tgfm_components = extract_valid_joint_susie_components(tgfm_obj.alpha_phi_no_weights, tgfm_obj.beta_phi_no_weights, tgfm_obj.ge_ld, tgfm_data['reference_ld'], .5)
	# Extract names of genetic elements
	genetic_element_names = np.hstack((tgfm_data['genes'], tgfm_data['variants']))
	# Extract dictionary list of genetic elements in the middel of this window
	middle_genetic_elements = extract_middle_genetic_elements(tgfm_data['genes'], trait_shared_data['middle_gene_indices'], tgfm_data['variants'], trait_shared_data['middle_variant_indices'])
	# Get number of eqtl components per gene
	n_eqtl_components_per_gene = get_n_eqtl_components_per_gene(trait_agnostic_data['valid_susie_components'])

	# loop through TGFM components for this window
	for tgfm_component in valid_tgfm_components:
		# Get probability component is mediated by gene expression in any cis tissue, gene
		mediated_probability = np.sum(tgfm_obj.alpha_phi[tgfm_component,:])
		# Get probability coming from each tissue
		tissue_mediated_probabilities = get_probability_coming_from_each_tissue(ordered_tissue_names, tissue_to_position_mapping, trait_shared_data['genes'], tgfm_obj.alpha_phi[tgfm_component,:])
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
		t_gene.write(window_name + '\t' + str(tgfm_component) + '\t' + ';'.join(trait_shared_data['genes']) + '\t' + ';'.join(n_eqtl_components_per_gene.astype(str)) + '\t' + ';'.join(tgfm_obj.alpha_phi[tgfm_component,:].astype(str)) + '\n')
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
	window_tgfm_output_file = tgfm_results_dir + trait_name + '_' + window_name + '_results.pkl'
	g = open(window_tgfm_output_file, "wb")
	pickle.dump(tgfm_results, g)
	g.close()


t_cs.close()
t_tiss.close()
t_gene.close()





