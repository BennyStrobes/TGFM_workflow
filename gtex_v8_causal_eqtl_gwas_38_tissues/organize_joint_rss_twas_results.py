import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os
import pdb
import pickle
import scipy.special
import scipy.stats

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



def extract_valid_susie_components(phi, ge_ld, ld_thresh):
	num_components = phi.shape[0]
	valid_components = []

	for component_num in range(num_components):
		cs_genes = get_credible_set_genes(phi[component_num,:], .95)

		# absolute ld among genes in credible set
		cs_gene_ld = np.abs(ge_ld[cs_genes,:][:,cs_genes])

		if np.min(cs_gene_ld) > ld_thresh:
			valid_components.append(component_num)

	return np.asarray(valid_components)

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

def extract_manhatten_plot_data_from_tgfm_joint_susie_tissue_prior(tgfm_const_prior_susie_twas_pickle_file, raw_twas_data, manhatten_plot_data_file):
	f = open(tgfm_const_prior_susie_twas_pickle_file, "rb")
	tgfm_twas_data = pickle.load(f)
	f.close()

	# Extract relevent fields from twas object
	gene_tissue_pairs = tgfm_twas_data.genes
	# Get gene names and tissue names
	gene_names = []
	tissue_names = []
	for gene_tissue_pair in gene_tissue_pairs:
		gene_name = gene_tissue_pair.split('_')[0]
		tissue_name = '_'.join(gene_tissue_pair.split('_')[1:])
		gene_names.append(gene_name)
		tissue_names.append(tissue_name)
	gene_names = np.asarray(gene_names)
	tissue_names = np.asarray(tissue_names)	

	valid_tgfm_components = extract_valid_joint_susie_components(tgfm_twas_data.alpha_phi_no_weights, tgfm_twas_data.beta_phi_no_weights, tgfm_twas_data.ge_ld, raw_twas_data['reference_ld'], .5)

	# Extract gwas z scores
	gwas_z_scores = raw_twas_data['gwas_beta']/raw_twas_data['gwas_beta_se']
	gwas_p_values = scipy.stats.norm.sf(np.abs(gwas_z_scores))*2
	variant_ids = raw_twas_data['variants']
	variant_chrom = raw_twas_data['bim'][:,0]
	variant_pos = raw_twas_data['bim'][:,3]

	component_names = []
	component_types = []
	variant_pis = []


	for tgfm_component in valid_tgfm_components:
		# Option 1 is component is mediated heavily by a gene, tissue pair
		if np.max(tgfm_twas_data.alpha_phi[tgfm_component,:]) > .2:
			top_gene_index = np.argmax(tgfm_twas_data.alpha_phi[tgfm_component,:])
			top_gene_name = gene_names[top_gene_index]
			top_tissue_name = tissue_names[top_gene_index]
			variant_pi = raw_twas_data['susie_alpha'][top_gene_index][0,:]
			variant_pi[variant_pi < .05] = 0.0

			variant_pis.append(variant_pi)
			component_names.append('Component_' + str(tgfm_component) + '-' + top_tissue_name + '-mediated')
		else: # non-mediation
			variant_pi = tgfm_twas_data.beta_phi[tgfm_component,:]
			variant_pi[variant_pi < .05] = 0.0

			variant_pis.append(variant_pi)
			component_names.append('Component_' + str(tgfm_component) + '-non-mediated')

	n_var = len(gwas_p_values)
	if len(valid_tgfm_components) == 0:
		fm_variant_anno = np.asarray(['NA']*n_var)
	else:
		fm_variant_anno = []
		variant_pis = np.transpose(np.asarray(variant_pis))
		for var_index in range(n_var):
			if np.max(variant_pis[var_index,:]) > 0.0:
				top_component = np.argmax(variant_pis[var_index,:])
				variant_anno = component_names[top_component]
				fm_variant_anno.append(variant_anno)
			else:
				fm_variant_anno.append('NA')
		fm_variant_anno = np.asarray(fm_variant_anno)

	t = open(manhatten_plot_data_file,'w')
	t.write('variant_id\tposition\tp_value\tfm_annotation\n')
	for var_index in range(n_var):
		t.write(variant_ids[var_index] + '\t' + variant_pos[var_index] + '\t' + str(gwas_p_values[var_index]) + '\t' + fm_variant_anno[var_index] + '\n')
	t.close()


def extract_data_from_tgfm_joint_susie_tissue_prior_twas_pickle_file(tgfm_const_prior_susie_twas_pickle_file, raw_twas_data, gwas_susie_alpha):
	f = open(tgfm_const_prior_susie_twas_pickle_file, "rb")
	tgfm_twas_data = pickle.load(f)
	f.close()

	# Extract relevent fields from twas object
	gene_tissue_pairs = tgfm_twas_data.genes
	# Get gene names and tissue names
	gene_names = []
	tissue_names = []
	for gene_tissue_pair in gene_tissue_pairs:
		gene_name = gene_tissue_pair.split('_')[0]
		tissue_name = '_'.join(gene_tissue_pair.split('_')[1:])
		gene_names.append(gene_name)
		tissue_names.append(tissue_name)
	gene_names = np.asarray(gene_names)
	tissue_names = np.asarray(tissue_names)


	#valid_susie_components = extract_valid_joint_susie_components(tgfm_twas_data.alpha_phi_no_weights, tgfm_twas_data.beta_phi_no_weights, tgfm_twas_data.ge_ld, raw_twas_data['reference_ld'], .5)
	

	corrz = np.corrcoef(gwas_susie_alpha, tgfm_twas_data.beta_phi)[0,1:]

	matched_component = np.argmax(np.abs(corrz))
	matched_corr = corrz[matched_component]
	if matched_corr < .2:
		print('matched corr: ' + str(matched_corr))

	#if np.sum(tgfm_twas_data.beta_phi[matched_component,:]) > .5:
		#print(np.sum(tgfm_twas_data.alpha_phi,axis=1))
	#valid_alpha_susie_components = extract_valid_susie_components(tgfm_twas_data.alpha_phi_no_weights, tgfm_twas_data.ge_ld, .5)
	#valid_beta_susie_components = extract_valid_susie_components(tgfm_twas_data.beta_phi_no_weights, raw_twas_data['reference_ld'], .5)

	'''

	if len(valid_susie_components) == 0:
		pips = np.zeros(len(gene_names))
	else:
		pips = 1.0-np.prod(1.0 - tgfm_twas_data.phi[valid_susie_components,:],axis=0)

	if np.min(np.abs(tgfm_twas_data.ge_ld)) > .5:
		print('situation 1')
		if np.max(np.abs(tgfm_twas_data.nominal_twas_z)) < 4.0:
			pips=pips*0.0
	elif len(valid_susie_components) > 8:
		print('situation 2')
	'''

	return gene_names, tissue_names, tgfm_twas_data.alpha_phi[matched_component,:], np.sum(tgfm_twas_data.beta_phi[matched_component,:])

def extract_data_from_tgfm_susie_tissue_prior_twas_pickle_file(tgfm_const_prior_susie_twas_pickle_file):
	f = open(tgfm_const_prior_susie_twas_pickle_file, "rb")
	tgfm_twas_data = pickle.load(f)
	f.close()

	# Extract relevent fields from twas object
	gene_tissue_pairs = tgfm_twas_data.genes
	# Get gene names and tissue names
	gene_names = []
	tissue_names = []
	for gene_tissue_pair in gene_tissue_pairs:
		gene_name = gene_tissue_pair.split('_')[0]
		tissue_name = '_'.join(gene_tissue_pair.split('_')[1:])
		gene_names.append(gene_name)
		tissue_names.append(tissue_name)
	gene_names = np.asarray(gene_names)
	tissue_names = np.asarray(tissue_names)

	valid_susie_components = extract_valid_susie_components(tgfm_twas_data.phi_no_weights, tgfm_twas_data.ge_ld, .5)

	if len(valid_susie_components) == 0:
		pips = np.zeros(len(gene_names))
	else:
		pips = 1.0-np.prod(1.0 - tgfm_twas_data.phi[valid_susie_components,:],axis=0)

	if np.min(np.abs(tgfm_twas_data.ge_ld)) > .5:
		print('situation 1')
		if np.max(np.abs(tgfm_twas_data.nominal_twas_z)) < 4.0:
			pips=pips*0.0
	elif len(valid_susie_components) > 8:
		print('situation 2')

	return gene_names, tissue_names, pips, len(valid_susie_components)

def extract_data_from_tgfm_susie_twas_pickle_file(tgfm_const_prior_susie_twas_pickle_file):
	f = open(tgfm_const_prior_susie_twas_pickle_file, "rb")
	tgfm_twas_data = pickle.load(f)
	f.close()

	# Extract relevent fields from twas object
	gene_tissue_pairs = tgfm_twas_data.genes
	# Get gene names and tissue names
	gene_names = []
	tissue_names = []
	for gene_tissue_pair in gene_tissue_pairs:
		gene_name = gene_tissue_pair.split('_')[0]
		tissue_name = '_'.join(gene_tissue_pair.split('_')[1:])
		gene_names.append(gene_name)
		tissue_names.append(tissue_name)
	gene_names = np.asarray(gene_names)
	tissue_names = np.asarray(tissue_names)

	valid_susie_components = extract_valid_susie_components(tgfm_twas_data.phi, tgfm_twas_data.ge_ld, .5)

	if len(valid_susie_components) == 0:
		pips = np.zeros(len(gene_names))
	else:
		pips = 1.0-np.prod(1.0 - tgfm_twas_data.phi[valid_susie_components,:],axis=0)

	if np.min(np.abs(tgfm_twas_data.ge_ld)) > .5:
		print('situation 1')
		if np.max(np.abs(tgfm_twas_data.nominal_twas_z)) < 4.0:
			pips=pips*0.0
	elif len(valid_susie_components) > 8:
		print('situation 2')

	return gene_names, tissue_names, pips, len(valid_susie_components)



def extract_data_from_tgfm_twas_pickle_file(tgfm_tiss_prior_twas_pickle_file):
	f = open(tgfm_tiss_prior_twas_pickle_file, "rb")
	tgfm_twas_data = pickle.load(f)
	f.close()

	# Extract relevent fields from twas object
	gene_tissue_pairs = tgfm_twas_data.genes
	twas_fusion_nominal_z_scores = tgfm_twas_data.nominal_twas_z
	twas_tgfm_nominal_z_scores = tgfm_twas_data.nominal_twas_rss_z
	twas_tgfm_multivariate_z_scores = tgfm_twas_data.alpha_mu/np.sqrt(tgfm_twas_data.alpha_var)
	
	# Get gene names and tissue names
	gene_names = []
	tissue_names = []
	for gene_tissue_pair in gene_tissue_pairs:
		gene_name = gene_tissue_pair.split('_')[0]
		tissue_name = '_'.join(gene_tissue_pair.split('_')[1:])
		gene_names.append(gene_name)
		tissue_names.append(tissue_name)
	gene_names = np.asarray(gene_names)
	tissue_names = np.asarray(tissue_names)
	
	return gene_names, tissue_names, twas_fusion_nominal_z_scores, twas_tgfm_multivariate_z_scores, tgfm_twas_data.alpha_mu, tgfm_twas_data.alpha_var

def extract_data_from_tgfm_fine_mapping_table(tgfm_tiss_prior_fine_mapping_file):
	data = np.loadtxt(tgfm_tiss_prior_fine_mapping_file, dtype=str, delimiter='\t')[1:,:]

	# Create mapping from tissue-gene to sum_posterior prob
	tg_to_sum_posterior_prob = {}
	tg_to_max_posterior_prob = {}
	nrows = data.shape[0]

	for row_num in range(nrows):
		tissue_gene_name = data[row_num, 0] + '_' + data[row_num, 1]
		posterior_prob = float(data[row_num, 3])
		if tissue_gene_name not in tg_to_sum_posterior_prob:
			tg_to_sum_posterior_prob[tissue_gene_name] = posterior_prob
		else:
			tg_to_sum_posterior_prob[tissue_gene_name] = tg_to_sum_posterior_prob[tissue_gene_name] + posterior_prob
		if tissue_gene_name not in tg_to_max_posterior_prob:
			tg_to_max_posterior_prob[tissue_gene_name] = posterior_prob
		else:
			tg_to_max_posterior_prob[tissue_gene_name] = np.max([posterior_prob, tg_to_max_posterior_prob[tissue_gene_name]])
	return tg_to_sum_posterior_prob

def extract_data_from_tgfm_fine_mapping_table_and_log_like(tgfm_tiss_prior_fine_mapping_file):
	data = np.loadtxt(tgfm_tiss_prior_fine_mapping_file, dtype=str, delimiter='\t')[1:,:]
	# Create mapping from tissue-gene to sum_posterior prob
	tg_to_sum_posterior_prob = {}
	tg_to_max_posterior_prob = {}
	tg_to_log_sum_exp = {}
	nrows = data.shape[0]

	for row_num in range(nrows):
		tissue_gene_name = data[row_num, 0] + '_' + data[row_num, 1]
		posterior_prob = float(data[row_num, 3])
		log_likelihood = float(data[row_num, 4])
		if tissue_gene_name not in tg_to_sum_posterior_prob:
			tg_to_sum_posterior_prob[tissue_gene_name] = posterior_prob
		else:
			tg_to_sum_posterior_prob[tissue_gene_name] = tg_to_sum_posterior_prob[tissue_gene_name] + posterior_prob
		if tissue_gene_name not in tg_to_log_sum_exp:
			tg_to_log_sum_exp[tissue_gene_name] = [log_likelihood]
		else:
			tg_to_log_sum_exp[tissue_gene_name].append(log_likelihood)
		if tissue_gene_name not in tg_to_max_posterior_prob:
			tg_to_max_posterior_prob[tissue_gene_name] = posterior_prob
		else:
			tg_to_max_posterior_prob[tissue_gene_name] = np.max([posterior_prob, tg_to_max_posterior_prob[tissue_gene_name]])
	# Now put log in log-sum-exp
	for tg in [*tg_to_log_sum_exp]:
		tg_to_log_sum_exp[tg] = scipy.special.logsumexp(tg_to_log_sum_exp[tg])
	return tg_to_sum_posterior_prob, tg_to_log_sum_exp

def get_var_dicti(dir_name, trait_name):
	file_name = dir_name + trait_name + '_component_organized_multivariate_twas_overlaps.txt'
	dicti = {}
	if os.path.exists(file_name) == False:
		return dicti
	f = open(file_name)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissue_gene_name = data[1] + '_' + data[2]
		sd = float(data[5])
		dicti[tissue_gene_name] = sd
	f.close()
	return dicti

def extract_gene_variances(twas_data_file):
	f = open(twas_data_file, "rb")
	twas_data = pickle.load(f)
	f.close()

	n_genes = len(twas_data['susie_mu'])

	varz = []
	for g_index in range(n_genes):
		var = np.sum((np.square(twas_data['susie_mu_sd'][g_index]) + np.square(twas_data['susie_mu'][g_index]))*twas_data['susie_alpha'][g_index])
		varz.append(var)
	return np.asarray(varz)

def get_coloc_prob_dicti_for_set_of_genes(coloc_results_dir, unique_genes, trait_name, adaptive=False):
	dicti = {}
	for gene_name in unique_genes:
		if adaptive:
			file_name = coloc_results_dir + trait_name + '_' + gene_name + '_adaptive_coloc_posterior_probabilities.txt'
		else:
			file_name = coloc_results_dir + trait_name + '_' + gene_name + '_coloc_posterior_probabilities.txt'
		f = open(file_name)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			tissue_name = data[0]
			posterior_prob = (data[-1])
			tissue_gene_name = tissue_name + '_' + gene_name
			dicti[tissue_gene_name] = posterior_prob
		f.close()
	return dicti

def extract_gwas_component_susie_pmces_in_same_order_as_twas_data(gwas_component_susie_pmces_file, twas_bim):
	num_twas_var = twas_bim.shape[0]
	# First create mapping from twas variants to position
	mapping = {}
	for var_iter in range(num_twas_var):
		var_name = twas_bim[var_iter, 1]
		mapping[var_name] = (var_iter, twas_bim[var_iter,4], twas_bim[var_iter,5])
		var_info = var_name.split('_')
		alt_var_name = var_info[0] + '_' + var_info[1] + '_' + var_info[3] + '_' + var_info[2] + '_' + var_info[4]
		mapping[alt_var_name] = (var_iter, twas_bim[var_iter,4], twas_bim[var_iter,5])

	# Load in gwas pmces data
	gwas_component_pmces_data = np.loadtxt(gwas_component_susie_pmces_file, dtype=str,delimiter='\t')
	gwas_component_pmces_data = gwas_component_pmces_data[1:,:]

	
	# Initialize output array
	gwas_pmces = np.zeros(num_twas_var)
	gwas_susie_mu = np.zeros(num_twas_var)
	gwas_susie_mu_sd = np.zeros(num_twas_var)
	gwas_susie_alpha = np.zeros(num_twas_var)


	used_positions = {}
	for var_iter in range(gwas_component_pmces_data.shape[0]):
		variant_id = gwas_component_pmces_data[var_iter,0] + '_b38'
		if variant_id in mapping:
			variant_info = variant_id.split('_')
			rev = 1.0
			if variant_info[2] != mapping[variant_id][1]:
				rev = -1.0
			gwas_pmces[mapping[variant_id][0]] = float(gwas_component_pmces_data[var_iter,1])*rev
			gwas_susie_mu[mapping[variant_id][0]] = float(gwas_component_pmces_data[var_iter,2])*rev
			gwas_susie_mu_sd[mapping[variant_id][0]] = float(gwas_component_pmces_data[var_iter,4])
			gwas_susie_alpha[mapping[variant_id][0]] = float(gwas_component_pmces_data[var_iter,3])

			used_positions[mapping[variant_id][0]] = 1
	if len(used_positions) != num_twas_var:
		print('assumption eroror')
		pdb.set_trace()
	return gwas_pmces, gwas_susie_mu, gwas_susie_mu_sd, gwas_susie_alpha

def get_variant_in_credible_set_prob_for_each_gene(gwas_susie_alpha, eqtl_susie_alphas):
	# Get gwas 95% cred set indices
	'''
	sorted_alpha = -np.sort(-gwas_susie_alpha)
	total = 0
	for ele in sorted_alpha:
		if total > .95:
			break
		total = total + ele
		min_val = ele
	gwas_95_cred_set_indices = gwas_susie_alpha >= min_val
	'''

	overlap_probs = []
	for eqtl_susie_alpha in eqtl_susie_alphas:
		valid_rows = np.max(eqtl_susie_alpha,axis=1) > .01
		if np.sum(valid_rows) == 0:
			overlap_prob = 0.0
		else:
			eqtl_susie_alpha_sub = eqtl_susie_alpha[valid_rows,:]
			for row in range(eqtl_susie_alpha_sub.shape[0]):
				ordered_alphas = -np.sort(-eqtl_susie_alpha_sub[row,:])
				total = 0.0
				for ele in ordered_alphas:
					if total > .95:
						break
					total = total + ele
					min_val = ele
				eqtl_susie_alpha_sub[row,:] = (eqtl_susie_alpha_sub[row,:] >= min_val)*1.0
			variant_in_a_eqtl_cs = np.max(eqtl_susie_alpha_sub,axis=0)
			overlap_prob = np.dot(gwas_susie_alpha, variant_in_a_eqtl_cs)
		overlap_probs.append(overlap_prob)
	return np.asarray(overlap_probs)



def get_variant_overlap_prob_for_each_gene(gwas_susie_alpha, eqtl_susie_alphas):
	# Get gwas 95% cred set indices
	'''
	sorted_alpha = -np.sort(-gwas_susie_alpha)
	total = 0
	for ele in sorted_alpha:
		if total > .95:
			break
		total = total + ele
		min_val = ele
	gwas_95_cred_set_indices = gwas_susie_alpha >= min_val
	'''

	overlap_probs = []
	for eqtl_susie_alpha in eqtl_susie_alphas:
		valid_rows = np.max(eqtl_susie_alpha,axis=1) > .01
		if np.sum(valid_rows) == 0:
			eqtl_susie_alpha_sub = eqtl_susie_alpha[:,:]
		else:
			eqtl_susie_alpha_sub = eqtl_susie_alpha[valid_rows,:]
		eqtl_pip = 1.0-np.prod(1 - eqtl_susie_alpha_sub,axis=0)
		overlap_prob = np.dot(gwas_susie_alpha,eqtl_pip)
		overlap_probs.append(overlap_prob)
	return np.asarray(overlap_probs)

def create_gene_tissue_to_pph4_mapping(gene_names, coloc_results_dir, trait_name):
	mapping = {}
	for gene_name in gene_names:
		gene_pph4_file = coloc_results_dir + trait_name + '_' + gene_name + '_coloc_posterior_probabilities.txt'
		f = open(gene_pph4_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			tissue_name = data[0]
			pph4_prob_str = data[-1]
			gene_tissue_name = gene_name + '_' + tissue_name
			if gene_tissue_name in mapping:
				print('assumption eroror')
				pdb.set_trace()
			mapping[gene_tissue_name] = pph4_prob_str

		f.close()
	return mapping


####################
# Command line args
####################

trait_name = sys.argv[1]
gtex_pseudotissue_file = sys.argv[2]
pseudotissue_gtex_rss_multivariate_twas_dir = sys.argv[3]
gene_version = sys.argv[4]
fusion_weights = sys.argv[5]
ukbb_genome_wide_susie_organized_results_dir = sys.argv[6]
coloc_results_dir = sys.argv[7]


# Version with tissue specific prior (NEED TO FIX ISSUE STILLLLLLL)
output_file = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + gene_version + '_fusion_weights_' + fusion_weights + '_component_organized_joint_tgfm_results.txt'
t = open(output_file,'w')
t.write('trait_component\ttissue\tgene\tnominal_fusion_twas_z_score\ttgfm_const_prior_twas_z_score\trobust_tgfm_tissue_prior_joint_susie_prob\ttgfm_tissue_prior_twas_z_score\trobust_tgfm_tissue_prior_twas_z_score\ttgfm_rss_regression_const_prior_posterior_prob\ttgfm_rss_regression_tissue_prior_posterior_prob\trobust_tgfm_rss_regression_tissue_prior_posterior_prob\trobust_tgfm_rss_regression_tissue_prior_log_likelihood\tvariant_in_a_cs_prob\tvariant_overlap_prob\tcoloc_pph4\tmanhatten_plot_data_file\n')

#susie_output_file = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + gene_version + '_fusion_weights_' + fusion_weights + '_number_of_joint_susie_components.txt'
#t2 = open(susie_output_file,'w')
#t2.write('trait_component\tnumber_susie_components_cp\tnumber_susie_components_robust_cp\tnumber_susie_components_tp\tnumber_susie_components_robust_tp\n')


aa = []
bb = []
# Loop through chromosomes
for chrom_num in range(1,23):
	print(chrom_num)
	chromosome_summary_file = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + str(chrom_num) + '_summary_tgfm_results_' + gene_version + '_const_1e-5_prior.txt'
	# Stream chromosome summary file
	f = open(chromosome_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent data
		component_name = data[0]
		num_genes = data[1]

		# Several different versions of analysis
		# Version 1: TGFM with constant prior
		if fusion_weights == 'False':
			tgfm_const_prior_twas_pickle_file = data[2]
		elif fusion_weights == 'True':
			tgfm_const_prior_twas_pickle_file = data[3]
		#tgfm_const_prior_fine_mapping_file = data[3]
		tgfm_const_prior_rss_regression_fine_mapping_file = data[4]
		#tgfm_const_prior_rss_prediction_fine_mapping_file = data[7]
		#tgfm_const_prior_rss_regression_nominal_fine_mapping_file = data[9]
		if tgfm_const_prior_twas_pickle_file == 'NA':
			continue


		twas_data_file = data[-1]
		gwas_component_susie_pmces_file = ukbb_genome_wide_susie_organized_results_dir + trait_name + '_' + component_name + '_component_posterior_mean_causal_effect_sizes.txt'
		f = open(twas_data_file, "rb")
		twas_data = pickle.load(f)
		f.close()
		gwas_component_susie_pmces, gwas_susie_mu, gwas_susie_mu_sd, gwas_susie_alpha = extract_gwas_component_susie_pmces_in_same_order_as_twas_data(gwas_component_susie_pmces_file, twas_data['bim'])
		variant_overlap_prob = get_variant_overlap_prob_for_each_gene(gwas_susie_alpha, twas_data['susie_alpha'])
		variant_in_credible_set_prob = get_variant_in_credible_set_prob_for_each_gene(gwas_susie_alpha, twas_data['susie_alpha'])
		

		#Version 2: TGFM susie with no tissue prior
		if fusion_weights == 'False':
			tgfm_const_prior_susie_twas_pickle_file = tgfm_const_prior_twas_pickle_file.split('tgfm_twas_results.pkl')[0]+ 'tgfm_susie_twas_results.pkl'
			tgfm_tiss_prior_susie_twas_pickle_file = tgfm_const_prior_twas_pickle_file.split('tgfm_twas_results.pkl')[0]+ 'tgfm_susie_twas_tissue_specific_prior_results.pkl'
			robust_tgfm_const_prior_susie_twas_pickle_file = tgfm_const_prior_twas_pickle_file.split('tgfm_twas_results.pkl')[0]+ 'tgfm_robust_joint_susie_twas_not_tissue_specific_prior_results.pkl'
			robust_tgfm_tiss_prior_susie_twas_pickle_file = tgfm_const_prior_twas_pickle_file.split('tgfm_twas_results.pkl')[0]+ 'tgfm_robust_joint_susie_twas_tissue_specific_prior_results.pkl'
		elif fusion_weights == 'True':
			print('currently not implemented')
			pdb.set_trace()

		if os.path.isfile(tgfm_const_prior_susie_twas_pickle_file) == False:
			continue 


		# Version 2: TGFM with tissue prior
		if fusion_weights == "False":
			tgfm_tiss_prior_twas_pickle_file = tgfm_const_prior_twas_pickle_file.split('tgfm_twas_results.pkl')[0] + 'fusion_' + fusion_weights + '_tgfm_twas_tissue_specific_prior_results.pkl'
		elif fusion_weights == "True":
			tgfm_tiss_prior_twas_pickle_file = tgfm_const_prior_twas_pickle_file.split('tgfm_fusion_twas_results.pkl')[0] + 'fusion_' + fusion_weights + '_tgfm_twas_tissue_specific_prior_results.pkl'
		tgfm_tiss_prior_rss_regression_fine_mapping_file = tgfm_const_prior_rss_regression_fine_mapping_file.split('tgfm_rss_regression_fine_mapping_table.txt')[0] + 'fusion_' + fusion_weights + '_tgfm_rss_regression_fine_mapping_tissue_specific_prior_table.txt'


		# Version 4: Robust TGFM with tissue prior
		if fusion_weights == "False":
			robust_tgfm_tissue_prior_twas_pickle_file = tgfm_const_prior_twas_pickle_file.split('tgfm_twas_results.pkl')[0] + 'fusion_' + fusion_weights + '_tgfm_robust_twas_tissue_specific_prior_results.pkl'
		elif fusion_weights == "True":
			robust_tgfm_tissue_prior_twas_pickle_file = tgfm_const_prior_twas_pickle_file.split('tgfm_fusion_twas_results.pkl')[0] + 'fusion_' + fusion_weights + '_tgfm_robust_twas_tissue_specific_prior_results.pkl'
		robust_tgfm_tissue_prior_rss_regression_fine_mapping_file = tgfm_const_prior_rss_regression_fine_mapping_file.split('tgfm_rss_regression_fine_mapping_table.txt')[0]+ 'fusion_' + fusion_weights + '_tgfm_robust_rss_regression_fine_mapping_tissue_specific_prior_table.txt'


		# Extract TGFM TWAS data for this disease component
		tp_gene_names, tp_tissue_names, tp_twas_fusion_nominal_z_scores, tp_twas_tgfm_multivariate_z_scores, tp_twas_tgfm_multivariate_alpha_mu, tp_twas_tgfm_multivariate_alpha_var  = extract_data_from_tgfm_twas_pickle_file(tgfm_tiss_prior_twas_pickle_file)
		cp_gene_names, cp_tissue_names, cp_twas_fusion_nominal_z_scores, cp_twas_tgfm_multivariate_z_scores, cp_twas_tgfm_multivariate_alpha_mu, cp_twas_tgfm_multivariate_alpha_var = extract_data_from_tgfm_twas_pickle_file(tgfm_const_prior_twas_pickle_file)
		r_tp_gene_names, r_tp_tissue_names, r_tp_twas_fusion_nominal_z_scores, r_tp_twas_tgfm_multivariate_z_scores, r_tp_twas_tgfm_multivariate_alpha_mu, r_tp_twas_tgfm_multivariate_alpha_var  = extract_data_from_tgfm_twas_pickle_file(robust_tgfm_tissue_prior_twas_pickle_file)


		# Extract TGFM-susie TWAS data from this disease component
		#cp_susie_gene_names, cp_susie_tissue_names, cp_susie_pip, num_cp_susie_components = extract_data_from_tgfm_susie_twas_pickle_file(tgfm_const_prior_susie_twas_pickle_file)
		#cp_robust_susie_gene_names, cp_robust_susie_tissue_names, cp_robust_susie_pip, num_cp_robust_susie_components = extract_data_from_tgfm_susie_twas_pickle_file(robust_tgfm_const_prior_susie_twas_pickle_file)
		#tp_susie_gene_names, tp_susie_tissue_names, tp_susie_pip, num_tp_susie_components = extract_data_from_tgfm_susie_tissue_prior_twas_pickle_file(tgfm_tiss_prior_susie_twas_pickle_file)
		#tp_robust_susie_gene_names, tp_robust_susie_tissue_names, tp_robust_susie_pip, num_tp_robust_susie_components = extract_data_from_tgfm_joint_susie_tissue_prior_twas_pickle_file(robust_tgfm_tiss_prior_susie_twas_pickle_file, twas_data)
		tp_robust_joint_susie_gene_names, tp_robust_joint_susie_tissue_names, tp_robust_joint_susie_gene_prob, tp_robust_joint_susie_non_med_prob = extract_data_from_tgfm_joint_susie_tissue_prior_twas_pickle_file(robust_tgfm_tiss_prior_susie_twas_pickle_file, twas_data, gwas_susie_alpha)

		manhatten_plot_data_file = pseudotissue_gtex_rss_multivariate_twas_dir + trait_name + '_' + gene_version + '_fusion_weights_' + fusion_weights + '_' + component_name + '_manhatten_plot_info.txt'
		extract_manhatten_plot_data_from_tgfm_joint_susie_tissue_prior(robust_tgfm_tiss_prior_susie_twas_pickle_file, twas_data, manhatten_plot_data_file)

		gene_tissue_to_coloc_pph4 = create_gene_tissue_to_pph4_mapping(np.unique(tp_robust_joint_susie_gene_names), coloc_results_dir, trait_name)

		#t2.write(component_name + '\t' + str(num_cp_susie_components) + '\t' + str(num_cp_robust_susie_components) + '\t' + str(num_tp_susie_components) + '\t' + str(num_tp_robust_susie_components) + '\n')

		# Quick error checking
		if np.array_equal(tp_gene_names, cp_gene_names) == False or np.array_equal(tp_tissue_names, cp_tissue_names) == False:
			print('assumption error')
			pdb.set_trace()

		# Extract TGFM fine mapping data for this disease component
		tp_sum_post_prob_rss_regression_dicti = extract_data_from_tgfm_fine_mapping_table(tgfm_tiss_prior_rss_regression_fine_mapping_file)


		cp_sum_post_prob_rss_regression_dicti = extract_data_from_tgfm_fine_mapping_table(tgfm_const_prior_rss_regression_fine_mapping_file)
	
		r_tp_sum_post_prob_rss_regression_dicti, r_tp_log_sum_exp_likelihood_dicti = extract_data_from_tgfm_fine_mapping_table_and_log_like(robust_tgfm_tissue_prior_rss_regression_fine_mapping_file)

		# Extract Unique genes
		unique_genes = np.unique(cp_gene_names)

		#adapive_coloc_prob_dicti = get_coloc_prob_dicti_for_set_of_genes(coloc_results_dir, unique_genes, trait_name, adaptive=True)
		#coloc_prob_dicti = get_coloc_prob_dicti_for_set_of_genes(coloc_results_dir, unique_genes, trait_name, adaptive=False)

		for g_index in range(len(tp_gene_names)):
			tissue_gene_name = cp_tissue_names[g_index] + '_' + cp_gene_names[g_index]
			gene_tissue_name = cp_gene_names[g_index] + '_' + cp_tissue_names[g_index]
			t.write(component_name + '\t' + cp_tissue_names[g_index] + '\t' + cp_gene_names[g_index] + '\t')
			t.write(str(cp_twas_fusion_nominal_z_scores[g_index]) + '\t')


			t.write(str(cp_twas_tgfm_multivariate_z_scores[g_index]) + '\t')
			t.write(str(tp_robust_joint_susie_gene_prob[g_index]) + '\t')

			t.write(str(tp_twas_tgfm_multivariate_z_scores[g_index]) + '\t')
			t.write(str(r_tp_twas_tgfm_multivariate_z_scores[g_index]) + '\t')


			t.write(str(cp_sum_post_prob_rss_regression_dicti[tissue_gene_name]) + '\t')
			t.write(str(tp_sum_post_prob_rss_regression_dicti[tissue_gene_name]) + '\t')
			t.write(str(r_tp_sum_post_prob_rss_regression_dicti[tissue_gene_name]) + '\t')
			t.write(str(r_tp_log_sum_exp_likelihood_dicti[tissue_gene_name]) + '\t')
			t.write(str(variant_in_credible_set_prob[g_index]) + '\t')
			t.write(str(variant_overlap_prob[g_index]) + '\t')
			t.write(gene_tissue_to_coloc_pph4[gene_tissue_name] + '\t')
			t.write(manhatten_plot_data_file + '\n')

	f.close()
t.close()
#t2.close()
