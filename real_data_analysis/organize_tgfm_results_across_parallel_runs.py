import numpy as np 
import os
import sys
import pdb
import pickle
import scipy.stats

def extract_trait_names(trait_names_file):
	arr = []
	f = open(trait_names_file)
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

def concatenate_results_across_parallel_jobs(file_stem, suffix, num_jobs, concat_output_file):
	#concat_output_file = file_stem + '_' +  suffix
	t = open(concat_output_file,'w')
	for job_number in range(num_jobs):
		job_input_file = file_stem + '_' + str(job_number) + '_' + str(num_jobs) + '_' + suffix
		f = open(job_input_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			if head_count == 0:
				head_count = head_count + 1
				if job_number == 0:
					t.write(line + '\n')
				continue
			t.write(line + '\n')
		f.close()
	t.close()

def extract_causal_effect_size(window_name, component_index, top_element, top_element_type, trait_name, tgfm_results_dir):
	window_pickle_file = tgfm_results_dir + trait_name + '_' + window_name + '_results.pkl'
	f = open(window_pickle_file, "rb")
	tgfm_res_window = pickle.load(f)
	f.close()	
	if top_element_type == 'variant':
		variant_index = np.where(tgfm_res_window['variants'] == top_element)[0]
		if len(variant_index) != 1:
			print('assumption eroror')
			pdb.set_trace()
		variant_index = variant_index[0]
		causal_effect_size = tgfm_res_window['beta_mu'][component_index,variant_index]
	elif top_element_type == 'gene':
		gene_index = np.where(tgfm_res_window['genes'] == top_element)[0]
		if len(gene_index) != 1:
			print('assumption eroror')
			pdb.set_trace()
		gene_index = gene_index[0]
		causal_effect_size = tgfm_res_window['alpha_mu'][component_index,gene_index]

		###########################
		'''
		eqtl_effect_file = '/n/scratch3/users/b/bes710/causal_eqtl_gwas/gtex/preprocessed_tgfm_data/' + 'cis_heritable_gene_' + window_name + '_tgfm_trait_agnostic_data_obj.pkl'
		f = open(eqtl_effect_file, "rb")
		eqtl_obj = pickle.load(f)
		f.close()
		eqtl_susie_mu = eqtl_obj['susie_mu'][gene_index]
		eqtl_susie_alpha = eqtl_obj['susie_alpha'][gene_index]
		eqtl_pmces = np.sum(eqtl_susie_alpha*eqtl_susie_mu,axis=0)
		eqtl_pmces2 = np.sum(eqtl_susie_alpha[0,:]*eqtl_susie_mu[0,:],axis=0)
		gene_var = np.dot(np.dot(eqtl_pmces, eqtl_obj['reference_ld']), eqtl_pmces)
		gene_std = np.sqrt(gene_var)
		pdb.set_trace()
		'''
	else:
		print('currently not implemented')
		pdb.set_trace()

	return causal_effect_size

def get_windows_with_a_large_number_of_components(cs_input_file, max_comp=7):
	f = open(cs_input_file)
	windows = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]
		if window_name not in windows:
			windows[window_name] = 0
		windows[window_name] = windows[window_name] + 1
	f.close()

	window_arr = np.asarray([*windows])
	bad_windows = {}
	for window in window_arr:
		if windows[window] > max_comp:
			bad_windows[window] = 1
	return bad_windows




def create_predicted_causal_effect_size_file(cs_input_file, pred_causal_effect_size_file, trait_name, tgfm_results_dir):
	# Extract windows with large numbr of components
	bad_windows = get_windows_with_a_large_number_of_components(cs_input_file)
	# Open output file handle
	t = open(pred_causal_effect_size_file,'w')
	# Write header to output file
	t.write('window_name\tcomponent_index\tgene_mediated_probability\ttop_element\ttop_element_prob\ttop_element_type\tfine_mappedcausal_effect_size\n')
	a_var = []
	a_gene = []
	probs = []
	# Loop through components
	f = open(cs_input_file)
	# Skip header
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Extract relevent fields from line
		window_name = data[0]
		if window_name in bad_windows:
			continue
		component_index = int(data[1])
		gene_mediated_prob = data[2]
		cs_elements_str = data[3]
		top_element = cs_elements_str.split(';')[0]
		cs_prob_str = data[4]
		top_element_prob = float(cs_prob_str.split(';')[0])
		if top_element.startswith('ENSG'):
			top_element_type = 'gene'
		else:
			top_element_type = 'variant'
		# Extract causal effect size of top element in the component
		causal_effect_size = extract_causal_effect_size(window_name, component_index, top_element, top_element_type, trait_name, tgfm_results_dir)
		t.write(window_name + '\t' + str(component_index) + '\t' + gene_mediated_prob + '\t' + top_element + '\t' + str(top_element_prob) + '\t' + top_element_type + '\t' + str(causal_effect_size) + '\n')
		if top_element_prob > .2:
			if top_element_type == 'gene':
				a_gene.append(causal_effect_size)
			elif top_element_type == 'variant':
				a_var.append(causal_effect_size)
		probs.append(float(gene_mediated_prob))
	f.close()
	t.close()
	a_gene = np.asarray(a_gene)
	a_var = np.asarray(a_var)
	print(np.mean(np.abs(a_gene)))
	print(np.mean(np.abs(a_var)))
	print(np.mean(probs))


def print_trait_pips_to_all_pip_file(ttt, concatenated_pip_summary_file, data_file_stem, file_stem, model_version, processed_tgfm_input_stem, trait_name):
	f = open(concatenated_pip_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]

		# Load in TGFM results pkl file for this window
		window_pkl_file = data_file_stem + '_' + window_name + '_results.pkl'
		if os.path.isfile(window_pkl_file) == False:
			continue
		# Load in tgfm results data
		g = open(window_pkl_file, "rb")
		tgfm_results = pickle.load(g)
		g.close()

		# Load in tgfm results data
		g = open(processed_tgfm_input_stem + '_' + window_name + '_tgfm_trait_agnostic_input_data_obj.pkl', "rb")
		tgfm_data = pickle.load(g)
		g.close()

		# GENE-TISSUE
		# Add middle gene indices to tgfm results
		tgfm_results['middle_gene_indices'] = np.copy(tgfm_data['middle_gene_indices'])
		# Middle gene-tissue pips
		middle_gene_pips = tgfm_results['expected_alpha_pips'][tgfm_results['middle_gene_indices']]
		middle_gene_names = tgfm_results['genes'][tgfm_results['middle_gene_indices']]
		for ii, gene_tissue_name in enumerate(middle_gene_names):
			gene_tissue_pip = middle_gene_pips[ii]
			ttt.write(trait_name + '\t' + 'gene_tissue' + '\t' + gene_tissue_name + '\t' + str(gene_tissue_pip) + '\n')

		# VARIANTS
		tgfm_results['middle_variant_indices'] = np.copy(tgfm_data['middle_variant_indices'])
		middle_variant_names = tgfm_results['variants'][tgfm_results['middle_variant_indices']]
		middle_variant_pips = tgfm_results['expected_beta_pips'][tgfm_results['middle_variant_indices']]
		for ii, variant_name in enumerate(middle_variant_names):
			variant_pip = middle_variant_pips[ii]
			ttt.write(trait_name + '\t' + 'variant' + '\t' + variant_name + '\t' + str(variant_pip) + '\n')


		# GENES
		middle_gene_pips = tgfm_results['expected_alpha_pips'][tgfm_results['middle_gene_indices']]
		middle_gene_names = tgfm_results['genes'][tgfm_results['middle_gene_indices']]

		alphas = tgfm_results['alpha_phis']
		middle_genes_indices = tgfm_results['middle_gene_indices']
		middle_gene_tissues = tgfm_results['genes'][middle_genes_indices]

		gene_to_indices = {}
		for ii, middle_gene_tissue in enumerate(middle_gene_tissues):
			gene_name = middle_gene_tissue.split('_')[0].split('.')[0]
			if gene_name not in gene_to_indices:
				gene_to_indices[gene_name] = []
			gene_to_indices[gene_name].append(ii)
		n_genes = len(gene_to_indices)
		ordered_genes = [*gene_to_indices]
		LL = len(alphas)
		gene_alphas = []
		for ll in range(LL):
			new_alpha = np.zeros((alphas[ll].shape[0], n_genes))
			for ii,gene_name in enumerate(ordered_genes):
				tmper = alphas[ll][:, middle_genes_indices]
				new_alpha[:, ii] = np.sum(tmper[:, gene_to_indices[gene_name]],axis=1)
			gene_alphas.append(new_alpha)
		expected_alpha_pips = compute_expected_pips_from_sampler_pis(gene_alphas)

		for ii, gene_name in enumerate(ordered_genes):
			agg_pip = expected_alpha_pips[ii]
			ttt.write(trait_name + '\t' + 'gene' + '\t' + gene_name + '\t' + str(agg_pip) + '\n')
	return ttt


def generate_per_gene_tissue_pip_full_summary_file(concatenated_pip_summary_file, per_gene_tissue_full_pip_summary_file, tissue_name_to_broad_category, data_file_stem, file_stem, model_version, processed_tgfm_input_stem, trait_name):
	f = open(concatenated_pip_summary_file)
	t = open(per_gene_tissue_full_pip_summary_file,'w')
	t.write('gene_tissue_name\tgene_name\ttissue_name\ttissue_visualization_category\twindow_name\tPIP\ttwas_p\tmin_window_gwas_p\n')
	head_count = 0

	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]

		# Load in TGFM results pkl file for this window
		window_pkl_file = data_file_stem + '_' + window_name + '_results.pkl'
		if os.path.isfile(window_pkl_file) == False:
			continue
		# Load in tgfm results data
		g = open(window_pkl_file, "rb")
		tgfm_results = pickle.load(g)
		g.close()

		# Load in tgfm results data
		g = open(processed_tgfm_input_stem + '_' + window_name + '_tgfm_trait_agnostic_input_data_obj.pkl', "rb")
		tgfm_data = pickle.load(g)
		g.close()
		# Load in tgfm trait data
		g = open(processed_tgfm_input_stem + '_' + window_name + '_tgfm_ukbb_data_obj.pkl', "rb")
		tgfm_trait_data = pickle.load(g)
		g.close()

		# Get gwas pvalues 
		trait_index = np.where(tgfm_trait_data['gwas_study_names'] == trait_name)[0][0]
		# Extract gwas betas and standard errrors
		gwas_beta = tgfm_trait_data['gwas_beta'][trait_index,:]
		gwas_beta_se = tgfm_trait_data['gwas_beta_se'][trait_index,:]
		gwas_z = gwas_beta/gwas_beta_se
		gwas_p = scipy.stats.norm.sf(abs(gwas_z))*2.0
		min_gwas_p = np.min(gwas_p)

		# Get TWAS pvalues
		twas_z_mat = np.asarray(tgfm_results['nominal_twas_z'])
		twas_p_mat = scipy.stats.norm.sf(abs(twas_z_mat))*2
		twas_pvalue = np.median(twas_p_mat,axis=0)


		# Add middle gene indices to tgfm results
		tgfm_results['middle_gene_indices'] = np.copy(tgfm_data['middle_gene_indices'])
		tgfm_results['middle_variant_indices'] = np.copy(tgfm_data['middle_variant_indices'])

		# Middle gene pips
		middle_gene_pips = tgfm_results['expected_alpha_pips'][tgfm_results['middle_gene_indices']]
		middle_gene_names = tgfm_results['genes'][tgfm_results['middle_gene_indices']]
		middle_gene_twas_p = twas_pvalue[tgfm_results['middle_gene_indices']]


		for ii, gene_tissue_name in enumerate(middle_gene_names):
			gene_tissue_pip = middle_gene_pips[ii]
			gene_name = gene_tissue_name.split('_')[0]
			tissue_name = '_'.join(gene_tissue_name.split('_')[1:])
			tissue_category = tissue_name_to_broad_category[tissue_name]
			t.write(gene_tissue_name + '\t' + gene_name + '\t' + tissue_name + '\t' + tissue_category + '\t' + window_name + '\t' + str(gene_tissue_pip) + '\t' + str(middle_gene_twas_p[ii]) + '\t' + str(min_gwas_p) + '\n')
	f.close()
	t.close()
	return


def generate_per_gene_tissue_pip_summary_file(concatenated_pip_summary_file, per_gene_pip_summary_file, tissue_name_to_broad_category):
	f = open(concatenated_pip_summary_file)
	t = open(per_gene_pip_summary_file,'w')
	t.write('gene_tissue_name\tgene_name\ttissue_name\ttissue_visualization_category\twindow_name\tPIP\n')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]
		if len(data) != 3:
			continue
		ele_names = np.asarray(data[1].split(';'))
		pips = np.asarray(data[2].split(';')).astype(float)
		for ii, ele_name in enumerate(ele_names):
			if ele_name.startswith('ENSG') == False:
				continue
			gene_name = ele_name.split('_')[0]
			tissue_name = '_'.join(ele_name.split('_')[1:])
			t.write(ele_name + '\t' + gene_name + '\t' + tissue_name + '\t' + tissue_name_to_broad_category[tissue_name] + '\t' + window_name + '\t' + str(pips[ii]) + '\n')
	f.close()
	t.close()
	return

def get_tissue_from_gene_tissue_pair(gene_tissue_names):
	tissue_names = []
	for gt in gene_tissue_names:
		tissue = '_'.join(gt.split('_')[1:])
		tissue_names.append(tissue)
	return np.asarray(tissue_names)

def extract_tissue_category_pips_from_tgfm_sampler_obj(tgfm_results, tissue_categories, tissue_category_to_tissue_names):
	middle_gene_tissue_names = tgfm_results['genes'][tgfm_results['middle_gene_indices']]
	middle_tissues = get_tissue_from_gene_tissue_pair(middle_gene_tissue_names)

	if len(middle_tissues) == 0:
		return np.zeros(len(tissue_categories))

	n_samples = tgfm_results['alpha_pips'].shape[0]
	n_components =len(tgfm_results['alpha_phis'])
	expected_tissue_category_pips = []
	for tissue_category	in tissue_categories:
		tissue_names = tissue_category_to_tissue_names[tissue_category]
		tissue_category_indices = []
		for tissue_name in tissue_names:
			tissue_indices = np.where(middle_tissues==tissue_name)[0]
			if len(tissue_indices) == 0.0:
				continue
			tissue_category_indices.append(tissue_indices)
		if len(tissue_category_indices) == 0.0:
			expected_tissue_category_pips.append(0.0)
			continue
		tissue_category_indices = np.hstack(tissue_category_indices)

		if len(tissue_category_indices) == 0.0:
			expected_tissue_category_pips.append(0.0)
			continue
		anti_pips = np.ones(n_samples)
		for component_iter in range(n_components):
			anti_pips = anti_pips*(1.0 - np.sum(tgfm_results['alpha_phis'][component_iter][:, tgfm_results['middle_gene_indices']][:, tissue_category_indices],axis=1))
		pips = 1.0 - anti_pips
		expected_tissue_category_pips.append(np.mean(pips))
	return np.asarray(expected_tissue_category_pips)

def extract_tissue_pips_from_tgfm_sampler_obj(tgfm_results, tissue_names):
	middle_gene_tissue_names = tgfm_results['genes'][tgfm_results['middle_gene_indices']]
	middle_tissues = get_tissue_from_gene_tissue_pair(middle_gene_tissue_names)

	if len(middle_tissues) == 0:
		return np.zeros(len(tissue_names))

	n_samples = tgfm_results['alpha_pips'].shape[0]
	n_components =len(tgfm_results['alpha_phis'])
	expected_tissue_pips = []
	for tissue_name in tissue_names:
		tissue_indices = np.where(middle_tissues==tissue_name)[0]

		if len(tissue_indices) == 0.0:
			expected_tissue_pips.append(0.0)
			continue

		anti_pips = np.ones(n_samples)
		for component_iter in range(n_components):
			anti_pips = anti_pips*(1.0 - np.sum(tgfm_results['alpha_phis'][component_iter][:, tgfm_results['middle_gene_indices']][:, tissue_indices],axis=1))
		pips = 1.0 - anti_pips
		expected_tissue_pips.append(np.mean(pips))
	return np.asarray(expected_tissue_pips)

def extract_tissue_category_pips_from_tgfm_pmces_obj(tgfm_results, tissue_categories, tissue_category_to_tissue_names):
	middle_gene_tissue_names = tgfm_results['genes'][tgfm_results['middle_gene_indices']]
	middle_tissues = get_tissue_from_gene_tissue_pair(middle_gene_tissue_names)
	middle_alpha_phi = tgfm_results['alpha_phi'][:, tgfm_results['middle_gene_indices']]	

	tissue_pips = []

	if len(middle_tissues) == 0:
		return np.zeros(len(tissue_categories))
	for tissue_category	in tissue_categories:
		tissue_names = tissue_category_to_tissue_names[tissue_category]
		tissue_category_indices = []
		for tissue_name in tissue_names:
			tissue_indices = np.where(middle_tissues==tissue_name)[0]
			if len(tissue_indices) == 0.0:
				continue
			tissue_category_indices.append(tissue_indices)
		if len(tissue_category_indices) == 0.0:
			tissue_pips.append(0.0)
			continue
		tissue_category_indices = np.hstack(tissue_category_indices)

		if len(tissue_category_indices) == 0.0:
			tissue_pips.append(0.0)
			continue

		component_pis = np.sum(middle_alpha_phi[:, tissue_category_indices],axis=1)

		anti_pip = 1.0
		for component_pi in component_pis:
			anti_pip = anti_pip*(1.0 - component_pi)
		pip = 1.0 - anti_pip
		tissue_pips.append(pip)
	return np.asarray(tissue_pips)


def extract_tissue_pips_from_tgfm_pmces_obj(tgfm_results, tissue_names):
	middle_gene_tissue_names = tgfm_results['genes'][tgfm_results['middle_gene_indices']]
	middle_tissues = get_tissue_from_gene_tissue_pair(middle_gene_tissue_names)
	middle_alpha_phi = tgfm_results['alpha_phi'][:, tgfm_results['middle_gene_indices']]

	tissue_pips = []

	if len(middle_tissues) == 0:
		return np.zeros(len(tissue_names))
	for tissue_name in tissue_names:
		tissue_indices = np.where(middle_tissues==tissue_name)[0]

		if len(tissue_indices) == 0.0:
			tissue_pips.append(0.0)
			continue

		component_pis = np.sum(middle_alpha_phi[:, tissue_indices],axis=1)

		anti_pip = 1.0
		for component_pi in component_pis:
			anti_pip = anti_pip*(1.0 - component_pi)
		pip = 1.0 - anti_pip
		tissue_pips.append(pip)
	return np.asarray(tissue_pips)

def generate_per_tissue_category_pip_summary_file(concatenated_pip_summary_file, per_tissue_category_pip_summary_file, tissue_names,tissue_categories, tissue_category_to_tissue_names, data_file_stem, file_stem, model_version, processed_tgfm_input_stem):
	f = open(concatenated_pip_summary_file)
	t = open(per_tissue_category_pip_summary_file,'w')
	t.write('tissue_category_name\twindow_name\tPIP\n')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]
		if len(data) != 3:
			continue
		# Load in TGFM results pkl file for this window
		window_pkl_file = data_file_stem + '_' + window_name + '_results.pkl'
		# Load in tgfm results data
		g = open(window_pkl_file, "rb")
		tgfm_results = pickle.load(g)
		g.close()

		# Load in tgfm results data
		g = open(processed_tgfm_input_stem + '_' + window_name + '_tgfm_trait_agnostic_input_data_obj.pkl', "rb")
		tgfm_data = pickle.load(g)
		g.close()
		# Add middle gene indices to tgfm results
		tgfm_results['middle_gene_indices'] = np.copy(tgfm_data['middle_gene_indices'])
		tgfm_results['middle_variant_indices'] = np.copy(tgfm_data['middle_variant_indices'])

		if model_version.startswith('susie_pmces_'):
			tissue_cat_pips = extract_tissue_category_pips_from_tgfm_pmces_obj(tgfm_results, tissue_categories, tissue_category_to_tissue_names)
		elif model_version.startswith('susie_sampler_'):
			tissue_cat_pips = extract_tissue_category_pips_from_tgfm_sampler_obj(tgfm_results, tissue_categories, tissue_category_to_tissue_names)
		else:
			print('assumption eroror: unknown model verison')
			pdb.set_trace()

		for tt, tissue_cat_name in enumerate(tissue_categories):
			t.write(tissue_cat_name + '\t' + window_name + '\t' + str(tissue_cat_pips[tt]) + '\n')
			'''
			if model_version.startswith('susie_sampler_'):
				if tissue_cat_pips[tt] > .4:
					print(tissue_cat_name + '\t' + str(tissue_cat_pips[tt]))
			'''
	f.close()
	t.close()
	return


def fill_in_causal_effect_size_matrix(null_mat, bs_eqtls_pmces_sparse):
	null_mat[bs_eqtls_pmces_sparse[:,0].astype(int), bs_eqtls_pmces_sparse[:,1].astype(int)] = bs_eqtls_pmces_sparse[:,2]
	return null_mat

def extract_extensive_per_gene_tissue_info_for_a_specific_high_sampler_tgfm_pip_gene(tgfm_results, gene_name, z_scores, sparse_eqtl_pmces,tgfm_data):
	# Extract index corresponding to this gene
	index_vec = np.where(tgfm_results['genes']==gene_name)[0]
	if len(index_vec) != 1:
		print('assumption eroorr')
		pdb.set_trace()
	gene_index = index_vec[0]

	avg_component_pis = []
	for component_iter in range(len(tgfm_results['alpha_phis'])):
		avg_component_pi = np.mean(tgfm_results['alpha_phis'][component_iter][:, gene_index])
		avg_component_pis.append(avg_component_pi)
	avg_component_pis = np.asarray(avg_component_pis)
	top_component = np.argmax(avg_component_pis)

	gene_eqtl_pmces = np.zeros((len(tgfm_data['genes']), len(tgfm_data['variants'])))
	twas_zs = []
	for itera in range(len(sparse_eqtl_pmces)):
		gene_eqtl_pmces = gene_eqtl_pmces*0.0
		gene_eqtl_pmces = fill_in_causal_effect_size_matrix(gene_eqtl_pmces, sparse_eqtl_pmces[itera])
		sampled_z = np.dot(z_scores, gene_eqtl_pmces[gene_index,:])
		twas_zs.append(sampled_z)
	twas_zs = np.asarray(twas_zs)

	pip_distr = tgfm_results['alpha_pips'][:,gene_index]
	pip_standard_dev = np.std(pip_distr)

	lb = np.sort(pip_distr)[2]
	ub = np.sort(pip_distr)[-3]

	return top_component, np.mean(twas_zs), pip_standard_dev, lb, ub

def extract_extensive_per_gene_tissue_info_for_a_specific_high_pmces_tgfm_pip_gene(tgfm_results, gene_name, z_scores, eqtl_pmces):
	# Extract index corresponding to this gene
	index_vec = np.where(tgfm_results['genes']==gene_name)[0]
	if len(index_vec) != 1:
		print('assumption eroorr')
		pdb.set_trace()
	gene_index = index_vec[0]

	top_component = np.argmax(tgfm_results['alpha_phi'][:,gene_index])

	twas_z = np.dot(eqtl_pmces[gene_index,:], z_scores)

	pip_standard_dev = 0

	gene_pip = tgfm_results['alpha_pip'][gene_index]

	return top_component, twas_z, pip_standard_dev, gene_pip, gene_pip


def generate_extensive_per_gene_tissue_pip_summary_file(per_gene_tissue_pip_summary_file, per_gene_tissue_high_pip_extensive_summary_file, tissue_names, file_stem, model_version, processed_tgfm_input_stem,pip_threshold, ukbb_preprocessed_for_genome_wide_susie_dir, trait_name, per_gene_tissue_high_pip_high_confidence_extensive_summary_file):

	f = open(per_gene_tissue_pip_summary_file)
	t2 = open(per_gene_tissue_high_pip_high_confidence_extensive_summary_file,'w')
	t = open(per_gene_tissue_high_pip_extensive_summary_file,'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\t' + 'component_mode\ttwas_z\tpip_standard_deviation\tpip_lb\tpip_ub\n')
			t2.write(line + '\t' + 'component_mode\ttwas_z\tpip_standard_deviation\tpip_lb\tpip_ub\n')
			continue
		full_gene_name = data[0]
		gene_name = data[1]
		tissue_name = data[2]
		tissue_category = data[3]
		window_name = data[4]
		pip = float(data[5])
		if pip < pip_threshold:
			continue
		# Load in TGFM results pkl file for this window
		window_pkl_file = file_stem + '_' + window_name + '_results.pkl'
		# Load in tgfm results data
		g = open(window_pkl_file, "rb")
		tgfm_results = pickle.load(g)
		g.close()

		# Load in tgfm results data
		g = open(processed_tgfm_input_stem + '_' + window_name + '_tgfm_trait_agnostic_input_data_obj.pkl', "rb")
		tgfm_data = pickle.load(g)
		g.close()
		g = open(processed_tgfm_input_stem + '_' + window_name + '_tgfm_ukbb_data_obj.pkl', "rb")
		tgfm_trait_data = pickle.load(g)
		g.close()

		# Get trait z-scores
		trait_index = np.where(tgfm_trait_data['gwas_study_names']==trait_name)[0][0]
		z_scores = tgfm_trait_data['gwas_beta'][trait_index,:]/tgfm_trait_data['gwas_beta_se'][trait_index,:]

		if model_version.startswith('susie_pmces_'):
			component_mode, twas_z, pip_standard_deviation, pip_lb, pip_ub = extract_extensive_per_gene_tissue_info_for_a_specific_high_pmces_tgfm_pip_gene(tgfm_results, full_gene_name, z_scores, tgfm_data['gene_eqtl_pmces'])
		elif model_version.startswith('susie_sampler_'):
			component_mode, twas_z, pip_standard_deviation, pip_lb, pip_ub = extract_extensive_per_gene_tissue_info_for_a_specific_high_sampler_tgfm_pip_gene(tgfm_results, full_gene_name, z_scores, tgfm_data['sparse_sampled_gene_eqtl_pmces'], tgfm_data)

		t.write(line + '\t' + str(component_mode) + '\t' + str(twas_z) + '\t' + str(pip_standard_deviation) + '\t' + str(pip_lb) + '\t' + str(pip_ub) + '\n')
		if pip_lb > .05:
			t2.write(line + '\t' + str(component_mode) + '\t' + str(twas_z) + '\t' + str(pip_standard_deviation) + '\t' + str(pip_lb) + '\t' + str(pip_ub) + '\n')

	f.close()
	t.close()
	t2.close()
	return

def generate_per_tissue_pip_summary_file(concatenated_pip_summary_file, per_tissue_pip_summary_file, tissue_names, file_stem, model_version, processed_tgfm_input_stem):
	f = open(concatenated_pip_summary_file)
	t = open(per_tissue_pip_summary_file,'w')
	t.write('tissue_name\twindow_name\tPIP\n')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]
		if len(data) != 3:
			continue
		# Load in TGFM results pkl file for this window
		window_pkl_file = file_stem + '_' + window_name + '_results.pkl'
		# Load in tgfm results data
		g = open(window_pkl_file, "rb")
		tgfm_results = pickle.load(g)
		g.close()

		# Load in tgfm results data
		g = open(processed_tgfm_input_stem + '_' + window_name + '_tgfm_trait_agnostic_input_data_obj.pkl', "rb")
		tgfm_data = pickle.load(g)
		g.close()
		# Add middle gene indices to tgfm results
		tgfm_results['middle_gene_indices'] = np.copy(tgfm_data['middle_gene_indices'])
		tgfm_results['middle_variant_indices'] = np.copy(tgfm_data['middle_variant_indices'])

		if model_version.startswith('susie_pmces_'):
			tissue_pips = extract_tissue_pips_from_tgfm_pmces_obj(tgfm_results, tissue_names)
		elif model_version.startswith('susie_sampler_'):
			tissue_pips = extract_tissue_pips_from_tgfm_sampler_obj(tgfm_results, tissue_names)
		else:
			print('assumption eroror: unknown model verison')
			pdb.set_trace()

		for tt, tissue_name in enumerate(tissue_names):
			t.write(tissue_name + '\t' + window_name + '\t' + str(tissue_pips[tt]) + '\n')
			#if tissue_pips[tt] > .5:
				#print(tissue_name + '\t' + str(tissue_pips[tt]))

	f.close()
	t.close()
	return

def extract_mediated_probs(alpha_phi_vec, n_tissues, tiss_to_position_mapping, tissue_gene_names):
	mediated_probs = np.zeros(n_tissues)
	for g_iter, tissue_gene_name in enumerate(tissue_gene_names):
		tissue_name = '_'.join(tissue_gene_name.split('_')[1:])
		mediated_probs[tiss_to_position_mapping[tissue_name]] = mediated_probs[tiss_to_position_mapping[tissue_name]] + alpha_phi_vec[g_iter]
	return mediated_probs

def component_in_middle_of_window(alpha_phi_vec, beta_phi_vec, middle_gene_indices_dicti, middle_variant_indices_dicti):
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

def component_in_middle_of_window_hack(alpha_phi_vec, middle_gene_indices_dicti, middle_variant_indices_dicti):
	booler = False
	# Gene wins
	if True:
		best_index = np.argmax(alpha_phi_vec)
		if best_index in middle_gene_indices_dicti:
			booler = True
	return booler


def generate_component_level_summary_data(concatenated_pip_summary_file, component_level_summary_file, tissue_names, data_file_stem, file_stem, model_version, processed_tgfm_input_stem):
	f = open(concatenated_pip_summary_file)
	t = open(component_level_summary_file,'w')
	t.write('window_name\tbootstrap_number\tcomponent_number\tnon_mediated_probability')
	tiss_to_position_mapping = {}
	for ii,tissue_name in enumerate(tissue_names):
		t.write('\t' + tissue_name)
		tiss_to_position_mapping[tissue_name] = ii
	t.write('\n')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]

		# Load in TGFM results pkl file for this window
		window_pkl_file = data_file_stem + '_' + window_name + '_results.pkl'
		# Load in tgfm results data
		g = open(window_pkl_file, "rb")
		tgfm_results = pickle.load(g)
		g.close()

		# Load in tgfm results data
		g = open(processed_tgfm_input_stem + '_' + window_name + '_tgfm_trait_agnostic_input_data_obj.pkl', "rb")
		tgfm_data = pickle.load(g)
		g.close()
		# Add middle gene indices to tgfm results
		tgfm_results['middle_gene_indices'] = np.copy(tgfm_data['middle_gene_indices'])
		tgfm_results['middle_variant_indices'] = np.copy(tgfm_data['middle_variant_indices'])

		middle_gene_indices_dicti = {}
		middle_variant_indices_dicti = {}
		for indexer in tgfm_results['middle_gene_indices']:
			middle_gene_indices_dicti[indexer] =1
		for indexer in tgfm_results['middle_variant_indices']:
			middle_variant_indices_dicti[indexer] = 1

		if 'valid_components' not in tgfm_results:
			continue

		if model_version.startswith('susie_pmces_'):
			for component_iter in tgfm_results['valid_components']:
				if component_in_middle_of_window(tgfm_results['alpha_phi'][component_iter, :], tgfm_results['beta_phi'][component_iter, :], middle_gene_indices_dicti, middle_variant_indices_dicti):
					non_mediated_prob = np.sum(tgfm_results['beta_phi'][component_iter,:])
					mediated_probs = extract_mediated_probs(tgfm_results['alpha_phi'][component_iter, :], len(tissue_names), tiss_to_position_mapping, tgfm_results['genes'])
					t.write(window_name + '\t' + 'NA' + '\t' + str(component_iter) + '\t' + str(non_mediated_prob) + '\t' + '\t'.join(mediated_probs.astype(str)) + '\n')
		elif model_version.startswith('susie_sampler_'):
			n_samples = len(tgfm_results['valid_components'])
			for sample_iter in range(n_samples):
				sample_valid_components = tgfm_results['valid_middle_components'][sample_iter]
				for component_iter in sample_valid_components:
					mediated_probs = extract_mediated_probs(tgfm_results['alpha_phis'][component_iter][sample_iter,:], len(tissue_names), tiss_to_position_mapping, tgfm_results['genes'])
					non_mediated_prob = 1.0 - np.sum(tgfm_results['alpha_phis'][component_iter][sample_iter,:])
					t.write(window_name + '\t' + str(sample_iter) + '\t' + str(component_iter) + '\t' + str(non_mediated_prob) + '\t' + '\t'.join(mediated_probs.astype(str)) + '\n')
					#if component_in_middle_of_window(tgfm_results['alpha_phis'][component_iter][sample_iter,:], tgfm_results['beta_phis'][component_iter][sample_iter,:], middle_gene_indices_dicti, middle_variant_indices_dicti):
					'''
					if component_in_middle_of_window_hack(tgfm_results['alpha_phis'][component_iter][sample_iter,:], middle_gene_indices_dicti, middle_variant_indices_dicti):
						#non_mediated_prob = np.sum(tgfm_results['beta_phis'][component_iter][sample_iter,:])
						mediated_probs = extract_mediated_probs(tgfm_results['alpha_phis'][component_iter][sample_iter,:], len(tissue_names), tiss_to_position_mapping, tgfm_results['genes'])
						non_mediated_prob = 1.0 - np.sum(tgfm_results['alpha_phis'][component_iter][sample_iter,:])
						t.write(window_name + '\t' + str(sample_iter) + '\t' + str(component_iter) + '\t' + str(non_mediated_prob) + '\t' + '\t'.join(mediated_probs.astype(str)) + '\n')
					'''
		else:
			print('assumption eroror: unknown model verison')
			pdb.set_trace()
	t.close()
	f.close()

	# Summarize results
	nm_probs = []
	mediated_probs = []
	f = open(component_level_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		nm_probs.append(float(data[3]))
		mediated_probs.append(np.asarray(data[4:]).astype(float))
	f.close()
	nm_probs = np.asarray(nm_probs)
	mediated_probs = np.asarray(mediated_probs)

	if model_version.startswith('susie_pmces_'):
		n_components = len(nm_probs)
	elif model_version.startswith('susie_sampler_'):
		n_components = len(nm_probs)/100.0

	avg_non_mediated_prob = np.mean(nm_probs)
	tissue_probs = np.sum(mediated_probs,axis=0)/np.sum(np.sum(mediated_probs,axis=0))

	print('N_components: ' + str(n_components))
	print('Avg non-mediated prob: ' + str(avg_non_mediated_prob))

	for tiss_iter in np.argsort(-tissue_probs):
		print(tissue_names[tiss_iter] + ': ' + str(tissue_probs[tiss_iter]))

	return


def extract_pseudotissue_names(gtex_pseudotissue_file,ignore_testis=False):
	f = open(gtex_pseudotissue_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0] == 'Testis' and ignore_testis:
			continue
		arr.append(data[0])
	f.close()
	return np.asarray(arr)


def generate_per_tissue_pip_summary_file_old(per_gene_pip_summary_file, pip_thresh, per_tissue_pip_summary_file, tissue_names):
	f = open(per_gene_pip_summary_file)
	t = open(per_tissue_pip_summary_file,'w')
	t.write('tissue_name\tn_expected_causal_genes\n')

	tiss_to_count = {}
	for tissue_name in tissue_names:
		tiss_to_count[tissue_name] = 0.0

	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissue_name = data[2]
		pip = float(data[4])
		if pip < pip_thresh:
			continue
		tiss_to_count[tissue_name] = tiss_to_count[tissue_name] + pip
	
	counts = []
	for tissue_name in tissue_names:
		t.write(tissue_name + '\t' + str(tiss_to_count[tissue_name]) + '\n')
		counts.append(tiss_to_count[tissue_name])
	counts = np.asarray(counts)
	f.close()
	t.close()

	#print(tissue_names[np.argsort(counts)])
	#print(counts[np.argsort(counts)]/np.sum(counts))


	return

def create_mapping_from_tissue_category_to_tissue_names(gtex_pseudotissue_category_file,ignore_testis=False):
	category_names = []
	mapping = {}
	mapping2 = {}
	mapping3 = {}
	alt_categories = []
	f = open(gtex_pseudotissue_category_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissue_name = data[0]
		category = data[-3]
		broad_category = data[-2]
		alt_category = data[-1]
		if tissue_name == 'Testis' and ignore_testis:
			continue
		mapping2[tissue_name] = broad_category
		mapping3[tissue_name] = alt_category
		category_names.append(category)
		if category not in mapping:
			mapping[category] = []
		mapping[category].append(tissue_name)
		alt_categories.append(alt_category)
	f.close()

	category_names = np.sort(np.unique(np.asarray(category_names)))

	for category_name in category_names:
		mapping[category_name] = np.asarray(mapping[category_name])

	return category_names, mapping, mapping2, mapping3, np.sort(np.unique(np.asarray(alt_categories)))


def correlate_per_tissue_h2_with_number_of_components_identified_in_each_tissue(component_level_summary_file, sldsc_mediated_h2_file, per_tissue_sldsc_comparison_summary_file, model_version, trait_name, tissue_names,tissue_name_to_broad_category):
	# Initialize data
	tissue_name_to_n_comp = {}
	tissue_name_to_n_comp_lb = {}
	tissue_name_to_n_comp_ub = {}
	tissue_name_to_h2 = {}
	tissue_name_to_h2_lb = {}
	tissue_name_to_h2_ub = {}
	for tissue in tissue_names:
		tissue_name_to_n_comp[tissue] = 0
		tissue_name_to_n_comp_lb[tissue] = 0
		tissue_name_to_n_comp_ub[tissue] = 0
		tissue_name_to_h2[tissue] = 0
		tissue_name_to_h2_lb[tissue] = 0
		tissue_name_to_h2_ub[tissue] = 0

	# Add SLDSC results
	f = open(sldsc_mediated_h2_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] not in tissue_name_to_h2:
			continue
		tissue_name = data[0]
		h2 = float(data[1])
		h2_se = float(data[2])
		h2_ub = h2 + 1.96*h2_se
		h2_lb = h2 - 1.96*h2_se
		tissue_name_to_h2[tissue_name] = h2
		tissue_name_to_h2_lb[tissue_name] = h2_lb
		tissue_name_to_h2_ub[tissue_name] = h2_ub
	f.close()

	# Add TGFM results
	tmp_data = np.loadtxt(component_level_summary_file, dtype=str,delimiter='\t')
	# Now compute average mediated probability (and standard error)
	tiss_probs = tmp_data[1:,4:].astype(float)
	tiss_frac = np.sum(tiss_probs,axis=0)/np.sum(tiss_probs)
	tiss_comp = np.sum(tiss_probs,axis=0)
	if model_version.startswith('susie_pmces_'):
		n_components = tiss_probs.shape[0]
	elif model_version.startswith('susie_sampler_'):
		n_components = (tiss_probs.shape[0])/100.0
	tiss_frac_se = np.sqrt((tiss_frac*(1.0-tiss_frac))/n_components)
	tiss_comp_se = tiss_frac_se*np.sum(tiss_probs)
	if model_version.startswith('susie_sampler_'):
		tiss_comp = tiss_comp/100.0
		tiss_comp_se = tiss_comp_se/100.0

	# Print to output
	t = open(per_tissue_sldsc_comparison_summary_file,'w')
	t.write('tissue\ttissue_category\tn_tgfm_components\tn_tgfm_components_lb\tn_tgfm_components_ub\ttglr_mediated_h2\ttglr_mediated_h2_lb\ttglr_mediated_h2_ub\n')
	comps = []
	h2s = []
	for ii,tissue in enumerate(tissue_names):
		t.write(tissue + '\t' + str(tissue_name_to_broad_category[tissue]) + '\t' + str(tiss_comp[ii]) + '\t' + str(tiss_comp[ii] - (1.96*tiss_comp_se[ii])) + '\t' + str(tiss_comp[ii] + (1.96*tiss_comp_se[ii])) + '\t' + str(tissue_name_to_h2[tissue]) + '\t' + str(tissue_name_to_h2_lb[tissue]) + '\t' + str(tissue_name_to_h2_ub[tissue]) + '\n')
		comps.append(tiss_comp[ii])
		h2s.append(tissue_name_to_h2[tissue])
	t.close()
	comps = np.asarray(comps)
	h2s = np.asarray(h2s)
	print(np.corrcoef(comps, h2s)[0,1])
	return

def correlate_number_of_high_pip_genes_with_per_tissue_sldsc_h2(tissue_names, per_gene_tissue_pip_summary_file, pip_threshold, sldsc_mediated_h2_file, per_tissue_sldsc_comparison_summary_file, tissue_name_to_broad_category):
	tissue_name_to_pip_count = {}
	tissue_name_to_h2 = {}
	tissue_name_to_h2_lb = {}
	tissue_name_to_h2_ub = {}
	for tissue in tissue_names:
		tissue_name_to_pip_count[tissue] = 0
		tissue_name_to_h2[tissue] = 0
		tissue_name_to_h2_lb[tissue] = 0
		tissue_name_to_h2_ub[tissue] = 0
	f = open(per_gene_tissue_pip_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissue_name = data[2]
		pip = float(data[5])
		if pip < pip_threshold:
			continue
		tissue_name_to_pip_count[tissue_name] = tissue_name_to_pip_count[tissue_name] + pip
	f.close()

	f = open(sldsc_mediated_h2_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] not in tissue_name_to_h2:
			continue
		tissue_name = data[0]
		h2 = float(data[1])
		h2_se = float(data[2])
		h2_ub = h2 + 1.96*h2_se
		h2_lb = h2 - 1.96*h2_se
		tissue_name_to_h2[tissue_name] = h2
		tissue_name_to_h2_lb[tissue_name] = h2_lb
		tissue_name_to_h2_ub[tissue_name] = h2_ub
	f.close()
	t = open(per_tissue_sldsc_comparison_summary_file,'w')
	t.write('tissue\ttissue_category\tn_tgfm_components\ttglr_mediated_h2\th2_ub\th2_lb\n')
	comps = []
	h2s = []
	for tissue in tissue_names:
		t.write(tissue + '\t' + str(tissue_name_to_broad_category[tissue]) + '\t' + str(tissue_name_to_pip_count[tissue]) + '\t' + str(tissue_name_to_h2[tissue]) + '\t' + str(tissue_name_to_h2_ub[tissue]) + '\t' + str(tissue_name_to_h2_lb[tissue]) + '\n')
		comps.append(tissue_name_to_pip_count[tissue])
		h2s.append(tissue_name_to_h2[tissue])
	t.close()
	comps = np.asarray(comps)
	h2s = np.asarray(h2s)
	print(np.corrcoef(comps, h2s)[0,1])
	return

def correlate_fraction_mediated_h2_with_fraction_of_components_mediated(component_level_summary_file, sldsc_fraction_mediated_h2_file, fraction_mediated_sldsc_comparison_summary_file, model_version, trait_name):
	# First extract non-mediated probabilities
	nm_probs = []
	head_count = 0
	f = open(component_level_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		nm_prob = float(data[3])
		nm_probs.append(nm_prob)
	f.close()
	nm_probs = np.asarray(nm_probs)
	med_probs = 1.0 - nm_probs

	# Now compute average mediated probability (and standard error)
	if model_version.startswith('susie_pmces_'):
		n_components = len(med_probs)
	elif model_version.startswith('susie_sampler_'):
		n_components = len(med_probs)/100.0
	avg_med_prob = np.mean(med_probs)
	se_avg_med_prob = np.sqrt((avg_med_prob*(1.0-avg_med_prob))/n_components)
	avg_med_prob_lb = avg_med_prob - 1.96*se_avg_med_prob
	avg_med_prob_ub = avg_med_prob + 1.96*se_avg_med_prob

	# Now extract fraction of heritability mediated by gene expression in any tissue
	tmp_data = np.loadtxt(sldsc_fraction_mediated_h2_file, dtype=str,delimiter='\t')
	tglr_fraction_med = float(tmp_data[1,0])
	tglr_fraction_med_se = float(tmp_data[1,1])
	tglr_fraction_med_lb = tglr_fraction_med - 1.96*tglr_fraction_med_se
	tglr_fraction_med_ub = tglr_fraction_med + 1.96*tglr_fraction_med_se

	# Write to output file
	t = open(fraction_mediated_sldsc_comparison_summary_file,'w')
	t.write('trait_name\ttgfm_fraction_med\ttgfm_fraction_med_lb\ttgfm_fraction_med_ub\ttglr_fraction_med\ttglr_fraction_med_lb\ttglr_fraction_med_ub\n')
	t.write(trait_name + '\t' + str(avg_med_prob) + '\t' + str(avg_med_prob_lb) + '\t' + str(avg_med_prob_ub) + '\t' + str(tglr_fraction_med) + '\t' + str(tglr_fraction_med_lb) + '\t' + str(tglr_fraction_med_ub) + '\n')
	t.close()
	return

def tally_number_of_causal_genetic_elements_sqrt_plot_input(concatenated_pip_summary_file, n_causal_genetic_elements_summary_cross_threshold_file, pip_windows):
	f = open(concatenated_pip_summary_file)
	gene_pips = []
	variant_pips = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 3:
			continue
		ele_names = data[1].split(';')
		ele_probs = np.asarray(data[2].split(';')).astype(float)
		n_ele = len(ele_names)
		for ele_iter in range(n_ele):
			ele_name = ele_names[ele_iter]
			ele_prob = ele_probs[ele_iter]
			if ele_name.startswith('ENSG'):
				gene_pips.append(ele_prob)
			else:
				variant_pips.append(ele_prob)
	f.close()
	gene_pips = np.asarray(gene_pips)
	total_genes = np.sum(gene_pips > .1)
	variant_pips = np.asarray(variant_pips)
	total_variants = np.sum(variant_pips > .1)
	max_gene_pip = np.max(gene_pips)
	max_variant_pip = np.max(variant_pips)
	prev_genes = 0.0
	t = open(n_causal_genetic_elements_summary_cross_threshold_file,'w')
	t.write('element_class\tPIP_threshold\tn_elements\n')

	gene_pips = gene_pips[gene_pips >=.1]
	total_count = 0
	for ii,gene_pip_upper in enumerate(pip_windows[:(len(pip_windows)-1)]):
		gene_pip_lower = pip_windows[ii+1]
		if gene_pip_upper == 1.0:
			gene_indices = (gene_pips <= gene_pip_upper) & (gene_pips >= gene_pip_lower)
			boundry = ']'
		else:
			gene_indices = (gene_pips < gene_pip_upper) & (gene_pips >= gene_pip_lower)
			boundry = ')'
		window_count = np.sum(gene_indices)
		prev_total_count = total_count
		total_count = total_count + window_count

		count = np.sqrt(total_count) - np.sqrt(prev_total_count)
		if count == 0:
			continue
		t.write('gene' + '\t' + '[' + str(gene_pip_lower) + ', ' + str(gene_pip_upper) + boundry + '\t' + str(count) + '\n')
	t.close()

	return



def tally_number_of_causal_genes_cross_pip_thresholds(per_gene_pip_summary_file, n_causal_genes_summary_cross_threshold_file):
	f = open(per_gene_pip_summary_file)
	gene_pips = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_name = data[0]
		gene_pip = float(data[2])
		gene_pips.append(gene_pip)
	f.close()
	gene_pips = np.asarray(gene_pips)
	total_genes = np.sum(gene_pips > .1)
	max_gene_pip = np.max(gene_pips)
	prev_genes = 0.0
	t = open(n_causal_genes_summary_cross_threshold_file,'w')
	t.write('element_class\tPIP_threshold\tn_elements\n')
	gene_pips = gene_pips[gene_pips >= .1]
	total_count = 0
	for gene_pip in np.sort(-gene_pips):
		prev_total_count = total_count
		total_count = total_count + 1

		count = np.sqrt(total_count) - np.sqrt(prev_total_count)

		t.write('gene' + '\t' + str(-gene_pip) + '\t' + str(count) + '\n')

	t.close()

	return

def tally_number_of_causal_variants_cross_pip_thresholds(concatenated_pip_summary_file, n_causal_variants_summary_cross_threshold_file):
	f = open(concatenated_pip_summary_file)
	gene_pips = []
	variant_pips = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 3:
			continue
		ele_names = data[1].split(';')
		ele_probs = np.asarray(data[2].split(';')).astype(float)
		n_ele = len(ele_names)
		for ele_iter in range(n_ele):
			ele_name = ele_names[ele_iter]
			ele_prob = ele_probs[ele_iter]
			if ele_name.startswith('ENSG'):
				gene_pips.append(ele_prob)
			else:
				variant_pips.append(ele_prob)
	f.close()
	gene_pips = np.asarray(gene_pips)
	total_genes = np.sum(gene_pips > .1)
	variant_pips = np.asarray(variant_pips)
	total_variants = np.sum(variant_pips > .1)
	max_gene_pip = np.max(gene_pips)
	max_variant_pip = np.max(variant_pips)
	prev_genes = 0.0
	t = open(n_causal_variants_summary_cross_threshold_file,'w')
	t.write('element_class\tPIP_threshold\tn_elements\n')
	variant_pips = variant_pips[variant_pips >= .1]
	total_count = 0
	for variant_pip in np.sort(-variant_pips):
		prev_total_count = total_count
		total_count = total_count + 1

		count = np.sqrt(total_count) - np.sqrt(prev_total_count)

		t.write('variant' + '\t' + str(-variant_pip) + '\t' + str(count) + '\n')

	t.close()

	return
def tally_number_of_causal_sc_gene_tissue_pairs_cross_pip_thresholds(concatenated_pip_summary_file, n_causal_gene_tissue_pairs_summary_cross_threshold_file, tissue_name_to_broad_category):
	f = open(concatenated_pip_summary_file)
	gene_pips = []
	variant_pips = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 3:
			continue
		ele_names = data[1].split(';')
		ele_probs = np.asarray(data[2].split(';')).astype(float)
		n_ele = len(ele_names)
		for ele_iter in range(n_ele):
			ele_name = ele_names[ele_iter]
			ele_prob = ele_probs[ele_iter]
			if ele_name.startswith('ENSG'):
				context_name = '_'.join(ele_name.split('_')[1:])
				broad_category = tissue_name_to_broad_category[context_name]
				if broad_category == 'sc_blood':
					gene_pips.append(ele_prob)
	f.close()
	gene_pips = np.asarray(gene_pips)
	total_genes = np.sum(gene_pips > .1)
	t = open(n_causal_gene_tissue_pairs_summary_cross_threshold_file,'w')
	t.write('element_class\tPIP_threshold\tn_elements\n')
	if total_genes == 0:
		t.write('gene\t' + str(0.2) + '\t' + str(0.0) + '\n')
	else:
		max_gene_pip = np.max(gene_pips)
		prev_genes = 0.0
		gene_pips = gene_pips[gene_pips >=.1]
		total_count = 0
		for gene_pip in np.sort(-gene_pips):
			prev_total_count = total_count
			total_count = total_count + 1
			count = np.sqrt(total_count) - np.sqrt(prev_total_count)
			t.write('gene' + '\t' + str(-gene_pip) + '\t' + str(count) + '\n')
	t.close()
	return

def tally_number_of_causal_sc_gene_tissue_pairs_cross_pip_thresholds_in_single_cell_type(concatenated_pip_summary_file, n_causal_gene_tissue_pairs_summary_cross_threshold_file, cell_type_name):
	f = open(concatenated_pip_summary_file)
	gene_pips = []
	variant_pips = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 3:
			continue
		ele_names = data[1].split(';')
		ele_probs = np.asarray(data[2].split(';')).astype(float)
		n_ele = len(ele_names)
		for ele_iter in range(n_ele):
			ele_name = ele_names[ele_iter]
			ele_prob = ele_probs[ele_iter]
			if ele_name.startswith('ENSG'):
				context_name = '_'.join(ele_name.split('_')[1:])
				if context_name == cell_type_name:
					gene_pips.append(ele_prob)
	f.close()
	gene_pips = np.asarray(gene_pips)
	total_genes = np.sum(gene_pips >= .2)
	t = open(n_causal_gene_tissue_pairs_summary_cross_threshold_file,'w')
	t.write('element_class\tPIP_threshold\tn_elements\n')
	if total_genes == 0:
		t.write('gene\t' + str(0.2) + '\t' + str(0.0) + '\n')
	else:
		max_gene_pip = np.max(gene_pips)
		prev_genes = 0.0
		gene_pips = gene_pips[gene_pips >=.2]
		total_count = 0
		for gene_pip in np.sort(-gene_pips):
			prev_total_count = total_count
			total_count = total_count + 1
			count = np.sqrt(total_count) - np.sqrt(prev_total_count)
			t.write('gene' + '\t' + str(-gene_pip) + '\t' + str(count) + '\n')
	t.close()
	return


def tally_number_of_causal_gene_tissue_pairs_cross_pip_thresholds(concatenated_pip_summary_file, n_causal_gene_tissue_pairs_summary_cross_threshold_file):
	f = open(concatenated_pip_summary_file)
	gene_pips = []
	variant_pips = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 3:
			continue
		ele_names = data[1].split(';')
		ele_probs = np.asarray(data[2].split(';')).astype(float)
		n_ele = len(ele_names)
		for ele_iter in range(n_ele):
			ele_name = ele_names[ele_iter]
			ele_prob = ele_probs[ele_iter]
			if ele_name.startswith('ENSG'):
				gene_pips.append(ele_prob)
			else:
				variant_pips.append(ele_prob)
	f.close()
	gene_pips = np.asarray(gene_pips)
	total_genes = np.sum(gene_pips > .1)
	variant_pips = np.asarray(variant_pips)
	total_variants = np.sum(variant_pips > .1)
	max_gene_pip = np.max(gene_pips)
	max_variant_pip = np.max(variant_pips)
	prev_genes = 0.0
	t = open(n_causal_gene_tissue_pairs_summary_cross_threshold_file,'w')
	t.write('element_class\tPIP_threshold\tn_elements\n')
	gene_pips = gene_pips[gene_pips >=.1]
	total_count = 0
	for gene_pip in np.sort(-gene_pips):
		prev_total_count = total_count
		total_count = total_count + 1

		count = np.sqrt(total_count) - np.sqrt(prev_total_count)

		t.write('gene' + '\t' + str(-gene_pip) + '\t' + str(count) + '\n')

	t.close()

	return

def tally_number_of_causal_genetic_elements(concatenated_pip_summary_file, n_causal_genetic_elements_summary_file, n_causal_genetic_elements_by_tissue_summary_file, pip_threshold, tissue_names):
	f = open(concatenated_pip_summary_file)
	counts_dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 3:
			continue
		ele_names = data[1].split(';')
		ele_probs = np.asarray(data[2].split(';')).astype(float)
		n_ele = len(ele_names)
		for ele_iter in range(n_ele):
			ele_name = ele_names[ele_iter]
			ele_prob = ele_probs[ele_iter]
			if ele_prob < pip_threshold:
				continue
			if ele_name.startswith('ENSG'):
				ele_class = '_'.join(ele_name.split('_')[1:])
			else:
				ele_class = 'variant'
			if ele_class not in counts_dicti:
				counts_dicti[ele_class] = 0
			counts_dicti[ele_class] = counts_dicti[ele_class] + 1
	f.close()
	if 'variant' not in counts_dicti:
		counts_dicti['variant'] = 0.0

	# Do tissue stratefied version
	t = open(n_causal_genetic_elements_by_tissue_summary_file, 'w')
	t.write('element_name\tcount\n')
	ele_name = 'variant'
	t.write(ele_name + '\t' + str(counts_dicti[ele_name]) + '\n')
	gene_counts = 0
	for tissue_name in tissue_names:
		if tissue_name in counts_dicti:
			gene_counts = gene_counts + counts_dicti[tissue_name]
			t.write(tissue_name + '\t' + str(counts_dicti[tissue_name]) + '\n')
		else:
			t.write(tissue_name + '\t' + '0.0\n')
	t.close()

	t = open(n_causal_genetic_elements_summary_file, 'w')
	t.write('element_name\tcount\n')
	ele_name = 'variant'
	t.write(ele_name + '\t' + str(counts_dicti[ele_name]) + '\n')
	t.write('gene\t' + str(gene_counts) + '\n')
	t.close()

	n_genes = gene_counts
	n_var = counts_dicti['variant']
	#print('total: ' + str(n_genes+n_var))
	#print('fraction:' + str(n_genes/(n_genes+n_var)))

	return n_genes, n_genes+n_var

def compute_tgfm_tissue_overlap_jaccard(concatenated_pip_summary_file, tissue_overlap_jaccard_file, tissue_names, pip_threshold=.1):
	# tissue info
	n_tissues = len(tissue_names)
	tissue_to_index = {}
	for tissue_index, tissue_name in enumerate(tissue_names):
		tissue_to_index[tissue_name] = tissue_index

	# Intersection
	intersection = np.zeros((n_tissues, n_tissues))
	# union
	marginal = np.zeros(n_tissues)

	f = open(concatenated_pip_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 3:
			continue
		ele_names = data[1].split(';')
		ele_probs = np.asarray(data[2].split(';')).astype(float)
		n_ele = len(ele_names)

		window_tissues = []
		for ele_iter in range(n_ele):
			ele_name = ele_names[ele_iter]
			ele_prob = ele_probs[ele_iter]
			if ele_prob < pip_threshold:
				continue
			if ele_name.startswith('ENSG'):
				tissue_name = '_'.join(ele_name.split('_')[1:])
				window_tissues.append(tissue_name)
		if len(window_tissues) == 0:
			continue

		window_tissues = np.sort(np.unique(window_tissues))
		n_window_tissues = len(window_tissues)
		for window_tissue in window_tissues:
			tissue_index = tissue_to_index[window_tissue]
			marginal[tissue_index] = marginal[tissue_index] + 1

		if len(window_tissues) <= 1:
			continue

		for tiss_iter in range(n_window_tissues):
			for tiss_iter2 in range(tiss_iter+1, n_window_tissues):

				tissue1_name = window_tissues[tiss_iter]
				tissue2_name = window_tissues[tiss_iter2]

				tissue1_index = tissue_to_index[tissue1_name]
				tissue2_index = tissue_to_index[tissue2_name]

				intersection[tissue1_index, tissue2_index] = intersection[tissue1_index, tissue2_index] + 1
				intersection[tissue2_index, tissue1_index] = intersection[tissue2_index, tissue1_index] + 1
	f.close()

	t = open(tissue_overlap_jaccard_file,'w')
	t.write('tissue1\ttissue2\tjaccard_index\n')
	for tiss_iter in range(n_tissues):
		for tiss_iter2 in range(tiss_iter+1, n_tissues):
			tissue1 = tissue_names[tiss_iter]
			tissue2 = tissue_names[tiss_iter2]
			if marginal[tiss_iter] == 0.0 or marginal[tiss_iter2] == 0.0:
				continue
			jacard_index = intersection[tiss_iter, tiss_iter2]/(marginal[tiss_iter] + marginal[tiss_iter2] - intersection[tiss_iter, tiss_iter2])
			t.write(tissue1 + '\t' + tissue2 + '\t' + str(jacard_index) + '\n')
	t.close()

	return


def compute_causal_tissue_pvalues(concatenated_pip_summary_file, causal_tissue_pvalue_file, tissue_names, pip_threshold=.1):
	f = open(concatenated_pip_summary_file)
	pvalue_dicti = {}
	for tissue_name in tissue_names:
		pvalue_dicti[tissue_name] = 1.0
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 3:
			continue
		ele_names = data[1].split(';')
		ele_probs = np.asarray(data[2].split(';')).astype(float)
		n_ele = len(ele_names)
		for ele_iter in range(n_ele):
			ele_name = ele_names[ele_iter]
			ele_prob = ele_probs[ele_iter]
			if ele_prob < pip_threshold:
				continue
			if ele_name.startswith('ENSG'):
				tissue_name = '_'.join(ele_name.split('_')[1:])
				pvalue_dicti[tissue_name] = pvalue_dicti[tissue_name]*(1.0 - ele_prob)
	f.close()		

	# Write to output
	t = open(causal_tissue_pvalue_file,'w')
	t.write('tissue_name\tpvalue\n')
	for tissue_name in tissue_names:
		t.write(tissue_name + '\t' + str(pvalue_dicti[tissue_name]) + '\n')
	t.close()
	return

def convert_from_genes_to_tissues(gene_names):
	tissue_names = []
	for gene_name in gene_names:
		tissue_name = '_'.join(gene_name.split('_')[1:])
		tissue_names.append(tissue_name)

	return np.asarray(tissue_names)

def tally_expected_number_of_causal_genetic_elements(concatenated_pip_summary_file, n_causal_genetic_elements_summary_file, n_causal_genetic_elements_by_tissue_summary_file, data_file_stem,file_stem, model_version, processed_tgfm_input_stem, tissue_names):
	f = open(concatenated_pip_summary_file)
	tiss_to_position_mapping = {}
	gene_tissue_counts = {}
	for ii,tissue_name in enumerate(tissue_names):
		#t.write('\t' + tissue_name)
		tiss_to_position_mapping[tissue_name] = ii
		gene_tissue_counts[tissue_name] = 0.0
	head_count = 0
	variant_counts = 0.0
	gene_counts = 0.0

	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]

		# Load in TGFM results pkl file for this window
		window_pkl_file = data_file_stem + '_' + window_name + '_results.pkl'
		# Load in tgfm results data
		g = open(window_pkl_file, "rb")
		tgfm_results = pickle.load(g)
		g.close()

		# Load in tgfm results data
		g = open(processed_tgfm_input_stem + '_' + window_name + '_tgfm_trait_agnostic_input_data_obj.pkl', "rb")
		tgfm_data = pickle.load(g)
		g.close()
		# Add middle gene indices to tgfm results
		tgfm_results['middle_gene_indices'] = np.copy(tgfm_data['middle_gene_indices'])
		tgfm_results['middle_variant_indices'] = np.copy(tgfm_data['middle_variant_indices'])

		# Beta PIPs
		beta_pips = tgfm_results['expected_beta_pips'][tgfm_results['middle_variant_indices']]
		# Alpha PIPs
		alpha_pips = tgfm_results['expected_alpha_pips'][tgfm_results['middle_gene_indices']]

		# Gene names
		middle_gene_names = tgfm_results['genes'][tgfm_results['middle_gene_indices']]
		middle_tissue_names = convert_from_genes_to_tissues(middle_gene_names)

		# Tally variant and gene counts
		variant_counts = variant_counts + np.sum(beta_pips)
		gene_counts = gene_counts + np.sum(alpha_pips)

		for ii, tissue_name in enumerate(middle_tissue_names):
			gene_tissue_counts[tissue_name] = gene_tissue_counts[tissue_name] + alpha_pips[ii]

	f.close()

	# Open output file handle and print results
	t = open(n_causal_genetic_elements_summary_file, 'w')
	t.write('element_name\tcount\n')
	t.write('variant' + '\t' + str(variant_counts) + '\n')
	t.write('gene\t' + str(gene_counts) + '\n')
	t.close()

	# Open tissue specific output file
	t = open(n_causal_genetic_elements_by_tissue_summary_file, 'w')
	t.write('element_name\tcount\n')
	t.write('variant' + '\t' + str(variant_counts) + '\n')
	for tissue_name in tissue_names:
		t.write(tissue_name + '\t' + str(gene_tissue_counts[tissue_name]) + '\n')
	t.close()

	print(gene_counts/(gene_counts+variant_counts))

	return

def compute_expected_pips_from_sampler_pis(gene_alphas):
	n_bs = gene_alphas[0].shape[0]
	n_genes = gene_alphas[0].shape[1]
	LL = len(gene_alphas)


	alpha_pips = np.ones((n_bs, n_genes))

	for component_iter in range(LL):
		alpha_pips = alpha_pips*(1.0 - gene_alphas[component_iter])

	alpha_pips = 1.0 - alpha_pips

	expected_alpha_pips = np.mean(alpha_pips,axis=0)
	return expected_alpha_pips


def generate_per_gene_tissue_group_pip_summary_file(concatenated_pip_summary_file, data_file_stem, file_stem, model_version, processed_tgfm_input_stem, tissue_group_arr):
	tissue_group_dicti = {}
	for tissue in tissue_group_arr:
		tissue_group_dicti[tissue] = 1

	f = open(concatenated_pip_summary_file)
	head_count = 0
	used_genes = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]
		if len(data) != 3:
			continue
		# Load in TGFM results pkl file for this window
		window_pkl_file = data_file_stem + '_' + window_name + '_results.pkl'
		# Load in tgfm results data
		g = open(window_pkl_file, "rb")
		tgfm_results = pickle.load(g)
		g.close()


		# Load in tgfm results data
		g = open(processed_tgfm_input_stem + '_' + window_name + '_tgfm_trait_agnostic_input_data_obj.pkl', "rb")
		tgfm_data = pickle.load(g)
		g.close()
		# Add middle gene indices to tgfm results
		tgfm_results['middle_gene_indices'] = np.copy(tgfm_data['middle_gene_indices'])
		tgfm_results['middle_variant_indices'] = np.copy(tgfm_data['middle_variant_indices'])

		# Middle gene pips
		middle_gene_pips = tgfm_results['expected_alpha_pips'][tgfm_results['middle_gene_indices']]
		middle_gene_names = tgfm_results['genes'][tgfm_results['middle_gene_indices']]

		alphas = tgfm_results['alpha_phis']
		middle_genes_indices = tgfm_results['middle_gene_indices']
		middle_gene_tissues = tgfm_results['genes'][middle_genes_indices]

		tissue_group_indices = []
		for ii, gene_name in enumerate(middle_gene_names):
			if '_'.join(gene_name.split('_')[1:]) in tissue_group_dicti:
				tissue_group_indices.append(ii)
		tissue_group_indices = np.asarray(tissue_group_indices)
		if len(tissue_group_indices) == 0:
			continue

		middle_gene_pips = middle_gene_pips[tissue_group_indices]
		middle_gene_names = middle_gene_names[tissue_group_indices]
		middle_gene_tissues = middle_gene_tissues[tissue_group_indices]


		gene_to_indices = {}
		for ii, middle_gene_tissue in enumerate(middle_gene_tissues):
			gene_name = middle_gene_tissue.split('_')[0].split('.')[0]
			if gene_name not in gene_to_indices:
				gene_to_indices[gene_name] = []
			gene_to_indices[gene_name].append(ii)
		n_genes = len(gene_to_indices)
		ordered_genes = [*gene_to_indices]
		LL = len(alphas)
		gene_alphas = []
		for ll in range(LL):
			new_alpha = np.zeros((alphas[ll].shape[0], n_genes))
			for ii,gene_name in enumerate(ordered_genes):
				tmper = alphas[ll][:, middle_genes_indices]
				new_alpha[:, ii] = np.sum(tmper[:, gene_to_indices[gene_name]],axis=1)
			gene_alphas.append(new_alpha)
		expected_alpha_pips = compute_expected_pips_from_sampler_pis(gene_alphas)
		for ii, gene_name in enumerate(ordered_genes):
			agg_pip = expected_alpha_pips[ii]
			if agg_pip >= .5:
				print(gene_name + '\t' + str(agg_pip))
	return


def generate_per_gene_pip_summary_file(concatenated_pip_summary_file, per_gene_pip_summary_file, data_file_stem, file_stem, model_version, processed_tgfm_input_stem):
	f = open(concatenated_pip_summary_file)
	t = open(per_gene_pip_summary_file,'w')
	t.write('gene_name\twindow_name\tPIP\n')
	head_count = 0
	used_genes = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]
		if len(data) != 3:
			continue
		# Load in TGFM results pkl file for this window
		window_pkl_file = data_file_stem + '_' + window_name + '_results.pkl'
		# Load in tgfm results data
		g = open(window_pkl_file, "rb")
		tgfm_results = pickle.load(g)
		g.close()


		# Load in tgfm results data
		g = open(processed_tgfm_input_stem + '_' + window_name + '_tgfm_trait_agnostic_input_data_obj.pkl', "rb")
		tgfm_data = pickle.load(g)
		g.close()
		# Add middle gene indices to tgfm results
		tgfm_results['middle_gene_indices'] = np.copy(tgfm_data['middle_gene_indices'])
		tgfm_results['middle_variant_indices'] = np.copy(tgfm_data['middle_variant_indices'])

		# Middle gene pips
		middle_gene_pips = tgfm_results['expected_alpha_pips'][tgfm_results['middle_gene_indices']]
		middle_gene_names = tgfm_results['genes'][tgfm_results['middle_gene_indices']]

		alphas = tgfm_results['alpha_phis']
		middle_genes_indices = tgfm_results['middle_gene_indices']
		middle_gene_tissues = tgfm_results['genes'][middle_genes_indices]

		gene_to_indices = {}
		for ii, middle_gene_tissue in enumerate(middle_gene_tissues):
			gene_name = middle_gene_tissue.split('_')[0].split('.')[0]
			if gene_name not in gene_to_indices:
				gene_to_indices[gene_name] = []
			gene_to_indices[gene_name].append(ii)
		n_genes = len(gene_to_indices)
		ordered_genes = [*gene_to_indices]
		LL = len(alphas)
		gene_alphas = []
		for ll in range(LL):
			new_alpha = np.zeros((alphas[ll].shape[0], n_genes))
			for ii,gene_name in enumerate(ordered_genes):
				tmper = alphas[ll][:, middle_genes_indices]
				new_alpha[:, ii] = np.sum(tmper[:, gene_to_indices[gene_name]],axis=1)
			gene_alphas.append(new_alpha)
		expected_alpha_pips = compute_expected_pips_from_sampler_pis(gene_alphas)
		for ii, gene_name in enumerate(ordered_genes):
			agg_pip = expected_alpha_pips[ii]
			if agg_pip >= .01:
				t.write(gene_name + '\t' + window_name + '\t' + str(agg_pip) + '\n')
	f.close()
	t.close()
	return


def generate_per_gene_alttissue_pip_summary_file(concatenated_pip_summary_file, per_gene_alttissue_pip_summary_file, tissue_name_to_alt_tissue, data_file_stem, file_stem, model_version, processed_tgfm_input_stem):
	f = open(concatenated_pip_summary_file)
	t = open(per_gene_alttissue_pip_summary_file,'w')
	t.write('gene_tissue_name\tgene_name\ttissue_name\ttissue_visualization_category\twindow_name\tPIP\n')
	head_count = 0
	used_genes = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]
		if len(data) != 3:
			continue
		# Load in TGFM results pkl file for this window
		window_pkl_file = data_file_stem + '_' + window_name + '_results.pkl'
		# Load in tgfm results data
		g = open(window_pkl_file, "rb")
		tgfm_results = pickle.load(g)
		g.close()

		# Load in tgfm results data
		g = open(processed_tgfm_input_stem + '_' + window_name + '_tgfm_trait_agnostic_input_data_obj.pkl', "rb")
		tgfm_data = pickle.load(g)
		g.close()
		# Add middle gene indices to tgfm results
		tgfm_results['middle_gene_indices'] = np.copy(tgfm_data['middle_gene_indices'])
		tgfm_results['middle_variant_indices'] = np.copy(tgfm_data['middle_variant_indices'])

		# Middle gene pips
		middle_gene_pips = tgfm_results['expected_alpha_pips'][tgfm_results['middle_gene_indices']]
		middle_gene_names = tgfm_results['genes'][tgfm_results['middle_gene_indices']]


		mapping = {}
		gene_alt_tissues = []
		genes = []
		for ii, gene_tissue_name in enumerate(middle_gene_names):
			gene_name = gene_tissue_name.split('_')[0]
			tissue_name = '_'.join(gene_tissue_name.split('_')[1:])
			alt_tissue_name = tissue_name_to_alt_tissue[tissue_name]
			gene_tissue_pip = middle_gene_pips[ii]
			genes.append(gene_name)

			gene_alt_tissue_name = gene_name + '_' + alt_tissue_name

			if gene_alt_tissue_name not in mapping:
				mapping[gene_alt_tissue_name] = gene_tissue_pip
				gene_alt_tissues.append(gene_alt_tissue_name)
			else:
				mapping[gene_alt_tissue_name] = mapping[gene_alt_tissue_name] + gene_tissue_pip

		for gene_alt_tissue in gene_alt_tissues:
			gene_alt_tissue_pip = mapping[gene_alt_tissue]
			if gene_alt_tissue_pip < .01:
				continue
			gene_name = gene_alt_tissue.split('_')[0]
			alt_tissue = '_'.join(gene_alt_tissue.split('_')[1:])

			t.write(gene_alt_tissue + '\t' + gene_name + '\t' + alt_tissue + '\t' + alt_tissue + '\t' + window_name + '\t' + str(gene_alt_tissue_pip) + '\n')
		unique_window_genes = np.unique(genes)
		for gene_name in unique_window_genes:
			if gene_name in used_genes:
				print('assumption eororor')
				pdb.set_trace()
			used_genes[gene_name] = 1

	t.close()
	f.close()
	return

def tally_number_of_causal_genes(per_gene_pip_summary_file, n_causal_genes_summary_file, pip_threshold):
	# Extract number of genes with pip above threshold
	f = open(per_gene_pip_summary_file)
	head_count = 0
	used = {}
	gene_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_name = data[0]
		pip = float(data[2])
		if gene_name in used:
			print('assumption eroroor')
			pdb.set_trace()
		used[gene_name] = 1
		if pip >= pip_threshold:
			gene_count = gene_count + 1
	f.close()
	
	# Print to output
	t = open(n_causal_genes_summary_file, 'w')
	t.write('element_name\tcount\n')
	t.write('gene\t' + str(gene_count) + '\n')
	t.close()
	return

def get_ensamble_id_and_gene_name(info_str):
	info = info_str.split('"')
	passed1 = False
	passed2 = False
	for ii, ele in enumerate(info):
		if ele == 'gene_id ':
			ensamble_id = info[(ii+1)]
			passed1 = True
		if ele == '; gene_name ':
			gene_name = info[(ii+1)]
			passed2 = True
	if passed1 == False or passed2 == False:
		print('assumption oerororr')
		pdb.set_trace()
	return ensamble_id, gene_name


def create_ensamble_id_to_gene_name_mapping(gene_annotation_file):
	f = open(gene_annotation_file)
	ensg_to_gene_name = {}
	gene_name_to_ensg = {}
	for line in f:
		if line.startswith('##'):
			continue
		line = line.rstrip()
		data = line.split('\t')
		if data[2] != 'gene':
			continue
		info_str = data[8]
		ensamble_id, gene_name = get_ensamble_id_and_gene_name(info_str)
		if ensamble_id in ensg_to_gene_name:
			print('repeat ensamble id')
			pdb.set_trace()
		ensg_to_gene_name[ensamble_id] = gene_name
		ensg_to_gene_name[ensamble_id.split('.')[0]] = gene_name
	f.close()
	return ensg_to_gene_name


def create_hit_summary_file(per_gene_tissue_pip_summary_file, per_gene_pip_summary_file, ensamble_id_to_gene_id, gene_hit_summary_file):
	# Create mapping from gene name to agg pip
	f = open(per_gene_pip_summary_file)
	gene_to_agg_pip = {}
	gene_to_window_name = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_name = data[0]
		pip = float(data[2])
		window_name = data[1]
		if gene_name in gene_to_agg_pip:
			print('assumption oeroro')
			pdb.set_trace()
		gene_to_agg_pip[gene_name] = pip
		gene_to_window_name[gene_name] = window_name
	f.close()


	# Open ouptut file handle
	t = open(gene_hit_summary_file,'w')
	# header
	t.write('gene_tissue\tensamble_id\tgene_id\ttissue_name\tgene_tissue_pip\tgene_pip\twindow_name\n')

	f = open(per_gene_tissue_pip_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_tissue_name = data[0]
		ensamble_id = data[1].split('.')[0]
		tissue_name = data[2]
		gene_tissue_pip = float(data[5])
		gene_pip = gene_to_agg_pip[ensamble_id]
		gene_id = ensamble_id_to_gene_id[ensamble_id]
		window_name = gene_to_window_name[ensamble_id]

		# Print to ouput
		if gene_tissue_pip > .5:
			t.write(gene_tissue_name + '\t' + ensamble_id + '\t' + gene_id + '\t' + tissue_name + '\t' + str(gene_tissue_pip) + '\t' + str(gene_pip) + '\t' + window_name + '\n')

	# Close file handles
	f.close()
	t.close()

	return


tgfm_results_dir = sys.argv[1]
gene_type = sys.argv[2]
num_jobs = int(sys.argv[3])
trait_names_file = sys.argv[4]
gtex_pseudotissue_file = sys.argv[5]
gtex_pseudotissue_category_file = sys.argv[6]
processed_tgfm_input_stem = sys.argv[7]
ukbb_preprocessed_for_genome_wide_susie_dir = sys.argv[8]
tgfm_organized_results_dir = sys.argv[9]
gene_annotation_file = sys.argv[10]



model_versions = ['susie_pmces_variant_gene', 'susie_sampler_variant_gene', 'susie_pmces_sparse_variant_gene_tissue', 'susie_sampler_sparse_variant_gene_tissue']
model_versions = ['susie_sampler_uniform', 'susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler', 'susie_sampler_tglr_bootstrapped_nonnegative_sampler']
model_versions = ['susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler']


#Extract tissue names
tissue_names = extract_pseudotissue_names(gtex_pseudotissue_file, ignore_testis=True)


# Create mapping from ensamble id to gene id
ensamble_id_to_gene_id = create_ensamble_id_to_gene_name_mapping(gene_annotation_file)


# Create mapping from category to tissue names
tissue_categories, tissue_category_to_tissue_names, tissue_name_to_broad_category, tissue_name_to_alt_tissue, alt_tissues = create_mapping_from_tissue_category_to_tissue_names(gtex_pseudotissue_category_file, ignore_testis=True)
'''
tissue_name_to_broad_category['B'] = 'sc_blood'
tissue_name_to_broad_category['NK'] = 'sc_blood'
tissue_name_to_broad_category['Prolif'] = 'sc_blood'
tissue_name_to_broad_category['T4'] = 'sc_blood'
tissue_name_to_broad_category['T8'] = 'sc_blood'
tissue_name_to_broad_category['cDC'] = 'sc_blood'
tissue_name_to_broad_category['cM'] = 'sc_blood'
tissue_name_to_broad_category['ncM'] = 'sc_blood'
tissue_name_to_broad_category['pDC'] = 'sc_blood'
'''

# Extract trait names
trait_names = extract_trait_names(trait_names_file)
#valid_trait_names = {'biochemistry_Cholesterol':1, 'blood_MEAN_PLATELET_VOL':1, 'blood_MONOCYTE_COUNT':1, 'body_BMIz':1, 'body_WHRadjBMIz':1, 'bp_DIASTOLICadjMEDz':1, 'biochemistry_VitaminD':1, 'blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT':1, 'lung_FEV1FVCzSMOKE':1}
print(trait_names)
print(len((trait_names)))
arr1 = []
arr2 = []











# Concatenate parrallelized results across runs for each trait
for trait_name in trait_names:
	print('###################################')
	print('###################################')
	print(trait_name)


	#if trait_name not in valid_trait_names:
		#continue
	for model_version in model_versions:
		###################################################
		# Concatenate PIP summary file across parallel runs (one line for each window)
		###################################################
		data_file_stem = tgfm_results_dir + 'tgfm_results_' + trait_name + '_' + gene_type + '_' + model_version
		file_stem = tgfm_organized_results_dir + 'tgfm_results_' + trait_name + '_' + gene_type + '_' + model_version
		suffix = 'tgfm_pip_summary.txt'
		concatenated_pip_summary_file = file_stem + '_' + suffix
		#concatenate_results_across_parallel_jobs(data_file_stem, suffix, num_jobs, concatenated_pip_summary_file)

		###################################################
		# Create full gene-Tissue pip summary file
		###################################################
		per_gene_tissue_full_pip_summary_file = file_stem + '_tgfm_per_gene_tissue_full_pip_summary.txt'
		#generate_per_gene_tissue_pip_full_summary_file(concatenated_pip_summary_file, per_gene_tissue_full_pip_summary_file, tissue_name_to_broad_category, data_file_stem, file_stem, model_version, processed_tgfm_input_stem, trait_name)

		###################################################
		# Update all PIP file
		###################################################
		# Open output file handle keeping track of all PIPs
		ttt = open(tgfm_organized_results_dir + 'GTEx_TGFM_PIPs_' + trait_name + '.txt','w')
		ttt.write('trait_name\tgenetic_element_class\tgenetic_element_name\tTGFM_PIP\n')
		ttt = print_trait_pips_to_all_pip_file(ttt, concatenated_pip_summary_file, data_file_stem, file_stem, model_version, processed_tgfm_input_stem, trait_name)
		ttt.close()
		# zip up file
		os.system('gzip ' + tgfm_organized_results_dir + 'GTEx_TGFM_PIPs_' + trait_name + '.txt')
		print(tgfm_organized_results_dir + 'GTEx_TGFM_PIPs_' + trait_name + '.txt')



		'''
		###################################################
		# Create gene-Tissue pip summary file
		###################################################
		per_gene_tissue_pip_summary_file = file_stem + '_tgfm_per_gene_tissue_pip_summary.txt'
		generate_per_gene_tissue_pip_summary_file(concatenated_pip_summary_file, per_gene_tissue_pip_summary_file, tissue_name_to_broad_category)

		###################################################
		# Create gene pip summary file
		###################################################
		per_gene_pip_summary_file = file_stem + '_tgfm_per_gene_pip_summary.txt'
		generate_per_gene_pip_summary_file(concatenated_pip_summary_file, per_gene_pip_summary_file, data_file_stem, file_stem, model_version, processed_tgfm_input_stem)


		###################################################
		# Create hit summary file
		###################################################
		gene_hit_summary_file = file_stem + '_tgfm_gene_tissue_hit_summary.txt'
		create_hit_summary_file(per_gene_tissue_pip_summary_file, per_gene_pip_summary_file, ensamble_id_to_gene_id, gene_hit_summary_file)



		###################################################
		# Tally up number of expected causal genetic elements
		###################################################
		n_causal_genetic_elements_summary_file = file_stem + '_tgfm_expected_n_causal_genetic_elements.txt'
		n_causal_genetic_elements_by_tissue_summary_file = file_stem + '_tgfm_expected_n_causal_genetic_elements_tissue_stratefied.txt'
		tally_expected_number_of_causal_genetic_elements(concatenated_pip_summary_file, n_causal_genetic_elements_summary_file, n_causal_genetic_elements_by_tissue_summary_file, data_file_stem,file_stem, model_version, processed_tgfm_input_stem, tissue_names)

		###################################################
		# Tally up number of causal genetic elements
		###################################################
		for pip_threshold in [.1, .2, .25, .3 , .5, .7, .75, .9, .95, .99]:
			n_causal_genetic_elements_summary_file = file_stem + '_tgfm_n_causal_genetic_elements_pip_' + str(pip_threshold) + '.txt'
			n_causal_genetic_elements_by_tissue_summary_file = file_stem + '_tgfm_n_causal_genetic_elements_tissue_stratefied_pip_' + str(pip_threshold) + '.txt'
			n_genes, n_tot = tally_number_of_causal_genetic_elements(concatenated_pip_summary_file, n_causal_genetic_elements_summary_file, n_causal_genetic_elements_by_tissue_summary_file, pip_threshold, tissue_names)
			#if pip_threshold == .5:
				#arr1.append(n_genes/n_tot)
				#arr2.append(n_tot)

		###################################################
		# Tally up number of causal genes (not gene-tissue pairs)
		###################################################
		for pip_threshold in [.1, .2, .25, .3 , .5, .7, .75, .9, .95, .99]:
			n_causal_genes_summary_file = file_stem + '_tgfm_n_causal_genes_pip_' + str(pip_threshold) + '.txt'
			tally_number_of_causal_genes(per_gene_pip_summary_file, n_causal_genes_summary_file, pip_threshold)

		###################################################
		# Tally up number of causal gene-tissue pairs across thresholds
		###################################################
		n_causal_gene_tissue_pairs_summary_cross_threshold_file = file_stem + '_tgfm_n_causal_gene_tissue_pairs_cross_pip_threshold_sqrt_plot_input.txt'
		tally_number_of_causal_gene_tissue_pairs_cross_pip_thresholds(concatenated_pip_summary_file, n_causal_gene_tissue_pairs_summary_cross_threshold_file)

		###################################################
		# Tally up number of variants across thresholds
		###################################################
		n_causal_variants_summary_cross_threshold_file = file_stem + '_tgfm_n_causal_variants_cross_pip_threshold_sqrt_plot_input.txt'
		tally_number_of_causal_variants_cross_pip_thresholds(concatenated_pip_summary_file, n_causal_variants_summary_cross_threshold_file)

		###################################################
		# Tally up number of genes across thresholds
		###################################################
		n_causal_genes_summary_cross_threshold_file = file_stem + '_tgfm_n_causal_genes_cross_pip_threshold_sqrt_plot_input.txt'
		tally_number_of_causal_genes_cross_pip_thresholds(per_gene_pip_summary_file, n_causal_genes_summary_cross_threshold_file)


		###################################################
		# Tally up number of causal sc gene-tissue pairs across thresholds
		###################################################
		n_causal_sc_gene_tissue_pairs_summary_cross_threshold_file = file_stem + '_tgfm_n_causal_sc_gene_tissue_pairs_cross_pip_threshold_sqrt_plot_input.txt'
		tally_number_of_causal_sc_gene_tissue_pairs_cross_pip_thresholds(concatenated_pip_summary_file, n_causal_sc_gene_tissue_pairs_summary_cross_threshold_file, tissue_name_to_broad_category)

		###################################################
		# Tally up number of causal sc gene-tissue pairs across thresholds in each cell type
		###################################################	
		if trait_name == 'blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT':
			continue
		blood_cell_types = np.asarray(['B', 'NK', 'Prolif', 'T4', 'T8', 'cDC', 'cM', 'ncM', 'pDC'])
		blood_cell_types = np.copy(tissue_names)
		for blood_cell_type in blood_cell_types:
			n_causal_sc_single_cell_type_gene_tissue_pairs_summary_cross_threshold_file = file_stem + '_tgfm_n_causal_tissue_gene_tissue_pairs_' + str(blood_cell_type) + '_cross_pip_threshold_sqrt_plot_input.txt'
			tally_number_of_causal_sc_gene_tissue_pairs_cross_pip_thresholds_in_single_cell_type(concatenated_pip_summary_file, n_causal_sc_single_cell_type_gene_tissue_pairs_summary_cross_threshold_file, blood_cell_type)

		'''



'''
# Tally up number of gene-{tissue-group} pairs
trait_name='body_WHRadjBMIz'
tissue_group = ['Adipose_Subcutaneous', 'Adipose_Visceral_Omentum']
data_file_stem = tgfm_results_dir + 'tgfm_results_' + trait_name + '_' + gene_type + '_' + model_versions[0]
file_stem = tgfm_organized_results_dir + 'tgfm_results_' + trait_name + '_' + gene_type + '_' + model_versions[0]
suffix = 'tgfm_pip_summary.txt'
concatenated_pip_summary_file = file_stem + '_' + suffix

generate_per_gene_tissue_group_pip_summary_file(concatenated_pip_summary_file, data_file_stem, file_stem, model_version, processed_tgfm_input_stem, tissue_group)
'''






