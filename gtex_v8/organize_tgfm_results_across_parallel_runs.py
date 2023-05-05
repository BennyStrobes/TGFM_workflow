import numpy as np 
import os
import sys
import pdb
import pickle


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

def generate_per_tissue_category_pip_summary_file(concatenated_pip_summary_file, per_tissue_category_pip_summary_file, tissue_names,tissue_categories, tissue_category_to_tissue_names, file_stem, model_version, processed_tgfm_input_stem):
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

def generate_component_level_summary_data(concatenated_pip_summary_file, component_level_summary_file, tissue_names, file_stem, model_version, processed_tgfm_input_stem):
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
				sample_valid_components = tgfm_results['valid_components'][sample_iter]
				for component_iter in sample_valid_components:
					if component_in_middle_of_window(tgfm_results['alpha_phis'][component_iter][sample_iter,:], tgfm_results['beta_phis'][component_iter][sample_iter,:], middle_gene_indices_dicti, middle_variant_indices_dicti):
						non_mediated_prob = np.sum(tgfm_results['beta_phis'][component_iter][sample_iter,:])
						mediated_probs = extract_mediated_probs(tgfm_results['alpha_phis'][component_iter][sample_iter,:], len(tissue_names), tiss_to_position_mapping, tgfm_results['genes'])
						t.write(window_name + '\t' + str(sample_iter) + '\t' + str(component_iter) + '\t' + str(non_mediated_prob) + '\t' + '\t'.join(mediated_probs.astype(str)) + '\n')
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

	'''
	print('N_components: ' + str(n_components))
	print('Avg non-mediated prob: ' + str(avg_non_mediated_prob))

	for tiss_iter in np.argsort(-tissue_probs):
		print(tissue_names[tiss_iter] + ': ' + str(tissue_probs[tiss_iter]))
	'''
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
	f = open(gtex_pseudotissue_category_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissue_name = data[0]
		category = data[-2]
		broad_category = data[-1]
		if tissue_name == 'Testis' and ignore_testis:
			continue
		mapping2[tissue_name] = broad_category
		category_names.append(category)
		if category not in mapping:
			mapping[category] = []
		mapping[category].append(tissue_name)
	f.close()

	category_names = np.sort(np.unique(np.asarray(category_names)))

	for category_name in category_names:
		mapping[category_name] = np.asarray(mapping[category_name])

	return category_names, mapping, mapping2


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


tgfm_results_dir = sys.argv[1]
gene_type = sys.argv[2]
num_jobs = int(sys.argv[3])
trait_names_file = sys.argv[4]
gtex_pseudotissue_file = sys.argv[5]
gtex_pseudotissue_category_file = sys.argv[6]
processed_tgfm_input_stem = sys.argv[7]
ukbb_preprocessed_for_genome_wide_susie_dir = sys.argv[8]
tgfm_sldsc_results_dir = sys.argv[9]


model_versions = ['susie_pmces_variant_gene', 'susie_sampler_variant_gene', 'susie_pmces_sparse_variant_gene_tissue', 'susie_sampler_sparse_variant_gene_tissue']


#Extract tissue names
tissue_names = extract_pseudotissue_names(gtex_pseudotissue_file, ignore_testis=True)

# Create mapping from category to tissue names
tissue_categories, tissue_category_to_tissue_names, tissue_name_to_broad_category = create_mapping_from_tissue_category_to_tissue_names(gtex_pseudotissue_category_file, ignore_testis=True)

# Extract trait names
trait_names = extract_trait_names(trait_names_file)
print(trait_names)

valid_trait_names = {'biochemistry_Cholesterol':1, 'blood_MEAN_PLATELET_VOL':1, 'blood_MONOCYTE_COUNT':1, 'body_BMIz':1, 'body_WHRadjBMIz':1, 'bp_DIASTOLICadjMEDz':1, 'biochemistry_VitaminD':1, 'blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT':1, 'lung_FEV1FVCzSMOKE':1}

valid_trait_names = {'blood_MEAN_PLATELET_VOL':1, 'blood_MONOCYTE_COUNT':1, 'body_BMIz':1, 'body_WHRadjBMIz':1, 'bp_DIASTOLICadjMEDz':1, 'lung_FEV1FVCzSMOKE':1}



# Concatenate parrallelized results across runs for each trait
for trait_name in trait_names:
	print('###################################')
	print('###################################')
	print(trait_name)
	if trait_name not in valid_trait_names:
		continue
	for model_version in model_versions:
		print(model_version)
		###################################################
		# Concatenate PIP summary file across parallel runs (one line for each window)
		###################################################
		file_stem = tgfm_results_dir + 'tgfm_results_' + trait_name + '_' + gene_type + '_' + model_version
		suffix = 'tgfm_pip_summary.txt'
		concatenated_pip_summary_file = file_stem + '_' + suffix
		concatenate_results_across_parallel_jobs(file_stem, suffix, num_jobs, concatenated_pip_summary_file)

		###################################################
		# Create component level summary data
		###################################################
		component_level_summary_file = file_stem + '_tgfm_component_level_summary.txt'
		generate_component_level_summary_data(concatenated_pip_summary_file, component_level_summary_file, tissue_names, file_stem, model_version, processed_tgfm_input_stem)

		##################################################
		# Correlate fraction mediated h2 with fraction of components mediated
		##################################################
		sldsc_fraction_mediated_h2_file = tgfm_sldsc_results_dir + trait_name + '_baseline_no_qtl_component_gene_no_testis_pmces_gene_adj_ld_scores_h2_5_50_med.txt'
		fraction_mediated_sldsc_comparison_summary_file = file_stem + '_tgfm_fraction_mediated_sldsc_comparison.txt'
		correlate_fraction_mediated_h2_with_fraction_of_components_mediated(component_level_summary_file, sldsc_fraction_mediated_h2_file, fraction_mediated_sldsc_comparison_summary_file, model_version, trait_name)

		##################################################
		# Correlate per tissue h2 with number of components identifiedin each tissue
		##################################################
		sldsc_mediated_h2_file = tgfm_sldsc_results_dir + trait_name + '_baseline_no_qtl_component_gene_no_testis_pmces_gene_adj_ld_scores_organized_mediated_h2_5_50.txt'
		per_tissue_sldsc_comparison_summary_file = file_stem + '_tgfm_per_tissue_sldsc_comparison_component_level.txt'
		correlate_per_tissue_h2_with_number_of_components_identified_in_each_tissue(component_level_summary_file, sldsc_mediated_h2_file, per_tissue_sldsc_comparison_summary_file, model_version, trait_name, tissue_names, tissue_name_to_broad_category)

		###################################################
		# Create gene-Tissue pip summary file
		###################################################
		per_gene_tissue_pip_summary_file = file_stem + '_tgfm_per_gene_tissue_pip_summary.txt'
		generate_per_gene_tissue_pip_summary_file(concatenated_pip_summary_file, per_gene_tissue_pip_summary_file, tissue_name_to_broad_category)

		###################################################
		# Correlate per tissue h2 with number of components identified in each tissue
		###################################################
		sldsc_mediated_h2_file = tgfm_sldsc_results_dir + trait_name + '_baseline_no_qtl_component_gene_no_testis_pmces_gene_adj_ld_scores_organized_mediated_h2_5_50.txt'
		per_tissue_sldsc_comparison_summary_file = file_stem + '_tgfm_per_tissue_sldsc_comparison_old.txt'
		pip_threshold = 0.2
		correlate_number_of_high_pip_genes_with_per_tissue_sldsc_h2(tissue_names, per_gene_tissue_pip_summary_file, pip_threshold, sldsc_mediated_h2_file, per_tissue_sldsc_comparison_summary_file, tissue_name_to_broad_category)
		
		###################################################
		# Extract nominal twas z and component-iter and PIP confidence intervals for high pip gene-tissue pair
		###################################################
		pip_threshold = .4
		per_gene_tissue_high_pip_extensive_summary_file = file_stem + '_tgfm_per_high_pip_gene_tissue_extensive_summary.txt'
		per_gene_tissue_high_pip_high_confidence_extensive_summary_file = file_stem + '_tgfm_per_high_pip_high_confidence_gene_tissue_extensive_summary.txt'
		generate_extensive_per_gene_tissue_pip_summary_file(per_gene_tissue_pip_summary_file, per_gene_tissue_high_pip_extensive_summary_file, tissue_names, file_stem, model_version, processed_tgfm_input_stem,pip_threshold, ukbb_preprocessed_for_genome_wide_susie_dir, trait_name, per_gene_tissue_high_pip_high_confidence_extensive_summary_file)

		###################################################
		# Create tissue pip summary file
		###################################################
		per_tissue_pip_summary_file = file_stem + '_tgfm_per_tissue_pip_summary.txt'
		generate_per_tissue_pip_summary_file(concatenated_pip_summary_file, per_tissue_pip_summary_file, tissue_names, file_stem, model_version, processed_tgfm_input_stem)

		###################################################
		# Create tissue-category pip summary file
		###################################################
		per_tissue_category_pip_summary_file = file_stem + '_tgfm_per_tissue_category_pip_summary.txt'
		generate_per_tissue_category_pip_summary_file(concatenated_pip_summary_file, per_tissue_category_pip_summary_file, tissue_names,tissue_categories, tissue_category_to_tissue_names, file_stem, model_version, processed_tgfm_input_stem)
