import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import pickle
import tgfm
import scipy.optimize
import rpy2
import rpy2.robjects.numpy2ri as numpy2ri
import rpy2.robjects as ro
ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
rcpp_num_R_pkg = importr('RcppNumerical')
base_R_pkg = importr('base')




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

def generate_component_level_abf_summary_data(concatenated_pip_summary_file, component_level_abf_summary_file, tissue_names, file_stem, model_version, processed_tgfm_input_stem, per_window_abf_output_stem, version='v1'):
	f = open(concatenated_pip_summary_file)
	t = open(component_level_abf_summary_file,'w')
	t.write('window_name\tbootstrap_number\tcomponent_number\telement_names\tmiddle_element_names\twindow_abf_file\twindow_component_samples_file\tannotation_mat_file\tcomponent_var_init\tmu_init\tmu_var_init\n')
	tiss_to_position_mapping = {}
	for ii,tissue_name in enumerate(tissue_names):
		tiss_to_position_mapping[tissue_name] = ii

	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]
		#print(window_name)


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

		per_window_abf_output_file = per_window_abf_output_stem + window_name + '.npy'
		per_window_sample_names_output_file = per_window_abf_output_stem + window_name + '_sample_names.npy'
		per_window_anno_output_file = per_window_abf_output_stem + window_name + '_anno.npy'
		per_window_mu_output_file = per_window_abf_output_stem + window_name + 'mu_init.npy'
		per_window_mu_var_output_file = per_window_abf_output_stem + window_name + 'mu_var_init.npy'
		per_window_abf = []
		per_window_mu = []
		per_window_mu_var = []
		window_has_component = False
		per_window_samples = []
		if model_version.startswith('pmces'):
			if len(tgfm_results['valid_components']) > 6:
				continue
			lines = []
			if version == 'v2' or version == 'v4':
				for component_iter in range(tgfm_results['alpha_lbf'].shape[0]):
					# V2
					log_abf = np.hstack((tgfm_results['alpha_lbf'][component_iter, :], tgfm_results['beta_lbf'][component_iter, :]))
					mu_init = np.hstack((tgfm_results['alpha_mu'][component_iter, :], tgfm_results['beta_mu'][component_iter, :]))
					mu_var_init = np.hstack((tgfm_results['alpha_var'][component_iter, :], tgfm_results['beta_var'][component_iter, :]))
					component_var_init = tgfm_results['component_variances'][component_iter]
					element_names = np.hstack((tgfm_results['genes'], tgfm_results['variants']))
					middle_element_names = np.hstack((tgfm_results['genes'][tgfm_results['middle_gene_indices']], tgfm_results['variants'][tgfm_results['middle_variant_indices']]))
					liner = window_name + '\t' + 'NA' + '\t' + str(component_iter) + '\t' + ';'.join(element_names) + '\t' + ';'.join(middle_element_names) + '\t' + per_window_abf_output_file + '\t' + per_window_sample_names_output_file + '\t' + per_window_anno_output_file + '\t' + str(component_var_init) + '\t' + per_window_mu_output_file + '\t' + per_window_mu_var_output_file + '\n'
					lines.append(liner)
					per_window_abf.append(log_abf)
					per_window_mu.append(mu_init)
					per_window_mu_var.append(mu_var_init)
					per_window_samples.append(0)
					if component_iter in tgfm_results['valid_components']:
						if component_in_middle_of_window(tgfm_results['alpha_phi'][component_iter, :], tgfm_results['beta_phi'][component_iter, :], middle_gene_indices_dicti, middle_variant_indices_dicti):
							window_has_component = True	
				if window_has_component:
					for liner in lines:
						t.write(liner)
		else:
			print('assumption eroror: unknown model verison')
			pdb.set_trace()
		# Convert to np array and save
		if window_has_component:
			per_window_abf = np.asarray(per_window_abf)
			np.save(per_window_abf_output_file, per_window_abf)
			per_window_samples = np.asarray(per_window_samples)
			np.save(per_window_sample_names_output_file, per_window_samples)
			np.save(per_window_anno_output_file, tgfm_data['annotation'])
			np.save(per_window_mu_output_file, np.asarray(per_window_mu))
			np.save(per_window_mu_var_output_file, np.asarray(per_window_mu_var))
	t.close()
	f.close()
	return




def update_prior_prob_emperical_distributions_for_variant_gene_tissue(component_level_abf_summary_file, tgfm_version, tissue_name_to_position, variant_prob_distr, tissue_probs_distr, tissue_names, window_to_class_to_indices, n_bootstraps):
	# General info
	n_tissues = len(tissue_names)
	window_names = np.asarray([*window_to_class_to_indices])
	n_windows = len(window_names)

	# Keep track of counts in each window
	variant_posterior_sum = np.zeros(n_windows)
	variant_counts = np.zeros(n_windows)
	tissue_posterior_sum = np.zeros((n_tissues, n_windows))
	tissue_counts = np.zeros((n_tissues, n_windows))

	# Create prior hash
	prior_hash = {}
	prior_hash['variant'] = variant_prob_distr
	for ii in range(n_tissues):
		prior_hash[tissue_names[ii]] = tissue_probs_distr[ii,:]

	# Create mapping from window name to expected prior probs
	window_to_prior_probs = {}
	for window_name in window_names:
		prior_prob_raw = np.zeros((window_to_class_to_indices[window_name]['n_elements'], n_bootstraps))
		class_name = 'variant'
		class_indices = window_to_class_to_indices[window_name][class_name]

		prior_prob_raw[class_indices, :] = prior_hash[class_name]
		for class_name in tissue_names:
			class_indices = window_to_class_to_indices[window_name][class_name]
			if len(class_indices) == 0:
				continue
			prior_prob_raw[class_indices, :] = prior_hash[class_name]
		# Normalize each bootstrapped sample
		prior_probs = prior_prob_raw/np.sum(prior_prob_raw,axis=0)

		E_ln_pi = np.mean(np.log(prior_probs),axis=1)
		E_pi = np.exp(E_ln_pi)
		window_to_prior_probs[window_name] = E_pi
	
	# Loop through components
	head_count = 0
	prev_window_name = 'NULL'
	f = open(component_level_abf_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		#print(data[0])
		window_name = data[0]

		if window_name != prev_window_name:
			if prev_window_name == 'NULL':
				window_counter = -1
			window_abfs = np.load(data[4])
			component_counter = 0
			prev_window_name = window_name
			window_counter = window_counter + 1

		ele_labfs = window_abfs[component_counter, :]
		component_counter = component_counter + 1

		# Normalize prior probabilities
		prior_probs = window_to_prior_probs[window_name]
		# Get alpha vec given lbf and prior probs
		alpha_vec = extract_alpha_vec_given_lbf_and_prior(ele_labfs, prior_probs)
		# update global counts
		indices = np.asarray(window_to_class_to_indices[window_name]['variant'])

		variant_posterior_sum[window_counter] = variant_posterior_sum[window_counter] + np.sum(alpha_vec[indices])
		variant_counts[window_counter] = variant_counts[window_counter] + len(indices)
		# tissues
		for tiss_iter, tissue_name in enumerate(tissue_names):
			indices = np.asarray(window_to_class_to_indices[window_name][tissue_name])
			if len(indices) == 0:
				continue
			tissue_posterior_sum[tiss_iter, window_counter] = tissue_posterior_sum[tiss_iter, window_counter] + np.sum(alpha_vec[indices])
			tissue_counts[tiss_iter, window_counter] = tissue_counts[tiss_iter, window_counter] + len(indices)
	f.close()

	# all window indices
	all_window_indices = np.arange(n_windows)
	# Loop through bootstraps
	for bs_iter in range(n_bootstraps):
		# For each bootstrap, calculate a probability
		bs_window_indices = np.random.choice(all_window_indices, size=n_windows, replace=True)
		variant_prob_distr[bs_iter] = np.sum(variant_posterior_sum[bs_window_indices])/np.sum(variant_counts[bs_window_indices])
		tissue_probs_distr[:, bs_iter]  = np.sum(tissue_posterior_sum[:, bs_window_indices],axis=1)/np.sum(tissue_counts[:, bs_window_indices],axis=1)

	return variant_prob_distr, tissue_probs_distr

def update_prior_probs_for_variant_gene_tissue(component_level_abf_summary_file, tgfm_version, tissue_name_to_position, variant_prob, tissue_probs, tissue_names, window_to_class_to_indices, window_to_class_to_middle_indices):
	variant_posterior_sum = 0.0
	variant_counts = 0.0
	tissue_posterior_sum = np.zeros(len(tissue_probs))
	tissue_counts = np.zeros(len(tissue_probs))

	# Create prior hash
	prior_hash = {}
	prior_hash['variant'] = variant_prob
	for ii, tissue_prob in enumerate(tissue_probs):
		prior_hash[tissue_names[ii]] = tissue_prob

	# Create mapping from window name to prior probs
	window_to_prior_probs = {}
	window_names = np.asarray([*window_to_class_to_indices])
	for window_name in window_names:
		prior_prob_raw = np.zeros(window_to_class_to_indices[window_name]['n_elements'])
		class_name = 'variant'
		class_indices = window_to_class_to_indices[window_name][class_name]
		prior_prob_raw[class_indices] = prior_hash[class_name]
		for class_name in tissue_names:
			class_indices = window_to_class_to_indices[window_name][class_name]
			if len(class_indices) == 0:
				continue
			prior_prob_raw[class_indices] = prior_hash[class_name]
		prior_probs = prior_prob_raw/np.sum(prior_prob_raw)
		window_to_prior_probs[window_name] = prior_probs

	# Loop through components
	head_count = 0
	prev_window_name = 'NULL'
	f = open(component_level_abf_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		#print(data[0])
		window_name = data[0]

		if window_name != prev_window_name:
			pdb.set_trace()
			window_abfs = np.load(data[4])
			component_counter = 0
			prev_window_name = window_name

		ele_labfs = window_abfs[component_counter, :]
		component_counter = component_counter + 1

		# Normalize prior probabilities
		prior_probs = window_to_prior_probs[window_name]
		# Get alpha vec given lbf and prior probs
		alpha_vec = extract_alpha_vec_given_lbf_and_prior(ele_labfs, prior_probs)
		# update global counts
		pdb.set_trace()
		indices = np.asarray(window_to_class_to_indices[window_name]['variant'])
		variant_posterior_sum = variant_posterior_sum + np.sum(alpha_vec[indices])
		variant_counts = variant_counts + len(indices)
		# tissues
		for tiss_iter, tissue_name in enumerate(tissue_names):
			indices = np.asarray(window_to_class_to_indices[window_name][tissue_name])
			if len(indices) == 0:
				continue
			tissue_posterior_sum[tiss_iter] = tissue_posterior_sum[tiss_iter] + np.sum(alpha_vec[indices])
			tissue_counts[tiss_iter] = tissue_counts[tiss_iter] + len(indices)

	f.close()

	variant_prob = variant_posterior_sum/variant_counts
	tissue_probs = tissue_posterior_sum/tissue_counts

	return variant_prob, tissue_probs

def extract_alpha_vec_given_lbf_and_prior(lbf, prior_probs):
	maxlbf = np.max(lbf)
	ww = np.exp(lbf - maxlbf)
	w_weighted = ww*prior_probs
	weighted_sum_w = sum(w_weighted)
	alpha = w_weighted / weighted_sum_w
	return alpha

def extract_alpha_vec_given_lbf_and_prior_matrix(lbf, prior_probs):
	maxlbf = np.max(lbf)
	ww = np.exp(lbf - maxlbf)
	w_weighted = (np.transpose(prior_probs)*ww)
	alpha = np.transpose(w_weighted)/np.sum(w_weighted,axis=1)
	return alpha

def extract_alpha_mat_given_lbf_and_prior_matrix(lbfs, prior_probs):
	max_lbfs = np.max(lbfs,axis=0)
	ww = np.exp(lbfs - max_lbfs)
	w_weighted = ww*prior_probs
	alpha = w_weighted/np.sum(w_weighted, axis=0)
	return alpha


def create_mapping_from_window_to_class_to_indices(component_level_abf_summary_file, tissue_name_to_position, tissue_names):
	f = open(component_level_abf_summary_file)
	dicti = {}
	dicti_middle = {}
	dicti_variant_middle = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]
		if window_name in dicti:
			continue
		# initialize sub-dictionary
		sub_dicti = {}
		sub_dicti_middle = {}
		arr_middle_variant = []
		sub_dicti['variant'] = []
		sub_dicti_middle['variant'] = []
		for tissue_name in tissue_names:
			sub_dicti[tissue_name] = []
			sub_dicti_middle[tissue_name] = []

		# Extract element names
		ele_names = np.asarray(data[3].split(';'))
		middle_ele_names = np.asarray(data[4].split(';'))
		middle_ele_dicti = {}
		for middle_ele_name in middle_ele_names:
			middle_ele_dicti[middle_ele_name] = 1


		# loop through elements and add to sub-dictionary
		variant_counter = -1
		for ele_iter, ele_name in enumerate(ele_names):
			if ele_name.startswith('ENSG'):
				ele_category = '_'.join(ele_name.split('_')[1:])
			else:
				ele_category = 'variant'
				variant_counter = variant_counter + 1
				if ele_name in middle_ele_dicti:
					arr_middle_variant.append(variant_counter)

			sub_dicti[ele_category].append(ele_iter)
			if ele_name in middle_ele_dicti:
				sub_dicti_middle[ele_category].append(ele_iter)
		
		# Change to np arrays
		sub_dicti['variant'] = np.asarray(sub_dicti['variant'])
		sub_dicti_middle['variant'] = np.asarray(sub_dicti_middle['variant'])
		for tissue_name in tissue_names:
			sub_dicti[tissue_name] = np.asarray(sub_dicti[tissue_name])
			sub_dicti_middle[tissue_name] = np.asarray(sub_dicti_middle[tissue_name])
		# Add n_elements to sub_dicti
		sub_dicti['n_elements'] = len(ele_names)

		# add to global dictionary
		dicti[window_name] = sub_dicti
		dicti_middle[window_name] = sub_dicti_middle
		dicti_variant_middle[window_name] = np.asarray(arr_middle_variant)

	f.close()
	return dicti, dicti_middle, dicti_variant_middle

def learn_iterative_variant_gene_tissue_prior(component_level_abf_summary_file, tgfm_version, tissue_names,per_window_abf_output_stem, max_iter=100):
	# Create mapping from tissue name to tissue position
	tissue_name_to_position = {}
	for ii, tissue_name in enumerate(tissue_names):
		tissue_name_to_position[tissue_name] = ii

	# Initialize prior probabilitie
	variant_prob = .1
	tissue_probs = np.ones(len(tissue_names))*.1

	# can probably precompute window to indices
	# Can probably precompute window to prior prob (in each iteration)
	# To do this
	## Get ordered list of window names
	## For each window name, create mapping from element class to indices

	# Create mapping from window to class to indices
	window_to_class_to_indices, window_to_class_to_middle_indices, window_to_variant_middle = create_mapping_from_window_to_class_to_indices(component_level_abf_summary_file, tissue_name_to_position, tissue_names)


	for itera in range(max_iter):
		variant_prob, tissue_probs = update_prior_probs_for_variant_gene_tissue(component_level_abf_summary_file, tgfm_version, tissue_name_to_position, variant_prob, tissue_probs, tissue_names, window_to_class_to_indices, window_to_class_to_middle_indices)
		print('################################')
		print('ITERATION: ' + str(itera))
		print('variant: ' + str(variant_prob))
		for ii, tissue_name in enumerate(tissue_names):
			print(tissue_name + ': ' + str(tissue_probs[ii]))

	return variant_prob, tissue_probs

def learn_iterative_variant_gene_tissue_emperical_distribution_prior(component_level_abf_summary_file, tgfm_version, tissue_names, per_window_abf_output_stem, max_iter=30, n_bootstraps=200):
	# Create mapping from tissue name to tissue position
	tissue_name_to_position = {}
	for ii, tissue_name in enumerate(tissue_names):
		tissue_name_to_position[tissue_name] = ii


	# Initialize prior probability emperical distribution
	variant_prob_distr = np.ones(n_bootstraps)*.1
	tissue_probs_distr = np.ones((len(tissue_names), n_bootstraps))*.1


	# Create mapping from window to class to indices
	window_to_class_to_indices = create_mapping_from_window_to_class_to_indices(component_level_abf_summary_file, tissue_name_to_position, tissue_names)

	# Iterative algorithm
	for itera in range(max_iter):
		variant_prob_distr, tissue_probs_distr = update_prior_prob_emperical_distributions_for_variant_gene_tissue(component_level_abf_summary_file, tgfm_version, tissue_name_to_position, variant_prob_distr, tissue_probs_distr, tissue_names, window_to_class_to_indices, n_bootstraps)
	return variant_prob_distr, tissue_probs_distr



def extract_iterative_variant_gene_tissue_log_priors(mapping, variants, genes):
	n_var = len(variants)
	n_genes = len(genes)
	n_bs = len(mapping['variant'])
	var_probs = np.zeros((n_var, n_bs))
	gene_probs = np.zeros((n_genes, n_bs))

	for var_iter in range(n_var):
		var_probs[var_iter, :] = mapping['variant']
	for gene_iter, gene_name in enumerate(genes):

		tissue_name = '_'.join(gene_name.split('_')[1:])
		gene_probs[gene_iter, :] = mapping[tissue_name]
	
	# Normalize rows
	normalizers = np.sum(gene_probs,axis=0) + np.sum(var_probs,axis=0)
	norm_var_probs = var_probs/normalizers
	norm_gene_probs = gene_probs/normalizers

	e_ln_pi_var = np.mean(np.log(norm_var_probs),axis=1)
	e_ln_pi_gene = np.mean(np.log(norm_gene_probs),axis=1)
	
	return e_ln_pi_var, e_ln_pi_gene


def print_log_prior_across_tgfm_windows(variant_gene_distr_prior_output_file, tgfm_input_summary_file, tmp_tgfm_stem):

	# Create mapping from element name to bs-probs
	mapping = {}
	f = open(variant_gene_distr_prior_output_file)
	head_count = 0
	for line in f:	
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		ele_name = data[0]
		bs_probs = np.asarray(data[3].split(';')).astype(float)
		mapping[ele_name] = bs_probs
	f.close()

	# Loop through windows
	f = open(tgfm_input_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]
		tgfm_input_pkl = data[2]


		# Load in tgfm input data
		g = open(tgfm_input_pkl, "rb")
		tgfm_data = pickle.load(g)
		g.close()

		# extract prior for this window
		variant_log_prior, gene_log_prior = extract_iterative_variant_gene_tissue_log_priors(mapping, tgfm_data['variants'], tgfm_data['genes'])

		# print to output
		output_file = tmp_tgfm_stem + '_' + window_name + '.txt'
		t = open(output_file,'w')
		t.write('element_name\tln_pi\n')
		for variant_iter, var_name in enumerate(tgfm_data['variants']):
			t.write(var_name + '\t' + str(variant_log_prior[variant_iter]) + '\n')
		for gene_iter, gene_name in enumerate(tgfm_data['genes']):
			t.write(gene_name + '\t' + str(gene_log_prior[gene_iter]) + '\n')
		t.close()
	f.close()

	return

def get_bs_indices_based_on_adjacent_window_groups(window_names, n_bootstraps):
	window_groups = []
	prev_start = 'null'
	window_indices = []
	for window_iter, window_name in enumerate(window_names):
		window_start = window_name.split(':')[1]
		window_end = window_name.split(':')[2]
		if prev_start != 'null':
			if int(window_start) - 1000000 != int(prev_start):
				window_groups.append(np.asarray(window_indices))
				window_indices = []
		window_indices.append(window_iter)
		prev_start = window_start
	window_groups.append(np.asarray(window_indices))

	n_groups = len(window_groups)

	bs_indices = []
	for bs_iter in range(n_bootstraps):
		group_indices = np.random.choice(np.arange(n_groups), size=n_groups, replace=True)
		tmp = []
		for group_index in group_indices:
			tmp.append(window_groups[group_index])
		indices = np.hstack(tmp)
		bs_indices.append(indices)
	return bs_indices


def create_mapping_from_window_component_to_prior_variance(component_level_abf_summary_file, n_bootstraps):
	dicti = {}
	f = open(component_level_abf_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]
		component_number = data[2]
		window_component_name = window_name + '_' + component_number
		prior_variance_init = float(data[-3])
		if window_component_name in dicti:
			print('assumption errorr')
			pdb.set_trace()
		dicti[window_component_name] = np.zeros(n_bootstraps) + prior_variance_init
	f.close()
	return dicti

def learn_iterative_variant_gene_tissue_prior_pip_level_bootstrapped(component_level_abf_summary_file, tgfm_version, tissue_names, per_window_abf_output_stem, max_iter=60, n_bootstraps=200, random_bootstrap=True, version='v1'):
	# Create mapping from tissue name to tissue position
	tissue_name_to_position = {}
	for ii, tissue_name in enumerate(tissue_names):
		tissue_name_to_position[tissue_name] = ii

	# Initialize prior probability emperical distribution
	variant_prob_distr = np.ones(n_bootstraps)*.1
	tissue_probs_distr = np.ones((len(tissue_names), n_bootstraps))*.1


	# Create mapping from window to class to indices
	window_to_class_to_indices, window_to_class_to_middle_indices, window_to_variant_middle = create_mapping_from_window_to_class_to_indices(component_level_abf_summary_file, tissue_name_to_position, tissue_names)

	# Get ordered window names
	window_names = np.asarray([*window_to_class_to_indices])
	n_windows = len(window_names)

	if random_bootstrap:
		# Get indices of windows corresponding to each bootstrap
		bs_indices = []
		# Generate bootstrap samples from full distribution
		for bs_iter in range(n_bootstraps):
			indices = np.random.choice(np.arange(n_windows), size=n_windows, replace=True)
			bs_indices.append(indices)
	else:
		# Get indices of windows corresponding to each bootstrap
		bs_indices = get_bs_indices_based_on_adjacent_window_groups(window_names, n_bootstraps)
	# bootstrapping mapping
	bs_mapping = []
	for bs_iter in range(n_bootstraps):
		indices = bs_indices[bs_iter]
		temper = {}
		for index in indices:
			if index not in temper:
				temper[index] = 1.0
			else:
				temper[index] = temper[index] + 1.0
		bs_mapping.append(temper)

	# Window to bootstrap vec
	window_to_bootstraps = {}
	for window_iter, window_name in enumerate(window_names):
		vec = np.zeros(n_bootstraps)
		for bs_iter in range(n_bootstraps):
			if window_iter in bs_mapping[bs_iter]:
				vec[bs_iter] = bs_mapping[bs_iter][window_iter]
		window_to_bootstraps[window_name] = vec

	window_component_to_prior_variance = create_mapping_from_window_component_to_prior_variance(component_level_abf_summary_file, n_bootstraps)

	# Iterative algorithm
	for itera in range(max_iter):
		print(itera)
		old_tissue_probs = np.copy(tissue_probs_distr)
		variant_prob_distr, tissue_probs_distr, window_component_to_prior_variance = update_prior_prob_for_variant_gene_tissue_bootstrapped(component_level_abf_summary_file, tgfm_version, tissue_name_to_position, variant_prob_distr, tissue_probs_distr, tissue_names, window_to_class_to_indices, window_to_class_to_middle_indices, n_bootstraps, bs_indices, bs_mapping, window_names, window_to_bootstraps, version, window_component_to_prior_variance)
		diff=tissue_probs_distr - old_tissue_probs
		print(np.sort(np.mean(diff,axis=1)))


	return variant_prob_distr, tissue_probs_distr

def get_number_of_annotations_from_component_level_summary_file(component_level_abf_summary_file):
	head_count = 0
	f = open(component_level_abf_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		anno_file = data[7]
		anno = np.load(anno_file)
		n_anno = anno.shape[1]
		break
	f.close()
	return n_anno

def learn_iterative_anno_variant_gene_tissue_prior_pip_level_bootstrapped(component_level_abf_summary_file, tgfm_version, tissue_names, per_window_abf_output_stem, max_iter=60, n_bootstraps=200):
	# Create mapping from tissue name to tissue position
	tissue_name_to_position = {}
	for ii, tissue_name in enumerate(tissue_names):
		tissue_name_to_position[tissue_name] = ii

	# Extract number of annotations from component summary file
	n_anno = get_number_of_annotations_from_component_level_summary_file(component_level_abf_summary_file)

	# Initialize prior probability emperical distribution
	# Initialize both probs to .1
	variant_prob_distr = np.zeros((n_anno, n_bootstraps))
	variant_prob_distr[0,:] = scipy.special.logit(.1)
	tissue_probs_distr = np.ones((len(tissue_names), n_bootstraps))*.1


	# Create mapping from window to class to indices
	window_to_class_to_indices, window_to_class_to_middle_indices, window_to_variant_middle = create_mapping_from_window_to_class_to_indices(component_level_abf_summary_file, tissue_name_to_position, tissue_names)

	# Get ordered window names
	window_names = np.asarray([*window_to_class_to_indices])
	n_windows = len(window_names)

	# Get indices of windows corresponding to each bootstrap
	bs_indices = []
	# Generate bootstrap samples from full distribution
	for bs_iter in range(n_bootstraps):
		indices = np.random.choice(np.arange(n_windows), size=n_windows, replace=True)
		bs_indices.append(indices)
	# bootstrapping mapping
	bs_mapping = []
	for bs_iter in range(n_bootstraps):
		indices = bs_indices[bs_iter]
		temper = {}
		for index in indices:
			if index not in temper:
				temper[index] = 1.0
			else:
				temper[index] = temper[index] + 1.0
		bs_mapping.append(temper)

	# Window to bootstrap vec
	window_to_bootstraps = {}
	for window_iter, window_name in enumerate(window_names):
		vec = np.zeros(n_bootstraps)
		for bs_iter in range(n_bootstraps):
			if window_iter in bs_mapping[bs_iter]:
				vec[bs_iter] = bs_mapping[bs_iter][window_iter]
		window_to_bootstraps[window_name] = vec



	# Iterative algorithm
	for itera in range(max_iter):
		print('###########################################')
		print('###########################################')
		print(itera)
		old_tissue_probs = np.copy(tissue_probs_distr)
		for tissue_iter, tissue_name in enumerate(tissue_names):
			print(tissue_name + '\t' + str(np.mean(old_tissue_probs[tissue_iter, :])) + '\t' + ';'.join(old_tissue_probs[tissue_iter, :].astype(str)))
		variant_prob_distr, tissue_probs_distr = update_prior_prob_for_anno_variant_gene_tissue_bootstrapped(component_level_abf_summary_file, tgfm_version, tissue_name_to_position, variant_prob_distr, tissue_probs_distr, tissue_names, window_to_class_to_indices, window_to_class_to_middle_indices, window_to_variant_middle, n_bootstraps, bs_indices, bs_mapping, window_names, window_to_bootstraps, n_anno)
		diff=tissue_probs_distr - old_tissue_probs
		print(np.mean(diff,axis=1))
		#print(np.mean(tissue_probs_distr,axis=1))


	return variant_prob_distr, tissue_probs_distr



def compute_expected_pips(alpha_mats):
	if len(alpha_mats) == 0:
		print('assumption errorr')
		pdb.set_trace()
	# Initialize expected pips
	pips = np.ones(alpha_mats[0].shape)
	for component_iter in range (len(alpha_mats)):
		pips = pips*(1.0-alpha_mats[component_iter])
	pips = 1.0 - pips

	return pips

def compute_expected_pips_from_sampled_distribution(alpha_mats, window_samples, n_bootstraps):
	if len(alpha_mats) == 0:
		print('assumption eroror')
		pdb.set_trace()
	if len(alpha_mats) != len(window_samples):
		print('assumption error')
		pdb.set_trace()
	unique_window_samples_to_pips = {}
	for window_sample in np.unique(window_samples):
		unique_window_samples_to_pips[window_sample] = np.ones(alpha_mats[0].shape)


	for component_iter, sample_name in enumerate(window_samples):
		alpha_mat = alpha_mats[component_iter]
		unique_window_samples_to_pips[sample_name] = (unique_window_samples_to_pips[sample_name])*(1.0-alpha_mat)


	expected_pips = np.zeros(alpha_mats[0].shape)
	for window_sample in np.unique(window_samples):
		expected_pips = expected_pips + (1.0 - unique_window_samples_to_pips[window_sample])
	expected_pips = expected_pips/len(np.unique(window_samples))

	weights = len(np.unique(window_samples))/n_bootstraps

	return expected_pips, weights

def update_prior_prob_for_variant_gene_tissue_bootstrapped(component_level_abf_summary_file, tgfm_version, tissue_name_to_position, variant_prob_distr, tissue_probs_distr, tissue_names, window_to_class_to_indices, window_to_class_to_middle_indices, n_bootstraps, bs_indices, bs_mapping, window_names, window_to_bootstraps, version, window_component_to_prior_variance):
	# General info
	n_tissues = len(tissue_names)	
	n_windows = len(window_names)

	# Keep track of counts in each window
	variant_posterior_sum = np.zeros(n_bootstraps)
	variant_counts = np.zeros(n_bootstraps)
	tissue_posterior_sum = np.zeros((n_tissues, n_bootstraps))
	tissue_counts = np.zeros((n_tissues, n_bootstraps))

	# Create prior hash
	prior_hash = {}
	prior_hash['variant'] = variant_prob_distr
	for ii in range(n_tissues):
		prior_hash[tissue_names[ii]] = tissue_probs_distr[ii,:]

	# Create mapping from window name to expected prior probs
	window_to_prior_probs = {}
	for window_name in window_names:
		prior_prob_raw = np.zeros((window_to_class_to_indices[window_name]['n_elements'], n_bootstraps))
		class_name = 'variant'
		class_indices = window_to_class_to_indices[window_name][class_name]

		prior_prob_raw[class_indices, :] = prior_hash[class_name]
		for class_name in tissue_names:
			class_indices = window_to_class_to_indices[window_name][class_name]
			if len(class_indices) == 0:
				continue
			prior_prob_raw[class_indices, :] = prior_hash[class_name]
		# Normalize each bootstrapped sample
		prior_probs = prior_prob_raw/np.sum(prior_prob_raw,axis=0)

		# Create mapping
		window_to_prior_probs[window_name] = prior_probs


	# Loop through components
	head_count = 0
	prev_window_name = 'NULL'
	f = open(component_level_abf_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		#print(data[0])
		window_name = data[0]
		component_number = data[2]

		if window_name != prev_window_name:
			if prev_window_name == 'NULL':
				window_counter = -1
			else:
				# need to tally up results
				if tgfm_version == 'pmces':
					expected_pips = compute_expected_pips(alpha_mats)
				elif tgfm_version == 'sampler':
					expected_pips, weights = compute_expected_pips_from_sampled_distribution(alpha_mats, window_samples, n_bootstraps)
				else:
					print('assumption error: tgfm version currently not implemented')
					pdb.set_trace()
				#if version == 'v4':
					#expected_pips[expected_pips < .5] = 0.0
				# update global counts
				# For variants
				indices = np.asarray(window_to_class_to_middle_indices[prev_window_name]['variant'])
				bs_scaling_factors = window_to_bootstraps[prev_window_name]
				variant_posterior_sum = variant_posterior_sum + (bs_scaling_factors*np.sum(expected_pips[indices,:],axis=0))
				variant_counts = variant_counts + (bs_scaling_factors*len(indices))
				# For genes
				for tiss_iter, tissue_name in enumerate(tissue_names):
					indices = np.asarray(window_to_class_to_middle_indices[prev_window_name][tissue_name])
					if len(indices) == 0:
						continue
					#if tissue_name == 'Whole_Blood' and np.sum(np.sum(expected_pips[indices,:],axis=0) > 0.1) > 0:
					#	pdb.set_trace()
					tissue_posterior_sum[tiss_iter, :] = tissue_posterior_sum[tiss_iter, :] + (bs_scaling_factors*np.sum(expected_pips[indices,:],axis=0))
					tissue_counts[tiss_iter, :] = tissue_counts[tiss_iter, :] + (bs_scaling_factors*len(indices))
			window_abfs = np.load(data[5])
			window_samples = np.load(data[6])
			window_mu_init = np.load(data[9])
			window_mu_var_init = np.load(data[10])
			component_counter = 0
			prev_window_name = window_name
			window_counter = window_counter + 1
			alpha_mats = []

		mu_init = window_mu_init[component_counter,:]
		mu_var_init = window_mu_var_init[component_counter, :]
		component_variance_init = float(data[8])
		component_variances_current = window_component_to_prior_variance[window_name + '_' + component_number]

		# Normalize prior probabilities
		prior_probs = window_to_prior_probs[window_name]


		# Reconstruct mu and mu_var given current component variances
		if component_variance_init == 0.0:
			component_variance_init = 1e-300
			ele_labfs = window_abfs[component_counter, :]
			alpha_mat = extract_alpha_vec_given_lbf_and_prior_matrix(ele_labfs, prior_probs)
			alpha_mats.append(alpha_mat)
			component_counter = component_counter + 1
		else:
			b_terms_init = mu_init/mu_var_init
			a_terms_init = -(1.0/(2.0*mu_var_init))
			a_terms_updated_tmp = a_terms_init + (1.0/(2.0*component_variance_init))

			alpha_mat = []
			'''
			import time
			t1 = time.time()
			for bs_iter in range(len(component_variances_current)):
				a_terms_updated = a_terms_updated_tmp - (1.0/(2.0*component_variances_current[bs_iter]))
				mu_var_updated = -1.0/(2.0*a_terms_updated)
				mu_updated = b_terms_init*mu_var_updated

				un_normalized_weights = (.5*np.log(component_variances_current[bs_iter])) + (.5*np.square(mu_updated)/mu_var_updated) + (.5*np.log(mu_var_updated))
				bs_alpha_vec = extract_alpha_vec_given_lbf_and_prior(un_normalized_weights, prior_probs[:,bs_iter])
				alpha_mat.append(bs_alpha_vec)
				# Update prior variance
				new_componenet_variance = np.sum((np.square(mu_updated) + mu_var_updated)*bs_alpha_vec)
				component_variances_current[bs_iter] = new_componenet_variance
			t2 = time.time()
			'''
			a_terms_updated = np.transpose(np.tile(a_terms_updated_tmp,(n_bootstraps, 1))) - (1.0/(2.0*component_variances_current))
			mu_var_updated = -1.0/(2.0*a_terms_updated)
			mu_updated = np.transpose(np.transpose(mu_var_updated)*b_terms_init)
			un_normalized_weights = (.5*np.log(component_variances_current)) + (.5*np.square(mu_updated)/mu_var_updated) + (.5*np.log(mu_var_updated))
			alpha_mat = extract_alpha_mat_given_lbf_and_prior_matrix(un_normalized_weights, prior_probs)

			component_variances_current = np.sum(alpha_mat*(np.square(mu_updated) + mu_var_updated),axis=0)


			#alpha_mat = np.transpose(np.asarray(alpha_mat))

			#ele_labfs = window_abfs[component_counter, :]
			component_counter = component_counter + 1

			# Normalize prior probabilities
			#prior_probs = window_to_prior_probs[window_name]
		
			# Get alpha vec given lbf and prior probs
			#alpha_mat2 = extract_alpha_vec_given_lbf_and_prior_matrix(ele_labfs, prior_probs)
			alpha_mats.append(alpha_mat)
			window_component_to_prior_variance[window_name + '_' + component_number] = component_variances_current

	f.close()

	# need to tally up results
	# need to tally up results
	if tgfm_version == 'pmces':
		expected_pips = compute_expected_pips(alpha_mats)
	elif tgfm_version == 'sampler':
		expected_pips, weights = compute_expected_pips_from_sampled_distribution(alpha_mats, window_samples, n_bootstraps)
	else:
		print('assumption error: tgfm version currently not implemented')
		pdb.set_trace()
	#if version == 'v4':
		#expected_pips[expected_pips < .5] = 0.0
	# update global counts
	# For variants
	indices = np.asarray(window_to_class_to_middle_indices[window_name]['variant'])
	bs_scaling_factors = window_to_bootstraps[window_name]
	variant_posterior_sum = variant_posterior_sum + (bs_scaling_factors*np.sum(expected_pips[indices,:],axis=0))
	variant_counts = variant_counts + (bs_scaling_factors*len(indices))
	# For genes
	for tiss_iter, tissue_name in enumerate(tissue_names):
		indices = np.asarray(window_to_class_to_middle_indices[window_name][tissue_name])
		if len(indices) == 0:
			continue
		tissue_posterior_sum[tiss_iter, :] = tissue_posterior_sum[tiss_iter, :] + (bs_scaling_factors*np.sum(expected_pips[indices,:],axis=0))
		tissue_counts[tiss_iter, :] = tissue_counts[tiss_iter, :] + (bs_scaling_factors*len(indices))


	variant_prob_distr = variant_posterior_sum/variant_counts
	tissue_counts[tissue_counts == 0.0] = .1
	tissue_probs_distr = tissue_posterior_sum/tissue_counts

	return variant_prob_distr, tissue_probs_distr, window_component_to_prior_variance


def logistic_regression_log_likelihood(beta, X, y):
	# get beta into correct format
	#beta_format = np.transpose(np.asmatrix(beta))
	#x_beta = np.dot(X, beta_format)
	# Compute probability each data point belongs to class 1
	#prob_1 = np.squeeze(np.asarray(scipy.special.expit(x_beta)))
	prob_1 = scipy.special.expit(np.dot(X,beta))
	# Compute probability each data point belongs to class 0
	prob_0 = 1.0 - prob_1
	# Return -log likelihood
	return -np.sum(np.log(prob_1)*y + np.log(prob_0)*(1-y))
	# Compute probability each data point belongs to class 1
	#log_prob_1 = np.squeeze(np.asarray(scipy.special.log_expit(x_beta)))
	# Compute probability each data point belongs to class 0
	#log_prob_0 = -np.squeeze(np.asarray(x_beta)) + log_prob_1

	# Old
	#prob_1 = np.squeeze(np.asarray(scipy.special.expit(x_beta)))
	#prob_0 = 1.0 - prob_1
	#old_log_prob_1 = np.log(prob_1)
	#old_log_prob_0 = np.log(prob_0)

	# Return -log likelihood
	#return -np.sum(log_prob_1*y + log_prob_0*(1-y))

# Take the gradient of beta
def logistic_regression_gradient(beta, X, y):
	# get beta into correct format (a matrix)
	#beta_format = np.transpose(np.asmatrix(beta))
	#x_beta = np.dot(X, beta_format)
	# Compute probability each data point belongs to class 1
	#prob_1 = np.squeeze(np.asarray(scipy.special.expit(x_beta)))
	prob_1 = scipy.special.expit(np.dot(X,beta))

	# Initialize gradient
	#gradient = np.zeros(X.shape[1])
	# Compute gradient at each element
	#for index in range(len(gradient)):
		#gradient[index] = np.sum(X[:,index]*y) - np.sum(X[:,index]*prob_1)
	gradient = np.dot(np.transpose(X),y-prob_1)
	return -gradient


def manual_logistic_regression(X,y, beta_init, fraction=1.0):
	beta = np.zeros(X.shape[1])  # Initialize betas to zero vec

	if fraction < 1.0:
		stochastic_indices = np.random.choice(np.arange(len(y)), size=int(np.floor(len(y)*fraction)), replace=False)
		X = X[stochastic_indices, :]
		y = y[stochastic_indices]

	# Run LBFGS
	## This uses 'logistic_regression_gradient' and 'logistic_regression_log_likelihood functions'
	#val = scipy.optimize.fmin_l_bfgs_b(logistic_regression_log_likelihood, beta,fprime=logistic_regression_gradient, args=(X,y))
	#t1 = time.time()
	#val = scipy.optimize.fmin_l_bfgs_b(logistic_regression_log_likelihood, beta,fprime=logistic_regression_gradient, args=(X,y))
	#t2 = time.time()
	res2=rcpp_num_R_pkg.fastLR(x=X,y=y, start=beta_init)
	#t3 = time.time()
	#pdb.set_trace()
	'''
	if val[2]['warnflag'] != 0:
		print('Scipy optimization did not converge successfully')
	else:
		print('Successful convergence')
		print('Negative log likelihood: ' + str(val[1]))
		print('Coefficient vector: ' + str('\t'.join(val[0].astype(str))))
	'''
	return res2.rx2('coefficients')


def update_prior_prob_for_anno_variant_gene_tissue_bootstrapped(component_level_abf_summary_file, tgfm_version, tissue_name_to_position, variant_prob_distr, tissue_probs_distr, tissue_names, window_to_class_to_indices, window_to_class_to_middle_indices, window_to_variant_middle, n_bootstraps, bs_indices, bs_mapping, window_names, window_to_bootstraps, n_anno):
	# General info
	n_tissues = len(tissue_names)	
	n_windows = len(window_names)

	# Keep track of counts in each window
	#variant_posterior_sum = np.zeros(n_bootstraps)
	#variant_counts = np.zeros(n_bootstraps)
	window_variant_pip_arr = []
	window_variant_anno_arr = []
	tissue_posterior_sum = np.zeros((n_tissues, n_bootstraps))
	tissue_counts = np.zeros((n_tissues, n_bootstraps))

	# Create prior hash
	prior_hash = {}
	prior_hash['variant'] = variant_prob_distr
	for ii in range(n_tissues):
		prior_hash[tissue_names[ii]] = tissue_probs_distr[ii,:]

	'''
	# Create mapping from window name to expected prior probs
	window_to_prior_probs = {}
	for window_name in window_names:
		prior_prob_raw = np.zeros((window_to_class_to_indices[window_name]['n_elements'], n_bootstraps))
		class_name = 'variant'
		class_indices = window_to_class_to_indices[window_name][class_name]

		prior_prob_raw[class_indices, :] = prior_hash[class_name]
		for class_name in tissue_names:
			class_indices = window_to_class_to_indices[window_name][class_name]
			if len(class_indices) == 0:
				continue
			prior_prob_raw[class_indices, :] = prior_hash[class_name]
		# Normalize each bootstrapped sample
		prior_probs = prior_prob_raw/np.sum(prior_prob_raw,axis=0)

		# Create mapping
		window_to_prior_probs[window_name] = prior_probs
	'''

	# Loop through components
	head_count = 0
	prev_window_name = 'NULL'
	f = open(component_level_abf_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		#print(data[0])
		window_name = data[0]

		if window_name != prev_window_name:
			if prev_window_name == 'NULL':
				window_counter = -1
			else:
				# need to tally up results
				if tgfm_version == 'pmces':
					expected_pips = compute_expected_pips(alpha_mats)
				elif tgfm_version == 'sampler':
					expected_pips, weights = compute_expected_pips_from_sampled_distribution(alpha_mats, window_samples, n_bootstraps)
				else:
					print('assumption error: tgfm version currently not implemented')
					pdb.set_trace()
				# update global counts
				# For variants
				indices = np.asarray(window_to_class_to_middle_indices[prev_window_name]['variant'])
				bs_scaling_factors = window_to_bootstraps[prev_window_name]
				window_variant_pip_arr.append(expected_pips[indices,:])
				window_variant_anno_arr.append(window_anno[window_to_variant_middle[prev_window_name], :])
				#variant_posterior_sum = variant_posterior_sum + (bs_scaling_factors*np.sum(expected_pips[indices,:],axis=0))
				#variant_counts = variant_counts + (bs_scaling_factors*len(indices))
				# For genes
				for tiss_iter, tissue_name in enumerate(tissue_names):
					indices = np.asarray(window_to_class_to_middle_indices[prev_window_name][tissue_name])
					if len(indices) == 0:
						continue
					tissue_posterior_sum[tiss_iter, :] = tissue_posterior_sum[tiss_iter, :] + (bs_scaling_factors*np.sum(expected_pips[indices,:],axis=0))
					tissue_counts[tiss_iter, :] = tissue_counts[tiss_iter, :] + (bs_scaling_factors*len(indices))
			window_abfs = np.load(data[5])
			window_samples = np.load(data[6])
			window_anno = np.load(data[7])
			window_genomic_anno_prob_pred = scipy.special.expit(np.dot(window_anno, variant_prob_distr))
			component_counter = 0
			prev_window_name = window_name
			window_counter = window_counter + 1
			alpha_mats = []

		ele_labfs = window_abfs[component_counter, :]
		component_counter = component_counter + 1

		# Extract prior probabilities
		#prior_probs = window_to_prior_probs[window_name]
		prior_prob_raw = np.zeros((window_to_class_to_indices[window_name]['n_elements'], n_bootstraps))
		class_name = 'variant'
		class_indices = window_to_class_to_indices[window_name][class_name]
		prior_prob_raw[class_indices, :] = window_genomic_anno_prob_pred
		for class_name in tissue_names:
			class_indices = window_to_class_to_indices[window_name][class_name]
			if len(class_indices) == 0:
				continue
			prior_prob_raw[class_indices, :] = prior_hash[class_name]
		# Normalize each bootstrapped sample
		prior_probs = prior_prob_raw/np.sum(prior_prob_raw,axis=0)

		
		# Get alpha vec given lbf and prior probs
		alpha_mat = extract_alpha_vec_given_lbf_and_prior_matrix(ele_labfs, prior_probs)
		alpha_mats.append(alpha_mat)


	f.close()

	# need to tally up results
	# need to tally up results
	if tgfm_version == 'pmces':
		expected_pips = compute_expected_pips(alpha_mats)
	elif tgfm_version == 'sampler':
		expected_pips, weights = compute_expected_pips_from_sampled_distribution(alpha_mats, window_samples, n_bootstraps)
	else:
		print('assumption error: tgfm version currently not implemented')
		pdb.set_trace()
	# update global counts
	# For variants
	indices = np.asarray(window_to_class_to_middle_indices[window_name]['variant'])
	bs_scaling_factors = window_to_bootstraps[window_name]
	window_variant_pip_arr.append(expected_pips[indices,:])
	window_variant_anno_arr.append(window_anno[window_to_variant_middle[window_name], :])
	#variant_posterior_sum = variant_posterior_sum + (bs_scaling_factors*np.sum(expected_pips[indices,:],axis=0))
	#variant_counts = variant_counts + (bs_scaling_factors*len(indices))
	# For genes
	for tiss_iter, tissue_name in enumerate(tissue_names):
		indices = np.asarray(window_to_class_to_middle_indices[window_name][tissue_name])
		if len(indices) == 0:
			continue
		tissue_posterior_sum[tiss_iter, :] = tissue_posterior_sum[tiss_iter, :] + (bs_scaling_factors*np.sum(expected_pips[indices,:],axis=0))
		tissue_counts[tiss_iter, :] = tissue_counts[tiss_iter, :] + (bs_scaling_factors*len(indices))

	#variant_prob_distr = variant_posterior_sum/variant_counts
	tissue_probs_distr = tissue_posterior_sum/tissue_counts

	# Now loop through bootstraps.
	# For each bootstap update variant_prob_distr
	for bs_iter in range(n_bootstraps):
		# Indices of windows corresponding to this bootstrap
		bs_window_indices = bs_indices[bs_iter]
		# Extract logistic regression input and output
		bs_anno = []
		bs_pip = []
		for bs_window_index in bs_window_indices:
			bs_anno.append(window_variant_anno_arr[bs_window_index])
			bs_pip.append(window_variant_pip_arr[bs_window_index][:,bs_iter])
		# Clean up data
		bs_anno = np.vstack(bs_anno)
		bs_pip = np.hstack(bs_pip)
		# Quick error check
		if bs_anno.shape[0] != len(bs_pip):
			print('assumptioner roror')
			pdb.set_trace()

		coefs = manual_logistic_regression(bs_anno,bs_pip, variant_prob_distr[:, bs_iter], fraction=1.0)
		variant_prob_distr[:, bs_iter] = coefs



	return variant_prob_distr, tissue_probs_distr







trait_name = sys.argv[1]
new_tgfm_stem = sys.argv[2]
tgfm_version = sys.argv[3]
processed_tgfm_input_stem = sys.argv[4]
gtex_pseudotissue_file = sys.argv[5]
tgfm_input_summary_file = sys.argv[6]
iterative_tgfm_prior_results_dir = sys.argv[7]

perm_iterative_prior_results = iterative_tgfm_prior_results_dir + new_tgfm_stem.split('/')[-1] 


#Extract tissue names
tissue_names = extract_pseudotissue_names(gtex_pseudotissue_file, ignore_testis=True)

###################################################
# Concatenate PIP summary file across parallel runs (one line for each window)
###################################################
suffix = 'tgfm_pip_summary.txt'
num_jobs=8
concatenated_pip_summary_file = new_tgfm_stem + '_iterative_prior_w_ard_' + suffix
concatenate_results_across_parallel_jobs(new_tgfm_stem, suffix, num_jobs, concatenated_pip_summary_file)


###################################################
# Create component level summary data
###################################################
component_level_abf_summary_file = new_tgfm_stem + '_iterative_prior_w_ard_v2' + '_tgfm_component_level_abf_summary.txt'
per_window_abf_output_stem = new_tgfm_stem + '_iterative_prior_w_ard_v2_per_window_abf_'
generate_component_level_abf_summary_data(concatenated_pip_summary_file, component_level_abf_summary_file, tissue_names, new_tgfm_stem, tgfm_version, processed_tgfm_input_stem, per_window_abf_output_stem, version='v2')


###################################################
# Learn iterative distribution variant-gene-tissue prior (bootstrapping ci intervals)
###################################################

n_bootstraps=100
variant_prob_emperical_distr, tissue_probs_emperical_distr = learn_iterative_variant_gene_tissue_prior_pip_level_bootstrapped(component_level_abf_summary_file, tgfm_version, tissue_names, per_window_abf_output_stem, max_iter=30, n_bootstraps=n_bootstraps)

# Print to output
variant_gene_distr_prior_output_file = perm_iterative_prior_results + '_iterative_variant_gene_prior_w_ard_v2_pip_level_bootstrapped.txt'
t = open(variant_gene_distr_prior_output_file,'w')
t.write('element_name\tprior\texp_E_ln_prior\tprior_distribution\n')
t.write('variant\t' + str(np.mean(variant_prob_emperical_distr)) + '\t' + str(np.exp(np.mean(np.log(variant_prob_emperical_distr)))) + '\t' + ';'.join(variant_prob_emperical_distr.astype(str)) + '\n')
for tiss_iter, tissue_name in enumerate(tissue_names):
	t.write(tissue_name + '\t' + str(np.mean(tissue_probs_emperical_distr[tiss_iter,:])) + '\t' + str(np.exp(np.mean(np.log(tissue_probs_emperical_distr[tiss_iter,:])))) + '\t' + ';'.join(tissue_probs_emperical_distr[tiss_iter,:].astype(str)) + '\n')
t.close()
print(variant_gene_distr_prior_output_file)
# Delete unneccessary files
os.system('rm ' + component_level_abf_summary_file)
os.system('rm ' + per_window_abf_output_stem + '*')




