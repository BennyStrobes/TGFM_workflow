import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import pickle
import tgfm




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

def update_prior_probs_for_variant_gene_tissue(component_level_abf_summary_file, tgfm_version, tissue_name_to_position, variant_prob, tissue_probs, tissue_names, window_to_class_to_indices):
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
	window_to_class_to_indices = create_mapping_from_window_to_class_to_indices(component_level_abf_summary_file, tissue_name_to_position, tissue_names)


	for itera in range(max_iter):
		variant_prob, tissue_probs = update_prior_probs_for_variant_gene_tissue(component_level_abf_summary_file, tgfm_version, tissue_name_to_position, variant_prob, tissue_probs, tissue_names, window_to_class_to_indices)
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



	# Iterative algorithm
	import time
	for itera in range(max_iter):
		print(itera)
		t1 = time.time()
		old_tissue_probs = np.copy(tissue_probs_distr)
		variant_prob_distr, tissue_probs_distr = update_prior_prob_for_variant_gene_tissue_bootstrapped(component_level_abf_summary_file, tgfm_version, tissue_name_to_position, variant_prob_distr, tissue_probs_distr, tissue_names, window_to_class_to_indices, window_to_class_to_middle_indices, n_bootstraps, bs_indices, bs_mapping, window_names, window_to_bootstraps, version)
		diff=tissue_probs_distr - old_tissue_probs
		print(np.sort(np.mean(diff,axis=1)))
		t2 = time.time()
		print(t2-t1)


	return variant_prob_distr, tissue_probs_distr

def extract_alpha_vec_given_lbf_and_prior_matrix(lbf, prior_probs):
	maxlbf = np.max(lbf)
	ww = np.exp(lbf - maxlbf)
	w_weighted = (np.transpose(prior_probs)*ww)
	alpha = np.transpose(w_weighted)/np.sum(w_weighted,axis=1)
	return alpha

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


def update_prior_prob_for_variant_gene_tissue_bootstrapped(component_level_abf_summary_file, tgfm_version, tissue_name_to_position, variant_prob_distr, tissue_probs_distr, tissue_names, window_to_class_to_indices, window_to_class_to_middle_indices, n_bootstraps, bs_indices, bs_mapping, window_names, window_to_bootstraps, version):
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
			component_counter = 0
			prev_window_name = window_name
			window_counter = window_counter + 1
			alpha_mats = []

		ele_labfs = window_abfs[component_counter, :]
		component_counter = component_counter + 1

		# Normalize prior probabilities
		prior_probs = window_to_prior_probs[window_name]
		
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

	return variant_prob_distr, tissue_probs_distr


def generate_component_level_abf_summary_data(concatenated_pip_summary_file, component_level_abf_summary_file, tissue_names, file_stem, model_version, per_window_abf_output_stem, version='v1'):
	f = open(concatenated_pip_summary_file)
	t = open(component_level_abf_summary_file,'w')
	t.write('window_name\tbootstrap_number\tcomponent_number\telement_names\tmiddle_element_names\twindow_abf_file\twindow_component_samples_file\tannotation_mat_file\n')
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
		if os.path.isfile(window_pkl_file) == False:
			continue
		# Load in tgfm results data
		g = open(window_pkl_file, "rb")
		tgfm_results = pickle.load(g)
		g.close()

		# Load in tgfm results data
		# Load in tgfm results data
		tgfm_input_data = data[2]
		g = open(tgfm_input_data, "rb")
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
		per_window_abf = []
		window_has_component = False
		per_window_samples = []
		if model_version.startswith('pmces'):
			lines = []
			if version == 'v1':
				for component_iter in tgfm_results['valid_components']:
					if component_in_middle_of_window(tgfm_results['alpha_phi'][component_iter, :], tgfm_results['beta_phi'][component_iter, :], middle_gene_indices_dicti, middle_variant_indices_dicti):
						# Extract log abf
						log_abf = np.hstack((tgfm_results['alpha_lbf'][component_iter, :], tgfm_results['beta_lbf'][component_iter, :]))
						element_names = np.hstack((tgfm_results['genes'], tgfm_results['variants']))
						middle_element_names = np.hstack((tgfm_results['genes'][tgfm_results['middle_gene_indices']], tgfm_results['variants'][tgfm_results['middle_variant_indices']]))
						t.write(window_name + '\t' + 'NA' + '\t' + str(component_iter) + '\t' + ';'.join(element_names) + '\t' + ';'.join(middle_element_names) + '\t' + per_window_abf_output_file + '\t' + per_window_sample_names_output_file + '\t' + per_window_anno_output_file + '\n')
						per_window_abf.append(log_abf)
						per_window_samples.append(0)
						window_has_component = True
			elif version == 'v2' or version == 'v4':
				for component_iter in range(tgfm_results['alpha_lbf'].shape[0]):
					# V2
					log_abf = np.hstack((tgfm_results['alpha_lbf'][component_iter, :], tgfm_results['beta_lbf'][component_iter, :]))
					element_names = np.hstack((tgfm_results['genes'], tgfm_results['variants']))
					middle_element_names = np.hstack((tgfm_results['genes'][tgfm_results['middle_gene_indices']], tgfm_results['variants'][tgfm_results['middle_variant_indices']]))
					liner = window_name + '\t' + 'NA' + '\t' + str(component_iter) + '\t' + ';'.join(element_names) + '\t' + ';'.join(middle_element_names) + '\t' + per_window_abf_output_file + '\t' + per_window_sample_names_output_file + '\t' + per_window_anno_output_file + '\n'
					lines.append(liner)
					per_window_abf.append(log_abf)
					per_window_samples.append(0)
					if component_iter in tgfm_results['valid_components']:
						if component_in_middle_of_window(tgfm_results['alpha_phi'][component_iter, :], tgfm_results['beta_phi'][component_iter, :], middle_gene_indices_dicti, middle_variant_indices_dicti):
							window_has_component = True	
				if window_has_component:
					for liner in lines:
						t.write(liner)
			elif version == 'v3':
				for component_iter in tgfm_results['valid_components']:
					# V2
					log_abf = np.hstack((tgfm_results['alpha_lbf'][component_iter, :], tgfm_results['beta_lbf'][component_iter, :]))
					element_names = np.hstack((tgfm_results['genes'], tgfm_results['variants']))
					middle_element_names = np.hstack((tgfm_results['genes'][tgfm_results['middle_gene_indices']], tgfm_results['variants'][tgfm_results['middle_variant_indices']]))
					liner = window_name + '\t' + 'NA' + '\t' + str(component_iter) + '\t' + ';'.join(element_names) + '\t' + ';'.join(middle_element_names) + '\t' + per_window_abf_output_file + '\t' + per_window_sample_names_output_file + '\t' + per_window_anno_output_file + '\n'
					lines.append(liner)
					per_window_abf.append(log_abf)
					per_window_samples.append(0)
					if component_iter in tgfm_results['valid_components']:
						if component_in_middle_of_window(tgfm_results['alpha_phi'][component_iter, :], tgfm_results['beta_phi'][component_iter, :], middle_gene_indices_dicti, middle_variant_indices_dicti):
							window_has_component = True	
				if window_has_component:
					for liner in lines:
						t.write(liner)
		elif model_version.startswith('sampler'):
			n_samples = len(tgfm_results['valid_components'])
			for sample_iter in range(n_samples):
				sample_valid_components = tgfm_results['valid_components'][sample_iter]
				for component_iter in sample_valid_components:
					if component_in_middle_of_window(tgfm_results['alpha_phis'][component_iter][sample_iter,:], tgfm_results['beta_phis'][component_iter][sample_iter,:], middle_gene_indices_dicti, middle_variant_indices_dicti):
						# Extract log abf
						log_abf = np.hstack((tgfm_results['alpha_lbfs'][component_iter][sample_iter, :], tgfm_results['beta_lbfs'][component_iter][sample_iter, :]))
						element_names = np.hstack((tgfm_results['genes'], tgfm_results['variants']))
						middle_element_names = np.hstack((tgfm_results['genes'][tgfm_results['middle_gene_indices']], tgfm_results['variants'][tgfm_results['middle_variant_indices']]))
						t.write(window_name + '\t' + str(sample_iter) + '\t' + str(component_iter) + '\t' + ';'.join(element_names) + '\t' + ';'.join(middle_element_names) + '\t' + per_window_abf_output_file + '\t' + per_window_sample_names_output_file + '\t' + per_window_anno_output_file+ '\n')
						per_window_abf.append(log_abf)
						window_has_component = True
						per_window_samples.append(sample_iter)
		else:
			print('assumption eroror: unknown model verison')
			pdb.set_trace()
		# Convert to np array and save
		if window_has_component:
			per_window_abf = np.asarray(per_window_abf)
			np.save(per_window_abf_output_file, per_window_abf)
			per_window_samples = np.asarray(per_window_samples)
			np.save(per_window_sample_names_output_file, per_window_samples)
			#np.save(per_window_anno_output_file, tgfm_data['annotation'])
	t.close()
	f.close()
	return




tgfm_input_summary_file = sys.argv[1]
new_tgfm_stem = sys.argv[2]
tgfm_version = sys.argv[3]



#Extract tissue names
tissue_names = []
for tiss_iter in range(10):
	tissue_names.append('tissue' + str(tiss_iter))
tissue_names = np.asarray(tissue_names)



###################################################
# Create component level summary data
###################################################
component_level_abf_summary_file = new_tgfm_stem + '_iterative_prior_pip_level' + '_tgfm_component_level_abf_summary_bs.txt'
per_window_abf_output_stem = new_tgfm_stem + '_iterative_prior_pip_level_per_window_abf_bs_'
generate_component_level_abf_summary_data(tgfm_input_summary_file, component_level_abf_summary_file, tissue_names, new_tgfm_stem, tgfm_version, per_window_abf_output_stem, version='v2')



###################################################
# Learn iterative distribution variant-gene-tissue prior (doing non-distribution based prior)
###################################################
n_bootstraps=100
variant_prob_emperical_distr, tissue_probs_emperical_distr = learn_iterative_variant_gene_tissue_prior_pip_level_bootstrapped(component_level_abf_summary_file, tgfm_version, tissue_names, per_window_abf_output_stem, max_iter=3, n_bootstraps=n_bootstraps)




# Print to output
variant_gene_distr_prior_output_file = new_tgfm_stem + '_iterative_variant_gene_prior_pip_level_bootstrapped.txt'
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

