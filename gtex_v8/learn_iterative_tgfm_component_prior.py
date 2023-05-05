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

def generate_component_level_abf_summary_data(concatenated_pip_summary_file, component_level_abf_summary_file, tissue_names, file_stem, model_version, processed_tgfm_input_stem, per_window_abf_output_stem):
	f = open(concatenated_pip_summary_file)
	t = open(component_level_abf_summary_file,'w')
	t.write('window_name\tbootstrap_number\tcomponent_number\telement_names\twindow_abf_file\n')
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
		per_window_abf = []
		window_has_component = False

		if model_version.startswith('pmces'):
			for component_iter in tgfm_results['valid_components']:
				if component_in_middle_of_window(tgfm_results['alpha_phi'][component_iter, :], tgfm_results['beta_phi'][component_iter, :], middle_gene_indices_dicti, middle_variant_indices_dicti):
					
					# Extract log abf
					log_abf = np.hstack((tgfm_results['alpha_lbf'][component_iter, :], tgfm_results['beta_lbf'][component_iter, :]))
					element_names = np.hstack((tgfm_results['genes'], tgfm_results['variants']))
					t.write(window_name + '\t' + 'NA' + '\t' + str(component_iter) + '\t' + ';'.join(element_names) + '\t' + per_window_abf_output_file + '\n')
					per_window_abf.append(log_abf)
					window_has_component = True
		elif model_version.startswith('sampler'):
			n_samples = len(tgfm_results['valid_components'])
			for sample_iter in range(n_samples):
				sample_valid_components = tgfm_results['valid_components'][sample_iter]
				for component_iter in sample_valid_components:
					if component_in_middle_of_window(tgfm_results['alpha_phis'][component_iter][sample_iter,:], tgfm_results['beta_phis'][component_iter][sample_iter,:], middle_gene_indices_dicti, middle_variant_indices_dicti):
						# Extract log abf
						log_abf = np.hstack((tgfm_results['alpha_lbfs'][component_iter][sample_iter, :], tgfm_results['beta_lbfs'][component_iter][sample_iter, :]))
						element_names = np.hstack((tgfm_results['genes'], tgfm_results['variants']))
						t.write(window_name + '\t' + str(sample_iter) + '\t' + str(component_iter) + '\t' + ';'.join(element_names) + '\t' + per_window_abf_output_file + '\n')
						per_window_abf.append(log_abf)
						window_has_component = True
		else:
			print('assumption eroror: unknown model verison')
			pdb.set_trace()
		# Convert to np array and save
		if window_has_component:
			per_window_abf = np.asarray(per_window_abf)
			np.save(per_window_abf_output_file, per_window_abf)
	t.close()
	f.close()
	return

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
		sub_dicti['variant'] = []
		for tissue_name in tissue_names:
			sub_dicti[tissue_name] = []

		# Extract element names
		ele_names = np.asarray(data[3].split(';'))

		# loop through elements and add to sub-dictionary
		for ele_iter, ele_name in enumerate(ele_names):
			if ele_name.startswith('ENSG'):
				ele_category = '_'.join(ele_name.split('_')[1:])
			else:
				ele_category = 'variant'
			sub_dicti[ele_category].append(ele_iter)
		
		# Change to np arrays
		sub_dicti['variant'] = np.asarray(sub_dicti['variant'])
		for tissue_name in tissue_names:
			sub_dicti[tissue_name] = np.asarray(sub_dicti[tissue_name])
		# Add n_elements to sub_dicti
		sub_dicti['n_elements'] = len(ele_names)

		# add to global dictionary
		dicti[window_name] = sub_dicti

	f.close()
	return dicti

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




trait_name = sys.argv[1]
new_tgfm_stem = sys.argv[2]
tgfm_version = sys.argv[3]
processed_tgfm_input_stem = sys.argv[4]
gtex_pseudotissue_file = sys.argv[5]



#Extract tissue names
tissue_names = extract_pseudotissue_names(gtex_pseudotissue_file, ignore_testis=True)

###################################################
# Concatenate PIP summary file across parallel runs (one line for each window)
###################################################
suffix = 'tgfm_pip_summary.txt'
num_jobs=50
concatenated_pip_summary_file = new_tgfm_stem + '_iterative_prior_' + suffix
concatenate_results_across_parallel_jobs(new_tgfm_stem, suffix, num_jobs, concatenated_pip_summary_file)


###################################################
# Create component level summary data
###################################################
component_level_abf_summary_file = new_tgfm_stem + '_iterative_prior' + '_tgfm_component_level_abf_summary.txt'
per_window_abf_output_stem = new_tgfm_stem + '_iterative_prior_per_window_abf_'
generate_component_level_abf_summary_data(concatenated_pip_summary_file, component_level_abf_summary_file, tissue_names, new_tgfm_stem, tgfm_version, processed_tgfm_input_stem, per_window_abf_output_stem)

###################################################
# Learn iterative variant-gene-tissue prior
###################################################
variant_prob, tissue_probs = learn_iterative_variant_gene_tissue_prior(component_level_abf_summary_file, tgfm_version, tissue_names, per_window_abf_output_stem, max_iter=30)

# Print to output
variant_gene_prior_output_file = new_tgfm_stem + '_iterative_variant_gene_prior.txt'
t = open(variant_gene_prior_output_file,'w')
t.write('element_name\tprior\n')
t.write('variant\t' + str(variant_prob) + '\n')
for tiss_iter, tissue_name in enumerate(tissue_names):
	t.write(tissue_name + '\t' + str(tissue_probs[tiss_iter]) + '\n')
t.close()



