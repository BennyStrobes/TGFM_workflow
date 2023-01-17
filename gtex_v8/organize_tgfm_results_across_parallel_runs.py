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

def concatenate_results_across_parallel_jobs(file_stem, num_jobs):
	concat_output_file = file_stem + '.txt'
	t = open(concat_output_file,'w')
	for job_number in range(num_jobs):
		job_input_file = file_stem + '_' + str(job_number) + '_' + str(num_jobs) + '.txt'
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


tgfm_results_dir = sys.argv[1]
gene_type = sys.argv[2]
num_jobs = int(sys.argv[3])
trait_names_file = sys.argv[4]



# Extract trait names
trait_names = extract_trait_names(trait_names_file)


# Concatenate parrallelized results across runs for each trait
for trait_name in trait_names:
	print(trait_name)
	# Concatenate cs-summary file
	file_stem = tgfm_results_dir + trait_name + '_tgfm_component_cs_summary'
	concatenate_results_across_parallel_jobs(file_stem, num_jobs)

	# Concatenate gene prob summary file
	file_stem = tgfm_results_dir + trait_name + '_tgfm_component_gene_prob_summary'
	concatenate_results_across_parallel_jobs(file_stem, num_jobs)

	# Concatenate tissue prob summary file
	file_stem = tgfm_results_dir + trait_name + '_tgfm_component_tissue_prob_summary'
	concatenate_results_across_parallel_jobs(file_stem, num_jobs)

	# get predicted causal effect sizes of fine-mapped snps
	cs_input_file = tgfm_results_dir + trait_name + '_tgfm_component_cs_summary.txt'
	pred_causal_effect_size_file = tgfm_results_dir + trait_name + '_tgfm_predicted_causal_effect_size.txt'
	create_predicted_causal_effect_size_file(cs_input_file, pred_causal_effect_size_file, trait_name, tgfm_results_dir)




