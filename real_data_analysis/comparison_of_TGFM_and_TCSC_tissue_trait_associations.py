import numpy as np
import os
import sys
import pdb




def get_list_of_trait_names(trait_names_file):
	f = open(trait_names_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0])
	f.close()
	return np.asarray(arr)



def extract_tgfm_trait_tissue_results(tgfm_results_dir, trait_names):
	dicti = {}
	all_tissue_names = []
	for trait_name in trait_names:
		# TGFM tissue specific prior summary
		tgfm_tissue_spec_prior_summary_file = tgfm_results_dir + 'tgfm_results_' + trait_name + '_component_gene_susie_pmces_uniform_iterative_variant_gene_prior_v2_pip_level_bootstrapped.txt'
		
		# Stream summary file
		f = open(tgfm_tissue_spec_prior_summary_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			ele_class = data[0]
			if ele_class == 'variant':
				continue

			tissue_name = data[0]
			avg_prob = float(data[1])
			prob_distr = np.asarray(data[3].split(';')).astype(float)

			mean_prob = np.mean(prob_distr)

			sorted_prob_distr = np.sort(prob_distr)
			prob_lb = (sorted_prob_distr[1] + sorted_prob_distr[2])/2.0
			prob_ub = (sorted_prob_distr[-2] + sorted_prob_distr[-3])/2.0

			dicti[trait_name + ':' + tissue_name] = (mean_prob, prob_lb, prob_ub)

			all_tissue_names.append(tissue_name)

		f.close()

	all_tissue_names = np.asarray(all_tissue_names)
	all_tissue_names = np.sort(np.unique(all_tissue_names))

	return dicti, all_tissue_names


def extract_tcsc_trait_tissue_results(tcsc_results_file, trait_names):
	valid_traits = {}
	for trait_name in trait_names:
		valid_traits['UKB_460K.' + trait_name] = trait_name
	dicti = {}
	all_tissue_names = []
	head_count = 0
	f = open(tcsc_results_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tcsc_trait_name = data[0]
		if tcsc_trait_name not in valid_traits:
			continue
		tgfm_trait_name = valid_traits[tcsc_trait_name]
		tissue = data[2]
		all_tissue_names.append(tissue)

		tcsc_h2 = float(data[3])
		tcsc_h2_se = float(data[4])

		tcsc_h2_lb = tcsc_h2 - 1.96*tcsc_h2_se
		tcsc_h2_ub = tcsc_h2 + 1.96*tcsc_h2_se

		dicti[tgfm_trait_name + ':' + tissue] = (tcsc_h2, tcsc_h2_lb, tcsc_h2_ub)

	f.close()

	all_tissue_names = np.asarray(all_tissue_names)
	all_tissue_names = np.sort(np.unique(all_tissue_names))

	return dicti, all_tissue_names




#######################
# Command line args
#######################
tgfm_results_dir = sys.argv[1]
trait_names_file = sys.argv[2]
tcsc_results_file = sys.argv[3]
tcsc_comparison_dir = sys.argv[4]


###############
# Get list of TGFM analyzed trait names
trait_names = get_list_of_trait_names(trait_names_file)

##############
# create mapping from trait_tissue_name to tgfm elements
tgfm_trait_tissue_results, tgfm_tissue_names = extract_tgfm_trait_tissue_results(tgfm_results_dir, trait_names)


# create mapping from trait_tissue_name to tcsc elements
tcsc_trait_tissue_results, tcsc_tissue_names = extract_tcsc_trait_tissue_results(tcsc_results_file, trait_names)



for trait_name in trait_names:
	if trait_name + ':' + 'Muscle_Skeletal' not in tcsc_trait_tissue_results:
		continue
	tcsc_vec = []
	tgfm_vec = []

	for tissue_name in tgfm_tissue_names:
		trait_tissue_pair = trait_name + ':' + tissue_name
		tcsc_vec.append(tcsc_trait_tissue_results[trait_tissue_pair][0])
		tgfm_vec.append(tgfm_trait_tissue_results[trait_tissue_pair][0])

		if tcsc_trait_tissue_results[trait_tissue_pair][1] > 0.0:
			print(trait_tissue_pair)
			print(tgfm_trait_tissue_results[trait_tissue_pair])

	tcsc_vec = np.asarray(tcsc_vec)
	tgfm_vec = np.asarray(tgfm_vec)

	print(trait_name)
	#print(np.corrcoef(tcsc_vec,tgfm_vec)[0,1])
	#if np.corrcoef(tcsc_vec,tgfm_vec)[0,1] < .1:
		#pdb.set_trace()








