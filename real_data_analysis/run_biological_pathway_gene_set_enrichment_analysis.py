import numpy as np 
import os
import sys
import pdb
import pickle
import scipy.stats



def create_file_containing_list_of_gene_models_tested_by_tgfm(gtex_susie_gene_models_dir, tgfm_gene_list_file, tgfm_gene_tissue_list_file):
	t1 = open(tgfm_gene_list_file,'w')
	t2 = open(tgfm_gene_tissue_list_file,'w')
	t1.write('ensamble_id\n')
	t2.write('ensambled_id\ttissue_name\n')

	used_genes = {}
	for tissue_name in os.listdir(gtex_susie_gene_models_dir):
		tissue_gene_list = gtex_susie_gene_models_dir + tissue_name + '/' + tissue_name + '_component_gene_pos_file.txt'
		head_count = 0
		f = open(tissue_gene_list)
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count = head_count + 1
				continue
			ensamble_id = data[1]
			t2.write(ensamble_id + '\t' + tissue_name + '\n')
			if ensamble_id not in used_genes:
				used_genes[ensamble_id] = 1
				t1.write(ensamble_id + '\n')
		f.close()
	t1.close()
	t2.close()
	return


def get_max_tgfm_sampler_z_v3(tgfm_res):
	alpha_zs = tgfm_res['alpha_mus']/np.sqrt(tgfm_res['alpha_vars'])
	combined_effects = np.zeros(alpha_zs[0].shape)
	for component_iter in range(len(alpha_zs)):
		component_e_alpha = tgfm_res['alpha_mus'][component_iter]*tgfm_res['alpha_phis'][component_iter]
		
		combined_effects =combined_effects + component_e_alpha

	return combined_effects

def get_gene_pip_dictionary(tgfm_results, middle_gene_indices):
	middle_gene_pips = tgfm_results['expected_alpha_pips'][middle_gene_indices]
	middle_gene_names = tgfm_results['genes'][middle_gene_indices]
		
	alphas = tgfm_results['alpha_phis']
	middle_gene_tissues = tgfm_results['genes'][middle_gene_indices]

	mapping = {}
	gene_to_indices = {}
	for ii, middle_gene_tissue in enumerate(middle_gene_tissues):
		gene_name = middle_gene_tissue.split('_')[0]
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
			tmper = alphas[ll][:, middle_gene_indices]
			new_alpha[:, ii] = np.sum(tmper[:, gene_to_indices[gene_name]],axis=1)
		gene_alphas.append(new_alpha)
	expected_alpha_pips = compute_expected_pips_from_sampler_pis(gene_alphas)
	for ii, gene_name in enumerate(ordered_genes):
		agg_pip = expected_alpha_pips[ii]
		mapping[gene_name.split('.')[0]] = agg_pip
	return mapping

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



def get_tgfm_genes_and_pips_for_this_trait(trait_name, tgfm_results_dir, tgfm_organized_results_dir, preprocessed_tgfm_data_dir):
	# Use this file just to get windows
	results_summary_file = tgfm_organized_results_dir + 'tgfm_results_' + trait_name + '_component_gene_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_tgfm_pip_summary.txt'

	# Use dictionary to keep track of gene-tissue pairs and max pips
	tgfm_genetissue_to_pip = {}
	tgfm_gene_to_max_pip = {}
	tgfm_gene_to_sum_pip = {}
	tgfm_gene_to_max_abs_twas_z = {}
	tgfm_gene_to_max_abs_tgfm_z = {}
	tgfm_gene_to_gene_pip = {}
	f = open(results_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]

		# Load in TGFM results dir
		tgfm_results_file = tgfm_results_dir + 'tgfm_results_' + trait_name + '_component_gene_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_' + window_name + '_results.pkl'
		g = open(tgfm_results_file, 'rb')
		tgfm_res = pickle.load(g)
		g.close()

		# TGFM PMCES (for twas z)
		'''
		tgfm_pmces_res_file = tgfm_results_dir + 'tgfm_results_' + trait_name + '_component_gene_susie_pmces_uniform_' + window_name + '_results.pkl'
		g = open(tgfm_pmces_res_file, 'rb')
		tgfm_pmces_res = pickle.load(g)
		g.close()
		'''

		# Load in TGFM input data file
		tgfm_input_file = preprocessed_tgfm_data_dir + 'component_gene_' + window_name + '_tgfm_trait_agnostic_input_data_obj.pkl'
		g = open(tgfm_input_file, 'rb')
		tgfm_input = pickle.load(g)
		g.close()

		# Extract relevent data
		middle_gene_indices = tgfm_input['middle_gene_indices']
		if len(middle_gene_indices) == 0:
			continue


		gene_max_tgfm_z = get_max_tgfm_sampler_z_v3(tgfm_res)[:, middle_gene_indices]


		gene_pips = tgfm_res['expected_alpha_pips'][middle_gene_indices]
		gene_names = tgfm_res['genes'][middle_gene_indices]
		gene_abs_twas_z = np.abs(np.median(np.asarray(tgfm_res['nominal_twas_z']),axis=0)[middle_gene_indices])
		tgfm_gene_to_gene_pip_tmp = get_gene_pip_dictionary(tgfm_res, middle_gene_indices)
		for gene_id in [*tgfm_gene_to_gene_pip_tmp]:
			tgfm_gene_to_gene_pip[gene_id] = tgfm_gene_to_gene_pip_tmp[gene_id]
		#gene_abs_twas_z = np.abs(tgfm_pmces_res['nominal_twas_z'][middle_gene_indices])
		#gene_twas_z = tgfm_res['nominal_twas_z'][middle_gene_indices]
		for ii,gene_name in enumerate(gene_names):
			if gene_name in tgfm_genetissue_to_pip:
				print('assumption eroror')
				pdb.set_trace()
			tgfm_genetissue_to_pip[gene_name] = gene_pips[ii]

			ensamble_id = gene_name.split('.')[0]
			if ensamble_id not in tgfm_gene_to_max_pip:
				tgfm_gene_to_max_pip[ensamble_id] = gene_pips[ii]
			else:
				tgfm_gene_to_max_pip[ensamble_id] = np.max([gene_pips[ii], tgfm_gene_to_max_pip[ensamble_id]])
			if ensamble_id not in tgfm_gene_to_sum_pip:
				tgfm_gene_to_sum_pip[ensamble_id] = gene_pips[ii]
			else:
				tgfm_gene_to_sum_pip[ensamble_id] = np.sum([gene_pips[ii], tgfm_gene_to_sum_pip[ensamble_id]])
			if ensamble_id not in tgfm_gene_to_max_abs_tgfm_z:
				tgfm_gene_to_max_abs_tgfm_z[ensamble_id] = gene_max_tgfm_z[:,ii]
			else:
				tgfm_gene_to_max_abs_tgfm_z[ensamble_id] = gene_max_tgfm_z[:,ii] + tgfm_gene_to_max_abs_tgfm_z[ensamble_id]
				#tgfm_gene_to_max_abs_tgfm_z[ensamble_id] = np.sum([gene_max_tgfm_z[:,ii], tgfm_gene_to_max_abs_tgfm_z[ensamble_id]])
				#tgfm_gene_to_max_abs_tgfm_z[ensamble_id] = np.max([np.abs(gene_max_tgfm_z[ii]), np.abs(tgfm_gene_to_max_abs_tgfm_z[ensamble_id])])
			if ensamble_id not in tgfm_gene_to_max_abs_twas_z:
				tgfm_gene_to_max_abs_twas_z[ensamble_id] = gene_abs_twas_z[ii]
			else:
				tgfm_gene_to_max_abs_twas_z[ensamble_id] = np.max([gene_abs_twas_z[ii], tgfm_gene_to_max_abs_twas_z[ensamble_id]])
	f.close()

	tgfm_gene_to_max_abs_tgfm_z2 = {}
	for gene in [*tgfm_gene_to_max_abs_tgfm_z]:
		tgfm_gene_to_max_abs_tgfm_z2[gene] = np.abs(np.mean(tgfm_gene_to_max_abs_tgfm_z[gene]))

	return tgfm_gene_to_gene_pip, tgfm_gene_to_max_abs_tgfm_z2, tgfm_gene_to_max_abs_twas_z

def compute_odds_ratio_quantities_for_this_pathways(gene_to_tgfm_gene_pip, pathway_ensamble_ids, pip_threshold=.5):
	# Pathway ensamble id dicti
	pathway_ensamble_id_dicti = {}
	for pathway_ensamble_id in pathway_ensamble_ids:
		pathway_ensamble_id_dicti[pathway_ensamble_id.split('.')[0]] = 1

	# Compute odds ratio quantities
	aa = 0  # high tgfm pip + in pathway
	bb = 0  # low tgfm pip + in pathway
	cc = 0  # high tgfm pip + not in pathway
	dd = 0  # Low tgfm pip but not in pathway
	ensamble_ids = np.asarray([*gene_to_tgfm_gene_pip])
	for ensamble_id in ensamble_ids:
		tgfm_pip = gene_to_tgfm_gene_pip[ensamble_id]
		if tgfm_pip >= pip_threshold and ensamble_id in pathway_ensamble_id_dicti:
			aa = aa + 1
		elif tgfm_pip < pip_threshold and ensamble_id in pathway_ensamble_id_dicti:
			bb = bb + 1
		elif tgfm_pip >= pip_threshold and ensamble_id not in pathway_ensamble_id_dicti:
			cc = cc + 1
		elif tgfm_pip < pip_threshold and ensamble_id not in pathway_ensamble_id_dicti:
			dd = dd + 1
	return aa, bb, cc, dd

def compute_probabilistic_odds_ratio_quantities_for_this_pathways(gene_to_tgfm_gene_pip, pathway_ensamble_ids, pip_threshold=.5):
	# Pathway ensamble id dicti
	pathway_ensamble_id_dicti = {}
	for pathway_ensamble_id in pathway_ensamble_ids:
		pathway_ensamble_id_dicti[pathway_ensamble_id.split('.')[0]] = 1

	# Compute odds ratio quantities
	aa = 0  # high tgfm pip + in pathway
	bb = 0  # low tgfm pip + in pathway
	cc = 0  # high tgfm pip + not in pathway
	dd = 0  # Low tgfm pip but not in pathway
	ensamble_ids = np.asarray([*gene_to_tgfm_gene_pip])
	for ensamble_id in ensamble_ids:
		tgfm_pip = gene_to_tgfm_gene_pip[ensamble_id]
		if tgfm_pip < pip_threshold:
			tgfm_pip = 0.0
		if ensamble_id in pathway_ensamble_id_dicti:
			aa = aa + tgfm_pip
			bb = bb + (1.0 - tgfm_pip)
		elif ensamble_id not in pathway_ensamble_id_dicti:
			cc = cc + tgfm_pip
			dd = dd + (1.0 - tgfm_pip)
	return aa, bb, cc, dd



###########################
# Command line args
###########################
tgfm_results_dir = sys.argv[1]
gene_type = sys.argv[2]
trait_name_file = sys.argv[3]
gene_annotation_file = sys.argv[4]
gtex_susie_gene_models_dir = sys.argv[5]
preprocessed_tgfm_data_dir = sys.argv[6]
biological_pathway_gene_set_file = sys.argv[7]
biological_pathway_enrichment_dir = sys.argv[8]
tgfm_organized_results_dir = sys.argv[9]


# Create file containing list of gene models tested by TGFM
tgfm_gene_tissue_list_file = biological_pathway_enrichment_dir + 'tgfm_gene_tissues_tested.txt'
tgfm_gene_list_file = biological_pathway_enrichment_dir + 'tgfm_genes_tested.txt'
#create_file_containing_list_of_gene_models_tested_by_tgfm(gtex_susie_gene_models_dir, tgfm_gene_list_file, tgfm_gene_tissue_list_file)


independent_traits= np.asarray(["body_HEIGHTz", "blood_MEAN_PLATELET_VOL", "bmd_HEEL_TSCOREz", "blood_MEAN_CORPUSCULAR_HEMOGLOBIN", "blood_MONOCYTE_COUNT", "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT", "pigment_HAIR", "lung_FEV1FVCzSMOKE", "body_BALDING1", "biochemistry_Cholesterol", "bp_DIASTOLICadjMEDz", "lung_FVCzSMOKE", "repro_MENARCHE_AGE", "disease_ALLERGY_ECZEMA_DIAGNOSED", "other_MORNINGPERSON", "repro_NumberChildrenEverBorn_Pooled"])

pip_thresholds = [0.0, .1, .25, .5]
for trait_name in independent_traits:
	for pip_threshold in pip_thresholds:
		print(trait_name)
		print(pip_threshold)
		# Get TGFM gene names, pips, and twas z scores
		gene_to_tgfm_gene_pip, gene_to_tgfm_pmces, gene_to_max_abs_twas_z = get_tgfm_genes_and_pips_for_this_trait(trait_name, tgfm_results_dir, tgfm_organized_results_dir, preprocessed_tgfm_data_dir)

		# Get TGFM esnamle ids
		tgfm_ensamble_ids = np.unique([*gene_to_tgfm_gene_pip])


		# Open output file handle for this trait
		t = open(biological_pathway_enrichment_dir + 'biologial_pathway_enrichment_analysis_' + trait_name + '_probabability_weighted_' + str(pip_threshold) + '.txt','w')
		t.write('pathway_index\tpathway_name\taa\tbb\tcc\tdd\torat\torat_pvalue\n')
		# Loop through pathways and do enrichment for each
		f = open(biological_pathway_gene_set_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			# Extract relevent fields
			pathway_index = data[0]
			pathway_name = data[1]
			pathway_ensamble_ids = np.asarray(data[5].split(','))
			aa, bb, cc, dd = compute_probabilistic_odds_ratio_quantities_for_this_pathways(gene_to_tgfm_gene_pip, pathway_ensamble_ids,pip_threshold=pip_threshold)
			fe_orat, fe_pvalue = scipy.stats.fisher_exact(np.asarray([[aa,bb],[cc,dd]]))		
			t.write(pathway_index + '\t' + pathway_name + '\t' + str(aa) + '\t' + str(bb) + '\t' + str(cc) + '\t' + str(dd) + '\t' + str(fe_orat) + '\t' + str(fe_pvalue) + '\n')
		f.close()
		t.close()


'''
pip_thresholds = [.25, .5, .7, .9]
for trait_name in independent_traits:
	for pip_threshold in pip_thresholds:
		print(trait_name)
		print(pip_threshold)
		# Get TGFM gene names, pips, and twas z scores
		gene_to_tgfm_gene_pip, gene_to_tgfm_pmces, gene_to_max_abs_twas_z = get_tgfm_genes_and_pips_for_this_trait(trait_name, tgfm_results_dir, tgfm_organized_results_dir, preprocessed_tgfm_data_dir)

		# Get TGFM esnamle ids
		tgfm_ensamble_ids = np.unique([*gene_to_tgfm_gene_pip])


		# Open output file handle for this trait
		t = open(biological_pathway_enrichment_dir + 'biologial_pathway_enrichment_analysis_' + trait_name + '_' + str(pip_threshold) + '.txt','w')
		t.write('pathway_index\tpathway_name\taa\tbb\tcc\tdd\torat\torat_pvalue\n')
		# Loop through pathways and do enrichment for each
		f = open(biological_pathway_gene_set_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			# Extract relevent fields
			pathway_index = data[0]
			pathway_name = data[1]
			pathway_ensamble_ids = np.asarray(data[5].split(','))
			aa, bb, cc, dd = compute_odds_ratio_quantities_for_this_pathways(gene_to_tgfm_gene_pip, pathway_ensamble_ids,pip_threshold=pip_threshold)

			fe_orat, fe_pvalue = scipy.stats.fisher_exact(np.asarray([[aa,bb],[cc,dd]]))
		
			t.write(pathway_index + '\t' + pathway_name + '\t' + str(aa) + '\t' + str(bb) + '\t' + str(cc) + '\t' + str(dd) + '\t' + str(fe_orat) + '\t' + str(fe_pvalue) + '\n')
		f.close()
		t.close()
'''


