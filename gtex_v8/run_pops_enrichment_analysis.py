import numpy as np 
import os
import sys
import pdb
import gzip
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



def extract_tgfm_tested_genes(tgfm_gene_list_file):
	f = open(tgfm_gene_list_file)
	arr = []
	dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0].split('.')[0])
		dicti[data[0].split('.')[0]] = 1
	f.close()
	return np.asarray(arr), dicti



def extract_trait_names(trait_name_file):
	f = open(trait_name_file)
	trait_names = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		trait_names.append(data[0])
	f.close()
	return np.asarray(trait_names)


def get_unique_pops_traits(pops_results_summary_file, dicti):
	f = gzip.open(pops_results_summary_file)
	head_count = 0
	traits = {}
	for line in f:
		line = line.rstrip().decode('utf-8')
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0] in dicti:
			traits[data[0]] = 1
	f.close()

	return traits

def create_dictionary_mapping_pops_trait_names_to_out_trait_names():
	dicti = {}
	dicti['CAD'] = 'CAD'
	dicti['TC'] = 'biochemistry_Cholesterol'
	dicti['Glucose'] = 'biochemistry_Glucose'
	dicti['HDLC'] = 'biochemistry_HDLcholesterol'
	dicti['HbA1c'] = 'biochemistry_HbA1c'
	dicti['LDLC'] = 'biochemistry_LDLdirect'
	dicti['TG'] = 'biochemistry_Triglycerides'
	dicti['VitD'] = 'biochemistry_VitaminD'
	dicti['Eosino'] = 'blood_EOSINOPHIL_COUNT'
	dicti['Lym'] = 'blood_LYMPHOCYTE_COUNT'
	dicti['MCH'] = 'blood_MEAN_CORPUSCULAR_HEMOGLOBIN'
	dicti['Mono'] = 'blood_MONOCYTE_COUNT'
	dicti['Plt'] = 'blood_PLATELET_COUNT'
	dicti['RBC'] = 'blood_RED_COUNT'
	dicti['WBC'] = 'blood_WHITE_COUNT'
	dicti['eBMD'] = 'bmd_HEEL_TSCOREz'
	dicti['BMI'] = 'body_BMIz'
	dicti['Height'] = 'body_HEIGHTz'
	dicti['WHRadjBMI'] = 'body_WHRadjBMIz'
	dicti['DBP'] = 'bp_DIASTOLICadjMEDz'
	dicti['SBP'] = 'bp_SYSTOLICadjMEDz'
	dicti['AID_Combined'] = 'disease_AID_ALL'
	dicti['Asthma'] = 'disease_ASTHMA_DIAGNOSED'
	dicti['Hypothyroidism'] = 'disease_HYPOTHYROIDISM_SELF_REP'
	dicti['T2D'] = 'disease_T2D'
	dicti['FEV1FVC'] = 'lung_FEV1FVCzSMOKE'
	dicti['Morning_Person'] = 'other_MORNINGPERSON'
	dicti['Age_at_Menarche'] = 'repro_MENARCHE_AGE'
	return dicti

def create_mapping_from_trait_name_to_gene_pops_scores(pops_results_summary_file, pops_trait_names_to_trait_names):
	dicti = {}
	f = gzip.open(pops_results_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip().decode('utf-8')
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pops_trait_name = data[0]
		if pops_trait_name not in pops_trait_names_to_trait_names:
			continue

		new_trait_name = pops_trait_names_to_trait_names[pops_trait_name]

		ensamble_id = data[2]
		pops_score = float(data[-1])

		if new_trait_name not in dicti:
			dicti[new_trait_name] = ([], [])

		dicti[new_trait_name][0].append(ensamble_id)
		dicti[new_trait_name][1].append(pops_score)

	trait_names = [*dicti]

	return dicti

def get_tgfm_genes_and_pips_for_this_trait(trait_name, tgfm_results_dir, preprocessed_tgfm_data_dir):
	# Use this file just to get windows
	results_summary_file = tgfm_results_dir + 'tgfm_results_' + trait_name + '_component_gene_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_tgfm_pip_summary.txt'

	# Use dictionary to keep track of gene-tissue pairs and max pips
	tgfm_genetissue_to_pip = {}
	tgfm_gene_to_max_pip = {}
	tgfm_gene_to_sum_pip = {}
	tgfm_gene_to_max_abs_twas_z = {}
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
		tgfm_pmces_res_file = tgfm_results_dir + 'tgfm_results_' + trait_name + '_component_gene_susie_pmces_uniform_' + window_name + '_results.pkl'
		g = open(tgfm_pmces_res_file, 'rb')
		tgfm_pmces_res = pickle.load(g)
		g.close()

		# Load in TGFM input data file
		tgfm_input_file = preprocessed_tgfm_data_dir + 'component_gene_' + window_name + '_tgfm_trait_agnostic_input_data_obj.pkl'
		g = open(tgfm_input_file, 'rb')
		tgfm_input = pickle.load(g)
		g.close()

		# Extract relevent data
		middle_gene_indices = tgfm_input['middle_gene_indices']
		if len(middle_gene_indices) == 0:
			continue

		gene_pips = tgfm_res['expected_alpha_pips'][middle_gene_indices]
		gene_names = tgfm_res['genes'][middle_gene_indices]
		gene_abs_twas_z = np.abs(tgfm_pmces_res['nominal_twas_z'][middle_gene_indices])
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

			if ensamble_id not in tgfm_gene_to_max_abs_twas_z:
				tgfm_gene_to_max_abs_twas_z[ensamble_id] = gene_abs_twas_z[ii]
			else:
				tgfm_gene_to_max_abs_twas_z[ensamble_id] = np.max([gene_abs_twas_z[ii], tgfm_gene_to_max_abs_twas_z[ensamble_id]])

	f.close()
	return tgfm_genetissue_to_pip, tgfm_gene_to_max_pip, tgfm_gene_to_max_abs_twas_z, tgfm_gene_to_sum_pip



######################
# Command line args
######################
tgfm_results_dir = sys.argv[1]
gene_type = sys.argv[2]
trait_names_file = sys.argv[3]
gtex_susie_gene_models_dir = sys.argv[4]
preprocessed_tgfm_data_dir = sys.argv[5]
pops_results_summary_file = sys.argv[6]
pops_enrichment_dir = sys.argv[7]

# Create file containing list of gene models tested by TGFM
tgfm_gene_tissue_list_file = pops_enrichment_dir + 'tgfm_gene_tissues_tested.txt'
tgfm_gene_list_file = pops_enrichment_dir + 'tgfm_genes_tested.txt'
create_file_containing_list_of_gene_models_tested_by_tgfm(gtex_susie_gene_models_dir, tgfm_gene_list_file, tgfm_gene_tissue_list_file)

# Extract list of TGFM genes
tgfm_tested_genes_arr, tgfm_tested_genes_dicti = extract_tgfm_tested_genes(tgfm_gene_list_file)

# Create mapping from pops trait names to trait names
pops_trait_names_to_trait_names = create_dictionary_mapping_pops_trait_names_to_out_trait_names()

# Extract all tgfm trait names
trait_names = extract_trait_names(trait_names_file)


# Create mapping from trait names to pops scores for all genes
trait_name_to_gene_pops_scores = create_mapping_from_trait_name_to_gene_pops_scores(pops_results_summary_file, pops_trait_names_to_trait_names)

# valid trait names (ie those that are found in our analysis and pops analysis)
valid_trait_names = np.asarray([*trait_name_to_gene_pops_scores])

# Open output file handle
output_file = pops_enrichment_dir + 'cross_traits_pops_tgfm_enrichment_summary.txt'
t = open(output_file,'w')
t.write('trait_name\tgene_name\tmax_tgfm_pip\tsum_tgfm_pip\tmax_abs_twas_z\tpops_score\n')

# Loop through trait names
for trait_name in valid_trait_names:

	# Gene names tested by pops for this trait
	pops_trait_gene_names = np.asarray(trait_name_to_gene_pops_scores[trait_name][0])
	pops_trait_gene_scores = np.asarray(trait_name_to_gene_pops_scores[trait_name][1])


	# Get TGFM gene names and pips
	tgfm_genetissue_to_pip, tgfm_gene_to_max_pip, tgfm_gene_to_max_abs_twas_z, tgfm_gene_to_sum_pip = get_tgfm_genes_and_pips_for_this_trait(trait_name, tgfm_results_dir, preprocessed_tgfm_data_dir)

	# Get pips for Pops scores
	pops_trait_gene_pips = []
	pops_trait_gene_sum_pips = []
	pops_trait_gene_twas_z = []
	for gene_name in pops_trait_gene_names:
		if gene_name not in tgfm_gene_to_max_pip:
			pops_trait_gene_pips.append(np.nan)
		else:
			pops_trait_gene_pips.append(tgfm_gene_to_max_pip[gene_name])
		if gene_name not in tgfm_gene_to_sum_pip:
			pops_trait_gene_sum_pips.append(np.nan)
		else:
			pops_trait_gene_sum_pips.append(tgfm_gene_to_sum_pip[gene_name])
		if gene_name not in tgfm_gene_to_max_abs_twas_z:
			pops_trait_gene_twas_z.append(np.nan)
		else:
			pops_trait_gene_twas_z.append(tgfm_gene_to_max_abs_twas_z[gene_name])
	pops_trait_gene_pips = np.asarray(pops_trait_gene_pips)
	pops_trait_gene_twas_z = np.asarray(pops_trait_gene_twas_z)
	pops_trait_gene_sum_pips = np.asarray(pops_trait_gene_sum_pips)

	# Filter to genes that we have TGFM AND POPS
	valid_indices = np.isnan(pops_trait_gene_pips) == False
	pops_scores = pops_trait_gene_scores[valid_indices]
	tgfm_pips = pops_trait_gene_pips[valid_indices]
	tgfm_sum_pips = pops_trait_gene_sum_pips[valid_indices]
	twas_abs_z = pops_trait_gene_twas_z[valid_indices]
	final_genes = pops_trait_gene_names[valid_indices]

	# Print to outputs
	# Trait name
	print('#############################')
	print(trait_name)
	# Correlation coefficient
	print(scipy.stats.spearmanr(pops_scores, tgfm_pips))
	print(scipy.stats.spearmanr(pops_scores, twas_abs_z))


	print(np.mean(pops_scores[tgfm_pips < .1]))
	print(np.mean(pops_scores[tgfm_pips > .4]))

	for itera, gene_name in enumerate(final_genes):
		t.write(trait_name + '\t' + gene_name + '\t' + str(tgfm_pips[itera]) + '\t' + str(tgfm_sum_pips[itera]) + '\t' + str(twas_abs_z[itera]) + '\t' + str(pops_scores[itera]) + '\n')

t.close()


