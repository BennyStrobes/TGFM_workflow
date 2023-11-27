import numpy as np 
import os 
import sys
import pdb
import scipy.stats
import pickle

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

def create_file_containing_mapping_from_ensamble_id_to_gene_symbol_id(gene_annotation_file, ensamble_to_gene_symbol_mapping_file):
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
	f.close()
	t = open(ensamble_to_gene_symbol_mapping_file,'w')
	t.write('ensamble_id\tgene_symbol\n')
	ensamble_names = np.sort([*ensg_to_gene_name])
	for ensamble_id in ensamble_names:
		t.write(ensamble_id + '\t' + ensg_to_gene_name[ensamble_id] + '\n')
	t.close()
	return


def create_mapping_from_ensamble_id_to_gene_symbol_id(ensamble_to_gene_symbol_mapping_file):
	mapping = {}

	f = open(ensamble_to_gene_symbol_mapping_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		ensamble_id = data[0].split('.')[0]
		gene_symbol = data[1]
		mapping[ensamble_id] = gene_symbol

	return mapping


def create_file_containing_list_of_gene_models_tested_by_tgfm(gtex_susie_gene_models_dir, ensg_to_gene_symbol, tgfm_gene_list_file, tgfm_gene_tissue_list_file):
	t1 = open(tgfm_gene_list_file,'w')
	t2 = open(tgfm_gene_tissue_list_file,'w')
	t1.write('ensamble_id\tgene_symbol_id\n')
	t2.write('ensambled_id\tgene_symbol_id\ttissue_name\n')

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
			gene_symbol_id = ensg_to_gene_symbol[ensamble_id]
			t2.write(ensamble_id + '\t' + gene_symbol_id + '\t' + tissue_name + '\n')
			if ensamble_id not in used_genes:
				used_genes[ensamble_id] = 1
				t1.write(ensamble_id + '\t' + gene_symbol_id + '\n')
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
		arr.append(data[1])
		dicti[data[1]] = 1
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

def create_mapping_from_trait_name_to_drug_target_gene_list(drug_target_gene_list_file):
	mapping = {}
	f = open(drug_target_gene_list_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		study_name = data[0]
		if study_name.startswith('UKB_460K.') == False:
			continue
		short_study_name = study_name.split('460K.')[1]
		if data[7] == '':
			continue
		drug_target_gene_names = np.asarray(data[7].split('"')[1].split(','))
		#drug_target_gene_names = np.asarray(data[4].split('"')[1].split(','))
		mapping[short_study_name] = drug_target_gene_names
	f.close()
	return mapping

def arr_to_dicti(drug_target_gene_arr):
	dicti = {}
	for ele in drug_target_gene_arr:
		dicti[ele] = 1
	return dicti

def get_drug_target_genes_overlapping_tgfm_genes(drug_target_gene_arr, tgfm_tested_genes_dicti):
	dicti = {}
	for ele in drug_target_gene_arr:
		if ele in tgfm_tested_genes_dicti:
			dicti[ele] = 1
	return dicti

def get_tgfm_prioritized_genes_based_on_hard_pip_threshold(tgfm_results_dir, trait_name, gene_type,ensg_to_gene_symbol, pip_thresh=.25):
	tgfm_summary_file = tgfm_results_dir + 'tgfm_results_' + trait_name + '_' + gene_type + '_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_tgfm_per_gene_tissue_pip_summary.txt'
	f = open(tgfm_summary_file)
	dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_name = data[1]
		pip = float(data[-1])
		if pip > pip_thresh:
			dicti[ensg_to_gene_symbol[gene_name]] = 1.0
	return dicti

def get_tgfm_prioritized_genes_based_on_weighting(tgfm_results_dir, trait_name, gene_type, ensg_to_gene_symbol,lower_thresh=.001):
	tgfm_summary_file = tgfm_results_dir + 'tgfm_results_' + trait_name + '_' + gene_type + '_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_tgfm_per_gene_tissue_pip_summary.txt'
	f = open(tgfm_summary_file)
	dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[1]
		pip = float(data[-1])
		if lower_thresh > pip:
			continue
		gene_name = ensg_to_gene_symbol[ensamble_id]
		if gene_name not in dicti:
			dicti[gene_name] = 0.0
		dicti[gene_name] = dicti[gene_name] + pip
		if dicti[gene_name] > 1.0:
			dicti[gene_name] = 1.0
	return dicti	

def compute_contingency_table(drug_target_gene_dicti, tgfm_prioritized_genes, tgfm_tested_genes_dicti):
	all_genes = np.asarray([*tgfm_tested_genes_dicti])
	aa = 0.0
	bb = 0.0
	cc = 0.0
	dd = 0.0

	for gene_name in all_genes:
		if gene_name in tgfm_prioritized_genes and gene_name in drug_target_gene_dicti:
			aa = aa + 1.0
		elif gene_name in tgfm_prioritized_genes and gene_name not in drug_target_gene_dicti:
			bb = bb + 1.0
		elif gene_name not in tgfm_prioritized_genes and gene_name in drug_target_gene_dicti:
			cc = cc + 1.0
		else:
			dd = dd + 1.0
	return aa, bb, cc,dd

def compute_weighted_contingency_table(drug_target_gene_dicti, tgfm_prioritized_genes, tgfm_tested_genes_dicti):
	all_genes = np.asarray([*tgfm_tested_genes_dicti])
	aa = 0.0
	bb = 0.0
	cc = 0.0
	dd = 0.0

	for gene_name in all_genes:
		if gene_name in tgfm_prioritized_genes and gene_name in drug_target_gene_dicti:
			aa = aa + tgfm_prioritized_genes[gene_name]
			cc = cc + (1.0 - tgfm_prioritized_genes[gene_name])
		elif gene_name in tgfm_prioritized_genes and gene_name not in drug_target_gene_dicti:
			bb = bb + tgfm_prioritized_genes[gene_name]
			dd = dd + (1.0 - tgfm_prioritized_genes[gene_name])
		elif gene_name not in tgfm_prioritized_genes and gene_name in drug_target_gene_dicti:
			cc = cc + 1.0
		else:
			dd = dd + 1.0
	return aa, bb, cc,dd

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


###########################
# Command line args
###########################
tgfm_results_dir = sys.argv[1]
gene_type = sys.argv[2]
trait_name_file = sys.argv[3]
gene_annotation_file = sys.argv[4]
gtex_susie_gene_models_dir = sys.argv[5]
preprocessed_tgfm_data_dir = sys.argv[6]
drug_target_gene_list_file = sys.argv[7]
drug_target_gene_set_enrichment_dir = sys.argv[8]
tgfm_organized_results_dir = sys.argv[9]


# Create file containing mapping from ensamble id to gene symbol id
ensamble_to_gene_symbol_mapping_file = drug_target_gene_set_enrichment_dir + 'ensamble_to_gene_symbol_mapping.txt'
#create_file_containing_mapping_from_ensamble_id_to_gene_symbol_id(gene_annotation_file, ensamble_to_gene_symbol_mapping_file)

# Extract mapping from ensg_to_gene_symbold
ensg_to_gene_symbol = create_mapping_from_ensamble_id_to_gene_symbol_id(ensamble_to_gene_symbol_mapping_file)

# Create file containing list of gene models tested by TGFM
tgfm_gene_tissue_list_file = drug_target_gene_set_enrichment_dir + 'tgfm_gene_tissues_tested.txt'
tgfm_gene_list_file = drug_target_gene_set_enrichment_dir + 'tgfm_genes_tested.txt'
#create_file_containing_list_of_gene_models_tested_by_tgfm(gtex_susie_gene_models_dir, ensg_to_gene_symbol, tgfm_gene_list_file, tgfm_gene_tissue_list_file)

# Extract list of TGFM genes
# Note: array and dictionary are off by 1 genes (this is due to the fact there is not a 1 to 1 mapping from ensamble ids to gene symbol ids)
tgfm_tested_genes_arr, tgfm_tested_genes_dicti = extract_tgfm_tested_genes(tgfm_gene_list_file)

# Extract all tgfm trait names
trait_names = extract_trait_names(trait_name_file)

# Create mapping from trait_name to gold standard gene list (if mapping exists)
trait_to_drug_target_gene_list = create_mapping_from_trait_name_to_drug_target_gene_list(drug_target_gene_list_file)


independent_traits= np.asarray(["body_HEIGHTz", "blood_MEAN_PLATELET_VOL", "bmd_HEEL_TSCOREz", "blood_MEAN_CORPUSCULAR_HEMOGLOBIN", "blood_MONOCYTE_COUNT", "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT", "pigment_HAIR", "lung_FEV1FVCzSMOKE", "body_BALDING1", "biochemistry_Cholesterol", "bp_DIASTOLICadjMEDz", "lung_FVCzSMOKE", "repro_MENARCHE_AGE", "disease_ALLERGY_ECZEMA_DIAGNOSED", "other_MORNINGPERSON", "repro_NumberChildrenEverBorn_Pooled"])



# Open up output file handle
output_file = drug_target_gene_set_enrichment_dir + 'mendelian_gene_set_overlap_summary_file.txt'
t = open(output_file,'w')
t.write('trait_name\tgene_name\tensamble_id\ttgfm_gene_pip\ttgfm_abs_gene_pmces\tmax_abs_twas_z\tmendelian_gene\n')


# Loop through trait names
for trait_name in trait_names:
	# only consider traits for which we have drug targets for
	if trait_name not in trait_to_drug_target_gene_list:
		continue

	if trait_name not in independent_traits:
		continue

	print(trait_name)

	# Get TGFM gene names, pips, and twas z scores
	gene_to_tgfm_gene_pip, gene_to_tgfm_pmces, gene_to_max_abs_twas_z = get_tgfm_genes_and_pips_for_this_trait(trait_name, tgfm_results_dir, tgfm_organized_results_dir, preprocessed_tgfm_data_dir)

	# Get TGFM esnamle ids
	tgfm_ensamble_ids = np.unique([*gene_to_tgfm_gene_pip])

	# Drug target gene symbols
	trait_drug_target_genes = {}
	for gene_name in trait_to_drug_target_gene_list[trait_name]:
		trait_drug_target_genes[gene_name] = 1

	# Loop through all genes and print to output
	for ensamble_id in tgfm_ensamble_ids:
		gene_symbol_id = ensg_to_gene_symbol[ensamble_id]

		boolean = 0
		if gene_symbol_id in trait_drug_target_genes:
			boolean = 1

		t.write(trait_name + '\t' + gene_symbol_id + '\t' + ensamble_id + '\t' + str(gene_to_tgfm_gene_pip[ensamble_id]) + '\t' + str(np.abs(gene_to_tgfm_pmces[ensamble_id])) + '\t' + str(gene_to_max_abs_twas_z[ensamble_id]) + '\t' +  str(boolean) + '\n')
	t.flush()
t.close()









'''
# Get gene sets
drug_target_gene_arr = trait_to_drug_target_gene_list[trait_name]
drug_target_gene_dicti = get_drug_target_genes_overlapping_tgfm_genes(drug_target_gene_arr, tgfm_tested_genes_dicti)


# Get tgfm prioritized genes
tgfm_prioritized_genes_pip_thresh_1 = get_tgfm_prioritized_genes_based_on_hard_pip_threshold(tgfm_results_dir, trait_name, gene_type, ensg_to_gene_symbol, pip_thresh=.1)
	
# Compute odds ratio contingency table
aa, bb, cc, dd = compute_contingency_table(drug_target_gene_dicti, tgfm_prioritized_genes_pip_thresh_1, tgfm_tested_genes_dicti)
print(trait_name)
print(str(aa) + '\t' + str(bb) + '\t' + str(cc) + '\t' + str(dd))
#print((aa/bb)/(cc/dd))
print(scipy.stats.fisher_exact(np.asarray([[aa,bb],[cc,dd]])))
table_1 = table_1 + np.asarray([[aa,bb],[cc,dd]])

# Get tgfm prioritized genes
tgfm_prioritized_genes_pip_thresh_25 = get_tgfm_prioritized_genes_based_on_hard_pip_threshold(tgfm_results_dir, trait_name, gene_type, ensg_to_gene_symbol, pip_thresh=.25)
	
# Compute odds ratio contingency table
aa, bb, cc, dd = compute_contingency_table(drug_target_gene_dicti, tgfm_prioritized_genes_pip_thresh_25, tgfm_tested_genes_dicti)
print(trait_name)
print(str(aa) + '\t' + str(bb) + '\t' + str(cc) + '\t' + str(dd))
print(scipy.stats.fisher_exact(np.asarray([[aa,bb],[cc,dd]])))
table_25 = table_25 + np.asarray([[aa,bb],[cc,dd]])

# Get tgfm prioritized genes
tgfm_prioritized_genes_pip_thresh_5 = get_tgfm_prioritized_genes_based_on_hard_pip_threshold(tgfm_results_dir, trait_name, gene_type, ensg_to_gene_symbol, pip_thresh=.5)
	
# Compute odds ratio contingency table
aa, bb, cc, dd = compute_contingency_table(drug_target_gene_dicti, tgfm_prioritized_genes_pip_thresh_5, tgfm_tested_genes_dicti)
print(trait_name)	
print(str(aa) + '\t' + str(bb) + '\t' + str(cc) + '\t' + str(dd))
print(scipy.stats.fisher_exact(np.asarray([[aa,bb],[cc,dd]])))
table_5 = table_5 + np.asarray([[aa,bb],[cc,dd]])
'''
