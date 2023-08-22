import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import pdb
import pickle
import scipy.stats
import statsmodels.api as sm




def create_file_containing_list_of_gene_models_tested_by_tgfm(gtex_susie_gene_models_dir, tgfm_gene_list_file, tgfm_gene_tissue_list_file):
	t1 = open(tgfm_gene_list_file,'w')
	t2 = open(tgfm_gene_tissue_list_file,'w')
	t1.write('ensamble_id\n')
	t2.write('ensambled_id\ttissue_name\n')

	used_genes = {}
	for tissue_name in os.listdir(gtex_susie_gene_models_dir):
		if tissue_name.endswith('.py'):
			continue
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


def extract_gene_set_info(non_disease_specific_gene_sets_file):
	f = open(non_disease_specific_gene_sets_file)
	head_count = 0
	universe_genes = {}
	for line in f:
		line = line.rstrip()
		data = line.split(',')
		if head_count == 0:
			head_count = head_count + 1
			total_columns = len(data)
			gene_set_names = np.asarray(data[1:60])
			gene_set_name_to_gene_dictionary = {}
			for gene_set_name in gene_set_names:
				gene_set_name_to_gene_dictionary[gene_set_name] = {}
			continue
		if len(data) != total_columns:
			print('assumption erroror')
			pdb.set_trace()
		# Add universe gene
		universe_genes[data[0]] = 1
		for gene_index, gene_name in enumerate(data[1:60]):
			gene_set_name = gene_set_names[gene_index]
			if gene_name == '':
				continue
			gene_set_name_to_gene_dictionary[gene_set_name][gene_name] = 1
	f.close()
	return universe_genes, gene_set_names, gene_set_name_to_gene_dictionary 


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

		tgfm_gene_to_gene_pip_tmp = get_gene_pip_dictionary(tgfm_res, middle_gene_indices)
		for gene_id in [*tgfm_gene_to_gene_pip_tmp]:
			tgfm_gene_to_gene_pip[gene_id] = tgfm_gene_to_gene_pip_tmp[gene_id]

	return tgfm_gene_to_gene_pip

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


tgfm_results_dir=sys.argv[1]
traits_file=sys.argv[2]
gtex_susie_gene_models_dir=sys.argv[3]
preprocessed_tgfm_data_dir=sys.argv[4]
tgfm_organized_results_dir=sys.argv[5]
non_disease_specific_gene_sets_file=sys.argv[6]
non_disease_specific_gene_set_enrichment_dir=sys.argv[7]

print(non_disease_specific_gene_set_enrichment_dir)


# Create file containing list of gene models tested by TGFM
tgfm_gene_tissue_list_file = non_disease_specific_gene_set_enrichment_dir + 'tgfm_gene_tissues_tested.txt'
tgfm_gene_list_file = non_disease_specific_gene_set_enrichment_dir + 'tgfm_genes_tested.txt'
create_file_containing_list_of_gene_models_tested_by_tgfm(gtex_susie_gene_models_dir, tgfm_gene_list_file, tgfm_gene_tissue_list_file)

# Extract list of TGFM genes
tgfm_tested_genes_arr, tgfm_tested_genes_dicti = extract_tgfm_tested_genes(tgfm_gene_list_file)


independent_traits= np.asarray(["body_HEIGHTz", "blood_MEAN_PLATELET_VOL", "bmd_HEEL_TSCOREz", "blood_MEAN_CORPUSCULAR_HEMOGLOBIN", "blood_MONOCYTE_COUNT", "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT", "pigment_HAIR", "lung_FEV1FVCzSMOKE", "body_BALDING1", "biochemistry_Cholesterol", "bp_DIASTOLICadjMEDz", "lung_FVCzSMOKE", "repro_MENARCHE_AGE", "disease_ALLERGY_ECZEMA_DIAGNOSED", "other_MORNINGPERSON", "repro_NumberChildrenEverBorn_Pooled"])

gene_set_universe, gene_set_names, gene_set_name_to_gene_dictionary = extract_gene_set_info(non_disease_specific_gene_sets_file)

# oepn output file handles
'''
tt = {}
for gene_set_name in gene_set_names:
	output_file = non_disease_specific_gene_set_enrichment_dir + gene_set_name + '_enrichment_summary.txt'
	tt[gene_set_name] = open(output_file,'w')
	tt[gene_set_name].write('trait_name\tgene_name\ttgfm_gene_pip\tgene_set_boolean\n')



# Loop through trait names
for trait_name in independent_traits:
	print(trait_name)

	# Get TGFM gene names and pips
	gene_to_tgfm_gene_pip = get_tgfm_genes_and_pips_for_this_trait(trait_name, tgfm_results_dir, tgfm_organized_results_dir, preprocessed_tgfm_data_dir)
	# Get tgfm genes
	tgfm_genes = np.asarray([*gene_to_tgfm_gene_pip])

	for tgfm_gene in tgfm_tested_genes_arr:
		if tgfm_gene not in gene_set_universe:
			continue
		if tgfm_gene not in gene_to_tgfm_gene_pip:
			tgfm_gene_pip = 0.0
		else:
			tgfm_gene_pip = gene_to_tgfm_gene_pip[tgfm_gene]
		for gene_set_name in gene_set_names:
			if tgfm_gene in gene_set_name_to_gene_dictionary[gene_set_name]:
				gene_set_boolean='1'
			else:
				gene_set_boolean='0'

			tt[gene_set_name].write(trait_name + '\t' + tgfm_gene + '\t' + str(tgfm_gene_pip) + '\t' + gene_set_boolean + '\n')


for gene_set_name in gene_set_names:
	output_file = non_disease_specific_gene_set_enrichment_dir + gene_set_name + '_enrichment_summary.txt'
	tt[gene_set_name].close()
'''

valid_gene_sets = {}
valid_gene_sets['MGI_essential'] = 'MGI essential'
valid_gene_sets['Haploinsufficient'] = 'Haploinsufficient'
valid_gene_sets['highPLI_Exac'] = 'High pLI'
valid_gene_sets['Olfactory'] = 'Olfactory receptors'
valid_gene_sets['high_EDS'] = 'High Enhancer Domain Score'
valid_gene_sets['highShet_Cassa'] = 'High Shet'
valid_gene_sets['Hart2017_nonessential'] = 'Non-essential'
valid_gene_sets['Genes_withSNPs_allSNPs_100kb'] = 'Most SNPs in 100kb'

output_file = non_disease_specific_gene_set_enrichment_dir +'global_enrichment_summary.txt'
print(output_file)
t = open(output_file,'w')
t.write('gene_set_name\todds_ratio\todds_ratio_lb\todds_ratio_ub\todds_ratio_pvalue\n')

for gene_set_name in gene_set_names:
	if gene_set_name not in valid_gene_sets:
		continue
	readable_gene_set_name = valid_gene_sets[gene_set_name]
	print(gene_set_name)
	output_file = non_disease_specific_gene_set_enrichment_dir + gene_set_name + '_enrichment_summary.txt'
	aa = np.loadtxt(output_file,dtype=str)
	pips = aa[1:,-2].astype(float)
	labels = aa[1:,-1].astype(float)
	denom = np.sum(labels)/len(labels)
	for pip in [.5]:
		'''
		indices = pips > pip
		numer = np.sum(labels[indices])/np.sum(indices)
		enrichment = numer/denom
		print(str(pip) + ':  ' + str(enrichment))
		'''
		hit_indices = pips > pip
		null_indices = pips <= pip
		aa = np.sum(labels[hit_indices])
		bb = np.sum(labels[null_indices])
		cc = np.sum((1.0-labels)[hit_indices])
		dd = np.sum((1.0-labels)[null_indices])
		res = scipy.stats.fisher_exact([[aa,bb],[cc,dd]])


		pips2 = np.copy(hit_indices*1.0)
		X2 = sm.add_constant(pips2)
		logit_res = sm.Logit(labels, X2*1.0).fit()

		pvalue = logit_res.pvalues[1]
		param = logit_res.params[1]
		param_se = logit_res.bse[1]

		param_ub = param + 1.96*param_se
		param_lb = param - 1.96*param_se

		odds_ratio = np.exp(param)
		odds_ratio_lb = np.exp(param_lb)
		odds_ratio_ub = np.exp(param_ub)
		t.write(readable_gene_set_name + '\t' + str(odds_ratio) + '\t' + str(odds_ratio_lb) + '\t' + str(odds_ratio_ub) + '\t' + str(pvalue) + '\n')
t.close()


