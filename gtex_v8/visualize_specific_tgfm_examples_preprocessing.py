import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import sys
import pdb
import pickle
import scipy.stats
import pyreadr


def extract_input_data_file_names_for_this_window(tgfm_input_summary_file, window_name):
	tgfm_input_file = 'none'
	trait_input_file = 'none'
	ld_input_file = 'none'
	f = open(tgfm_input_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		line_window_name = data[0]
		if line_window_name != window_name:
			continue
		tgfm_input_file = data[2]
		trait_input_file = data[3]
		ld_input_file = data[1]
	f.close()
	if tgfm_input_file == 'none':
		print('assumption error: input file doesnt exist for window')
		pdb.set_trace()

	return tgfm_input_file, trait_input_file, ld_input_file



def generate_snp_df_input_data(trait_name, window_name, tgfm_input_summary_file, tgfm_results_dir, snp_df_output_file, snp_id_to_rsid):
	# Get file names of input data corresponding to this specific window
	window_tgfm_input_data_file, window_trait_input_data_file, ld_input_file = extract_input_data_file_names_for_this_window(tgfm_input_summary_file, window_name)

	# Load in tgfm input data
	f = open(window_tgfm_input_data_file, "rb")
	window_tgfm_input_data = pickle.load(f)
	f.close()

	# Load in trait input data
	f = open(window_trait_input_data_file, "rb")
	window_trait_input_data = pickle.load(f)
	f.close()

	# Load in TGFM results data
	tgfm_res_file = tgfm_results_dir + 'tgfm_results_' + trait_name + '_component_gene_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_' + window_name + '_results.pkl'
	f = open(tgfm_res_file, "rb")
	tgfm_res = pickle.load(f)
	f.close()

	# Extract trait data
	trait_index = np.where(window_trait_input_data['gwas_study_names']==trait_name)[0][0]
	gwas_beta = window_trait_input_data['gwas_beta'][trait_index, :]
	gwas_beta_se = window_trait_input_data['gwas_beta_se'][trait_index, :]
	gwas_z = gwas_beta/gwas_beta_se
	gwas_p = scipy.stats.norm.sf(abs(gwas_z))*2
	neg_log_10_p = -np.log10(gwas_p)
	# Variant information
	variant_names = window_tgfm_input_data['variants']
	variant_positions = window_tgfm_input_data['varriant_positions']
	middle_variants = np.asarray([False]*len(variant_positions))
	middle_variants[window_tgfm_input_data['middle_variant_indices']] =True

	# Snp expected pips
	snp_expected_tgfm_pips = tgfm_res['expected_beta_pips']

	pdb.set_trace()

	# Open output file and print to outpu
	t = open(snp_df_output_file,'w')
	# Print header
	t.write('snp_name\tsnp_position\tgwas_z\tgwas_p\tgwas_neg_log10_p\tTGFM_PIP\tmiddle_variant_boolean\trs_id\n')
	# loop through snps
	n_snps = len(snp_expected_tgfm_pips)
	for snp_iter in range(n_snps):
		t.write(variant_names[snp_iter] + '\t' + str(variant_positions[snp_iter]) + '\t' + str(gwas_z[snp_iter]) + '\t' + str(gwas_p[snp_iter]) + '\t')
		t.write(str(neg_log_10_p[snp_iter]) + '\t' + str(snp_expected_tgfm_pips[snp_iter]) + '\t' + str(middle_variants[snp_iter]) + '\t' + snp_id_to_rsid[variant_names[snp_iter]] + '\n')
	t.close()

	return

def generate_gene_df_input_data(trait_name, window_name, tgfm_input_summary_file, tgfm_results_dir, gene_df_output_file, ensamble_id_to_gene_id):
	# Get file names of input data corresponding to this specific window
	window_tgfm_input_data_file, window_trait_input_data_file, ld_input_file = extract_input_data_file_names_for_this_window(tgfm_input_summary_file, window_name)

	# Load in tgfm input data
	f = open(window_tgfm_input_data_file, "rb")
	window_tgfm_input_data = pickle.load(f)
	f.close()

	# Load in trait input data
	f = open(window_trait_input_data_file, "rb")
	window_trait_input_data = pickle.load(f)
	f.close()

	# Load in TGFM results data
	tgfm_res_file = tgfm_results_dir + 'tgfm_results_' + trait_name + '_component_gene_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_' + window_name + '_results.pkl'
	f = open(tgfm_res_file, "rb")
	tgfm_res = pickle.load(f)
	f.close()


	# QUick error check
	if np.array_equal(window_tgfm_input_data['genes'], tgfm_res['genes']) == False:
		print('assumption eroror')
		pdb.set_trace()

	# Extract relevent fields
	gene_tissue_names = tgfm_res['genes']
	gene_tss = window_tgfm_input_data['tss']
	middle_genes = np.asarray([False]*len(gene_tss))
	middle_genes[window_tgfm_input_data['middle_gene_indices']] =True

	# gene expected pips
	gene_expected_tgfm_pips = tgfm_res['expected_alpha_pips']

	# Distributions of TWAS z-scors
	twas_z_mat = np.asarray(tgfm_res['nominal_twas_z'])
	neg_log10_p_mat = -np.log10(scipy.stats.norm.sf(abs(twas_z_mat))*2)

	n_genes = twas_z_mat.shape[1]
	mean_neg_log10_p = []
	neg_log_10_p_lb = []
	neg_log_10_p_ub = []
	for gene_iter in range(n_genes):
		mean_neg_log10_p.append(np.mean(neg_log10_p_mat[:, gene_iter]))
		sorted_p = np.sort(neg_log10_p_mat[:, gene_iter])
		neg_log_10_p_lb.append(sorted_p[4])
		neg_log_10_p_ub.append(sorted_p[95])
	neg_log_10_p_lb = np.asarray(neg_log_10_p_lb)
	neg_log_10_p_ub = np.asarray(neg_log_10_p_ub)
	mean_neg_log10_p = np.asarray(mean_neg_log10_p)

	# Print to output file
	t = open(gene_df_output_file,'w')
	# Print header
	t.write('gene_name\tgene_tss\tgwas_neg_log10_p_mean\tgwas_neg_log10_p_lb\tgwas_neg_log10_p_ub\tTGFM_PIP\tmiddle_gene_boolean\tgene_id_tissue_name\ttissue_name\tgene_id\n')	
	# Loop through genes
	for gene_iter in range(n_genes):
		gene_tissue_name = gene_tissue_names[gene_iter]
		ensamble_id = gene_tissue_name.split('_')[0]
		tissue_name = '_'.join(gene_tissue_name.split('_')[1:])
		gene_id = ensamble_id_to_gene_id[ensamble_id]
		gene_id_tissue_name = gene_id + ' (' + tissue_name + ')'

		t.write(gene_tissue_names[gene_iter] + '\t' + str(gene_tss[gene_iter]) + '\t' + str(mean_neg_log10_p[gene_iter]) + '\t')
		t.write(str(neg_log_10_p_lb[gene_iter]) + '\t' + str(neg_log_10_p_ub[gene_iter]) + '\t' + str(gene_expected_tgfm_pips[gene_iter]) + '\t' + str(middle_genes[gene_iter]) + '\t' + gene_id_tissue_name + '\t' + tissue_name + '\t' + gene_id + '\n')
	# Close file handle
	t.close()
	return


def generate_link_df_input_data(snp_df_output_file, gene_df_output_file, gtex_susie_gene_models_dir, link_df_output_file):
	# Create mapping from variant to (position, neg_log_10_p)
	variant_mapping = {}
	f = open(snp_df_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_name = data[0]
		variant_pos = data[1]
		var_neg_log_10_p = data[4]
		variant_mapping[variant_name] = (variant_pos, var_neg_log_10_p)
	f.close()

	t = open(link_df_output_file,'w')
	t.write('gene_name\tvariant_name\tgene_tss\tvariant_position\tgene_neg_log_10_p\tvariant_neg_log_10_p\n')
	f = open(gene_df_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_tissue_name = data[0]
		tissue_name = '_'.join(gene_tissue_name.split('_')[1:])
		gene_name = gene_tissue_name.split('_')[0]
		gene_tss = data[1]
		gene_neg_log10_p = data[2]

		# Gene weight file
		gene_model_file = gtex_susie_gene_models_dir + tissue_name + '/' + tissue_name + '_' + gene_name + '_1KG_only_fusion_output.wgt.RDat'
		# Load gene-tissue model from gene-tissue weight file
		gene_tissue_model = pyreadr.read_r(gene_model_file)

		# Loop through components
		gene_components = np.asarray(gene_tissue_model['susie_cs'])[:,0] - 1
		gene_snp_names = np.asarray(gene_tissue_model['variant_names'])[:,0]

		for component_iter in gene_components:
			best_snp_index = np.argmax(np.asarray(gene_tissue_model['susie_alpha'])[0,:])
			best_snp = gene_snp_names[best_snp_index]
			best_snp = best_snp.split('_b3')[0]
			snp_position = variant_mapping[best_snp][0]
			variant_neg_log_10_p = variant_mapping[best_snp][1]

			t.write(gene_tissue_name + '\t' + best_snp + '\t' + gene_tss + '\t' + snp_position + '\t' + gene_neg_log10_p + '\t' + variant_neg_log_10_p + '\n')
	f.close()


	return




def generate_visualization_input_data_for_specific_example(trait_name, window_name, tgfm_input_summary_file, tgfm_results_dir, tgfm_organized_results_dir, gtex_susie_gene_models_dir, ensamble_id_to_gene_id,snp_id_to_rsid, example_output_root):
	# Generate snp_df input data
	snp_df_output_file = example_output_root + '_snp_df.txt'
	generate_snp_df_input_data(trait_name, window_name, tgfm_input_summary_file, tgfm_results_dir, snp_df_output_file, snp_id_to_rsid)

	# Generate gene input data
	gene_df_output_file = example_output_root + '_gene_df.txt'
	generate_gene_df_input_data(trait_name, window_name, tgfm_input_summary_file, tgfm_results_dir, gene_df_output_file, ensamble_id_to_gene_id)


	# Generate link df input data
	link_df_output_file = example_output_root + '_link_df.txt'
	#generate_link_df_input_data(snp_df_output_file, gene_df_output_file, gtex_susie_gene_models_dir, link_df_output_file)




	return


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


def create_ensamble_id_to_gene_name_mapping(gene_annotation_file):
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
	return ensg_to_gene_name


def create_mapping_from_snp_id_to_rsid(sumstat_file):
	mapping = {}
	f = open(sumstat_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rs_id = data[0]
		snp_id = 'chr' + data[1] + '_' + data[2] + '_' + data[4] + '_' + data[5]
		snp_id2 = 'chr' + data[1] + '_' + data[2] + '_' + data[5] + '_' + data[4]
		mapping[snp_id] = rs_id
		mapping[snp_id2] = rs_id
	f.close()
	return mapping





#######################
# Command line args
#######################
specific_examples_input_file = sys.argv[1]
tgfm_input_summary_file = sys.argv[2]
tgfm_results_dir = sys.argv[3]
tgfm_organized_results_dir = sys.argv[4]
gtex_susie_gene_models_dir = sys.argv[5]
gene_annotation_file = sys.argv[6]
visualize_specific_tgfm_examples_dir = sys.argv[7]
ukbb_sumstats_dir = sys.argv[8]  # Used to create mapping from snp id to rs id


# Create mapping from snp_id to rsid using gwas summary stat file
snp_id_to_rsid = create_mapping_from_snp_id_to_rsid(ukbb_sumstats_dir + 'disease_THYROID_ANY_SELF_REP_hg38_liftover.bgen.stats')


# Create mapping from ensamble id to gene id
ensamble_id_to_gene_id = create_ensamble_id_to_gene_name_mapping(gene_annotation_file)


# Loop through examples to visualize
head_count = 0
f = open(specific_examples_input_file)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	# Skip header
	if head_count == 0:
		head_count = head_count + 1
		continue

	# Extract relevent fields
	trait_name = data[0]
	gene_tissue_name = data[1]
	gene_id = data[2]
	ensamble_id = data[3]
	tissue_name = data[4]
	window_name = data[5]
	if trait_name != 'disease_ALLERGY_ECZEMA_DIAGNOSED':
		continue

	# Output root
	example_output_root = visualize_specific_tgfm_examples_dir + trait_name + '_' + window_name

	# Generate input data for this window
	generate_visualization_input_data_for_specific_example(trait_name, window_name, tgfm_input_summary_file, tgfm_results_dir, tgfm_organized_results_dir, gtex_susie_gene_models_dir, ensamble_id_to_gene_id,snp_id_to_rsid, example_output_root)


f.close()