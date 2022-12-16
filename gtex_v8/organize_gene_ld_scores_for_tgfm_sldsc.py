import numpy as np 
import os
import sys
import pdb
import gzip




def get_pseudotissue_names(gtex_pseudotissue_file):
	f = open(gtex_pseudotissue_file)
	arr = []
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

def get_n_genes_per_chromosome_per_tissue(pseudotissue_names, gtex_susie_gene_models_dir):
	n_genes = np.zeros((22, len(pseudotissue_names)))
	for tissue_iter, pseudotissue_name in enumerate(pseudotissue_names):
		pos_file = gtex_susie_gene_models_dir + pseudotissue_name + '/' + pseudotissue_name + '_cis_heritable_gene_pos_file.txt'
		tmp = np.loadtxt(pos_file, dtype=str,delimiter='\t')
		for chrom_num in range(1,23):
			gene_count = np.sum(tmp[1:,2] == str(chrom_num))
			n_genes[(chrom_num-1), tissue_iter] = gene_count
	return n_genes

def load_in_gene_ld_scores(gene_ld_score_files):
	arr = []
	for gene_ld_score_file in gene_ld_score_files:
		gene_ld_scores = np.loadtxt(gene_ld_score_file)
		arr.append(gene_ld_scores)
	return np.transpose(np.asarray(arr))

def merge_variant_and_gene_ld_score_files(variant_ld_score_file, gene_ld_score_files, pseudotissue_names, merged_ld_score_file):
	gene_ld_scores = load_in_gene_ld_scores(gene_ld_score_files)
	gene_ld_scores[np.isnan(gene_ld_scores)] = 0.0

	f = gzip.open(variant_ld_score_file)
	t = open(merged_ld_score_file,'w')
	head_count = 0
	line_counter = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(np.asarray(data)) + '\t' + '\t'.join(pseudotissue_names) + '\n')
			continue
		t.write('\t'.join(np.asarray(data)) + '\t' + '\t'.join(gene_ld_scores[line_counter,:].astype(str)) + '\n')
		line_counter = line_counter + 1
	f.close()
	t.close()

	# Quick error checks before quiting
	if line_counter != gene_ld_scores.shape[0]:
		print('fundamental assumption eroror')
		pdb.set_trace()
	return 

def merge_variant_and_gene_ld_score_files_for_genotype_intercept(variant_ld_score_file, gene_ld_score_files, pseudotissue_names, merged_ld_score_file):
	gene_ld_scores = load_in_gene_ld_scores(gene_ld_score_files)
	gene_ld_scores[np.isnan(gene_ld_scores)] = 0.0

	f = gzip.open(variant_ld_score_file)
	t = open(merged_ld_score_file,'w')
	head_count = 0
	line_counter = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(np.asarray(data[:4])) + '\t' + '\t'.join(pseudotissue_names) + '\n')
			continue
		t.write('\t'.join(np.asarray(data[:4])) + '\t' + '\t'.join(gene_ld_scores[line_counter,:].astype(str)) + '\n')
		line_counter = line_counter + 1
	f.close()
	t.close()

	# Quick error checks before quiting
	if line_counter != gene_ld_scores.shape[0]:
		print('fundamental assumption eroror')
		pdb.set_trace()
	return 

def merge_m_files(variant_m_file, n_genes_per_pseudotissue, merged_m_file):
	t = open(merged_m_file,'w')
	m_vec = np.loadtxt(variant_m_file)
	t.write('\t'.join(m_vec.astype(str)) + '\t' + '\t'.join(n_genes_per_pseudotissue.astype(str)) + '\n')
	t.close()

def merge_m_files_for_genotype_intercept(variant_m_file, n_genes_per_pseudotissue, merged_m_file):
	t = open(merged_m_file,'w')
	m_vec = np.loadtxt(variant_m_file)
	t.write(str(m_vec[0]) + '\t' + '\t'.join(n_genes_per_pseudotissue.astype(str)) + '\n')
	t.close()

def make_ld_score_input_shell(chrom_num,  n_genes_per_pseudotissue, pseudotissue_names, variant_models, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir):
	for variant_model in variant_models:
		variant_ld_score_file = preprocessed_tgfm_sldsc_data_dir + variant_model + '.' + str(chrom_num) + '.l2.ldscore.gz'
		variant_m_file = preprocessed_tgfm_sldsc_data_dir + variant_model + '.' + str(chrom_num) + '.l2.M'
		variant_m_5_50_file = preprocessed_tgfm_sldsc_data_dir + variant_model + '.' + str(chrom_num) + '.l2.M_5_50'
		for gene_model_suffix in gene_model_suffixes:
			gene_ld_score_files = []
			for pseudotissue_name in pseudotissue_names:
				filer = preprocessed_tgfm_sldsc_data_dir + 'tissue_eqtl.' + str(chrom_num) + '.' + pseudotissue_name + '_' + gene_model_suffix
				gene_ld_score_files.append(filer)
			gene_ld_score_files = np.asarray(gene_ld_score_files)

			# Create output root
			output_root = preprocessed_tgfm_sldsc_data_dir + variant_model + '_' + gene_model_suffix + '.' + str(chrom_num) + '.l2'

			print(variant_model + '  ' + gene_model_suffix)
			# Merge ld scores files
			merged_ld_score_file = output_root + '.ldscore'
			merge_variant_and_gene_ld_score_files(variant_ld_score_file, gene_ld_score_files, pseudotissue_names, merged_ld_score_file)

			merged_m_file = output_root + '.M'
			merge_m_files(variant_m_file, n_genes_per_pseudotissue, merged_m_file)

			merged_m_5_50_file = output_root + '.M_5_50'
			merge_m_files(variant_m_5_50_file, n_genes_per_pseudotissue, merged_m_5_50_file)

def make_ld_score_input_shell_for_genotype_intercept(chrom_num, n_genes_per_pseudotissue, pseudotissue_names, reference_variant_model, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir):
	variant_ld_score_file = preprocessed_tgfm_sldsc_data_dir + reference_variant_model + '.' + str(chrom_num) + '.l2.ldscore.gz'
	variant_m_file = preprocessed_tgfm_sldsc_data_dir + reference_variant_model + '.' + str(chrom_num) + '.l2.M'
	variant_m_5_50_file = preprocessed_tgfm_sldsc_data_dir + reference_variant_model + '.' + str(chrom_num) + '.l2.M_5_50'
	for gene_model_suffix in gene_model_suffixes:
		gene_ld_score_files = []
		for pseudotissue_name in pseudotissue_names:
			filer = preprocessed_tgfm_sldsc_data_dir + 'tissue_eqtl.' + str(chrom_num) + '.' + pseudotissue_name + '_' + gene_model_suffix
			gene_ld_score_files.append(filer)
		gene_ld_score_files = np.asarray(gene_ld_score_files)

		# Create output root
		output_root = preprocessed_tgfm_sldsc_data_dir + 'genotype_intercept' + '_' + gene_model_suffix + '.' + str(chrom_num) + '.l2'

		print('genotype_intercept' + '  ' + gene_model_suffix)
		# Merge ld scores files
		merged_ld_score_file = output_root + '.ldscore'
		merge_variant_and_gene_ld_score_files_for_genotype_intercept(variant_ld_score_file, gene_ld_score_files, pseudotissue_names, merged_ld_score_file)

		merged_m_file = output_root + '.M'
		merge_m_files_for_genotype_intercept(variant_m_file, n_genes_per_pseudotissue, merged_m_file)

		merged_m_5_50_file = output_root + '.M_5_50'
		merge_m_files_for_genotype_intercept(variant_m_5_50_file, n_genes_per_pseudotissue, merged_m_5_50_file)


def create_annotation_sdev_file(ld_score_file_stem):
	output_file = ld_score_file_stem + '_annotation_sdev.txt'
	scores_arr = []
	for chrom_num in range(1,23):
		print(chrom_num)
		chrom_ldscore_file = ld_score_file_stem + '.' + str(chrom_num) + '.l2.ldscore'
		f = open(chrom_ldscore_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				if chrom_num == 1:
					anno_names = np.asarray(data[3:])
				continue
			scores_arr.append(np.asarray(data[3:]).astype(float))
		f.close()
	scores_arr = np.asarray(scores_arr)
	sdevs = np.std(scores_arr,axis=0)
	if len(sdevs) != len(anno_names):
		print('assumption eroror')
		pdb.set_trace()

	t = open(output_file,'w')
	t.write('annotation_name\tannotation_sdev\n')
	for ii in range(len(sdevs)):
		t.write(anno_names[ii] + '\t' + str(sdevs[ii]) + '\n')
	t.close()



def generate_annotation_sdev_files_shell(variant_models, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir):
	for variant_model in variant_models:
		for gene_model in gene_model_suffixes:
			ld_score_file_stem = preprocessed_tgfm_sldsc_data_dir + variant_model + '_' + gene_model #.11.l2.ldscore
			create_annotation_sdev_file(ld_score_file_stem)

def generate_annotation_sdev_file_for_genotype_intercept_shell(reference_variant_model, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir):
	for gene_model in gene_model_suffixes:
		ld_score_file_stem = preprocessed_tgfm_sldsc_data_dir + 'genotype_intercept' + '_' + gene_model
		create_annotation_sdev_file(ld_score_file_stem)


gtex_pseudotissue_file = sys.argv[1]
preprocessed_tgfm_sldsc_data_dir = sys.argv[2]
gene_type = sys.argv[3]
gtex_susie_gene_models_dir = sys.argv[4]


# Extract names of pseudotissues
pseudotissue_names = get_pseudotissue_names(gtex_pseudotissue_file)

# Get number of genes per chromosome per tissue
#n_genes_per_chromosome_per_tissue = get_n_genes_per_chromosome_per_tissue(pseudotissue_names, gtex_susie_gene_models_dir)

# Various iterations to run over
variant_models = ['baselineLD_no_qtl', 'baseline_no_qtl']

# Gene modedls
gene_model_suffixes = ['gene_ld_scores', 'gene_adj_ld_scores', 'pmces_gene_ld_scores', 'pmces_gene_adj_ld_scores']

'''
for chrom_num in range(1,23):
	print(chrom_num)
	make_ld_score_input_shell(chrom_num,  n_genes_per_chromosome_per_tissue[(chrom_num-1),:], pseudotissue_names, variant_models, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir)
'''

generate_annotation_sdev_files_shell(variant_models, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir)


# Make version without genomic annotaitons (call it genotype intercept)
reference_variant_model = 'baselineLD_no_qtl'
'''
for chrom_num in range(1,23):
	make_ld_score_input_shell_for_genotype_intercept(chrom_num,  n_genes_per_chromosome_per_tissue[(chrom_num-1),:], pseudotissue_names, reference_variant_model, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir)
'''
generate_annotation_sdev_file_for_genotype_intercept_shell(reference_variant_model, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir)

