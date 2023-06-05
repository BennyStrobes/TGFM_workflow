import numpy as np 
import os
import sys
import pdb
import gzip




def get_pseudotissue_names(gtex_pseudotissue_file, tissue_version):
	sex_tissues_dicti = {}
	sex_tissues_dicti['Ovary'] = 1
	sex_tissues_dicti['Prostate'] = 1
	sex_tissues_dicti['Testis'] = 1
	sex_tissues_dicti['Uterus'] = 1
	sex_tissues_dicti['Vagina'] = 1
	f = open(gtex_pseudotissue_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if tissue_version == 'non_sex_tissues' and data[0] in sex_tissues_dicti:
			continue
		if tissue_version == 'no_testis' and data[0] == 'Testis':
			continue
		arr.append(data[0])
	f.close()

	return np.asarray(arr)

def get_n_genes_per_chromosome_per_tissue(pseudotissue_names, gtex_susie_gene_models_dir, gene_type):
	n_genes = np.zeros((22, len(pseudotissue_names)))
	for tissue_iter, pseudotissue_name in enumerate(pseudotissue_names):
		pos_file = gtex_susie_gene_models_dir + pseudotissue_name + '/' + pseudotissue_name + '_' + gene_type + '_pos_file.txt'
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

def print_weight_output_file(variant_ld_score_file, gene_weight_files, pseudotissue_names, weight_output_file):
	gene_weights = load_in_gene_ld_scores(gene_weight_files)
	gene_weights[np.isnan(gene_weights)] = 0.0
	agg_gene_weights = np.sum(gene_weights,axis=1)

	f = gzip.open(variant_ld_score_file)
	t = open(weight_output_file,'w')
	head_count = 0
	line_counter = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\tgene_weight\n')
			continue
		t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + str(agg_gene_weights[line_counter]) + '\n')
		line_counter = line_counter + 1
	f.close()
	t.close()

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

def make_gene_weight_input_file(chrom_num, pseudotissue_names, variant_model, gene_model_suffix, preprocessed_tgfm_sldsc_data_dir, gene_type, tissue_version):
	variant_ld_score_file = preprocessed_tgfm_sldsc_data_dir + variant_model + '.' + str(chrom_num) + '.l2.ldscore.gz'
	gene_weight_files = []
	for pseudotissue_name in pseudotissue_names:

		filer = preprocessed_tgfm_sldsc_data_dir + 'tissue_eqtl.' + str(chrom_num) + '.' + pseudotissue_name + '_' + gene_type + '_' + gene_model_suffix
		gene_weight_files.append(filer)
	gene_weight_files = np.asarray(gene_weight_files)
	
	# Create output root
	weight_output_file = preprocessed_tgfm_sldsc_data_dir + gene_type + '_' + tissue_version + '_' + gene_model_suffix + '.' + str(chrom_num) + '.gene_weights'

	# Merge ld scores files
	print_weight_output_file(variant_ld_score_file, gene_weight_files, pseudotissue_names, weight_output_file)


def make_ld_score_input_shell(chrom_num,  n_genes_per_pseudotissue, pseudotissue_names, variant_models, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir, gene_type, tissue_version):
	for variant_model in variant_models:
		variant_ld_score_file = preprocessed_tgfm_sldsc_data_dir + variant_model + '.' + str(chrom_num) + '.l2.ldscore.gz'
		variant_m_file = preprocessed_tgfm_sldsc_data_dir + variant_model + '.' + str(chrom_num) + '.l2.M'
		variant_m_5_50_file = preprocessed_tgfm_sldsc_data_dir + variant_model + '.' + str(chrom_num) + '.l2.M_5_50'
		for gene_model_suffix in gene_model_suffixes:
			gene_ld_score_files = []
			for pseudotissue_name in pseudotissue_names:
				filer = preprocessed_tgfm_sldsc_data_dir + 'tissue_eqtl.' + str(chrom_num) + '.' + pseudotissue_name + '_' + gene_type + '_' + gene_model_suffix
				gene_ld_score_files.append(filer)
			gene_ld_score_files = np.asarray(gene_ld_score_files)

			# Create output root
			output_root = preprocessed_tgfm_sldsc_data_dir + variant_model + '_' + gene_type + '_' + tissue_version + '_' + gene_model_suffix + '.' + str(chrom_num) + '.l2'

			print(variant_model + '  ' + gene_model_suffix)
			# Merge ld scores files
			merged_ld_score_file = output_root + '.ldscore'
			merge_variant_and_gene_ld_score_files(variant_ld_score_file, gene_ld_score_files, pseudotissue_names, merged_ld_score_file)

			merged_m_file = output_root + '.M'
			merge_m_files(variant_m_file, n_genes_per_pseudotissue, merged_m_file)

			merged_m_5_50_file = output_root + '.M_5_50'
			merge_m_files(variant_m_5_50_file, n_genes_per_pseudotissue, merged_m_5_50_file)

def make_ld_score_input_shell_for_genotype_intercept(chrom_num, n_genes_per_pseudotissue, pseudotissue_names, reference_variant_model, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir, gene_type, tissue_version):
	variant_ld_score_file = preprocessed_tgfm_sldsc_data_dir + reference_variant_model + '.' + str(chrom_num) + '.l2.ldscore.gz'
	variant_m_file = preprocessed_tgfm_sldsc_data_dir + reference_variant_model + '.' + str(chrom_num) + '.l2.M'
	variant_m_5_50_file = preprocessed_tgfm_sldsc_data_dir + reference_variant_model + '.' + str(chrom_num) + '.l2.M_5_50'
	for gene_model_suffix in gene_model_suffixes:
		gene_ld_score_files = []
		for pseudotissue_name in pseudotissue_names:
			filer = preprocessed_tgfm_sldsc_data_dir + 'tissue_eqtl.' + str(chrom_num) + '.' + pseudotissue_name + '_' + gene_type + '_' + gene_model_suffix
			gene_ld_score_files.append(filer)
		gene_ld_score_files = np.asarray(gene_ld_score_files)

		# Create output root
		output_root = preprocessed_tgfm_sldsc_data_dir + 'genotype_intercept' + '_' + gene_type + '_' + tissue_version + '_' + gene_model_suffix + '.' + str(chrom_num) + '.l2'

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



def generate_annotation_sdev_files_shell(variant_models, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir, gene_type, tissue_version):
	for variant_model in variant_models:
		for gene_model in gene_model_suffixes:
			ld_score_file_stem = preprocessed_tgfm_sldsc_data_dir + variant_model + '_' + gene_type + '_' + tissue_version + '_' + gene_model #.11.l2.ldscore
			create_annotation_sdev_file(ld_score_file_stem)

def generate_annotation_sdev_file_for_genotype_intercept_shell(reference_variant_model, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir, gene_type, tissue_version):
	for gene_model in gene_model_suffixes:
		ld_score_file_stem = preprocessed_tgfm_sldsc_data_dir + 'genotype_intercept' + '_' + gene_type + '_' + tissue_version + '_' + gene_model
		create_annotation_sdev_file(ld_score_file_stem)

def get_ld_score_correlation_heatmap(variant_model, gene_model_suffix, preprocessed_tgfm_sldsc_data_dir):
	ld_scores = []
	for chrom_num in range(1,10):
		print(chrom_num)
		output_root = preprocessed_tgfm_sldsc_data_dir + variant_model + '_' + gene_model_suffix + '.' + str(chrom_num) + '.l2'
		# Merge ld scores files
		merged_ld_score_file = output_root + '.ldscore'
		head_count = 0
		f = open(merged_ld_score_file)
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				if chrom_num == 1:
					header = np.asarray(data[3:])
				continue
			ld_scores.append(np.asarray(data[3:]).astype(float))
		f.close()
	ld_scores = np.asarray(ld_scores)

	correlation_mat_output_root = preprocessed_tgfm_sldsc_data_dir + variant_model + '_' + gene_model_suffix + '_ld_score_annotation_correlation_matrix'
	full_anno_corr = np.corrcoef(np.transpose(ld_scores))
	full_anno_row_names = header
	full_anno_col_names = header
	np.savetxt(correlation_mat_output_root + '.txt', full_anno_corr, fmt="%s", delimiter='\t')
	np.savetxt(correlation_mat_output_root + '_row_names.txt', full_anno_row_names, fmt="%s", delimiter='\t')
	np.savetxt(correlation_mat_output_root + '_col_names.txt', full_anno_col_names, fmt="%s", delimiter='\t')


	subset_anno_corr = full_anno_corr[:93, 93:]
	subset_anno_row_names = header[:93]
	subset_anno_col_names = header[93:]
	np.savetxt(correlation_mat_output_root + '_subset.txt', subset_anno_corr, fmt="%s", delimiter='\t')
	np.savetxt(correlation_mat_output_root + '_subset_row_names.txt', subset_anno_row_names, fmt="%s", delimiter='\t')
	np.savetxt(correlation_mat_output_root + '_subset_col_names.txt', subset_anno_col_names, fmt="%s", delimiter='\t')


def get_number_of_variant_annotations(variant_ld_score_file):
	f = gzip.open(variant_ld_score_file)
	head_count = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			n_anno = len(data[3:])
			break
	f.close()


	return n_anno


gtex_pseudotissue_file = sys.argv[1]
preprocessed_tgfm_sldsc_data_dir = sys.argv[2]
gene_type = sys.argv[3]
gtex_susie_gene_models_dir = sys.argv[4]
tissue_version = sys.argv[5]


# Extract names of pseudotissues
pseudotissue_names = get_pseudotissue_names(gtex_pseudotissue_file, tissue_version)


# Get number of genes per chromosome per tissue
n_genes_per_chromosome_per_tissue = get_n_genes_per_chromosome_per_tissue(pseudotissue_names, gtex_susie_gene_models_dir, gene_type)


# Various iterations to run over
#variant_models = ['baselineLD_no_qtl']
variant_models = ['baselineLD_no_qtl', 'baseline_no_qtl']

# Gene models
gene_model_suffixes = ['pmces_gene_adj_ld_scores']
#gene_model_suffixes = ['gene_ld_scores', 'gene_adj_ld_scores', 'pmces_gene_ld_scores', 'pmces_gene_adj_ld_scores']

for variant_model in variant_models:
	chrom_num = 21 # Actual chromosome choice doesn't matter because all chromosomes have the same number of annotations
	variant_ld_score_file = preprocessed_tgfm_sldsc_data_dir + variant_model + '.' + str(chrom_num) + '.l2.ldscore.gz'
	n_variant_annotations = get_number_of_variant_annotations(variant_ld_score_file)
	n_gene_annotations = len(pseudotissue_names)

	# Print non-negative coefficient file
	non_neg_coef_file = preprocessed_tgfm_sldsc_data_dir + variant_model + '_' + tissue_version + '_nonnegative_coefficients.txt'
	t = open(non_neg_coef_file,'w')
	for variant_iter in range(n_variant_annotations):
		t.write('False\n')
	for gene_iter in range(n_gene_annotations):
		t.write('True\n')
	t.close()
# Do the same for genotype intercept model
non_neg_coef_file = preprocessed_tgfm_sldsc_data_dir + 'genotype_intercept' + '_' + tissue_version + '_nonnegative_coefficients.txt'
t = open(non_neg_coef_file,'w')
n_gene_annotations = len(pseudotissue_names)
t.write('False\n')
for gene_iter in range(n_gene_annotations):
	t.write('True\n')
t.close()



for chrom_num in range(1,23):
	print(chrom_num)
	make_ld_score_input_shell(chrom_num,  n_genes_per_chromosome_per_tissue[(chrom_num-1),:], pseudotissue_names, variant_models, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir, gene_type, tissue_version)
generate_annotation_sdev_files_shell(variant_models, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir, gene_type, tissue_version)


# Make version without genomic annotaitons (call it genotype intercept)
reference_variant_model = 'baselineLD_no_qtl'
for chrom_num in range(1,23):
	make_ld_score_input_shell_for_genotype_intercept(chrom_num,  n_genes_per_chromosome_per_tissue[(chrom_num-1),:], pseudotissue_names, reference_variant_model, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir, gene_type, tissue_version)
generate_annotation_sdev_file_for_genotype_intercept_shell(reference_variant_model, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir, gene_type, tissue_version)





















# Generate gene weights files
'''
variant_model = 'baselineLD_no_qtl'  # Simply using this to get rs-ids (could also use baseline or intercept)
gene_model ='pmces_gene_weights'
for chrom_num in range(1,23):
	print(chrom_num)
	make_gene_weight_input_file(chrom_num, pseudotissue_names, variant_model, gene_model, preprocessed_tgfm_sldsc_data_dir, gene_type, tissue_version)
'''

'''
variant_model = 'baselineLD_no_qtl'
gene_model_suffix = 'pmces_gene_adj_ld_scores'

get_ld_score_correlation_heatmap(variant_model, gene_model_suffix, preprocessed_tgfm_sldsc_data_dir)
'''
