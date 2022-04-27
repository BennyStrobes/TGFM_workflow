import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import os
import pdb
import numpy as np
from pandas_plink import read_plink1_bin
import pickle


def get_gtex_tissues(gtex_tissue_file):
	f = open(gtex_tissue_file)
	head_count = 0
	arr = []
	arr2 = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0])
		arr2.append(data[2])
	f.close()
	return np.asarray(arr), np.asarray(arr2)


def get_genes_on_chromosome(cafeh_gene_list, chrom_num):
	f = open(cafeh_gene_list)
	counter = 0
	head_count = 0
	genes = []
	gene_dicti = {}
	gene_to_tissues = {}
	gene_to_tss = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		counter = counter + 1
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_name = data[0]
		line_chrom_num = data[1]
		if line_chrom_num != 'chr' + chrom_num:
			continue
		genes.append(gene_name)
		gene_dicti[gene_name] = []
		gene_tissues = np.sort(data[2].split(';'))
		gene_to_tissues[gene_name] = gene_tissues
		gene_to_tss[gene_name] = int(data[3])
	f.close()
	return np.asarray(genes), gene_to_tissues, gene_to_tss

def array_to_dictionary_of_lists(ordered_genes):
	dicti = {}
	for gene in ordered_genes:
		dicti[gene] = []
	return dicti

def fill_in_gene_dicti_from_single_tissue_summary_stats(sum_stats_file, gene_dicti, gtex_tissue):
	f = open(sum_stats_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[1]
		if gene_id not in gene_dicti:
			continue
		variant_id = data[0]
		beta = float(data[2])
		t_stat = float(data[3])
		beta_std_err = beta/t_stat
		stringer = gtex_tissue + ',' + gene_id + ',' + variant_id + ',' + str(beta) + ',' + str(beta_std_err)
		gene_dicti[gene_id].append(stringer)
	f.close()
	return gene_dicti

def fill_in_gene_dicti(gtex_tissues, meta_analyzed_booleans, chrom_num, eqtl_summary_stats_dir, gene_dicti):
	for i,gtex_tissue in enumerate(gtex_tissues):
		print(gtex_tissue)
		meta_analyzed_boolean = meta_analyzed_booleans[i]
		if meta_analyzed_boolean == 'False':
			sum_stats_file = eqtl_summary_stats_dir + gtex_tissue + '_chr' + chrom_num + '_matrix_eqtl_results.txt'
		elif meta_analyzed_boolean == 'True':
			sum_stats_file = eqtl_summary_stats_dir + gtex_tissue + '_meta_analyis_chr' + chrom_num + '_matrix_eqtl_results.txt'
		else:
			print('assumption eroror')
			pdb.set_trace()

		gene_dicti = fill_in_gene_dicti_from_single_tissue_summary_stats(sum_stats_file, gene_dicti, gtex_tissue)
	return gene_dicti

def load_in_ukbb_variants_and_their_flips(ukbb_variant_file, chr_num):
	dicti = {}
	chrom_string = 'chr' + chr_num
	f = open(ukbb_variant_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		snp_id1 = data[0]
		snp_info = snp_id1.split('_')
		if snp_info[0] != chrom_string:
			continue
		if len(snp_info) != 4:
			continue
		a1 = snp_info[2]
		a2 = snp_info[3]
		snp_id2 = snp_info[0] + '_' + snp_info[1] + '_' + a2 + '_' + a1
		dicti[snp_id1] = 1
		dicti[snp_id2] = 1
	f.close()
	return dicti

def organize_gene_test_arr(gene_test_arr):
	# This will take two pass
	# In first pass get list of tissues and get list of variants
	tissues = {}
	variants = {}
	for gene_test in gene_test_arr:
		gene_test_info = gene_test.split(',')
		tissues[gene_test_info[0]] = 1
		variants[gene_test_info[2]] = 1
	# Convert tissue and variant dictionaries to arrays
	tissue_vec = np.sort([*tissues])
	variant_vec = np.sort([*variants])

	# Quick error checking
	if len(tissue_vec)*len(variant_vec) != len(gene_test_arr):
		print('assumption eroror')
		pdb.set_trace()

	# Second pass
	# Here we extract beta and standard error matrices
	beta_mat = np.ones((len(tissue_vec), len(variant_vec)))*-800
	std_err_mat = np.ones((len(tissue_vec), len(variant_vec)))*-800

	# Create mapping from variant to position
	variant_to_col_index = {}
	for itera, variant_id in enumerate(variant_vec):
		variant_to_col_index[variant_id] = itera
	# Create mapping from tissue to position
	tissue_to_row_index = {}
	for itera, tissue_id in enumerate(tissue_vec):
		tissue_to_row_index[tissue_id] = itera

	# Now loop through gene_test_arr again
	for gene_test in gene_test_arr:
		gene_test_info = gene_test.split(',')
		tissue_name = gene_test_info[0]
		variant_id = gene_test_info[2]
		beta = float(gene_test_info[3])
		std_err_beta = float(gene_test_info[4])

		row_index = tissue_to_row_index[tissue_name]
		col_index = variant_to_col_index[variant_id]

		beta_mat[row_index, col_index] = beta
		std_err_mat[row_index, col_index] = std_err_beta

	# Quick erorr checking
	if np.sum(beta_mat== -800) != 0.0 or np.sum(std_err_mat== -800) != 0.0:
		print('assumption error')
		pdb.set_trace()

	new_variant_vec = []
	for variant_id in variant_vec:
		variant_info = variant_id.split('_')
		if len(variant_info) != 5:
			print('assumption eroror')
			pdb.set_trace()
		new_variant_id = '_'.join(variant_info[:4])
		new_variant_vec.append(new_variant_id)

	new_variant_vec = np.asarray(new_variant_vec)
	if len(new_variant_vec) != len(variant_vec):
		print('assumption eroror)')
		pdb.set_trace()

	return beta_mat, std_err_mat, new_variant_vec, tissue_vec

def get_indices_of_variants_in_ukbb(gene_variant_arr, ukbb_variants):
	indices = []
	for index, variant in enumerate(gene_variant_arr):
		if variant in ukbb_variants:
			indices.append(index)
	return np.asarray(indices)

def flip_effect_and_reference_allele_in_snp_arr(gene_variant_arr):
	new_variant_arr = []
	for variant_id in gene_variant_arr:
		info = variant_id.split('_')
		new_variant_id = info[0] + '_' + info[1] + '_' + info[3] + '_' + info[2]
		new_variant_arr.append(new_variant_id)
	return np.asarray(new_variant_arr)

######################
# Command line args
######################
chrom_num = int(sys.argv[1])  # Chrom num
gtex_tissue_file = sys.argv[2]  # GTEx tissue files
xt_gene_list_file = sys.argv[3]  #  File containing list of genes to be analyzed
eqtl_summary_stats_dir = sys.argv[4]  # eQTL summary stats (should be limited to variants in ref_1kg_genotype_dir) and should be the same in all tissues
ref_1kg_genotype_dir = sys.argv[5]  # Reference genotype data in plink format
ukbb_sumstats_hg38_dir = sys.argv[6]  # Filter to variants in UKBB. Should be a file in this directory called 'ukbb_hg38_liftover_bgen_snps.txt' that has names of all snps used in their analysis
gtex_preprocessed_for_susie_dir = sys.argv[7]  # Output dir


# Extract ordered list of gtex tissues
gtex_tissues, meta_analyzed_booleans = get_gtex_tissues(gtex_tissue_file)

# Extract ordered list of genes on this chromosome
ordered_genes, gene_to_tissues, gene_to_tss = get_genes_on_chromosome(xt_gene_list_file, str(chrom_num))


# Create dictionary where each key is a gene in ordered_genes and each value is an empty list
gene_dicti = array_to_dictionary_of_lists(ordered_genes)
# Fill in the dictionary with each elent in list is a string corresponding to info on a cis snp
gene_dicti = fill_in_gene_dicti(gtex_tissues, meta_analyzed_booleans, str(chrom_num), eqtl_summary_stats_dir, gene_dicti)
# Temp saving (for code dev.)
#f = open('gene_dicti_temp.pkl', 'wb')
#pickle.dump(gene_dicti, f)
#f.close()
# Temp loading
#gene_dicti = pickle.load(open('gene_dicti_temp.pkl', 'rb'))


# Load in dictionary containing all variants used in ukbb_sumstats_hg38_dir
ukbb_variant_file = ukbb_sumstats_hg38_dir + 'ukbb_hg38_liftover_bgen_snps.txt'
ukbb_variants = load_in_ukbb_variants_and_their_flips(ukbb_variant_file, str(chrom_num))


# Load in Reference Genotype data
geno_stem = ref_1kg_genotype_dir + '1000G.EUR.hg38.' + str(chrom_num) + '.'
G_obj = read_plink1_bin(geno_stem + 'bed', geno_stem + 'bim', geno_stem + 'fam', verbose=False)
G = G_obj.values # Numpy 2d array of dimension num samples X num snps
ref_chrom = np.asarray(G_obj.chrom)
ref_pos = np.asarray(G_obj.pos)
# For our purposes, a0 is the effect allele
# For case of plink package, a0 is the first column in the plink bim file
ref_a0 = np.asarray(G_obj.a0)
ref_a1 = np.asarray(G_obj.a1)
# Extract reference snps names
reference_snp_names = []
reference_alt_snp_names = []
for var_iter in range(len(G_obj.a1)):
	snp_name = 'chr' + ref_chrom[var_iter] + '_' + str(ref_pos[var_iter]) + '_' + ref_a0[var_iter] + '_' + ref_a1[var_iter]
	snp_name_alt = 'chr' + ref_chrom[var_iter] + '_' + str(ref_pos[var_iter]) + '_' + ref_a1[var_iter] + '_' + ref_a0[var_iter]
	reference_snp_names.append(snp_name)
	reference_alt_snp_names.append(snp_name_alt)
reference_snp_names = np.asarray(reference_snp_names)
reference_alt_snp_names = np.asarray(reference_alt_snp_names)
# Create mapping from snp_name to (reference position, refValt)
ref_snp_mapping = {}
for itera in range(len(reference_snp_names)):
	ref_snp_name = reference_snp_names[itera]
	ref_alt_snp_name = reference_alt_snp_names[itera]

	if ref_snp_name in ref_snp_mapping:
		print('extra')
	if ref_alt_snp_name in ref_snp_mapping:
		print('extra')

	ref_snp_mapping[ref_snp_name] = (itera, 1.0)
	ref_snp_mapping[ref_alt_snp_name] = (itera, -1.0)



# Open file handle to keep track of genes and data
organization_output_file = gtex_preprocessed_for_susie_dir + 'chr' + str(chrom_num) + '_susie_input_gene_organization_file.txt'
t = open(organization_output_file,'w')
t.write('gene_id\tchrom_num\tbeta_file\tstd_err_file\tvariant_file\ttissue_file\tref_genotype_file\tsample_size_file\n')

# Loop through genes
for gene_index, test_gene in enumerate(ordered_genes):
	print(test_gene)
	# Extract cis-snps for this gene
	gene_test_arr = gene_dicti[test_gene]

	# Convert from gene_test_arr to beta-matrix, standard-errr matrix, and ordered-variant list
	gene_beta_mat, gene_std_err_mat, gene_variant_arr, gene_tissue_arr = organize_gene_test_arr(gene_test_arr)

	# Currently gene_variant_arr is of format chrom_pos_ref_eff, whereas we will be using nomenclature of: chrom_pos_eff_ref
	gene_variant_arr = flip_effect_and_reference_allele_in_snp_arr(gene_variant_arr)

	# Filter to variants in UKBB
	variants_in_ukbb_indices = get_indices_of_variants_in_ukbb(gene_variant_arr, ukbb_variants)
	if len(variants_in_ukbb_indices) == 0:
		continue
	gene_beta_mat = gene_beta_mat[:, variants_in_ukbb_indices]
	gene_std_err_mat = gene_std_err_mat[:, variants_in_ukbb_indices]
	gene_variant_arr = gene_variant_arr[variants_in_ukbb_indices]

	# Extract indices in reference data that the above snps correspond to
	# Also extract whether above snp is flipped relative to reference data
	reference_indices = []
	flip_values = []
	for variant_id in gene_variant_arr:
		info = ref_snp_mapping[variant_id]
		reference_indices.append(info[0])
		flip_values.append(info[1])
	reference_indices = np.asarray(reference_indices)
	flip_values = np.asarray(flip_values)
	
	# Get reference genotype for this gene
	gene_reference_genotype = G[:, reference_indices]

	# Correct for flips in GTEx data (assuming 1KG is reference)
	for snp_index, flip_value in enumerate(flip_values):
		if flip_value == -1.0:
			gene_beta_mat[:, snp_index] = gene_beta_mat[:, snp_index]*-1.0
			old_snp = gene_variant_arr[snp_index]
			snp_info = old_snp.split('_')
			new_snp = snp_info[0] + '_' + snp_info[1] + '_' + snp_info[3] + '_' + snp_info[2]
			gene_variant_arr[snp_index] = new_snp
	
	# Quick error check
	if len(gene_variant_arr) != len(np.unique(gene_variant_arr)):
		print('assumptione roror')
		pdb.set_trace()

	# Save data to output file
	# Beta file
	beta_file = gtex_preprocessed_for_susie_dir + test_gene + '_beta.txt'
	np.savetxt(beta_file, gene_beta_mat, fmt="%s", delimiter='\t')
	# stderr file
	stderr_file = gtex_preprocessed_for_susie_dir + test_gene + '_beta_std_err.txt'
	np.savetxt(stderr_file, gene_std_err_mat, fmt="%s", delimiter='\t')
	# Variant file
	variant_file = gtex_preprocessed_for_susie_dir + test_gene + '_variant_ids.txt'
	np.savetxt(variant_file, gene_variant_arr, fmt="%s", delimiter='\t')
	# Tissue file
	tissue_file = gtex_preprocessed_for_susie_dir + test_gene + '_tissues.txt'
	np.savetxt(tissue_file, gene_tissue_arr, fmt="%s", delimiter='\t')
	# Reference Genotype
	ref_geno_file = gtex_preprocessed_for_susie_dir + test_gene + '_ref_1kg_genotype.txt'
	np.savetxt(ref_geno_file, gene_reference_genotype, fmt="%s", delimiter='\t')
	# Sample sizes
	sample_size_file = gtex_preprocessed_for_susie_dir + test_gene + '_study_sample_sizes.txt'
	samp_size_vec = np.ones(len(gene_tissue_arr))*320
	np.savetxt(sample_size_file, samp_size_vec, fmt="%s", delimiter='\t')


	t.write(test_gene + '\t' + str(chrom_num) + '\t' + beta_file + '\t' + stderr_file + '\t' + variant_file + '\t' + tissue_file + '\t' + ref_geno_file + '\t' + sample_size_file + '\n')

t.close()


