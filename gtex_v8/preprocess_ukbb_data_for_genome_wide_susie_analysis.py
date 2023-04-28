import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import os
import pdb
import numpy as np
from pandas_plink import read_plink1_bin
import pickle
import scipy.sparse



def extract_ukbb_trait_names_and_file(ukbb_trait_file):
	f = open(ukbb_trait_file)
	head_count = 0
	trait_names = []
	trait_files = []
	trait_sample_sizes = []
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		trait_names.append(data[0])
		trait_files.append(data[1])
		trait_sample_sizes.append(data[2])
	f.close()
	return np.asarray(trait_names), np.asarray(trait_files), np.asarray(trait_sample_sizes)



def get_gtex_variants_on_this_chromosome(gtex_chrom_bim_file):
	gtex_variants = {}
	f = open(gtex_chrom_bim_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 6:
			print('assumption eroror')
			pdb.set_trace()
		variant_id = data[1]
		if variant_id in gtex_variants:
			print('assumption eroror')
			pdb.set_trace()
		gtex_variants[variant_id] = 1
	f.close()
	return gtex_variants

def get_window_names_on_this_chromosome(genome_wide_window_file, chrom_num):
	f = open(genome_wide_window_file)
	window_dictionary = {}
	chromosome = ['NULL']*300000000
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		line_chrom_num = data[0]
		# Limit to windows on this chromosome
		if line_chrom_num != chrom_num:
			continue
		window_name = data[0] + ':' + data[1] + ':' + data[2]
		if window_name in window_dictionary:
			print('assumption eroror')
			pdb.set_trace()
		window_dictionary[window_name] = []
		start = int(data[1])
		end = int(data[2])
		for bp_pos in range(start,end):
			# Position has no elements
			if chromosome[bp_pos] == 'NULL':
				chromosome[bp_pos] = window_name
			else:
				chromosome[bp_pos] = chromosome[bp_pos] + ',' + window_name
	f.close()
	return chromosome, window_dictionary

def identify_ukbb_repeat_variants(trait_file, chrom_num):
	f = open(trait_file)
	all_variants = {}
	repeat_variants = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Skip variants not on this chromosome
		line_chrom_num = data[1]
		if line_chrom_num != chrom_num:
			continue
		line_variant_id = 'chr' + data[1] + '_' + data[2] + '_' + data[4] + '_' + data[5] + '_b38'
		line_variant_id_alt = 'chr' + data[1] + '_' + data[2] + '_' + data[5] + '_' + data[4] + '_b38'

		if line_variant_id in all_variants or line_variant_id_alt in all_variants:
			repeat_variants[line_variant_id] = 1
			repeat_variants[line_variant_id_alt] = 1
		all_variants[line_variant_id] = 1
		all_variants[line_variant_id_alt] = 1
	f.close()
	return repeat_variants		

def fill_in_window_dictionary_for_single_trait(trait_name, trait_file, chrom_num, window_names_chromosome, gtex_variants, window_dictionary):
	repeat_variants = identify_ukbb_repeat_variants(trait_file, chrom_num)
	f = open(trait_file)
	head_count = 0
	used_variants = {}
	afs = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Skip variants not on this chromosome
		line_chrom_num = data[1]
		if line_chrom_num != chrom_num:
			continue
		# Skip variants with MAF < .01
		#af = float(data[6])
		#if af > .99 or af < .01:
		#	continue
		# Get variant id for this line
		line_variant_id = 'chr' + data[1] + '_' + data[2] + '_' + data[4] + '_' + data[5] + '_b38'
		line_variant_id_alt = 'chr' + data[1] + '_' + data[2] + '_' + data[5] + '_' + data[4] + '_b38'
		variant_pos = int(data[2])
		rs_id = data[0]

		# Ignore strand ambiguous variants from analysis
		if data[4] == 'A' and data[5] == 'T':
			continue
		if data[4] == 'T' and data[5] == 'A':
			continue
		if data[4] == 'C' and data[5] == 'G':
			continue
		if data[4] == 'G' and data[5] == 'C':
			continue

		# Some repeats in this file..
		# Skip variants seen twice in UKBB file
		if line_variant_id in repeat_variants or line_variant_id_alt in repeat_variants:
			continue
		# This just double checks that we have correctly filtered out repear variants
		if line_variant_id in used_variants or line_variant_id_alt in used_variants:
			print('variant repeat assumption eroror')
			pdb.set_trace()
		used_variants[line_variant_id] =1
		used_variants[line_variant_id_alt] = 1

		'''
		if rs_id.startswith('rs') == False:
			print('non rsid variant')
			continue
		beta = float(data[10])
		std_err = float(data[11])
		stringer = trait_name + ',' + line_variant_id + ',' + str(beta) + ',' + str(std_err) + ',' + rs_id
		window_string = window_names_chromosome[variant_pos]
		windows = window_string.split(',')
		for window in windows:
			if window == 'NULL':  # This is necessary because we threw out some long range LD windows
				continue
			window_dictionary[window].append(stringer)
		'''
		# Variant id in gtex variants
		if line_variant_id in gtex_variants:
			if line_variant_id_alt in gtex_variants:
				print('assumption error')
				pdb.set_trace()
			if rs_id.startswith('rs') == False:
				print('non rsid variant')
				continue
			beta = float(data[10])
			std_err = float(data[11])
			stringer = trait_name + ',' + line_variant_id + ',' + str(beta) + ',' + str(std_err) + ',' + rs_id
			window_string = window_names_chromosome[variant_pos]
			windows = window_string.split(',')
			for window in windows:
				if window == 'NULL':  # This is necessary because we threw out some long range LD windows
					continue
				window_dictionary[window].append(stringer)
		# Flipped variant id in gtex variants 
		elif line_variant_id_alt in gtex_variants:
			if rs_id.startswith('rs') == False:
				print('non rsid variant')
				continue
			beta = float(data[10])
			std_err = float(data[11])
			stringer = trait_name + ',' + line_variant_id + ',' + str(beta) + ',' + str(std_err) + ',' + rs_id
			window_string = window_names_chromosome[variant_pos]
			windows = window_string.split(',')
			for window in windows:
				if window == 'NULL':  # This is necessary because we threw out some long range LD windows
					continue
				window_dictionary[window].append(stringer)
	f.close()

	return window_dictionary

def fill_in_window_dictionary(trait_names, trait_files, chrom_num, window_names_chromosome, gtex_variants, window_dictionary):
	for trait_index, trait_name in enumerate(trait_names):
		print(trait_name)
		trait_file = trait_files[trait_index]
		# need to first line valid variants
		window_dictionary = fill_in_window_dictionary_for_single_trait(trait_name, trait_file, chrom_num, window_names_chromosome, gtex_variants, window_dictionary)
	return window_dictionary

def organize_window_test_arr(gene_test_arr):
	# This will take two pass
	# In first pass get list of tissues and get list of variants
	tissues = {}
	variants = {}
	variant_to_rs_id = {}
	for gene_test in gene_test_arr:
		gene_test_info = gene_test.split(',')
		tissues[gene_test_info[0]] = 1
		variant_id = gene_test_info[1]
		rs_id = gene_test_info[4]
		variants[variant_id] = 1

		if variant_id in variant_to_rs_id:
			if variant_to_rs_id[variant_id] != rs_id:
				print('assumption erroror')
				pdb.set_trace()
		else:
			variant_to_rs_id[variant_id] = rs_id

	# Convert tissue and variant dictionaries to arrays
	tissue_vec = np.sort([*tissues])
	variant_vec = np.sort([*variants])

	rs_id_vec = []
	for variant in variant_vec:
		rs_id = variant_to_rs_id[variant]
		rs_id_vec.append(rs_id)
	rs_id_vec = np.asarray(rs_id_vec)
	if len(rs_id_vec) != len(np.unique(rs_id_vec)):
		print('rs id assumptino eororro')
		pdb.set_trace()




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
		variant_id = gene_test_info[1]
		beta = float(gene_test_info[2])
		std_err_beta = float(gene_test_info[3])

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

	return beta_mat, std_err_mat, new_variant_vec, tissue_vec, rs_id_vec


def create_mapping_from_rsid_to_in_sample_variant_index(chrom_pvar_file):
	f = open(chrom_pvar_file)
	dicti = {}
	rs_id_to_alleles = {}
	head_count = 0
	indexer = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[2]
		alleles = data[4] + '_' + data[3]
		if rsid in rs_id_to_alleles:
			print('assumption eroror')
			pdb.set_trace()
		rs_id_to_alleles[rsid] = alleles
		if rsid in dicti:
			print('assumption eroror')
			pdb.set_trace()
		dicti[rsid] = indexer
		indexer = indexer + 1
	f.close()
	return dicti, rs_id_to_alleles

def quick_error_check_to_make_sure_ukbb_in_sample_variants_match_ukbb_sumstat_variants(window_variant_arr, window_rs_id, rs_id_to_in_sample_alleles):
	missing_count = 0
	for variant_index, variant_name in enumerate(window_variant_arr):
		rs_id = window_rs_id[variant_index]
		if rs_id not in rs_id_to_in_sample_alleles:
			missing_count = missing_count + 1
			continue
		in_sample_alleles = rs_id_to_in_sample_alleles[rs_id]
		variant_alleles = '_'.join(variant_name.split('_')[2:4])
		if variant_alleles != in_sample_alleles:
			print('in sample ld variant id assumption erororor')
			pdb.set_trace()
	pdb.set_trace()
	return

def extract_overlapping_variants(window_rsids, window_variant_arr, rs_id_to_in_sample_variant, rs_id_to_in_sample_alleles):
	valid_window_indices = []
	in_sample_variant_indices = []
	flip_values = []

	for ii, window_rsid in enumerate(window_rsids):
		if window_rsid in rs_id_to_in_sample_variant:
			#valid_window_indices.append(ii)
			#in_sample_variant_indices.append(rs_id_to_in_sample_variant[window_rsid])

			window_variant_id = window_variant_arr[ii]
			variant_alleles = '_'.join(window_variant_id.split('_')[2:4])
			alt_variant_alleles = window_variant_id.split('_')[3] + '_' + window_variant_id.split('_')[2]
			in_sample_alleles = rs_id_to_in_sample_alleles[window_rsid]
			if variant_alleles == in_sample_alleles:
				flip_values.append(1.0)
				valid_window_indices.append(ii)
				in_sample_variant_indices.append(rs_id_to_in_sample_variant[window_rsid])
			elif alt_variant_alleles == in_sample_alleles:
				flip_values.append(-1.0)
				valid_window_indices.append(ii)
				in_sample_variant_indices.append(rs_id_to_in_sample_variant[window_rsid])

	
	if len(valid_window_indices) != len(flip_values):
		print('assumptino eroror')
		pdb.set_trace()

	if len(valid_window_indices) != len(in_sample_variant_indices):
		print('assumptino eroror')
		pdb.set_trace()

	if len(in_sample_variant_indices) != len(np.unique(in_sample_variant_indices)):
		print('assumptione rororo')
		pdb.set_trace()

	return np.asarray(valid_window_indices), np.asarray(in_sample_variant_indices), np.asarray(flip_values)


def parse_snp_range(snp_range):
	"""Get range of SNPs to analyze. 
	
	Parameters
	----------
	snp_range: str
		Example: 'c1_s0_e2000_r1'
	Returns
	-------
	dic_range : dict
		
		- dic_range['chr'] : chromosome
		- dic_range['start'] : start position
		- dic_range['end'] : end position
		- dic_range['chr_ref'] : reference chromosome list (List)
	"""

	dic_range = {x: None for x in ["chr", "start", "end", "chr_ref"]}

	for x in snp_range.split("_"):

		if x[0] == "c":
			dic_range["chr"] = int(x.replace("c", "").strip())

		if x[0] == "s":
			dic_range["start"] = int(x.replace("s", "").strip())

		if x[0] == "e":
			dic_range["end"] = int(x.replace("e", "").strip())

		if x[0] == "r":
			temp_str = x.replace("r", "").strip()
			if temp_str == "all":
				dic_range["chr_ref"] = list(np.arange(1, 23))
			else:
				dic_range["chr_ref"] = [int(x) for x in temp_str.split(",")]

	return dic_range

def read_ld(fpath):
	"""
	Read LD files.
		- `_fullld.npy` : full_ld matrix, np.array(dtype=np.float32)
		- `_ld.npz` : ld matrix with SNPs in 10MB window, sp.sparse.csc_matrix(dtype=np.float32)
	Parameters
	----------
	fpath: str
		LD file path.
	Returns
	-------
	mat_ld : np.array(dtype=np.float32) or sp.sparse.csc_matrix(dtype=np.float32)
		LD matrix of dimension (n_ref_snp, n_snp)
	dic_range : dict
		
		- dic_range['chr'] : chromosome
		- dic_range['start'] : start position
		- dic_range['end'] : end position
		- dic_range['chr_ref'] : reference chromosome list (List)      
	"""
	
	# Check fpath
	err_msg = "fpath should end with one of ['_fullld.npy', '_ld.npz'] : %s" % fpath
	assert fpath.endswith("_fullld.npy") | fpath.endswith("_ld.npz"), err_msg
	
	if fpath.endswith("_fullld.npy"):
		mat_ld = np.load(fpath)
		temp_str = [x for x in fpath.split('.') if x.endswith('_fullld')][0]
		dic_range = parse_snp_range(temp_str)
		
	if fpath.endswith("_ld.npz"):
		mat_ld = scipy.sparse.load_npz(fpath)
		temp_str = [x for x in fpath.split('.') if x.endswith('_ld')][0]
		dic_range = parse_snp_range(temp_str)

	return mat_ld,dic_range


def indices_dont_lie_in_file(sample_ld_variant_indices, file_index_start, file_index_end):
	booler = True
	for indexer in sample_ld_variant_indices:
		if indexer >= file_index_start and indexer < file_index_end:
			booler = False
	return booler

def extract_ld_mat_from_in_sample_ld(sample_ld_variant_indices, ukbb_in_sample_ld_dir, chrom_num):
	min_index = np.min(sample_ld_variant_indices)
	max_index = np.max(sample_ld_variant_indices)

	num_var = len(sample_ld_variant_indices)

	ld_mat = np.zeros((num_var, num_var)) + -2000.0

	for file_name in os.listdir(ukbb_in_sample_ld_dir):
		if file_name.startswith('ukb_imp_v3.c' + chrom_num + '_') == False:
			continue
		if file_name.endswith('.compute_ld.sbatch.log'):
			continue
		full_file_name = ukbb_in_sample_ld_dir + file_name

		file_info = file_name.split('_')
		file_index_start = int(file_info[3].split('s')[1])
		file_index_end = int(file_info[4].split('e')[1])

		if indices_dont_lie_in_file(sample_ld_variant_indices, file_index_start, file_index_end):
			continue

		file_indices = []
		col_names = []
		for ii, sample_index in enumerate(sample_ld_variant_indices):
			if sample_index >= file_index_start and sample_index < file_index_end:
				file_indices.append(ii)
				col_names.append(sample_index - file_index_start)
		file_indices = np.asarray(file_indices)
		col_names = np.asarray(col_names)

		sparse_mat, sparse_mat_info = read_ld(full_file_name)
		ld_mat[:, file_indices] = (sparse_mat[sample_ld_variant_indices,:][:,col_names]).toarray()

	if np.sum(ld_mat == -2000.0) > 0:
		print('assumption eroror')
		pdb.set_trace()
	return ld_mat


def get_ukbb_variants_on_this_chromosome(bim_file):
	variants = {}
	f = open(bim_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		snp_id = 'chr' + data[0] + '_' + data[3] + '_' + data[4] + '_' + data[5] + '_b38'
		variants[snp_id] = 1
	f.close()
	return variants

def correct_ld_mat_for_af_standardization(ld_mat):
	n_snps = ld_mat.shape[0]
	correction = 1.0/np.diag(ld_mat)
	for snp_iter in range(n_snps):
		ld_mat[:,snp_iter] = ld_mat[:,snp_iter]*np.sqrt(correction[snp_iter])
		ld_mat[snp_iter,:] = ld_mat[snp_iter,:]*np.sqrt(correction[snp_iter])
	return ld_mat


#######################
# Command line args
chrom_num = sys.argv[1]
genome_wide_window_file = sys.argv[2]
ukbb_sumstats_hg38_dir = sys.argv[3]
gtex_genotype_dir = sys.argv[4]
ref_1kg_genotype_dir = sys.argv[5]
ukbb_preprocessed_for_genome_wide_susie_dir = sys.argv[6]
ukbb_in_sample_ld_dir = sys.argv[7]
ukbb_in_sample_genotype_dir = sys.argv[8]
##########################





# create mapping from RS_ID to UKBB in_sample LD variant INDEX
chrom_pvar_file = ukbb_in_sample_genotype_dir + 'ukb_imp_chr' + chrom_num + '_v3.pvar'
rs_id_to_in_sample_variant, rs_id_to_in_sample_alleles = create_mapping_from_rsid_to_in_sample_variant_index(chrom_pvar_file)



# Extract names of UKBB studies for this analysis
ukbb_trait_file = ukbb_sumstats_hg38_dir + 'ukbb_hg38_sumstat_files_with_samp_size_and_h2.txt'
#ukbb_trait_file = ukbb_sumstats_hg38_dir + 'ukbb_hg38_sumstat_files_with_samp_size_and_h2_temp_small.txt'  # Just used for debugging
trait_names, trait_files, trait_sample_sizes = extract_ukbb_trait_names_and_file(ukbb_trait_file)
num_traits = len(trait_files)

# create mapping from trait name to sample size
trait_name_to_sample_size = {}
for trait_index, trait_name in enumerate(trait_names):
	trait_sample_size = trait_sample_sizes[trait_index]
	trait_name_to_sample_size[trait_name] = trait_sample_size


# First extract list of gtex variants in UKBB [we will be using gtex variant orientation]
# Also note that these gtex variants are also found in UKBB
gtex_variants = get_gtex_variants_on_this_chromosome(gtex_genotype_dir + 'Whole_Blood_GTEx_v8_genotype_EUR_' + chrom_num + '.bim')  #note: whole blood is randomly choosen but really doesn't matter b/c all tissues have the same variants
ukbb_variants = get_ukbb_variants_on_this_chromosome(ref_1kg_genotype_dir + '1000G.EUR.hg38.' + str(chrom_num) + '.bim')

# Get names of windows on this chromosome (in data structure where of array of length chromosome)
window_names_chromosome, window_dictionary = get_window_names_on_this_chromosome(genome_wide_window_file, chrom_num)



# Fill in the dictionary with each elent in list is a string corresponding to info on a cis snp
window_dictionary = fill_in_window_dictionary(trait_names, trait_files, chrom_num, window_names_chromosome, ukbb_variants, window_dictionary)



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


output_file = ukbb_preprocessed_for_genome_wide_susie_dir + 'genome_wide_susie_windows_and_processed_data_chrom_' + chrom_num + '.txt'
t = open(output_file,'w')
f = open(genome_wide_window_file)



head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\tbeta_file\tstd_err_file\tvariant_file\ttissue_file\tref_genotype_file\tukbb_in_sample_ld_file\tsample_size_file\n')
		continue
	line_chrom_num = data[0]
	# Limit to windows on this chromosome
	if line_chrom_num != chrom_num:
		continue
	# Info about this window
	window_name = data[0] + ':' + data[1] + ':' + data[2]

	# Extract SNPs and gwas effects in this window
	window_test_arr = window_dictionary[window_name]

	# Convert from window_test_arr to beta-matrix, standard-errr matrix, and ordered-variant list
	window_beta_mat, window_std_err_mat, window_variant_arr, window_study_arr, window_rs_id = organize_window_test_arr(window_test_arr)

	# Throw out windows with fewer than 50 variants
	if len(window_variant_arr) < 50:
		continue

	# Get sample sizes
	window_sample_sizes = []
	for study_name in window_study_arr:
		window_sample_sizes.append(trait_name_to_sample_size[study_name])
	window_sample_sizes = np.asarray(window_sample_sizes)

	# Extract indices in reference data that the above snps correspond to
	# Also extract whether above snp is flipped relative to reference data
	reference_indices = []
	flip_values = []
	for variant_id in window_variant_arr:
		info = ref_snp_mapping[variant_id]
		reference_indices.append(info[0])
		flip_values.append(info[1])
	reference_indices = np.asarray(reference_indices)
	flip_values = np.asarray(flip_values)
	
	# Get reference genotype for this gene
	gene_reference_genotype = G[:, reference_indices]


	# Correct for flips in 1KG data (assuming UKBB is reference)
	for snp_index, flip_value in enumerate(flip_values):
		if flip_value == -1.0:
			# correct for flips (assuming UKBB is reference)
			gene_reference_genotype[:, snp_index] = 2.0 - gene_reference_genotype[:, snp_index]
			# correct for flips (assuming 1KG is reference)
			#window_beta_mat[:, snp_index] = window_beta_mat[:, snp_index]*-1.0
			#old_snp = window_variant_arr[snp_index]
			#snp_info = old_snp.split('_')
			#new_snp = snp_info[0] + '_' + snp_info[1] + '_' + snp_info[3] + '_' + snp_info[2]
			#window_variant_arr[snp_index] = new_snp

	# Quick error check
	if len(window_variant_arr) != len(np.unique(window_variant_arr)):
		print('assumptione roror')
		pdb.set_trace()


	# Now deal with UKBB in-sample LD-matrix
	# Make sure our variants line up with UKBB in sample LD (otherwise code will fail)
	#quick_error_check_to_make_sure_ukbb_in_sample_variants_match_ukbb_sumstat_variants(window_variant_arr, window_rs_id, rs_id_to_in_sample_alleles)
	# Extract valid window indices, as well as sample ld variant indices corresponding to those
	valid_window_indices, sample_ld_variant_indices, sample_ld_flips = extract_overlapping_variants(window_rs_id, window_variant_arr, rs_id_to_in_sample_variant, rs_id_to_in_sample_alleles)

	# Throw out windows with fewer than 50 variants
	if len(sample_ld_variant_indices) < 50:
		continue

	# NOW GET LD matrix
	ukbb_in_sample_ld_mat = extract_ld_mat_from_in_sample_ld(sample_ld_variant_indices, ukbb_in_sample_ld_dir, chrom_num)
	# Flip alleles to make sure it matches summary statistics
	for var_index, sample_ld_flip_value in enumerate(sample_ld_flips):
		if sample_ld_flip_value == -1.0:
			ukbb_in_sample_ld_mat[var_index, :] = ukbb_in_sample_ld_mat[var_index, :]*-1.0
			ukbb_in_sample_ld_mat[:, var_index] = ukbb_in_sample_ld_mat[:, var_index]*-1.0

	#ukbb_in_sample_ld_mat = correct_ld_mat_for_af_standardization(ukbb_in_sample_ld_mat)

	# Save data to output file
	# Beta file
	beta_file = ukbb_preprocessed_for_genome_wide_susie_dir + window_name + '_beta.txt'
	np.savetxt(beta_file, window_beta_mat[:,valid_window_indices], fmt="%s", delimiter='\t')
	# stderr file
	stderr_file = ukbb_preprocessed_for_genome_wide_susie_dir + window_name + '_beta_std_err.txt'
	np.savetxt(stderr_file, window_std_err_mat[:,valid_window_indices], fmt="%s", delimiter='\t')
	# Variant file
	variant_file = ukbb_preprocessed_for_genome_wide_susie_dir + window_name + '_variant_ids.txt'
	np.savetxt(variant_file, window_variant_arr[valid_window_indices], fmt="%s", delimiter='\t')
	# Tissue file
	study_file = ukbb_preprocessed_for_genome_wide_susie_dir + window_name + '_studies.txt'
	np.savetxt(study_file, window_study_arr, fmt="%s", delimiter='\t')
	# Reference Genotype
	ref_geno_file = ukbb_preprocessed_for_genome_wide_susie_dir + window_name + '_ref_1kg_genotype.txt'
	np.savetxt(ref_geno_file, gene_reference_genotype[:, valid_window_indices], fmt="%s", delimiter='\t')
	# In sample LD
	in_sample_ld_file = ukbb_preprocessed_for_genome_wide_susie_dir + window_name + '_ukbb_in_sample_ld.npy'
	np.save(in_sample_ld_file, ukbb_in_sample_ld_mat)
	# Sample sizes
	sample_size_file = ukbb_preprocessed_for_genome_wide_susie_dir + window_name + '_study_sample_sizes.txt'
	np.savetxt(sample_size_file, window_sample_sizes, fmt="%s", delimiter='\t')

	t.write(line + '\t' + beta_file + '\t' + stderr_file + '\t' + variant_file + '\t' + study_file + '\t' + ref_geno_file + '\t' + in_sample_ld_file + '\t' + sample_size_file + '\n')

f.close()
