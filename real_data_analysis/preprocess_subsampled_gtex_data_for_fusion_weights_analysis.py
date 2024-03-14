import numpy as np 
import os
import sys
import pdb

def create_gene_to_chrom_num_and_tss_mapping(xt_pc_gene_list_file):
	f = open(xt_pc_gene_list_file)
	head_count = 0
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_name = data[0]
		chrom_num = data[1]
		tss = data[3]
		if gene_name in dicti:
			print('assumption eroor')
			pdb.set_trace()
		dicti[gene_name] = (chrom_num, tss)
	f.close()
	return dicti


def get_sample_names_from_covariate_file(covariate_file):
	cov = np.loadtxt(covariate_file, dtype=str,delimiter='\t')
	return cov[0,1:]

def convert_from_sample_name_to_family_and_within_family_ids(sample_names):
	family_ids = []
	within_family_ids = []
	for sample_name in sample_names:
		family_ids.append('0')
		within_family_ids.append(sample_name)
	return np.asarray(family_ids), np.asarray(within_family_ids)

def make_covariate_file(study_covariate_file, family_ids, within_family_ids, pb_covariate_file):
	valid_indis = {}
	for indi in within_family_ids:
		valid_indis[indi] = 1

	# Open file handle for gene pheno file
	t3 = open(study_covariate_file,'w')

	# Load in covariates
	cov_data = np.loadtxt(pb_covariate_file, dtype=str, delimiter='\t')

	valid_columns = []
	valid_columns.append(0)
	for ii, indi_id in enumerate(cov_data[0,:]):
		if indi_id in valid_indis:
			valid_columns.append(ii)
	valid_columns = np.asarray(valid_columns)
	cov_data = cov_data[:, valid_columns]

	if np.array_equal(cov_data[0,1:], within_family_ids) == False:
		print('assumption eroorro')
		pdb.set_trace()

	cov_data = np.transpose(cov_data[1:, 1:])

	num_samples = cov_data.shape[0]

	if num_samples != len(family_ids):
		print('assumption errror')
		pdb.set_trace()

	for sample_num in range(num_samples):
		t3.write(family_ids[sample_num] + '\t' + within_family_ids[sample_num] + '\t' + '\t'.join(cov_data[sample_num,:]) + '\n')

	# Close filehandle
	t3.close()
	return


def make_gene_pheno_file(gene_pheno_file_name, family_ids, within_family_ids, gene_expression_vector):
	# Open file handle for gene pheno file
	t3 = open(gene_pheno_file_name,'w')

	# QUick erorr checking
	if len(family_ids) != len(gene_expression_vector):
		print('assumption eroror')
		pdb.set_trace()

	# Get number of smaples we have measurement for this gene
	num_samples = len(gene_expression_vector)

	# Loop through samples and print sample measurement to pheno file
	for sample_num in range(num_samples):
		t3.write(family_ids[sample_num] + '\t' + within_family_ids[sample_num] + '\t' + gene_expression_vector[sample_num] + '\n')
	# CLose filehandle
	t3.close()
	return

tissue_name = sys.argv[1]
xt_pc_gene_list_file = sys.argv[2]
gtex_expression_dir = sys.argv[3]
gtex_covariate_dir = sys.argv[4]
gtex_fusion_weights_data_dir = sys.argv[5]
subsample_size = int(sys.argv[6])
original_tissue_name = sys.argv[7]


# Create mapping from gene to (chrom_num, tss)
gene_to_chrom_num_and_tss = create_gene_to_chrom_num_and_tss_mapping(xt_pc_gene_list_file)

# Covariate files
covariate_file = gtex_covariate_dir + original_tissue_name + '_covariates.txt'

# Get sample names from covariate file
full_sample_names = get_sample_names_from_covariate_file(covariate_file)

# Now subsample sample names
sample_names_tmp = np.sort(np.random.choice(full_sample_names, size=subsample_size,replace=False))
sample_names = []
for sample_name in full_sample_names:
	if sample_name in sample_names_tmp:
		sample_names.append(sample_name)
sample_names = np.asarray(sample_names)


# Convert from sample names to family_id and within_family_id
family_ids, within_family_ids = convert_from_sample_name_to_family_and_within_family_ids(sample_names)

# Save subsampled individual ids
subsampled_indi_file = gtex_fusion_weights_data_dir + tissue_name + '_subsampled_individuals.txt'
t = open(subsampled_indi_file,'w')
for sample_name in sample_names:
	t.write('0' + '\t' + sample_name + '\n')
t.close()

# Create covariate file for plink
study_covariate_file = gtex_fusion_weights_data_dir + tissue_name + '_cov'
make_covariate_file(study_covariate_file, family_ids, within_family_ids, covariate_file)


# Open file handle for plink summary file
plink_gene_summary_file = gtex_fusion_weights_data_dir + tissue_name + '_gene_summary.txt'
t = open(plink_gene_summary_file,'w')
# Header
t.write('Gene_id\tchrom_num\tTSS\tgene_pheno_file\tcovariate_file\n')



for chrom_num in range(1,23):
	print(chrom_num)
	chrom_expression_file = gtex_expression_dir + original_tissue_name + '_normalized_expression_matrix_eqtl_ready_chr' + str(chrom_num) + '.txt'
	f = open(chrom_expression_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1

			valid_subsampled_indices = []
			for ii, indi_id in enumerate(data[1:]):
				if indi_id in sample_names:
					valid_subsampled_indices.append(ii)
			valid_subsampled_indices = np.asarray(valid_subsampled_indices)

			# Error checking
			if np.array_equal(np.asarray(data[1:])[valid_subsampled_indices], within_family_ids) == False:
				print('assumption eroro')
				pdb.set_trace()
			continue
		gene_name = data[0]
		# Ingnore non-protein coding genes
		if gene_name not in gene_to_chrom_num_and_tss:
			continue
		raw_gene_expression = (np.asarray(data[1:])[valid_subsampled_indices]).astype(float)
		gene_expression = (raw_gene_expression - np.mean(raw_gene_expression))/np.std(raw_gene_expression)
		gene_expression = gene_expression.astype(str)
		# Chrom num and TSS of this gene
		string_chrom_num, gene_tss = gene_to_chrom_num_and_tss[gene_name]

		# Create Pheno file for this gene
		gene_pheno_file_name = gtex_fusion_weights_data_dir + tissue_name + '_' + gene_name + '_pheno'
		make_gene_pheno_file(gene_pheno_file_name, family_ids, within_family_ids, gene_expression)

		# Update plink gene summary file with new gene
		t.write(gene_name + '\t' + str(chrom_num) + '\t' + str(gene_tss) + '\t' + gene_pheno_file_name + '\t' + study_covariate_file + '\n')
	f.close()
t.close()
