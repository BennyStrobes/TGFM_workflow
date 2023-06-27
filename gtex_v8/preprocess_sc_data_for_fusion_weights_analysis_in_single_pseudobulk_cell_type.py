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
		tss = data[2]
		if gene_name in dicti:
			print('assumption eroor')
			pdb.set_trace()
		dicti[gene_name] = (chrom_num, tss)
	f.close()
	return dicti



def get_sample_names_from_covariate_file(covariate_file):
	cov = np.loadtxt(covariate_file, dtype=str,delimiter='\t')
	return cov[1:,0]


def convert_from_sample_name_to_family_and_within_family_ids(sample_names):
	family_ids = []
	within_family_ids = []
	for sample_name in sample_names:
		family_ids.append('0')
		within_family_ids.append(sample_name)
	return np.asarray(family_ids), np.asarray(within_family_ids)

def make_covariate_file(study_covariate_file, family_ids, within_family_ids, pb_covariate_file):
	# Open file handle for gene pheno file
	t3 = open(study_covariate_file,'w')

	# Load in covariates
	cov_data = np.loadtxt(pb_covariate_file, dtype=str, delimiter='\t')
	cov_data = cov_data[1:, 1:]

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


######################
# Command line args
######################
pb_cell_type_name = sys.argv[1]
sc_pseudobulk_expression_dir = sys.argv[2]
output_dir = sys.argv[3]
num_jobs = int(sys.argv[4])


# Create mapping from gene to tss
gene_tss_file = sc_pseudobulk_expression_dir + 'cell_types_' + pb_cell_type_name + '_gene_tss_hg38.txt'
gene_to_gene_pos_mapping = create_gene_to_chrom_num_and_tss_mapping(gene_tss_file)

# Covariate files
covariate_file = sc_pseudobulk_expression_dir + 'cell_types_' + pb_cell_type_name + '_pseudobulk_covariates.txt'

# Get sample names from covariate file
sample_names = get_sample_names_from_covariate_file(covariate_file)


# Convert from sample names to family_id and within_family_id
family_ids, within_family_ids = convert_from_sample_name_to_family_and_within_family_ids(sample_names)



# Create covariate file for plink
study_covariate_file = output_dir + 'cell_types_' + pb_cell_type_name + '_fusion_ready_cov'
make_covariate_file(study_covariate_file, family_ids, within_family_ids, covariate_file)


# Open file handle for plink summary file
plink_gene_summary_file = output_dir + 'cell_types_' + pb_cell_type_name + '_fusion_ready_gene_summary.txt'
t = open(plink_gene_summary_file,'w')
# Header
t.write('Gene_id\tchrom_num\tTSS\tgene_pheno_file\tcovariate_file\tcomposit_tissues\n')


for chrom_num in range(1,23):
	print(chrom_num)
	chrom_expression_file = sc_pseudobulk_expression_dir + 'cell_types_' + pb_cell_type_name + '_pseudobulk_expression.txt'
	f = open(chrom_expression_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			# Error checking
			if np.array_equal(np.asarray(data[1:]), within_family_ids) == False:
				print('assumption eroro')
				pdb.set_trace()
			continue
		gene_name = data[0]
		# Ingnore non-protein coding genes
		if gene_name not in gene_to_gene_pos_mapping:
			continue
		gene_expression = data[1:]
		# Chrom num and TSS of this gene
		string_chrom_num, gene_tss = gene_to_gene_pos_mapping[gene_name]

		if string_chrom_num != 'chr' + str(chrom_num):
			continue

		# Create Pheno file for this gene
		gene_pheno_file_name = output_dir + 'cell_types_' + pb_cell_type_name + '_fusion_ready_' + gene_name + '_pheno'
		make_gene_pheno_file(gene_pheno_file_name, family_ids, within_family_ids, gene_expression)

		# Update plink gene summary file with new gene
		t.write(gene_name + '\t' + str(chrom_num) + '\t' + str(gene_tss) + '\t' + gene_pheno_file_name + '\t' + study_covariate_file + '\t' + pb_cell_type_name + '\n')
	f.close()
t.close()


# Now parallelize the jobs
# Parallelize
data_raw = np.loadtxt(plink_gene_summary_file, dtype=str,delimiter='\t')[1:,:]
nrows = data_raw.shape[0]
splits = np.array_split(np.arange(nrows),num_jobs)

for job_num, split_mat in enumerate(splits):
	split_output_file = output_dir + 'cell_types_' + pb_cell_type_name + '_fusion_ready_gene_summary_' + str(job_num) + '.txt'
	valid_rows = {}
	for row_num in split_mat:
		valid_rows[row_num] = 1
	t = open(split_output_file,'w')
	t.write('Gene_id\tchrom_num\tTSS\tgene_pheno_file\tcovariate_file\tcomposit_tissues\n')
	f = open(plink_gene_summary_file)
	head_count = 0
	line_num = 0
	for line in f:
		line = line.rstrip()
		if head_count == 0:
			head_count = head_count + 1
			continue
		if line_num in valid_rows:
			t.write(line + '\n')
		line_num = line_num + 1
	f.close()
	t.close()

