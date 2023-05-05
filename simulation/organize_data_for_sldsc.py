import numpy as np 
import os
import sys
import pdb
import gzip




def reformat_gwas_summary_statistics_for_sldsc(genotype_bim_file, standard_gwas_results_file, ldsc_gwas_results_file, n_gwas_individuals):
	# First create mapping from rsid to alleles
	rsid_to_alleles = {}
	f = open(genotype_bim_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		rs_id = data[1]
		alleles = '\t'.join(data[4:])
		if rs_id in rsid_to_alleles:
			print('assumptino eroror')
			pdb.set_trace()
		rsid_to_alleles[rs_id] = alleles
	f.close()


	# Stream standard gwas results and print to sldsc gwas results
	f = open(standard_gwas_results_file)
	t = open(ldsc_gwas_results_file,'w')
	head_count = 0
	# Write new header
	t.write('SNP\tA1\tA2\tN\tCHISQ\tZ\n')

	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rs_id = data[0]
		z_score = float(data[3])
		alleles = rsid_to_alleles[rs_id]
		t.write(rs_id + '\t' + alleles + '\t' + n_gwas_individuals + '\t' + str(np.square(z_score)) + '\t' + str(z_score) + '\n')
	f.close()
	t.close()

	return

def load_in_string_matrix(file_name, delimiter='\t'):
	arr = []
	f = open(file_name)
	for line in f:
		line = line.rstrip()
		data = line.split(delimiter)
		arr.append(data)
	f.close()
	return np.asarray(arr)

def create_joint_ld_score_file(variant_ld_score_file, gene_ld_score_file, joint_ld_score_file):
	# Load in gene-ld scores into memory
	gene_ld_scores_raw = load_in_string_matrix(gene_ld_score_file)

	n_tiss = gene_ld_scores_raw.shape[1] - 1
	# create names of annotation defining eqtls
	tissue_anno_names = []
	for tiss_iter in range(n_tiss):
		tissue_anno_names.append('tissue' + str(tiss_iter))
	tissue_anno_names = np.asarray(tissue_anno_names)

	# Open input and output files
	f = gzip.open(variant_ld_score_file)
	t = open(joint_ld_score_file,'w')

	# Stream input file
	head_count = 0
	variant_counter = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		# QUick error check
		if len(data) != 56:
			print('assumption eroror')
			pdb.set_trace()

		# Print header
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\t' + '\t'.join(tissue_anno_names) + '\n')
			continue

		# Standard line

		# Quick error check
		if gene_ld_scores_raw[variant_counter, 0] != data[1]:
			print('assumption eroror')
			pdb.set_trace()

		# Print joint line
		variant_geno_ld_scores = gene_ld_scores_raw[variant_counter, 1:]
		t.write(line + '\t' + '\t'.join(variant_geno_ld_scores) + '\n')
		variant_counter = variant_counter + 1

	f.close()
	t.close()

	return


def create_joint_M_file(variant_m_file, gene_m_file, joint_m_file):
	variant_m = np.loadtxt(variant_m_file)
	gene_m = np.loadtxt(gene_m_file)
	t = open(joint_m_file,'w')
	t.write('\t'.join(variant_m.astype(int).astype(str)) + '\t' + '\t'.join(gene_m.astype(int).astype(str)) + '\n')
	t.close()

def reformat_variant_weights_for_ldsc(regression_snps_file, existing_variant_weights_file, new_variant_weights_file):
	# Create mapping from rsid to weight string
	f = gzip.open(existing_variant_weights_file)
	head_count = 0
	dicti = {}
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			weight_header = line
			continue
		rsid = data[1]
		if rsid in dicti:
			print('assumtpion eroror')
			pdb.set_trace()
		dicti[rsid] = line
	f.close()

	f = open(regression_snps_file)
	t = open(new_variant_weights_file,'w')
	t.write(weight_header + '\n')
	for line in f:
		rsid = line.rstrip()
		t.write(dicti[rsid] + '\n')
	f.close()
	t.close()


	return



############################
# Command line args
############################
simulation_number = sys.argv[1]
chrom_num = sys.argv[2]
simulation_name_string = sys.argv[3]
n_gwas_individuals = sys.argv[4]
processed_genotype_data_dir = sys.argv[5]
simulated_gwas_dir = sys.argv[6]
ldsc_weights_dir = sys.argv[7]
simulated_ld_scores_dir = sys.argv[8]


############################
#  Reformat gwas summary statistics for sldsc
############################
genotype_bim_file = processed_genotype_data_dir + 'simulated_gwas_data_' + chrom_num + '.bim'  # Used to create mapping from rsid to A1, A2
standard_gwas_results_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results_hm3_noHMC_snps_only.txt'
ldsc_gwas_results_file = simulated_ld_scores_dir + simulation_name_string + '_ldsc_ready_summary_statistics.txt'
reformat_gwas_summary_statistics_for_sldsc(genotype_bim_file, standard_gwas_results_file, ldsc_gwas_results_file, n_gwas_individuals)


############################
#  Generate Joint-variant-gene LD-score files for each of the eqtl sample sizes
# Using eqtl pmces
############################
variant_ld_score_file = simulated_ld_scores_dir + simulation_name_string + '_baseline.' + chrom_num + '.l2.ldscore.gz'  # Variant level ld-score (doesn't change as we vary expression data)
# Loop through eqtl sample sizes
eqtl_sample_sizes = np.asarray([300,500,1000, 'inf'])  # Various eqtl data sets

for eqtl_sample_size in eqtl_sample_sizes:
	# File containing gene ld scores for this data set
	gene_ld_score_file = simulated_ld_scores_dir + simulation_name_string + '_gene_weighted_ld_scores_eqtlss_' + str(eqtl_sample_size) + '_ld_scores.txt'
	# Output file
	joint_ld_score_file = simulated_ld_scores_dir + simulation_name_string + '_joint_baseline_variant_' + str(eqtl_sample_size) + '_gene_ld_scores.l2.ldscore'
	create_joint_ld_score_file(variant_ld_score_file, gene_ld_score_file, joint_ld_score_file)


############################
#  Generate Joint-variant-gene LD-score files for each of the eqtl sample sizes
# Using susie distribution eqtls
############################
variant_ld_score_file = simulated_ld_scores_dir + simulation_name_string + '_baseline.' + chrom_num + '.l2.ldscore.gz'  # Variant level ld-score (doesn't change as we vary expression data)
# Loop through eqtl sample sizes
eqtl_sample_sizes = np.asarray([300,500,1000])  # Various eqtl data sets

for eqtl_sample_size in eqtl_sample_sizes:
	# File containing gene ld scores for this data set
	gene_ld_score_file = simulated_ld_scores_dir + simulation_name_string + '_susie_distr_gene_weighted_ld_scores_eqtlss_' + str(eqtl_sample_size) + '_ld_scores.txt'
	# Output file
	joint_ld_score_file = simulated_ld_scores_dir + simulation_name_string + '_joint_baseline_variant_' + str(eqtl_sample_size) + '_susie_distr_gene_ld_scores.l2.ldscore'
	create_joint_ld_score_file(variant_ld_score_file, gene_ld_score_file, joint_ld_score_file)



############################
#  Generate Joint-variant-gene M files for each of the eqtl sample sizes
# Using eqtl pmces
############################
variant_m_file = simulated_ld_scores_dir + simulation_name_string + '_baseline.' + chrom_num + '.l2.M'  # Variant level m file (doesn't change as we vary expression data)
variant_m_5_50_file = simulated_ld_scores_dir + simulation_name_string + '_baseline.' + chrom_num + '.l2.M_5_50'  # Variant level m file (doesn't change as we vary expression data)

# Loop through eqtl sample sizes
eqtl_sample_sizes = np.asarray([300,500,1000, 'inf'])  # Various eqtl data sets
for eqtl_sample_size in eqtl_sample_sizes:
	# File containing gene ld scores for this data set
	gene_m_file = simulated_ld_scores_dir + simulation_name_string + '_gene_weighted_ld_scores_eqtlss_' + str(eqtl_sample_size) + '_M.txt'
	# Output file
	joint_m_5_50_file = simulated_ld_scores_dir + simulation_name_string + '_joint_baseline_variant_' + str(eqtl_sample_size) + '_gene_ld_scores.l2.M_5_50'
	create_joint_M_file(variant_m_5_50_file, gene_m_file, joint_m_5_50_file)
	# Output file
	joint_m_file = simulated_ld_scores_dir + simulation_name_string + '_joint_baseline_variant_' + str(eqtl_sample_size) + '_gene_ld_scores.l2.M'
	create_joint_M_file(variant_m_file, gene_m_file, joint_m_file)


############################
#  Generate Joint-variant-gene M files for each of the eqtl sample sizes
# Using susie distribution eqtls
############################
variant_m_file = simulated_ld_scores_dir + simulation_name_string + '_baseline.' + chrom_num + '.l2.M'  # Variant level m file (doesn't change as we vary expression data)
variant_m_5_50_file = simulated_ld_scores_dir + simulation_name_string + '_baseline.' + chrom_num + '.l2.M_5_50'  # Variant level m file (doesn't change as we vary expression data)

# Loop through eqtl sample sizes
eqtl_sample_sizes = np.asarray([300,500,1000])  # Various eqtl data sets
for eqtl_sample_size in eqtl_sample_sizes:
	# File containing gene ld scores for this data set
	gene_m_file = simulated_ld_scores_dir + simulation_name_string + '_susie_distr_gene_weighted_ld_scores_eqtlss_' + str(eqtl_sample_size) + '_M.txt'
	# Output file
	joint_m_5_50_file = simulated_ld_scores_dir + simulation_name_string + '_joint_baseline_variant_' + str(eqtl_sample_size) + '_susie_distr_gene_ld_scores.l2.M_5_50'
	create_joint_M_file(variant_m_5_50_file, gene_m_file, joint_m_5_50_file)
	# Output file
	joint_m_file = simulated_ld_scores_dir + simulation_name_string + '_joint_baseline_variant_' + str(eqtl_sample_size) + '_susie_distr_gene_ld_scores.l2.M'
	create_joint_M_file(variant_m_file, gene_m_file, joint_m_file)





'''
############################
#  Generate Joint-variant-gene LD-score files for each of the eqtl sample sizes
# Using unbiased marginal
############################
variant_ld_score_file = simulated_ld_scores_dir + simulation_name_string + '_baseline.' + chrom_num + '.l2.ldscore.gz'  # Variant level ld-score (doesn't change as we vary expression data)
# Loop through eqtl sample sizes
eqtl_sample_sizes = np.asarray([100,300,500, 1000])  # Various eqtl data sets

for eqtl_sample_size in eqtl_sample_sizes:
	# File containing gene ld scores for this data set
	gene_ld_score_file = simulated_ld_scores_dir + simulation_name_string + '_unbiased_marginal_gene_weighted_ld_scores_eqtlss_' + str(eqtl_sample_size) + '_ld_scores.txt'
	# Output file
	joint_ld_score_file = simulated_ld_scores_dir + simulation_name_string + '_joint_baseline_variant_' + str(eqtl_sample_size) + '_unbiased_marginal_gene_ld_scores.l2.ldscore'
	create_joint_ld_score_file(variant_ld_score_file, gene_ld_score_file, joint_ld_score_file)
'''



'''
############################
#  Generate Joint-variant-gene M files for each of the eqtl sample sizes
# Using unbiased marginal
############################
variant_m_file = simulated_ld_scores_dir + simulation_name_string + '_baseline.' + chrom_num + '.l2.M'  # Variant level m file (doesn't change as we vary expression data)
variant_m_5_50_file = simulated_ld_scores_dir + simulation_name_string + '_baseline.' + chrom_num + '.l2.M_5_50'  # Variant level m file (doesn't change as we vary expression data)

# Loop through eqtl sample sizes
eqtl_sample_sizes = np.asarray([100,300,500, 1000])  # Various eqtl data sets
for eqtl_sample_size in eqtl_sample_sizes:
	# File containing gene ld scores for this data set
	gene_m_file = simulated_ld_scores_dir + simulation_name_string + '_unbiased_marginal_gene_weighted_ld_scores_eqtlss_' + str(eqtl_sample_size) + '_M.txt'
	# Output file
	joint_m_5_50_file = simulated_ld_scores_dir + simulation_name_string + '_joint_baseline_variant_' + str(eqtl_sample_size) + '_unbiased_marginal_gene_ld_scores.l2.M_5_50'
	create_joint_M_file(variant_m_5_50_file, gene_m_file, joint_m_5_50_file)
	# Output file
	joint_m_file = simulated_ld_scores_dir + simulation_name_string + '_joint_baseline_variant_' + str(eqtl_sample_size) + '_unbiased_marginal_gene_ld_scores.l2.M'
	create_joint_M_file(variant_m_file, gene_m_file, joint_m_file)
'''











