import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import pdb
import gzip
import pandas as pd



def extract_gwas_summary_statistics(sumstat_file):
	dicti = {}
	f = open(sumstat_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent fiels
		variant_id = data[0]
		a1 = data[1]
		a2 = data[2]
		variant_info = variant_id.split('_')
		# Quick error check
		if variant_info[2] != a1:
			print('assumption eroror')
			pdb.set_trace()
		z_score = float(data[3])
		beta = float(data[4])
		beta_se = float(data[5])
		# Two versions of variant
		variant_id_ref = variant_info[0] + '_' + variant_info[1] + '_' + variant_info[2] + '_' + variant_info[3]
		variant_id_alt = variant_info[0] + '_' + variant_info[1] + '_' + variant_info[3] + '_' + variant_info[2]

		if variant_id_ref in dicti or variant_id_alt in dicti:
			print('assumption eroror')
			pdb.set_trace()

		dicti[variant_id_ref] = (z_score, beta, beta_se)
		dicti[variant_id_alt] = (-z_score, -beta, beta_se)
	f.close()
	return dicti

def extract_gwas_data_for_vector_of_variants(variants, gwas_sumstat_data):
	z_scores = []
	betas = []
	betas_se = []
	for variant in variants:
		variant_tuple = gwas_sumstat_data[variant]
		z_scores.append(variant_tuple[0])
		betas.append(variant_tuple[1])
		betas_se.append(variant_tuple[2])
	return np.asarray(z_scores), np.asarray(betas), np.asarray(betas_se)

def make_sample_size_df(full_studies, gtex_variants, gtex_sample_sizes, trait_sample_size):
	# Initialize matrix
	n_mat = np.zeros((len(full_studies), len(gtex_variants)))
	# Add gtex rows
	for row_num, gtex_sample_size in enumerate(gtex_sample_sizes):
		n_mat[row_num,:] = n_mat[row_num,:] + float(gtex_sample_size)
	# Add gwas row
	n_mat[-1,:] = n_mat[-1,:] + trait_sample_size
	
	# Put in pandas df
	n_df = pd.DataFrame(n_mat.astype(int), index=full_studies, columns=gtex_variants)
	return n_df


gtex_susie_input_data_file = sys.argv[1]
trait_name = sys.argv[2]
sumstat_file = sys.argv[3]
trait_sample_size = int(sys.argv[4])
coloc_data_dir = sys.argv[5]  # output file


# Extract summary statistics from gwas into dictionary
gwas_sumstat_data = extract_gwas_summary_statistics(sumstat_file)


# Create output file handle
output_file = coloc_data_dir + trait_name + '_processed_gene_list.txt'
t = open(output_file,'w')
t.write('gene_id\tchrom_num\tgenotype_file\tz_file\tn_file\tbeta_file\tstd_err_file\n')

f = open(gtex_susie_input_data_file)
head_count = 0
counter = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	# Skip header
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Extract relevant fields from line
	counter = counter + 1
	gene_id = data[0]
	chrom_num = data[1]
	gtex_beta_file = data[2]
	gtex_beta_std_err_file = data[3]
	gtex_variant_file = data[4]
	gtex_tissue_file = data[5]
	gtex_ref_genotype_file = data[6]
	gtex_sample_size_file = data[7]
	
	# Extract gtex variants
	gtex_variants = np.loadtxt(gtex_variant_file, dtype=str, delimiter='\t')

	# Ignore genes with only 1 variant
	if len(gtex_variants.shape) == 0:
		continue

	# Extract gwas z-scores, gwas betas, gwas beta_se
	gwas_z_scores, gwas_betas, gwas_beta_se = extract_gwas_data_for_vector_of_variants(gtex_variants, gwas_sumstat_data)

	# Extract gtex beta and gtex beta_se
	gtex_beta = np.loadtxt(gtex_beta_file)
	gtex_beta_se = np.loadtxt(gtex_beta_std_err_file)
	gtex_z = gtex_beta/gtex_beta_se

	# Extract gtex study names and sample sizes
	gtex_studies = np.loadtxt(gtex_tissue_file, dtype=str, delimiter='\t')
	gtex_sample_sizes = np.loadtxt(gtex_sample_size_file, dtype=str, delimiter='\t')

	# Only one tissue present, need to fix
	if len(gtex_z.shape) == 1:
		gtex_beta = np.reshape(gtex_beta, (1, len(gtex_beta)))
		gtex_beta_se = np.reshape(gtex_beta_se, (1, len(gtex_beta_se)))
		gtex_z = np.reshape(gtex_z, (1, len(gtex_z)))
		gtex_studies = np.reshape(gtex_studies,1)
		gtex_sample_sizes = np.reshape(gtex_sample_sizes, 1)

	# Quick error checking
	if gtex_beta.shape[1] != len(gwas_z_scores):
		print('assumption eroorro')
		pdb.set_trace()

	# Append studie names together
	full_studies = np.hstack((gtex_studies, [trait_name]))

	# Put everything in pandas dfs
	z_df = pd.DataFrame(np.vstack((gtex_z, gwas_z_scores)), index=full_studies, columns=gtex_variants)
	beta_df = pd.DataFrame(np.vstack((gtex_beta, gwas_betas)), index=full_studies, columns=gtex_variants)
	std_err_df = pd.DataFrame(np.vstack((gtex_beta_se, gwas_beta_se)), index=full_studies, columns=gtex_variants)
	
	n_df = make_sample_size_df(full_studies, gtex_variants, gtex_sample_sizes, trait_sample_size)

	# Save summary stats to pickles
	# Save z_df to pkl
	z_pkl_output_file = coloc_data_dir + trait_name + '_' + gene_id + '_z_df.txt'
	z_df.to_csv(z_pkl_output_file, sep="\t")
	# Save beta_df to pkl
	beta_pkl_output_file = coloc_data_dir + trait_name + '_' + gene_id + '_beta_df.txt'
	beta_df.to_csv(beta_pkl_output_file, sep="\t")
	# Save std_err_df to pkl
	std_err_pkl_output_file = coloc_data_dir + trait_name + '_' + gene_id + '_std_err_df.txt'
	std_err_df.to_csv(std_err_pkl_output_file, sep="\t")
	# Save n_df to pkl
	n_pkl_output_file = coloc_data_dir + trait_name + '_' + gene_id + '_n_df.txt'
	n_df.to_csv(n_pkl_output_file, sep="\t")

	# Save locations of pickled file
	t.write(gene_id + '\t' + str(chrom_num) + '\t' + gtex_ref_genotype_file + '\t' + z_pkl_output_file + '\t' + n_pkl_output_file + '\t' + beta_pkl_output_file + '\t' + std_err_pkl_output_file + '\n')

f.close()
t.close()

