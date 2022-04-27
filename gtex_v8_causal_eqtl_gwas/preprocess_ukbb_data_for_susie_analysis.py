import numpy as np 
import os
import sys
import pdb




def extract_ukbb_trait_names_and_file(ukbb_trait_file):
	f = open(ukbb_trait_file)
	head_count = 0
	trait_names = []
	trait_files = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		trait_names.append(data[0])
		trait_files.append(data[1])
	f.close()
	return np.asarray(trait_names), np.asarray(trait_files)


def fill_in_snp_dictionaries_with_gwas_summary_stats(trait_file, trait_num, num_traits, beta_dictionary, std_err_dictionary):
	f = open(trait_file)
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		counter = counter + 1
		if np.mod(counter, 10000) == 0:
			print(counter)
		chrom_num = 'chr' + data[1]
		pos = data[2]
		a1 = data[4]
		a2 = data[5]

		snp_id1 = chrom_num + '_' + pos + '_' + a1 + '_' + a2
		snp_id2 = chrom_num + '_' + pos + '_' + a2 + '_' + a1

		if snp_id1 not in beta_dictionary:
			beta_dictionary[snp_id1] = np.zeros(num_traits)
			std_err_dictionary[snp_id1] = np.zeros(num_traits)
		if snp_id2 not in beta_dictionary:
			beta_dictionary[snp_id2] = np.zeros(num_traits)
			std_err_dictionary[snp_id2] = np.zeros(num_traits)

		beta = float(data[10])
		std_err = float(data[11])


		beta_dictionary[snp_id1][trait_num] = beta
		beta_dictionary[snp_id2][trait_num] = -beta
		std_err_dictionary[snp_id1][trait_num] = std_err
		std_err_dictionary[snp_id2][trait_num] = std_err

	f.close()
	return beta_dictionary, std_err_dictionary




gtex_gene_file = sys.argv[1]
ukbb_sumstats_hg38_dir = sys.argv[2]
ukbb_preprocessed_for_susie_dir = sys.argv[3]


ukbb_trait_file = ukbb_sumstats_hg38_dir + 'ukbb_hg38_sumstat_files_small.txt'

trait_names, trait_files = extract_ukbb_trait_names_and_file(ukbb_trait_file)
num_traits = len(trait_files)
print(num_traits)

# Map snps to vectors of length number of traits 
beta_dictionary = {}
std_err_dictionary = {}

for trait_num, trait_name in enumerate(trait_names):
	print(trait_name)
	trait_file = trait_files[trait_num]

	# Fill in snp dictionary with gwas summary stats
	beta_dictionary, std_err_dictionary = fill_in_snp_dictionaries_with_gwas_summary_stats(trait_file, trait_num, num_traits, beta_dictionary, std_err_dictionary)

f = open(gtex_gene_file)
ukbb_gene_file = ukbb_preprocessed_for_susie_dir + 'susie_input_ukkbb_gene_organization_file.txt'
t = open(ukbb_gene_file,'w')
head_count = 0

for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	# Parse line
	test_gene = data[0]
	print(test_gene)
	chrom_num = data[1]
	variant_file = data[4]
	ref_genotype_file = data[6]
	# Missing: beta_file, std_err_file, tissue_file, sample_size_file
	gtex_variants = np.atleast_1d(np.loadtxt(variant_file,dtype=str,delimiter='\t'))


	# Create beta and standard error matrix
	beta_arr = []
	std_err_arr = []
	for gtex_variant in gtex_variants:
		beta_arr.append(beta_dictionary[gtex_variant])
		std_err_arr.append(std_err_dictionary[gtex_variant])
	beta_arr = np.transpose(np.asarray(beta_arr))
	std_err_arr = np.transpose(np.asarray(std_err_arr))

	# Save data to output file
	# Beta file
	beta_file = ukbb_preprocessed_for_susie_dir + test_gene + '_beta.txt'
	np.savetxt(beta_file, beta_arr, fmt="%s", delimiter='\t')
	# stderr file
	stderr_file = ukbb_preprocessed_for_susie_dir + test_gene + '_beta_std_err.txt'
	np.savetxt(stderr_file, std_err_arr, fmt="%s", delimiter='\t')
	# Tissue file
	tissue_file = ukbb_preprocessed_for_susie_dir + test_gene + '_studies.txt'
	np.savetxt(tissue_file, trait_names, fmt="%s", delimiter='\t')
	# Sample sizes
	sample_size_file = ukbb_preprocessed_for_susie_dir + test_gene + '_study_sample_sizes.txt'
	samp_size_vec = np.ones(len(trait_names))*381496
	np.savetxt(sample_size_file, samp_size_vec, fmt="%s", delimiter='\t')
	
	t.write(test_gene + '\t' + str(chrom_num) + '\t' + beta_file + '\t' + stderr_file + '\t' + variant_file + '\t' + tissue_file + '\t' + ref_genotype_file + '\t' + sample_size_file + '\n')

t.close()

