import numpy as np 
import os
import sys
import pdb




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



input_file = ukbb_preprocessed_for_genome_wide_susie_dir + 'genome_wide_susie_windows_and_processed_data_chrom_' + chrom_num + '.txt'
f = open(input_file)
output_file = ukbb_preprocessed_for_genome_wide_susie_dir + 'genome_wide_susie_windows_and_processed_data_chrom_fixed_ld' + chrom_num + '.txt'
t = open(output_file,'w')

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	print('_'.join(data[:3]))
	old_in_sample_ld_file = data[10]
	ld_mat = np.load(old_in_sample_ld_file)
	fixed_ld_mat = correct_ld_mat_for_af_standardization(ld_mat)
	new_in_sample_ld_file = old_in_sample_ld_file.split('.np')[0] + '_af_corrected.npy'
	np.save(new_in_sample_ld_file, fixed_ld_mat)

	string_to_print = '\t'.join(data[:10]) + '\t' + new_in_sample_ld_file + '\t' + data[11] + '\n'
	t.write(string_to_print)

f.close()
t.close()