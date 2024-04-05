import numpy as np 
import os
import sys
import pdb
import scipy.stats


def create_mapping_from_rsid_to_gwas_sum_stats(simulated_gwas_file):
	mapping = {}

	f = open(simulated_gwas_file)
	head_count = 0
	for line in f:
		line =line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[1]
		beta = float(data[6])
		se = np.sqrt(float(data[7]))
		z_score = beta/se
		pvalue = scipy.stats.norm.sf(abs(z_score))*2
		if rsid in mapping:
			print('assumption erororor')
			pdb.set_trace()
		mapping[rsid] = (beta, pvalue)
	f.close()
	return mapping

#######################################
# Command line args
#######################################
global_window_file = sys.argv[1]
simulation_number = sys.argv[2]
chrom_num = sys.argv[3]
simulation_name_string = sys.argv[4]
eqtl_sample_size = sys.argv[5]
simulated_gwas_dir = sys.argv[6] # input dir
processed_genotype_data_dir = sys.argv[7]
n_gwas_individuals = sys.argv[8]
focus_gwas_summary_stat_file = sys.argv[9]  # output file

# Extract window names
window_names = np.loadtxt(global_window_file, dtype=str, delimiter='\t')[1:,0]

# Create mapping from rsid to beta and pvalue
simulated_gwas_file = simulated_gwas_dir + simulation_name_string + '_merged_gwas_summary_stats.txt'
rsid_to_gwas_sum_stats = create_mapping_from_rsid_to_gwas_sum_stats(simulated_gwas_file)

# Genotype bim file
genotype_bim_file = processed_genotype_data_dir + 'simulated_gwas_data_' + str(chrom_num) + '.bim'
genotype_df = np.loadtxt(genotype_bim_file,dtype=str, delimiter='\t')

# Open output file handle
t = open(focus_gwas_summary_stat_file,'w')
t.write('CHR\tSNP\tBP\tA1\tA2\tN\tBETA\tP\n')

n_snps = genotype_df.shape[0]

for snp_iter in range(n_snps):
	rsid = genotype_df[snp_iter,1]
	line_chrom_num = genotype_df[snp_iter,0]
	snp_pos = genotype_df[snp_iter, 3]
	snp_a1 = genotype_df[snp_iter, 4]
	snp_a2 = genotype_df[snp_iter, 5]
	t.write(line_chrom_num + '\t' + rsid + '\t' + snp_pos + '\t' + snp_a1 + '\t' + snp_a2 + '\t' + n_gwas_individuals + '\t' + str(rsid_to_gwas_sum_stats[rsid][0]) + '\t' + str(rsid_to_gwas_sum_stats[rsid][1]) + '\n')

t.close()

