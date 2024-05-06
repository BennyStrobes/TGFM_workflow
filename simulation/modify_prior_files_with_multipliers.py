import numpy as np
import os
import sys
import pdb








#######################
# Command line args
#######################
tgfm_output_stem = sys.argv[1]
nm_var_prior_multiplier_str = sys.argv[2]
gene_tissue_prior_multiplier_str = sys.argv[3]

# Convert multipliers to floats
nm_var_prior_multiplier = float(nm_var_prior_multiplier_str)
gene_tissue_prior_multiplier = float(gene_tissue_prior_multiplier_str)


original_prior_file = tgfm_output_stem + '_iterative_variant_gene_prior_pip_level_bootstrapped.txt'
output_prior_file = tgfm_output_stem + '_iterative_variant_gene_prior_pip_level_bootstrapped_scaled_nm_' + str(nm_var_prior_multiplier_str) + '_gt_' + gene_tissue_prior_multiplier_str +'.txt'

causal_tissues = {}
causal_tissues['tissue0'] = 1
causal_tissues['tissue3'] = 1


# Open input file handle
f = open(original_prior_file)
t = open(output_prior_file,'w')

head_count = 0

for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	if data[0] == 'variant':
		new_mean = float(data[1])*nm_var_prior_multiplier
		new_bs_values = np.asarray(data[3].split(';')).astype(float)*nm_var_prior_multiplier
		new_exp_log_mean = np.exp(np.mean(np.log(new_bs_values)))
		t.write(data[0] + '\t' + str(new_mean) + '\t' + str(new_exp_log_mean) + '\t' + ';'.join(new_bs_values.astype(str)) + '\n')
	else:
		# Tissues
		if data[0] not in causal_tissues:
			new_mean = float(data[1])*gene_tissue_prior_multiplier
			new_bs_values = np.asarray(data[3].split(';')).astype(float)*gene_tissue_prior_multiplier
			new_exp_log_mean = np.exp(np.mean(np.log(new_bs_values)))
			t.write(data[0] + '\t' + str(new_mean) + '\t' + str(new_exp_log_mean) + '\t' + ';'.join(new_bs_values.astype(str)) + '\n')
		else:
			# Causal tissues (no change)
			t.write(line + '\n')

f.close()
t.close()


