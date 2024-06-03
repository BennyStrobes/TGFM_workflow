import numpy as np
import os
import sys
import pdb
import scipy.stats







tgfm_iterative_prior_file = sys.argv[1]
two_step_fine_mapping_tissues_file = sys.argv[2]
gene_trait_architecture = sys.argv[3]

pvalues = []
tissue_names = []
f = open(tgfm_iterative_prior_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	if data[0] == 'variant':
		continue

	# Get tissue name
	tissue_name = data[0]
	# Get bootstrapping probabilities
	bootstrapped_coefficients = np.asarray(data[3].split(';')).astype(float)

	# Get pvalue
	if np.std(bootstrapped_coefficients, ddof=1) == 0.0:
		pvalue = 1.0
	else:
		zz = np.mean(bootstrapped_coefficients)/np.std(bootstrapped_coefficients, ddof=1)
		pvalue = scipy.stats.norm.sf(abs(zz))  # 1-sided pvalue

	pvalues.append(pvalue)
	tissue_names.append(tissue_name)
f.close()

# Put into organized vector
pvalues = np.asarray(pvalues)
tissue_names = np.asarray(tissue_names)

ordered_indices = np.argsort(pvalues)
ordered_pvalues = pvalues[ordered_indices]
ordered_tissue_names = tissue_names[ordered_indices]


# Print to output file
t = open(two_step_fine_mapping_tissues_file,'w')
t.write('tissue_name\tpvalue\n')


for ii, tissue_name in enumerate(ordered_tissue_names):
	pvalue = ordered_pvalues[ii]
	# Always print if best
	if ii == 0:
		t.write(tissue_name + '\t' + str(pvalue) + '\n')
		best_tissue_name = tissue_name


if gene_trait_architecture == '1_caus_t':
	true_causal_tissues = np.asarray([tissue_names[0]])
	true_causal_tissue_pvalues = np.asarray([pvalues[0]])
elif gene_trait_architecture == '2_caus_t':
	true_causal_tissues = np.asarray([tissue_names[0], tissue_names[3]])
	true_causal_tissue_pvalues = np.asarray([pvalues[0], pvalues[3]])


for ii, tissue_name in enumerate(true_causal_tissues):
	pvalue = true_causal_tissue_pvalues[ii]
	if tissue_name == best_tissue_name:
		continue
	t.write(tissue_name + '\t' + str(pvalue) + '\n')
t.close()

print(pvalues)
print(two_step_fine_mapping_tissues_file)
