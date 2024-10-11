import numpy as np
import os
import sys
import pdb



def extract_trait_names(trait_names_file):
	f = open(trait_names_file)
	trait_names = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		trait_names.append(data[0])
	f.close()
	return np.asarray(trait_names)




####################
# Command line args
####################
tgfm_results_dir = sys.argv[1]
gene_type = sys.argv[2]
num_jobs = int(sys.argv[3])
trait_names_file = sys.argv[4]
tgfm_organized_results_dir = sys.argv[5]


# Extract trait names
trait_names = extract_trait_names(trait_names_file)

# Open output file handle
t = open(tgfm_organized_results_dir + 'component_matching_organized_res_16_independent_traits.txt','w')
# Print header
t.write('trait_name\twindow_name\tcomponent_num\tfraction_of_correlated_samples\tgene_fraction\n')


# Loop through traits
for trait_name in trait_names:
	# Loop through parallele jobs
	for job_num in range(num_jobs):
		res_file = tgfm_results_dir + 'tgfm_results_' + trait_name + '_' + gene_type + '_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_' + str(job_num) + '_8_tgfm_component_matching_summary.txt'
		f = open(res_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			t.write(trait_name + '\t' + data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\n')
		f.close()


# Close output file handle
t.close()
