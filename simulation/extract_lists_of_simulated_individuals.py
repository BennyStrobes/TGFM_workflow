import numpy as np 
import os
import sys
import pdb



def print_individual_list_to_output(individual_list, start_index, end_index, output_file):
	new_list = individual_list[start_index:end_index]
	t = open(output_file,'w')
	for ele in new_list:
		t.write(ele + '\n')
	t.close()





plink_file_stem = sys.argv[1]
n_gwas_individuals = int(sys.argv[2])
gwas_individual_file = sys.argv[3]
eqtl_individual_stem = sys.argv[4]
ref_genotype_individual_file = sys.argv[5]

eqtl_sample_sizes = np.asarray([100, 200, 300, 500, 1000])


plink_fam_file = plink_file_stem + '.fam'

individual_list = []
f = open(plink_fam_file)
for line in f:
	line = line.rstrip()
	data = line.split()
	indi_id = data[0] + '\t' + data[1]
	individual_list.append(indi_id)
f.close()
individual_list = np.asarray(individual_list)



print_individual_list_to_output(individual_list, 0, n_gwas_individuals, gwas_individual_file)

start_index = n_gwas_individuals


for eqtl_sample_size in eqtl_sample_sizes:
	print_individual_list_to_output(individual_list, start_index, (start_index+eqtl_sample_size), eqtl_individual_stem + str(eqtl_sample_size) + '.txt')
	start_index = start_index + eqtl_sample_size


ref_genotype_ss=500
print_individual_list_to_output(individual_list, start_index, (start_index+ref_genotype_ss), ref_genotype_individual_file)
start_index = start_index + ref_genotype_ss

