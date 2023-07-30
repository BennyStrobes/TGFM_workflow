import numpy as np 
import os
import sys
import pdb


def filter_pos_file_to_single_tissue(all_tissue_gene_model_pos_file, tissue_specific_pos_file, tissue_name_string):
	f = open(all_tissue_gene_model_pos_file)
	t = open(tissue_specific_pos_file,'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		if data[1].endswith(tissue_name_string):
			t.write(line + '\n')
	f.close()
	t.close()
	return








########################
# Command line args
########################
all_tissue_gene_model_pos_file = sys.argv[1]
tissue_specific_pos_file_stem = sys.argv[2]

for tissue_iter in range(10):
	tissue_specific_pos_file = tissue_specific_pos_file_stem + 'tissue' + str(tissue_iter) + '_gene_models_python_summary.pos'
	filter_pos_file_to_single_tissue(all_tissue_gene_model_pos_file, tissue_specific_pos_file, 'tissue' + str(tissue_iter))
