import numpy as np 
import os
import sys
import pdb






gtex_preprocessed_for_susie_dir = sys.argv[1]


output_file = gtex_preprocessed_for_susie_dir + 'susie_input_gene_organization_file.txt'
t = open(output_file,'w')


for chrom_num in range(1,23):
	input_file = gtex_preprocessed_for_susie_dir + 'chr' + str(chrom_num) + '_susie_input_gene_organization_file.txt'
	head_count = 0
	f = open(input_file)

	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			if chrom_num == 1:  # Only print header first go around
				t.write(line + '\n')
			continue
		t.write(line + '\n')
	f.close()
t.close()
