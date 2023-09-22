import numpy as np 
import os
import sys
import pdb






ukbb_preprocessed_for_genome_wide_susie_dir = sys.argv[1]

output_file = ukbb_preprocessed_for_genome_wide_susie_dir + 'genome_wide_susie_windows_and_processed_data.txt'
t = open(output_file,'w')


for chrom_num in range(1,23):
	input_file = ukbb_preprocessed_for_genome_wide_susie_dir + 'genome_wide_susie_windows_and_processed_data_chrom_' + str(chrom_num) + '.txt'

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
