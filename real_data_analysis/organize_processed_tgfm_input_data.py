import numpy as np 
import os
import sys
import pdb







num_jobs = sys.argv[1]
gene_type = sys.argv[2]
preprocessed_tgfm_data_dir = sys.argv[3]


agg_output_file = preprocessed_tgfm_data_dir + gene_type + '_tgfm_input_data_summary.txt'
t = open(agg_output_file,'w')

for job_number in range(int(num_jobs)):
	file_name = preprocessed_tgfm_data_dir + gene_type + '_tgfm_input_data_summary_' + str(job_number) + '_' + num_jobs + '.txt'
	head_count = 0
	f = open(file_name)
	for line in f:
		line = line.rstrip()
		if head_count == 0:
			head_count = head_count + 1
			if job_number == 0:
				t.write(line + '\n')
			continue
		t.write(line + '\n')
	f.close()
t.close()

print(agg_output_file)