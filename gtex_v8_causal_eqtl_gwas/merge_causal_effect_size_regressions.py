import numpy as np 
import os
import sys
import pdb












trait_name = sys.argv[1]
diry = sys.argv[2]
total_jobs = int(sys.argv[3])




output_file = diry +'causal_effect_regression_merged.txt'
t = open(output_file,'w')

for job_number in range(total_jobs):
	input_file = diry + 'causal_effect_regression_' + str(job_number) + '_' + str(total_jobs) + '.txt'
	f = open(input_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 7:
			print('assumption eroror')
			pdb.set_trace()
		if head_count == 0:
			head_count = 1
			if job_number == 0:
				t.write(line + '\n')
			continue
		t.write(line + '\n')
	f.close()
t.close()