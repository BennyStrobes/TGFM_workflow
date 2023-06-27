import numpy as np 
import os
import sys
import pdb




filtered_sample_info_file = sys.argv[1]
plink_format_filter_sample_info_file = sys.argv[2]


f = open(filtered_sample_info_file)
t = open(plink_format_filter_sample_info_file,'w')
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count ==0:
		head_count = head_count + 1
		continue
	t.write(data[1] + '\n')

f.close()
t.close()