import numpy as np 
import os
import sys
import pdb





sldsc_input_file = sys.argv[1] + '.l2.ldscore'
nonnegative_coef_output_file = sys.argv[2]


boolean_arr = []
head_count = 0
f = open(sldsc_input_file)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		for ele in data[3:]:
			if ele.startswith('tissue'):
				boolean_arr.append('True')
			else:
				boolean_arr.append('False')
		continue
	break
f.close()

t = open(nonnegative_coef_output_file,'w')

for ele in boolean_arr:
	t.write(ele + '\n')
t.close()