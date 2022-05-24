import numpy as np 
import os
import sys
import pdb







input_dir = sys.argv[1]
output_file = sys.argv[2]
chrom_num = sys.argv[3]


input_file = input_dir + 'blood_WHITE_COUNT_hg38_liftover.bgen.stats'


f = open(input_file)
t = open(output_file, 'w')

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	line_chrom_num = data[1]
	chrom_pos = data[2]
	a1 = data[4]
	a2 = data[5]
	if line_chrom_num != chrom_num:
		continue
	gtex_id1 = 'chr' + line_chrom_num + '_' + chrom_pos + '_' + a1 + '_' + a2 + '_b38'
	gtex_id2 = 'chr' + line_chrom_num + '_' + chrom_pos + '_' + a2 + '_' + a1 + '_b38'
	t.write(gtex_id1 + '\n')
	t.write(gtex_id2 + '\n')
f.close()
t.close()