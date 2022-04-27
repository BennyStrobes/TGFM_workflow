import numpy as np 
import os
import sys
import pdb








gene_summary_file = sys.argv[1]
gtex_fusion_weights_dir = sys.argv[2]
tissue_name = sys.argv[3]


pos_file = gtex_fusion_weights_dir + tissue_name + '.pos'
t = open(pos_file,'w')
t.write('WGT\tID\tCHR\tP0\tP1\n')
counting = 0

f = open(gene_summary_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	gene_id = data[0]
	chrom_num = data[1]
	tss = data[2]
	gene_weight_file = gtex_fusion_weights_dir + tissue_name + '_' + gene_id + '_1KG_only_fusion_output.wgt.RDat'
	if os.path.exists(gene_weight_file) == False:
		continue
	# Weight file exists so print to output
	t.write(gene_weight_file + '\t' + gene_id + '\t' + chrom_num + '\t' + tss + '\t' + tss + '\n')
	counting = counting + 1

print(tissue_name)
print(counting)

t.close()
f.close()

