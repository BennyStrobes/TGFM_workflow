import numpy as np
import os 
import sys
import pdb








###################
# Command line args
####################
tmp_pos_file = sys.argv[1]
pos_file = sys.argv[2]
gene_model_base_root = sys.argv[3]


f = open(tmp_pos_file)
t = open(pos_file,'w')

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	# get RDat file
	#Rdat_file = gene_model_base_root + data[1] + '.wgt.RDat'
	Rdat_file = data[1] + '.wgt.RDat'
	t.write(data[0] + '\t' + Rdat_file + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\t' + data[5] + '\t' + data[6] + '\n')

	# Remove unnecessary files
	bim_file = gene_model_base_root + data[1] + '_snp_bim.txt'
	os.system('rm ' + bim_file)
	weights_file = gene_model_base_root + data[1] + '_snp_gene_effects.txt'
	os.system('rm ' + weights_file)
f.close()
t.close()
os.system('rm ' + tmp_pos_file)
