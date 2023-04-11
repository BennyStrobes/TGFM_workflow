import numpy as np 
import os
import sys
import pdb
import gzip


def extract_list_of_1kg_rsids(bim_file):
	f = open(bim_file)
	rs_ids = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		rsid = data[1]
		rs_ids[rsid] =1
	f.close()
	return rs_ids




ldsc_baseline_annotation_dir = sys.argv[1] # input file
chrom_num = sys.argv[2]
kg_genotype_dir = sys.argv[3]
ldsc_annotation_rs_id_file = sys.argv[4] # output file

# This is currently a hack and should be fixed.
# These are 3 variants on chrom 1 that have no variation across the 100 eqtl samples
variants_to_exclude = {}
variants_to_exclude['rs116631404'] = 1
variants_to_exclude['rs17103322'] = 1
variants_to_exclude['rs12033962'] = 1



kg_rs_ids = extract_list_of_1kg_rsids(kg_genotype_dir + '1000G.EUR.QC.' + chrom_num + '.bim')


annot_file = ldsc_baseline_annotation_dir + 'baseline.' + str(chrom_num) + '.annot.gz'

t = open(ldsc_annotation_rs_id_file,'w')

f = gzip.open(annot_file)
used = {}
head_count = 0
skipped = 0
for line in f:
	line = line.decode('utf-8').rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		continue
	rs_id = data[2]
	if rs_id.startswith('rs') == False:
		print('assumption erororo')
		pdb.set_trace()
	if rs_id not in kg_rs_ids:
		skipped = skipped + 1
		continue
	if rs_id in variants_to_exclude:
		continue

	t.write(rs_id + '\n')
	if rs_id in used:
		print('assumption eroror')
		pdb.set_trace()
	used[rs_id] = 1
t.close()
f.close()

