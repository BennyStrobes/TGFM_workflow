import numpy as np 
import os
import sys
import pdb



def extract_all_cell_types(abc_file):
	f = open(abc_file)
	cell_types = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 24:
			print('assumption erroror')
			pdb.set_trace()
		cell_types[data[-1]] = 1
	f.close()
	return np.sort([*cell_types])






abc_dir = sys.argv[1]



abc_file = abc_dir + 'abc_links_by_cell_type_hg38_liftover.tsv'


cell_types = extract_all_cell_types(abc_file)


# Filter cell types to roadmap and encode cell types
filtered_cell_types = []
for ct in cell_types:
	if 'Roadmap' in ct or 'ENCODE' in ct:
		filtered_cell_types.append(ct)


t_dicti = {}
for ct in filtered_cell_types:
	t_dicti[ct] = open(abc_dir + ct + '_abc_links_by_single_cell_type_hg38_liftover.tsv','w')
	t_dicti[ct].write('gene\tname\tscore\tgroup\tchr\tstart\tend\tname\n')

head_count = 0
f = open(abc_file)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	gene = data[6]
	name = data[4]
	abc_score = data[-4]
	ct = data[-1]
	chrom_num = data[0]
	start = data[1]
	end = data[2]
	if ct in t_dicti:
		t_dicti[ct].write(gene + '\t' + name + '\t' + abc_score + '\t' + ct + '\t' + chrom_num + '\t' + start + '\t' + end + '\t' + name + '\n')
f.close()


for ct in filtered_cell_types:
	t_dicti[ct].close()


