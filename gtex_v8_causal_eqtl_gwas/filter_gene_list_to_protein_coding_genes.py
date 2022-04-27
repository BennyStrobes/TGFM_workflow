import numpy as np 
import os
import sys
import pdb


def create_mapping_from_gene_name_to_gene_info(gene_annotation_file):
	f = open(gene_annotation_file)
	mapping = {}
	for line in f:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		data = line.split('\t')
		if len(data) != 9:
			print('assumption eroror')
			pdb.set_trace()
		if data[2] != 'gene':
			continue
		ensamble_id = 'null'
		gene_type = 'null'
		gene_info = data[8].split(';')
		for info in gene_info:
			if info.startswith('gene_id'):
				ensamble_id = info.split('"')[1]
			elif info.startswith(' gene_type'):
				gene_type = info.split('"')[1]
		if ensamble_id == 'null' or gene_type == 'null':
			print('assumption eroror')
			pdb.set_trace()
		gene_chrom_num = data[0]
		gene_strand = data[6]
		if float(data[3]) > float(data[4]):
			print('assumption erroror')
			pdb.set_trace()
		if gene_strand == '+':
			tss = data[3]
		elif gene_strand == '-':
			tss = data[4]
		else:
			print('assumption error')


		# Add to info
		if ensamble_id not in mapping:
			mapping[ensamble_id] = (gene_type, gene_chrom_num, gene_strand, tss)
		else:
			if mapping[ensamble_id][0] != gene_type:
				print('assumption eroror')
				pdb.set_trace()
			if mapping[ensamble_id][1] != gene_chrom_num:
				print('assumption eroror')
				pdb.set_trace()
			if mapping[ensamble_id][2] != gene_strand:
				print('assumption eroror')
				pdb.set_trace()
			if mapping[ensamble_id][3] != tss:
				print('assumption eroror')
				pdb.set_trace()
	f.close()
	return mapping


original_gene_file = sys.argv[1]
gene_file = sys.argv[2]
gene_annotation_file = sys.argv[3]


mapping = create_mapping_from_gene_name_to_gene_info(gene_annotation_file)



f = open(original_gene_file)
t = open(gene_file,'w')

head_count = 0

for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	ensamble_id = data[0]
	gene_info = mapping[ensamble_id]
	if gene_info[0] != 'protein_coding':
		continue
	t.write(line + '\n')
t.close()
f.close()