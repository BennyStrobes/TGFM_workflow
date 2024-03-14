import numpy as np 
import os
import sys
import pdb
import gzip


def create_mapping_from_gene_name_to_tss_info_based_on_gtex_data(gtex_gene_list):
	f = open(gtex_gene_list)
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[0].split('.')[0]
		chrom_num = data[1]
		tss = data[-1]
		if ensamble_id in mapping:
			print('assumption oerororor')
			pdb.set_trace()
		mapping[ensamble_id] = (chrom_num, tss)
	f.close()

	return mapping


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
		if gene_chrom_num == 'chrX' or gene_chrom_num == 'chrY':
			continue
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
		ensamble_id = ensamble_id.split('.')[0]
		if ensamble_id not in mapping:
			mapping[ensamble_id] = (gene_type, gene_chrom_num, gene_strand, tss)
		else:
			if mapping[ensamble_id][0] != gene_type:
				print('assumption eroror')
				pdb.set_trace()
			if mapping[ensamble_id][1] != gene_chrom_num:
				print('assumption eroror')
				pdb.set_trace()
			'''
			if mapping[ensamble_id][2] != gene_strand:
				print('assumption eroror')
				pdb.set_trace()
			'''
			if mapping[ensamble_id][3] != tss:
				print('assumption eroror')
				pdb.set_trace()
	f.close()
	return mapping


def get_list_of_pseudobulk_cell_types(ct_summmary_file):
	f = open(ct_summmary_file)
	cell_types = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		cell_types.append(data[2])
	f.close()
	return np.asarray(cell_types)

def extract_cell_type_gene_names_from_expression_file(pb_cell_type_expression_file):
	gene_names = []
	f = open(pb_cell_type_expression_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_names.append(data[0])
	f.close()
	gene_names = np.asarray(gene_names)
	if len(np.unique(gene_names)) != len(gene_names):
		print('assumption erororo')
		pdb.set_trace()
	return gene_names

######################################
# Command line args
######################################
pseudobulk_expression_dir = sys.argv[1]
hg38_gene_annotation_file = sys.argv[2]
gtex_gene_list = sys.argv[3]


######################
# Create mapping from gene name to chrom_num and tss based on gtex data
gtex_gene_mapping = create_mapping_from_gene_name_to_tss_info_based_on_gtex_data(gtex_gene_list)


######################
# Create mapping from gene name to chrom_num and tss based on gene annotation file
mapping = create_mapping_from_gene_name_to_gene_info(hg38_gene_annotation_file)


'''
#######################
# Get list of pseudobulk cell types
ct_summmary_file = pseudobulk_expression_dir + 'pseudobulk_data_set_summary_filtered.txt'
pb_cell_types = get_list_of_pseudobulk_cell_types(ct_summmary_file)

# Loop through cell types
for pb_cell_type in pb_cell_types:
	print(pb_cell_type)
	# Extract ordered gene names for this pb cell type
	pb_cell_type_expression_file = pseudobulk_expression_dir + 'cell_types_' + pb_cell_type + '_pseudobulk_expression.txt'
	cell_type_gene_names = extract_cell_type_gene_names_from_expression_file(pb_cell_type_expression_file)

	# Open output file handle
	pb_cell_type_gene_tss_file = pseudobulk_expression_dir + 'cell_types_' + pb_cell_type + '_gene_tss_hg38.txt'
	t = open(pb_cell_type_gene_tss_file,'w')
	t.write('gene_name\tchrom_num\ttss\n')
	for gene_name in cell_type_gene_names:
		if gene_name not in mapping:
			print('skipped')
			#pdb.set_trace()
			if gene_name in gtex_gene_mapping:
				print('assumption eroror')
				pdb.set_trace()
			continue

		gene_info = mapping[gene_name]
		chrom_num = gene_info[1]
		tss = gene_info[3]

		# Error checking
		if gene_name in gtex_gene_mapping:
			if chrom_num != gtex_gene_mapping[gene_name][0]:
				print('assumption eroroor')
				pdb.set_trace()
			if tss != gtex_gene_mapping[gene_name][1]:
				print('assumption eroror')
				pdb.set_trace()
		t.write(gene_name + '\t' + chrom_num + '\t' + tss + '\n')
	t.close()
'''


# Also do for bulk_PBMC
# Extract ordered gene names for this pb cell type
pb_cell_type_expression_file = pseudobulk_expression_dir + 'cell_types_PBMC' + '_pseudobulk_expression.txt'
cell_type_gene_names = extract_cell_type_gene_names_from_expression_file(pb_cell_type_expression_file)

# Open output file handle
pb_cell_type_gene_tss_file = pseudobulk_expression_dir + 'cell_types_PBMC'+ '_gene_tss_hg38.txt'
t = open(pb_cell_type_gene_tss_file,'w')
t.write('gene_name\tchrom_num\ttss\n')
for gene_name in cell_type_gene_names:
	if gene_name not in mapping:
		print('skipped')
		#pdb.set_trace()
		if gene_name in gtex_gene_mapping:
			print('assumption eroror')
			pdb.set_trace()
		continue

	gene_info = mapping[gene_name]
	chrom_num = gene_info[1]
	tss = gene_info[3]

	# Error checking
	if gene_name in gtex_gene_mapping:
		if chrom_num != gtex_gene_mapping[gene_name][0]:
			print('assumption eroroor')
			pdb.set_trace()
		if tss != gtex_gene_mapping[gene_name][1]:
			print('assumption eroror')
			pdb.set_trace()
	t.write(gene_name + '\t' + chrom_num + '\t' + tss + '\n')
t.close()

print(pb_cell_type_gene_tss_file)
