import numpy as np 
import os
import sys
import pdb
import gzip

def get_transcript_and_gene_id_from_gencode_gene_info(gene_info):
	transcript_id = 'NA'
	gene_id = 'NA'
	fields = gene_info.split(';')

	for field in fields:
		if field.startswith('gene_id'):
			gene_id = field.split('"')[1]
		elif field.startswith(' transcript_id'):
			transcript_id = field.split('"')[1]
	return transcript_id, gene_id



def create_mapping_from_enst_transcript_ids_to_ensamble_ids(hg38_gene_annotation_file):
	f = gzip.open(hg38_gene_annotation_file)
	dicti = {}

	for line in f:
		line = line.decode('utf-8').rstrip()
		if line.startswith('##'):
			continue
		data = line.split('\t')
		if len(data) != 9:
			print('assumption eroror')
			pdb.set_trace()
		gene_info = data[8]
		transcript_id, gene_id = get_transcript_and_gene_id_from_gencode_gene_info(gene_info)
		if transcript_id == 'NA' or gene_id == 'NA':
			continue

		if transcript_id not in dicti:
			dicti[transcript_id] = gene_id
		else:
			if dicti[transcript_id] != gene_id:
				print('assumptin eororr')
				pdb.set_trace()
	f.close()
	return dicti


def create_mapping_from_entrez_ids_to_ensamble_ids(hg38_ensamble_trascript_to_entrez_id_file, enst_ids_to_ensamble_ids):
	dicti = {}
	f = open(hg38_ensamble_trascript_to_entrez_id_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 2:
			print('assumption eororor')
		enst_id = data[0]
		entrez_id = data[1]
		if enst_id not in enst_ids_to_ensamble_ids:
			continue
		ensamble_id = enst_ids_to_ensamble_ids[enst_id]

		if entrez_id not in dicti:
			dicti[entrez_id] = [ensamble_id]
		else:
			old_arr = dicti[entrez_id]
			old_arr.append(ensamble_id)
			dicti[entrez_id] = old_arr
	f.close()

	new_dicti = {}
	arr = []
	for entrez_id in [*dicti]:
		unique_ensambles = np.unique(dicti[entrez_id])
		if len(unique_ensambles) != 1:
			continue
		new_dicti[entrez_id] = unique_ensambles[0]
	return new_dicti



hg38_gene_annotation_file = sys.argv[1]
hg38_ensamble_trascript_to_entrez_id_file = sys.argv[2]
biological_pathway_gene_set_file = sys.argv[3]
ensamble_biological_pathway_gene_set_file = sys.argv[4]


# First create mapping from enst transcript ids to ensamble ids
enst_ids_to_ensamble_ids = create_mapping_from_enst_transcript_ids_to_ensamble_ids(hg38_gene_annotation_file)

# Next create mapping from entrez id to ensamble ids
entrez_to_to_ensamble_id = create_mapping_from_entrez_ids_to_ensamble_ids(hg38_ensamble_trascript_to_entrez_id_file, enst_ids_to_ensamble_ids)


# Convert and print to output
f = open(biological_pathway_gene_set_file, encoding = "ISO-8859-1")
t = open(ensamble_biological_pathway_gene_set_file,'w')

head_count = 0
diffs = []
counter = 0
for line in f:
	counter = counter + 1

	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	if head_count == 1:
		t.write(line + '\n')
		head_count = head_count + 1
		continue
	entrez_id_string = data[5]
	entrez_ids = np.asarray(entrez_id_string.split('[')[1].split(']')[0].split(', '))
	ensamble_ids = []
	for entrez_id in entrez_ids:
		if entrez_id in entrez_to_to_ensamble_id:
			ensamble_ids.append(entrez_to_to_ensamble_id[entrez_id])
	ensamble_ids = np.asarray(ensamble_ids)
	ensamble_ids = np.unique(ensamble_ids)
	diff = len(entrez_ids) - len(ensamble_ids)
	diffs.append(diff)

	num_genes = len(ensamble_ids)

	if num_genes > 20:
		t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + data[4])
		t.write('\t' + ','.join(ensamble_ids) + '\t' + str(num_genes) + '\n')

f.close()
t.close()










