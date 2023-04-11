import numpy as np 
import os
import sys
import gzip
import pdb



def create_mapping_from_rsid_to_cm_position(kg_genotype_bim_file):
	f = open(kg_genotype_bim_file)
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		rsid = data[1]
		cm_pos = data[2]
		if rsid in dicti:
			print('assumption errorr')
			pdb.set_trace()
		dicti[rsid] = cm_pos
	f.close()
	return dicti


def create_mapping_from_rsid_to_genomic_annotation(baseline_annotation_file):
	f = gzip.open(baseline_annotation_file)
	head_count = 0
	dicti = {}
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			header = np.asarray(data)
			head_count = head_count + 1
			continue
		rsid = data[2]
		if rsid in dicti:
			print('assumption eroror')
			pdb.set_trace()
		dicti[rsid] = np.asarray(data)
	f.close()
	return dicti, header






genotype_bim_file = sys.argv[1]
kg_genotype_bim_file = sys.argv[2] # For CM distance
baseline_annotation_file = sys.argv[3]
output_annotation_file = sys.argv[4]


# Create mapping from rsid to cm position
rsid_to_cm_position = create_mapping_from_rsid_to_cm_position(kg_genotype_bim_file)


# Create mapping from rsid to genomic annotation
rsid_to_genomic_annotation, annotation_header = create_mapping_from_rsid_to_genomic_annotation(baseline_annotation_file)


# Stream simulation bim file and for each variant create annotaiotn file
t = open(output_annotation_file,'w')
t.write('\t'.join(annotation_header) + '\n')
f = open(genotype_bim_file)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	line_rsid = data[1]
	line_pos = data[3]

	kg_cm = rsid_to_cm_position[line_rsid]
	annotation_vec = rsid_to_genomic_annotation[line_rsid]

	# Quick error checks
	if kg_cm != annotation_vec[3]:
		print('assumption erroro')
		pdb.set_trace()
	if line_pos != annotation_vec[1]:
		print('assumption eroror')
		pdb.set_trace()
	# Print to output
	t.write('\t'.join(annotation_vec) + '\n')
f.close()
t.close()



