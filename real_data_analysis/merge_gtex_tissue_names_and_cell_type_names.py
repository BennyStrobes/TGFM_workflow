import numpy as np 
import os
import sys
import pdb


def get_number_of_genes_with_a_gene_model_in_gtex(tissue_name, gtex_gene_models_dir):
	tissue_summary_file = gtex_gene_models_dir + tissue_name + '/' + tissue_name + '_component_gene_pos_file.txt'
	aa = np.loadtxt(tissue_summary_file,dtype=str,delimiter='\t')
	n_genes = aa.shape[0] - 1
	return n_genes

def get_mapping_from_ct_to_n_cells_per_individual(filer):
	mapping = {}
	head_count = 0
	f = open(filer)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		ct = data[2]
		n_cells_per_indi_file = data[7]
		tmp = np.loadtxt(n_cells_per_indi_file,dtype=str,delimiter ='\t')
		avg_n_cells_per_indi = np.mean(tmp[1:,1].astype(float))
		if ct in mapping:
			print('assumption eroror')
			pdb.set_trace()
		mapping[ct] = avg_n_cells_per_indi
	f.close()
	return mapping


gtex_pseudotissue_file = sys.argv[1]
pb_cell_type_file = sys.argv[2]
merged_tissue_cell_type_file = sys.argv[3]
gtex_gene_models_dir = sys.argv[4]
sc_pbmc_gene_models_dir = sys.argv[5]
sc_pseudobulk_expression_dir = sys.argv[6]

# Get number of cells per individual from the pseudobulk data
ct_to_n_cells_per_indi = get_mapping_from_ct_to_n_cells_per_individual(sc_pseudobulk_expression_dir + 'pseudobulk_data_set_summary.txt')

t = open(merged_tissue_cell_type_file,'w')
t.write('context\tcontext_sample_size\tdata_set\tn_gene_model_genes\tavg_n_cells_per_individual\n')

f = open(gtex_pseudotissue_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	tissue_name = data[0]
	n_genes = get_number_of_genes_with_a_gene_model_in_gtex(tissue_name, gtex_gene_models_dir)
	t.write(data[0] + '\t' + data[1] + '\t' + 'GTEx' + '\t' + str(n_genes) + '\tNA' + '\n')
f.close()

f = open(pb_cell_type_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	tissue_name = data[2]
	n_genes = get_number_of_genes_with_a_gene_model_in_gtex(tissue_name, sc_pbmc_gene_models_dir)
	t.write(data[2] + '\t'  + data[5] + '\t' + 'Perez_sc_pbmc\t' + str(n_genes) + '\t' + str(ct_to_n_cells_per_indi[tissue_name]) + '\n')
f.close()


t.close()

print(merged_tissue_cell_type_file)