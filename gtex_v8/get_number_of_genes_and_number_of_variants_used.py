import numpy as np 
import os
import sys
import pdb



def compute_num_reference_variants(ukkbb_window_summary_file, num_variants_output_file):
	f = open(ukkbb_window_summary_file)
	head_count = 0
	varz = {}
	count = 0
	for line in f:
		count = count + 1
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_var_file = data[7]
		var_list = np.loadtxt(window_var_file, dtype=str)
		for var in var_list:
			varz[var] = 1
	f.close()
	num_variants = len(varz)
	t = open(num_variants_output_file,'w')
	t.write(str(num_variants) + '\n')
	t.close()


def compute_num_reference_genes(gtex_susie_gene_models_dir, tissues, num_genes_output_file):
	t = open(num_genes_output_file,'w')
	t.write('tissue\tnum_genes\n')
	for tissue in tissues:
		tissue_pos_file = gtex_susie_gene_models_dir + tissue + '/' + tissue + '_cis_heritable_gene_pos_file.txt'
		tmp = np.loadtxt(tissue_pos_file, dtype=str)
		num_genes = tmp.shape[0] - 1
		t.write(tissue + '\t' + str(num_genes) + '\n')
	t.close()






ukkbb_window_summary_file = sys.argv[1]
gtex_susie_gene_models_dir = sys.argv[2]
gtex_tissue_file = sys.argv[3]
num_genes_and_variants_dir = sys.argv[4]

tissues = np.loadtxt(gtex_tissue_file,dtype=str)[1:,0]



num_variants_output_file = num_genes_and_variants_dir + 'number_of_reference_variants.txt'
#compute_num_reference_variants(ukkbb_window_summary_file, num_variants_output_file)

num_genes_output_file = num_genes_and_variants_dir + 'number_of_reference_genes.txt'
compute_num_reference_genes(gtex_susie_gene_models_dir, tissues, num_genes_output_file)
