import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pyreadr
import numpy as np 
import os
import sys
import pdb








pseudotissue_name = sys.argv[1]
composit_tissue_string = sys.argv[2]
gtex_fusion_weights_data_dir = sys.argv[3]
gtex_susie_pmces_fusion_weights_dir = sys.argv[4]
pseudotissue_gtex_susie_pmces_fusion_weights_dir = sys.argv[5]  # Output file

print(pseudotissue_name)


# Extract composit tissues from composit tissue string
composit_tissues = composit_tissue_string.split(',')

genes = {}
# Loop through composit tissues
for composit_tissue in composit_tissues:
	# Get gene summary file for composit tissue
	gene_summary_file = gtex_fusion_weights_data_dir + composit_tissue + '/' + composit_tissue + '_gene_summary.txt'
	head_count = 0
	counter = 0
	f = open(gene_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		chrom_num = data[1]
		tss = data[2]
		gene_weight_file = gtex_susie_pmces_fusion_weights_dir + composit_tissue + '/' + composit_tissue + '_' + gene_id + '_1KG_only_fusion_output.wgt.RDat'
		if os.path.exists(gene_weight_file) == False:
			continue
		result = pyreadr.read_r(gene_weight_file)
		if np.var(np.asarray(result['susie_mu'])) <= 0.0:
			continue
		if np.sum(np.sum(np.asarray(result['susie_V']))) <= 0.0:
			continue
		if np.asarray(result['susie_converged']['susie_converged'])[0] == False:
			continue
		hsq = np.asarray(result['hsq'])[0][0]
		hsq_p = np.asarray(result['hsq.pv'])[0,0]
		if hsq < 0.0 or hsq_p > .01:
			continue
		output_line = gene_weight_file + '\t' + gene_id + '\t' + chrom_num + '\t' + tss + '\t' + tss + '\t' + str(hsq) + '\t' + str(hsq_p) + '\n'
		if gene_id not in genes:
			genes[gene_id] = output_line
	f.close()



pos_file = pseudotissue_gtex_susie_pmces_fusion_weights_dir + pseudotissue_name + '_cis_heritable_genes.pos'
t = open(pos_file,'w')
t.write('WGT\tID\tCHR\tP0\tP1\thsq\thsq_p\n')

gene_names = np.asarray([*genes])


for gene_name in gene_names:
	t.write(genes[gene_name])
t.close()








# Extract composit tissues from composit tissue string
composit_tissues = composit_tissue_string.split(',')

genes = {}
# Loop through composit tissues
for composit_tissue in composit_tissues:
	# Get gene summary file for composit tissue
	gene_summary_file = gtex_fusion_weights_data_dir + composit_tissue + '/' + composit_tissue + '_gene_summary.txt'
	head_count = 0
	counter = 0
	f = open(gene_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		chrom_num = data[1]
		tss = data[2]
		gene_weight_file = gtex_susie_pmces_fusion_weights_dir + composit_tissue + '/' + composit_tissue + '_' + gene_id + '_1KG_only_fusion_output.wgt.RDat'
		if os.path.exists(gene_weight_file) == False:
			continue
		result = pyreadr.read_r(gene_weight_file)
		if np.var(np.asarray(result['susie_mu'])) <= 0.0:
			continue
		if np.sum(np.sum(np.asarray(result['susie_V']))) <= 0.0:
			continue
		if np.asarray(result['susie_converged']['susie_converged'])[0] == False:
			continue
		hsq = np.asarray(result['hsq'])[0][0]
		hsq_p = np.asarray(result['hsq.pv'])[0,0]
		output_line = gene_weight_file + '\t' + gene_id + '\t' + chrom_num + '\t' + tss + '\t' + tss + '\t' + str(hsq) + '\t' + str(hsq_p) + '\n'
		if gene_id not in genes:
			genes[gene_id] = output_line
	f.close()



pos_file = pseudotissue_gtex_susie_pmces_fusion_weights_dir + pseudotissue_name + '_all_genes.pos'
t = open(pos_file,'w')
t.write('WGT\tID\tCHR\tP0\tP1\thsq\thsq_p\n')

gene_names = np.asarray([*genes])


for gene_name in gene_names:
	t.write(genes[gene_name])
t.close()
