import numpy as np 
import os
import sys
import pdb



def fill_in_chrom_genes(chrom_genes, chrom_num, composit_tissue_gene_summary_file):
	f = open(composit_tissue_gene_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_name = data[0]
		if data[1] != str(chrom_num):
			continue
		if gene_name not in chrom_genes:
			chrom_genes[gene_name] = []
		chrom_genes[gene_name].append(line)
	f.close()
	return chrom_genes

def extract_relevent_fields_from_unparsed_array(arr):
	tss = -100
	phen_files = "NA"
	cov_files = "NA"
	tissues = []
	for stringer in arr:
		data = stringer.split('\t')
		if tss == -100:
			tss = int(data[2])
		else:
			if int(data[2]) != tss:
				print('assumption eroror')
				pdb.set_trace()

		if phen_files == 'NA':
			phen_files = data[3]
		else:
			phen_files = phen_files + ',' + data[3]

		if cov_files == 'NA':
			cov_files = data[4]
		else:
			cov_files = cov_files + ',' + data[4]

		tiss_name = data[3].split('/')[-2]
		tissues.append(tiss_name)

	return tss, phen_files, cov_files, tissues


pseudotissue_name = sys.argv[1]
composit_tissue_string = sys.argv[2]
gtex_processed_expression_dir = sys.argv[3]
gtex_pseudotissue_gene_model_input_dir = sys.argv[4]
num_jobs = int(sys.argv[5])



composit_tissues = composit_tissue_string.split(',')

output_file = gtex_pseudotissue_gene_model_input_dir + pseudotissue_name + '_gene_summary.txt'

t = open(output_file,'w')

t.write('Gene_id\tchrom_num\tTSS\tgene_pheno_file\tcovariate_file\tcomposit_tissues\n')


for chrom_num in range(1,23):
	chrom_genes = {}

	for composit_tissue in composit_tissues:
		composit_tissue_gene_summary_file = gtex_processed_expression_dir + composit_tissue + '/' + composit_tissue + '_gene_summary.txt'
		chrom_genes = fill_in_chrom_genes(chrom_genes, chrom_num, composit_tissue_gene_summary_file)

	# loop through genes on this chrosome
	for gene_name in np.asarray([*chrom_genes]):

		tss, pheno_files, cov_files, gene_composit_tissues = extract_relevent_fields_from_unparsed_array(chrom_genes[gene_name])

		if len(gene_composit_tissues) > len(composit_tissues):
			print('assumption eroror')
			pdb.set_trace()
		if len(gene_composit_tissues) == 0:
			print('assusmptione roror')
			pdb.set_trace()
		t.write(gene_name + '\t' + str(chrom_num) + '\t' + str(tss) + '\t' + pheno_files + '\t' + cov_files + '\t' + ','.join(np.asarray(gene_composit_tissues)) + '\n')

t.close()





# Parallelize
data_raw = np.loadtxt(output_file, dtype=str,delimiter='\t')[1:,:]
nrows = data_raw.shape[0]
splits = np.array_split(np.arange(nrows),num_jobs)

for job_num, split_mat in enumerate(splits):
	split_output_file = gtex_pseudotissue_gene_model_input_dir + pseudotissue_name + '_gene_summary_' + str(job_num) + '.txt'
	valid_rows = {}
	for row_num in split_mat:
		valid_rows[row_num] = 1
	t = open(split_output_file,'w')
	t.write('Gene_id\tchrom_num\tTSS\tgene_pheno_file\tcovariate_file\tcomposit_tissues\n')
	f = open(output_file)
	head_count = 0
	line_num = 0
	for line in f:
		line = line.rstrip()
		if head_count == 0:
			head_count = head_count + 1
			continue
		if line_num in valid_rows:
			t.write(line + '\n')
		line_num = line_num + 1
	f.close()
	t.close()

