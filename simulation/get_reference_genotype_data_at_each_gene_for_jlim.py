import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np
import os
import pdb
from pandas_plink import read_plink1_bin


def convert_genotype_data_to_jlim_format(plink_window_stem, gene_name, processed_jlim_genotype_data_dir):
	# Extract genotype data
	G_obj = read_plink1_bin(plink_window_stem + '.bed', plink_window_stem + '.bim', plink_window_stem + '.fam', verbose=False)

	G_obj_geno = G_obj.values # Numpy 2d array of dimension num samples X num snps
	G_obj_chrom = np.asarray(G_obj.chrom)
	G_obj_pos = np.asarray(G_obj.pos)
	# For our purposes, a0 is the effect allele
	# For case of plink package, a0 is the first column in the plink bim file
	G_obj_a0 = np.asarray(G_obj.a0)
	G_obj_a1 = np.asarray(G_obj.a1)
	# RSids
	G_obj_rsids = np.asarray(G_obj.snp)
	G_obj_sample_names = np.asarray(G_obj.sample)
	# Snp ids
	G_obj_snp_ids = 'chr' + G_obj_chrom + '_' + (G_obj_pos.astype(str)) + '_' + G_obj_a0 + '_' + G_obj_a1

	# Output file
	jlim_geno_output_file = processed_jlim_genotype_data_dir + 'chr' + G_obj_chrom[0] + '.' + str(np.min(G_obj_pos)) + '_' + str(np.max(G_obj_pos)) + '.txt'
	t = open(jlim_geno_output_file,'w')

	for ii, snp_pos in enumerate(G_obj_pos):
		snp_a0 = G_obj_a0[ii]
		snp_a1 = G_obj_a1[ii]
		t.write(G_obj_chrom[0] + '\t' + str(snp_pos) + '\t' + snp_a0 + '\t' + snp_a1 + '\t' + 'GT')

		for genotype_val in G_obj_geno[:, ii]:
			if np.isnan(genotype_val):
				af = np.nanmean(G_obj_geno[:, ii])/2
				random_a0 = np.random.choice([0,1],p=[1-af,af])
				random_a1 = np.random.choice([0,1],p=[1-af,af])

				if random_a0 == 0:
					tmp_0 = snp_a0
				else:
					tmp_0 = snp_a1
				if random_a1 == 0:
					tmp_1 = snp_a0
				else:
					tmp_1 = snp_a1
				t.write('\t' + tmp_0 + '\t' + tmp_1)
			elif genotype_val == 0.0:
				t.write('\t' + snp_a0 + '\t' + snp_a0)
			elif genotype_val == 1.0:
				t.write('\t' + snp_a0 + '\t' + snp_a1)
			elif genotype_val == 2.0:
				t.write('\t' + snp_a1 + '\t' + snp_a1)
			else:
				print('genotype assumptino eroror')
				pdb.set_trace()
		t.write('\n')
	t.close()

	return jlim_geno_output_file, G_obj_rsids



reference_genotype_stem = sys.argv[1]
simulated_gene_position_file = sys.argv[2]
processed_jlim_genotype_data_dir = sys.argv[3]  # output dir
jlim_window_summary_file = sys.argv[4]


np.random.seed(0)


t_summary = open(jlim_window_summary_file,'w')
t_summary.write('gene_name\tgene_tss\tjlim_geno_file_name\twindow_snps\n')
cis_window=100000

# Loop through genes
f = open(simulated_gene_position_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	gene_name = data[1]
	chrom = data[0]
	gene_tss = int(data[2])
	if chrom != '1':
		print('assumption error')
		pdb.set_trace()


	start_pos = gene_tss - int(cis_window)
	end_pos = gene_tss + int(cis_window)

	if start_pos < 0:
		start_pos=0


	plink_window_stem = processed_jlim_genotype_data_dir + gene_name + '_tmp'
	command_string = 'plink --bfile ' + reference_genotype_stem + ' --keep-allele-order --threads 1 --make-bed --out ' + plink_window_stem + ' --chr ' + chrom + ' --from-bp ' + str(start_pos) + ' --to-bp ' + str(end_pos) +' --allow-no-sex'
	os.system(command_string)

	# Check if genotype file is made (not made due to no variants in window)
	if os.path.isfile(plink_window_stem + '.bed') == False:
		os.system('rm ' + plink_window_stem + '.log')
		continue

	# Convert genotype file into JLIM format
	jlim_geno_file_name, window_rsids = convert_genotype_data_to_jlim_format(plink_window_stem, gene_name, processed_jlim_genotype_data_dir)

	# Remove unncessary files
	os.system('rm ' + plink_window_stem + '.bed')
	os.system('rm ' + plink_window_stem + '.bim')
	os.system('rm ' + plink_window_stem + '.fam')
	os.system('rm ' + plink_window_stem + '.log')

	t_summary.write(gene_name + '\t' + str(gene_tss) + '\t' + jlim_geno_file_name + '\t' + ';'.join(window_rsids) + '\n')
	t_summary.flush()

f.close()
t_summary.close()