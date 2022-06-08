import numpy as np 
import os
import sys
import pdb


def get_dictionary_list_of_gtex_variants(gtex_genotype_dir):
	gtex_variants = {}
	for chrom_num in range(1,23):
		chrom_bim_file = gtex_genotype_dir + 'Adipose_Subcutaneous_GTEx_v8_genotype_EUR_overlap_1kg_and_ukbb_' + str(chrom_num) + '.bim'
		f = open(chrom_bim_file)
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			line_variant_id = data[1]
			if line_variant_id in gtex_variants:
				print('assumption eroror')
				pdb.set_trace()
			gtex_variants[line_variant_id] = 0
		f.close()
	return gtex_variants


def reformat_ukbb_sumstats_to_be_ammendable_with_gtex(study_name, full_file_name, gtex_variants, output_file):
	used_variants = {}
	t = open(output_file,'w')
	t.write('SNP\tA1\tA2\tZ\tbeta\tbeta_se\n')

	head_count = 0
	f = open(full_file_name)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_id1 = 'chr' + data[1] + '_' + data[2] + '_' + data[4] + '_' + data[5] + '_b38'
		variant_id2 = 'chr' + data[1] + '_' + data[2] + '_' + data[5] + '_' + data[4] + '_b38'
		beta = float(data[10])
		std_err = float(data[11])
		zscore = beta/std_err
		if variant_id1 in gtex_variants and variant_id2 in gtex_variants:
			print('assumption eroror')
			pdb.set_trace()
		elif variant_id1 in gtex_variants and variant_id1 not in used_variants:
			t.write(variant_id1 + '\t' + data[4] + '\t' + data[5] + '\t' + str(zscore) + '\t' + str(beta) + '\t' + str(std_err) + '\n')
			used_variants[variant_id1] = 1
		elif variant_id2 in gtex_variants and variant_id2 not in used_variants:
			t.write(variant_id2 + '\t' + data[5] + '\t' + data[4] + '\t' + str(-zscore) + '\t' + str(-beta) + '\t' + str(std_err) + '\n')
			used_variants[variant_id2] = 1
	f.close()
	t.close()


gtex_genotype_dir = sys.argv[1]
ukbb_sumstats_hg38_dir = sys.argv[2]
gtex_fusion_processed_intermediate_data = sys.argv[3]


gtex_variants = get_dictionary_list_of_gtex_variants(gtex_genotype_dir)

study_organization_file = ukbb_sumstats_hg38_dir + 'ukbb_hg38_sumstat_files_with_samp_size_and_h2.txt'
f = open(study_organization_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	study_name = data[0]
	full_file_name = data[1]

	print(study_name)
	output_file = gtex_fusion_processed_intermediate_data + study_name + '_fusion_processed_sumstats_hg38.txt'
	reformat_ukbb_sumstats_to_be_ammendable_with_gtex(study_name, full_file_name, gtex_variants, output_file)
f.close()
