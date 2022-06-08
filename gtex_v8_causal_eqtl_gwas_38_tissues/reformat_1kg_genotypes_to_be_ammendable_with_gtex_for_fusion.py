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





gtex_genotype_dir = sys.argv[1]
ref_1kg_genotype_dir = sys.argv[2]
gtex_fusion_processed_intermediate_data = sys.argv[3]


gtex_variants = get_dictionary_list_of_gtex_variants(gtex_genotype_dir)


for chrom_num in range(1,23):
	used_variants = {}
	input_bim = ref_1kg_genotype_dir + '1000G.EUR.hg38.' + str(chrom_num) + '.bim'
	output_bim = gtex_fusion_processed_intermediate_data + '1000G.EUR.gtex_formatted.hg38.' + str(chrom_num) + '.bim'
	f = open(input_bim)
	t = open(output_bim,'w')

	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		variant_id1 = 'chr' + data[0] + '_' + data[3] + '_' + data[4] + '_' + data[5] + '_b38'
		variant_id2 = 'chr' + data[0] + '_' + data[3] + '_' + data[5] + '_' + data[4] + '_b38'
		#if variant_id1 in used_variants or variant_id2 in used_variants:
		#	print('assumptino eroror')
		#	pdb.set_trace()
		if variant_id1 in gtex_variants and variant_id2 in gtex_variants:
			print('assumption erororo')
			pdb.set_trace()
		elif variant_id1 in gtex_variants:
			t.write(data[0] + '\t' + variant_id1 + '\t' + '\t'.join(data[2:]) + '\n')
			used_variants[variant_id1] = 1
		elif variant_id2 in gtex_variants:
			t.write(data[0] + '\t' + variant_id2 + '\t' + '\t'.join(data[2:]) + '\n')
			used_variants[variant_id2] = 1
		else:
			t.write(line + '\n')
	f.close()
	t.close()

	# Cp bed and fam files
	command1 = 'cp ' + ref_1kg_genotype_dir + '1000G.EUR.hg38.' + str(chrom_num) + '.fam ' + gtex_fusion_processed_intermediate_data + '1000G.EUR.gtex_formatted.hg38.' + str(chrom_num) + '.fam'
	command2 = 'cp ' + ref_1kg_genotype_dir + '1000G.EUR.hg38.' + str(chrom_num) + '.bed ' + gtex_fusion_processed_intermediate_data + '1000G.EUR.gtex_formatted.hg38.' + str(chrom_num) + '.bed'
	os.system(command1)
	os.system(command2)

