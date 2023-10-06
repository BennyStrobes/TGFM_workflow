import numpy as np 
import os
import sys
import pdb




def extract_chromosome_variant_list_from_winow_file(window_file):
	f = open(window_file)
	variant_list = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_variant_file = data[7]
		variants = np.loadtxt(window_variant_file,dtype=str)
		for variant in variants:
			variant_list[variant] = 1
	f.close()

	# Quick error checking
	new_variant_list = {}
	for variant in np.asarray([*variant_list]):
		variant_info = variant.split('_')
		a1 = variant_info[2]
		a2 = variant_info[3]
		alt_variant = variant_info[0] + '_' + variant_info[1] + '_' + variant_info[3] + '_' + variant_info[2]
		if alt_variant in variant_list:
			print('assumption errorororo')
			pdb.set_trace()
		if a1 == 'A' and a2 == 'T':
			continue
		if a1 == 'T' and a2 == 'A':
			continue
		if a1 == 'C' and a2 == 'G':
			continue
		if a1 == 'G' and a2 == 'C':
			continue
		new_variant_list[variant] = 1

	return new_variant_list



ukbb_preprocessed_for_genome_wide_susie_dir = sys.argv[1]
ukbb_gtex_formatted_bim = sys.argv[2]


t = open(ukbb_gtex_formatted_bim, 'w')

for chrom_num_int in range(1,23):
	chrom_num = str(chrom_num_int)


	# Extract dictionary list of all variants used by UKBB in this analysis
	window_file = ukbb_preprocessed_for_genome_wide_susie_dir + 'genome_wide_susie_windows_and_processed_data_chrom_' + chrom_num + '.txt'
	variant_list = extract_chromosome_variant_list_from_winow_file(window_file)

	# Print those variants to Bim formatted file
	for variant in np.asarray([*variant_list]):
		variant_info = variant.split('_')
		alt_variant = variant_info[0] + '_' + variant_info[1] + '_' + variant_info[3] + '_' + variant_info[2] + '_b38'
		std_variant = variant + '_b38'
		t.write(std_variant + '\n')
		t.write(alt_variant + '\n')
t.close()
