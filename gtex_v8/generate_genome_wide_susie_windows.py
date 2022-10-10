import numpy as np 
import os
import sys
import pdb



def create_mapping_from_chrom_num_to_min_max_variant_position(sumstats_file):
	# Initialize dictionary
	dicti = {}
	for chrom_num in range(1,23):
		chrom_string = str(chrom_num)
		dicti[chrom_string] = (1000000000000, -10)

	f = open(sumstats_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		chrom_string = data[1]
		pos = float(data[2])
		af = float(data[6])
		if af < .01 or af > .99:
			continue
		if chrom_string not in dicti:
			continue
		old_tuple = dicti[chrom_string]
		dicti[chrom_string] = (min(old_tuple[0], pos), max(old_tuple[1], pos))
	f.close()
	return dicti

def window_does_not_overlap_long_ld_region(long_ld_regions, chrom_string, window_start, window_end):
	boolean = True
	for long_ld_region in long_ld_regions:
		long_ld_region_info = long_ld_region.split('_')
		long_ld_region_chrom = long_ld_region_info[0]
		long_ld_region_start = int(long_ld_region_info[1])
		long_ld_region_end = int(long_ld_region_info[2])
		if 'chr' + chrom_string != long_ld_region_chrom:
			continue
		# window start in the long ld region
		if window_start >= long_ld_region_start and window_start <= long_ld_region_end:
			boolean = False
		# window end in the long ld region
		if window_end >= long_ld_region_start and window_end <= long_ld_region_end:
			boolean = False
	return boolean

sumstats_dir = sys.argv[1]
output_file = sys.argv[2]

# Regions identified in methods of Weissbrod et al 2020 Nature Genet (hg19) and lifted over to hg38
# We will throw out windows in these regions
long_ld_regions = ['chr6_25499772_33532223', 'chr8_8142478_12142491', 'chr11_45978449_57232526']

# we are simply using this file to assess which snps are in it
# as all traits in this dir have the same snps, it doesn't matter which trait we use
sumstats_file = sumstats_dir + 'biochemistry_Cholesterol_hg38_liftover.bgen.stats'

# Create mapping from chrom num to min and max variant position
chrom_to_min_max_variant_pos = create_mapping_from_chrom_num_to_min_max_variant_position(sumstats_file)

# Print windows to output file
t = open(output_file, 'w')
t.write('chrom_num\tstart_pos_inclusive\tend_position_exclusive\tchrom_first_window_boolean\tchrom_last_window_boolean\n')

for chrom_num in range(1,23):
	print(chrom_num)
	chrom_string = str(chrom_num)
	# Starting position and ending position of variants on this chromosome
	start_pos = int(chrom_to_min_max_variant_pos[chrom_string][0])
	end_pos = int(chrom_to_min_max_variant_pos[chrom_string][1])

	cur_pos = start_pos - 1

	converged_boolean = False
	first_itera = 'True'
	last_itera = 'False'

	while converged_boolean == False:
		window_start = cur_pos
		window_end = cur_pos + 3000000

		if window_end > end_pos:
			converged_boolean = True
			last_itera = 'True'

		if window_does_not_overlap_long_ld_region(long_ld_regions, chrom_string, window_start, window_end):
			t.write(chrom_string + '\t' + str(window_start) + '\t' + str(window_end) + '\t' + first_itera + '\t' + last_itera + '\n')

		cur_pos = cur_pos + 1000000
		first_itera = 'False'





