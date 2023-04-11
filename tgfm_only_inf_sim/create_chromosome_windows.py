import numpy as np 
import os
import sys
import pdb




def extract_variant_positions(reference_bim):
	f = open(reference_bim)
	var_pos_arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		var_pos = int(data[3])
		var_pos_arr.append(var_pos)
	f.close()

	var_pos_arr = np.asarray(var_pos_arr)

	# Quick error checking
	if np.array_equal(var_pos_arr, np.sort(var_pos_arr)) == False:
		print('assumption eororor')
		pdb.set_trace()
	return var_pos_arr






#####################
# Command line args
#####################
window_size_mb = float(sys.argv[1])
reference_bim = sys.argv[2]
window_file = sys.argv[3]  # Output file
chrom_num = sys.argv[4]


# Extract variant positions from bim file
variant_positions = extract_variant_positions(reference_bim)

# Start and stop positions
start_pos = np.min(variant_positions)
end_pos = np.max(variant_positions)

mb_size = 1000000

# Initialize window file
t = open(window_file,'w')
t.write('window_name' + '\t' + 'window_start_inclusive\twindow_middle_start_inclusive\twindow_middle_end_exclusive\twindow_end_exclusive\tnum_variants_in_middle\n')

# Loop through windows and print to output
cur_window_start = start_pos
arr = []
while end_pos >= (cur_window_start + mb_size):
	cur_window_middle_start = cur_window_start + (1*mb_size)
	cur_window_middle_end = cur_window_start + (2*mb_size)
	cur_window_end = cur_window_start + (3*mb_size)

	window_name = 'chr' + chrom_num + '_' + str(cur_window_start) + '_' + str(cur_window_end)

	num_middle_snps = np.sum((variant_positions >= cur_window_middle_start) & (variant_positions < cur_window_middle_end))

	t.write(window_name + '\t' + str(cur_window_start) + '\t' + str(cur_window_middle_start) + '\t' + str(cur_window_middle_end) + '\t' + str(cur_window_end) + '\t' + str(num_middle_snps) + '\n')

	cur_window_start = cur_window_start + mb_size
	arr.append(num_middle_snps)

t.close()


