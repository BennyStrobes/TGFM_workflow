import numpy as np 
import os
import sys
import pdb


def extract_positions_of_start_and_of_simulation(simulation_window_list_file):
	start_pos = 'NULL'
	end_pos = 'NULL'
	head_count = 0
	start_count = 0
	f = open(simulation_window_list_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_start = float(data[1])
		window_end = float(data[4])
		if start_count == 0:
			start_count = start_count + 1
			start_pos = window_start
			end_pos = window_end
			prev_end_pos = window_end
		else:
			end_pos = window_end
			if window_start > prev_end_pos:
				print('assumptino eroror')
			prev_end_pos = window_end

	f.close()
	if start_pos == 'NULL' or end_pos == 'NULL':
		print('assumption eroror')
		pdb.set_trace()
	return start_pos, end_pos






simulation_window_list_file = sys.argv[1]
simulated_gene_position_file = sys.argv[2]
simulated_revised_gene_position_file = sys.argv[3]


global_window_start_pos, global_window_end_pos = extract_positions_of_start_and_of_simulation(simulation_window_list_file)

f = open(simulated_gene_position_file)
t = open(simulated_revised_gene_position_file,'w')

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	tss = float(data[2])
	if tss >= global_window_start_pos and tss < global_window_end_pos:
		t.write(line + '\n')
f.close()
t.close()


