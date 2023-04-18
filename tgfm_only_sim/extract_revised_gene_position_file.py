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



def extract_window_starts_and_ends(simulation_window_list_file):
	window_starts = []
	window_ends = []
	f = open(simulation_window_list_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_start = int(data[1])
		window_end = int(data[4])
		window_starts.append(window_start)
		window_ends.append(window_end)
	return np.asarray(window_starts), np.asarray(window_ends)

def gene_tss_lies_in_a_window(tss, window_starts, window_ends):
	booler = False
	n_windows = len(window_starts)
	if n_windows != len(window_ends):
		print('assumption errror')
		pdb.set_trace()

	for window_iter in range(n_windows):
		if tss >= window_starts[window_iter] and tss < window_ends[window_iter]:
			booler = True
	return booler



simulation_window_list_file = sys.argv[1]
simulated_gene_position_file = sys.argv[2]
simulated_revised_gene_position_file = sys.argv[3]


window_starts, window_ends = extract_window_starts_and_ends(simulation_window_list_file)
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
	if gene_tss_lies_in_a_window(tss, window_starts, window_ends):
		t.write(line + '\n')
f.close()
t.close()


