import numpy as np
import os
import sys
import pdb












simulated_gwas_dir = sys.argv[1]
simulation_name_string = sys.argv[2]
global_window_file = sys.argv[3]

# Loop through global windows
f = open(global_window_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	window_name = data[0]

	# Get name of gwas sum stat file
	gwas_sumstat_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results_window_' + window_name + '.txt'
	# Remove gwas sumstat
	os.system('rm ' + gwas_sumstat_file)
f.close()