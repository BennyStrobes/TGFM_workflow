import numpy as np 
import os
import sys
import pdb











########################
# Command line args
########################
global_window_file = sys.argv[1]
simulation_window_list_file = sys.argv[2]
simulation_number = int(sys.argv[3])

# Set set
np.random.seed(simulation_number)

# Load in all windows
all_windows_data_raw = np.loadtxt(global_window_file, dtype=str, delimiter='\t')
header = all_windows_data_raw[0,:]
all_windows_data = all_windows_data_raw[1:,:]
total_windows = all_windows_data.shape[0]

# Randomly select 1 window
randomly_selected_window = np.random.choice(total_windows)

# Print to output file
t = open(simulation_window_list_file,'w')
t.write('\t'.join(header) + '\n')
t.write('\t'.join(all_windows_data[randomly_selected_window,:]) + '\n')
t.close()