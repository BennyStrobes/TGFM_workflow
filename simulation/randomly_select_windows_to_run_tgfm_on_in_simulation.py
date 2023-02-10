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

# Filter windows to those with significant number of middle var
num_middle_var = all_windows_data[:,-1].astype(float)
subset_windows_data = all_windows_data[num_middle_var > 500,:]

# Number of window subsets
total_windows = subset_windows_data.shape[0]

# Randomly select 1 window
randomly_selected_windows = np.random.choice(total_windows, size=20, replace=False)

# Print to output file
t = open(simulation_window_list_file,'w')
t.write('\t'.join(header) + '\n')
for randomly_selected_window in np.sort(randomly_selected_windows):
	t.write('\t'.join(subset_windows_data[randomly_selected_window,:]) + '\n')
t.close()

