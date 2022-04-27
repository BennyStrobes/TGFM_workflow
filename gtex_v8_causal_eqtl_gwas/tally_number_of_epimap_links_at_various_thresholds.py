import numpy as np 
import os
import sys
import pdb

def get_epimap_cell_types_and_file_names(epimap_input_dir):
	cell_types = []
	for file_name in os.listdir(epimap_input_dir):
		if file_name.endswith('links_by_group_hg38_liftover.tsv') == False:
			continue
		cell_type = file_name.split('_links_by_group_hg38')[0]
		cell_types.append(cell_type)
	cell_types = np.sort(cell_types)
	file_names = []
	for cell_type in cell_types:
		file_name = epimap_input_dir + cell_type + '_links_by_group_hg38_liftover.tsv'
		file_names.append(file_name)
	return cell_types, np.asarray(file_names)



def get_num_links(epimap_file_name, threshold):
	f = open(epimap_file_name)
	counter = 0
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if float(data[2]) > threshold:
			counter = counter + 1
	f.close()
	return counter




epimap_input_dir = sys.argv[1]
epimap_ukbb_overlap_dir = sys.argv[2]


epimap_cell_types, epimap_file_names = get_epimap_cell_types_and_file_names(epimap_ukbb_overlap_dir)


output_file = epimap_ukbb_overlap_dir + 'number_epimap_links.txt'
t = open(output_file,'w')
t.write('cell_type\tthreshold\tnum_links\n')

thresholds = [.1, .3, .5, .7, .9]

for threshold in thresholds:
	print(threshold)
	for cell_type_iter, epimap_cell_type in enumerate(epimap_cell_types):
		epimap_file_name = epimap_file_names[cell_type_iter]

		num_links = get_num_links(epimap_file_name, threshold)

		t.write(epimap_cell_type + '\t' + str(threshold) + '\t' + str(num_links) + '\n')

t.close()