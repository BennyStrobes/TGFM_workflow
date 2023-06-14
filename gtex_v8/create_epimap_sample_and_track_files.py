import numpy as np 
import os
import sys
import pdb




def create_mapping_from_bs_id_to_tissue_group(file_name):
	f = open(file_name)
	mapping = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		bs_id = data[0]
		tissue_type = data[1]
		mapping[bs_id] = tissue_type
		if len(data) != 2:
			print('assumption eroror')
			pdb.set_trace()
	f.close()
	return mapping

def get_track_names(dir_name):
	arr = []
	for file_name in os.listdir(dir_name):
		info = file_name.split('.')
		track_type = info[1]
		sample_num = info[2]
		track_name = track_type + '.' + sample_num
		arr.append(track_name)
	arr = np.sort(np.asarray(arr))
	return arr

def create_mapping_from_bs_id_to_tissue_class(file_name):
	f = open(file_name)
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		mapping[data[0]] = data[1]
	f.close()
	return mapping



epimap_data_dir = sys.argv[1]
epimap_sample_file = sys.argv[2]
epimap_track_file = sys.argv[3]

# Get ordered list of bs_ids
ordered_bs_ids = np.loadtxt(epimap_data_dir + 'samplelist_833_final.tsv', dtype=str)

# create mapping from bs_id to tissue group
bs_id_to_tissue_group = create_mapping_from_bs_id_to_tissue_group(epimap_data_dir + 'mnemonic_mapping.tsv')
bs_id_to_tissue_class = create_mapping_from_bs_id_to_tissue_class(epimap_data_dir + 'main_metadata_table.tsv')

# Print Epimap sample file
t = open(epimap_sample_file, 'w')
t.write('standard_sample_name\tbs_id\ttissue_class\ttissue_group\n')
for itera, bs_id in enumerate(ordered_bs_ids):
	if bs_id_to_tissue_class[bs_id] == 'Cancer':
		continue
	t.write(str(itera+1) + '\t' + bs_id + '\t' + bs_id_to_tissue_class[bs_id] + '\t' + bs_id_to_tissue_group[bs_id] + '\n')
t.close()


# Extract track names
track_names = get_track_names(epimap_data_dir + 'bed_from_signal/')

# Print track names
t = open(epimap_track_file,'w')
t.write('track_name\tstandard_sample_name\tbs_id\ttissue_class\ttissue_group\ttrack_type\n')
for track_name in track_names:
	info = track_name.split('.')
	track_type = info[0]
	sample_name = info[1]
	bs_id = ordered_bs_ids[int(sample_name) -1]
	tissue_group = bs_id_to_tissue_group[bs_id]
	tissue_class = bs_id_to_tissue_class[bs_id]
	if tissue_class == 'Cancer':
		continue
	t.write(track_name + '\t' + sample_name + '\t' + bs_id + '\t' + tissue_class + '\t' + tissue_group + '\t' + track_type + '\n')
t.close()









