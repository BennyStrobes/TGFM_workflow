import numpy as np 
import os
import sys
import pdb
import gzip



def extract_track_names(epimap_track_file):
	track_names = []
	f = open(epimap_track_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		track_names.append(data[0])
	f.close()
	track_names = np.asarray(track_names)

	if len(track_names) != len(np.unique(track_names)):
		print('assumption eroror')
		pdb.set_trace()

	return track_names



#Run liftOver with parameters specified by Chris Wilks (https://github.com/ChristopherWilks/snaptron/blob/master/scripts/liftover_snaptron_coords.pl)
def run_liftover(input_file, output_file, missing_file, liftover_directory):
    stringer = liftover_directory + 'liftOver -minMatch=1.00 -ends=2 ' + input_file + ' ' + liftover_directory + 'hg19ToHg38.over.chain.gz ' + output_file + ' ' + missing_file
    os.system(stringer)


epimap_track_file = sys.argv[1]
epimap_data_dir = sys.argv[2] # input dir
epimap_enrichment_raw_data_dir = sys.argv[3] # Ouput dir
liftover_directory = sys.argv[4]


# Extrack array of names of all epimap tracks
track_names = extract_track_names(epimap_track_file)

# Loop through epimap tracks (and for each track liftover)
for track_name in track_names:

	# track hg19 bed file
	track_hg19_bed_file = epimap_data_dir + 'bed_from_signal/' + 'Epimap.' + track_name + '.bed'

	# Run liftover
	liftover_output_file = epimap_enrichment_raw_data_dir + 'Epimap.' + track_name + '.hg38.bed'
	liftover_missing_file = epimap_enrichment_raw_data_dir + track_name + '_liftover_missing.txt'
	run_liftover(track_hg19_bed_file, liftover_output_file, liftover_missing_file, liftover_directory)


