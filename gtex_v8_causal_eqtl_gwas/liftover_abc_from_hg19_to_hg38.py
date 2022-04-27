import numpy as np 
import os
import pdb
import sys
import gzip


def get_epimap_cell_types_and_file_names(epimap_input_dir):
	cell_types = []
	for file_name in os.listdir(epimap_input_dir):
		if file_name.endswith('.tsv.gz') == False:
			continue
		if file_name.startswith('links_by_group') == False:
			continue
		cell_type = file_name.split('by_group.')[1].split('.tsv')[0]
		cell_types.append(cell_type)
	cell_types = np.sort(cell_types)
	file_names = []
	for cell_type in cell_types:
		file_name = epimap_input_dir + 'links_by_group.' + cell_type + '.tsv.gz'
		file_names.append(file_name)
	return cell_types, np.asarray(file_names)

#Run liftOver with parameters specified by Chris Wilks (https://github.com/ChristopherWilks/snaptron/blob/master/scripts/liftover_snaptron_coords.pl)
def run_liftover(input_file, output_file, missing_file, liftover_directory):
    stringer = liftover_directory + 'liftOver -minMatch=1.00 -ends=2 ' + input_file + ' ' + liftover_directory + 'hg19ToHg38.over.chain.gz ' + output_file + ' ' + missing_file
    os.system(stringer)


def abc_file_to_bed_file(abc_file, temp_hg19_bed):
	f = gzip.open(abc_file)
	t = open(temp_hg19_bed,'w')
	head_count = 0
	for line in f:
		line =line.rstrip().decode('utf-8')
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\n')
	f.close()
	t.close()

def recreate_abc_file_in_hg38_format(hg19_association_file, hg38_association_file, liftover_output_file, liftover_missing_file):
	missing = {}
	f = open(liftover_missing_file)
	for line in f:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		data = line.split('\t')
		missing[data[0] + '_' + data[1] + '_' + data[2]] = 1
	f.close()
	f = gzip.open(hg19_association_file)
	t = open(hg38_association_file,'w')
	g = open(liftover_output_file)
	head_count = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		chrom_num = data[0]
		pos_start = data[1]
		pos_end = data[2]
		snp_key = chrom_num + '_' + pos_start + '_' + pos_end
		if snp_key in missing:
			continue
		hg38_info = next(g).rstrip().split('\t')
		hg38_pos_start = hg38_info[1]
		hg38_pos_end = hg38_info[2]	
		hg38_chrom = hg38_info[0]
		t.write(hg38_chrom + '\t' + hg38_pos_start + '\t' + hg38_pos_end + '\t' + '\t'.join(data[3:]) + '\n')
	f.close()
	g.close()
	t.close()


abc_input_dir = sys.argv[1]
abc_ukbb_genome_wide_susie_overlap_dir = sys.argv[2]
liftover_directory = sys.argv[3]


abc_input_file = abc_input_dir + 'AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz'

# First convert abc file to bed file
temp_hg19_bed = abc_ukbb_genome_wide_susie_overlap_dir + 'temp_bed_hg19.txt'
abc_file_to_bed_file(abc_input_file, temp_hg19_bed)

# Run liftover
liftover_output_file = abc_ukbb_genome_wide_susie_overlap_dir + 'liftover_output.txt'
liftover_missing_file = abc_ukbb_genome_wide_susie_overlap_dir + 'liftover_missing.txt'
run_liftover(temp_hg19_bed, liftover_output_file, liftover_missing_file, liftover_directory)


# Recreate abc file in hg38 format
hg38_abc_file = abc_ukbb_genome_wide_susie_overlap_dir + 'abc_links_by_cell_type_hg38_liftover.tsv'
recreate_abc_file_in_hg38_format(abc_input_file, hg38_abc_file, liftover_output_file, liftover_missing_file)

# Remove temporary files
os.system('rm ' + liftover_output_file)
os.system('rm ' + liftover_missing_file)
os.system('rm ' + temp_hg19_bed)

