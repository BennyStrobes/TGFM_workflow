import numpy as np 
import os
import sys
import pdb
import gzip




def make_tmp_hg19_bed_file_from_bim(input_bim, temp_hg19_bed_file):
	f = open(input_bim)
	t = open(temp_hg19_bed_file,'w')
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		chrom_num = data[0]
		var_id = data[1]
		pos = data[3]
		t.write('chr' + str(chrom_num) + '\t' + pos + '\t' + str(int(pos) + 1) + '\t' + var_id + '\n')
	f.close()
	t.close()
	return



#Run liftOver with parameters specified by Chris Wilks (https://github.com/ChristopherWilks/snaptron/blob/master/scripts/liftover_snaptron_coords.pl)
def run_liftover(input_file, output_file, missing_file, liftover_directory):
    stringer = liftover_directory + 'liftOver -minMatch=1.00 -ends=2 ' + input_file + ' ' + liftover_directory + 'hg19ToHg38.over.chain.gz ' + output_file + ' ' + missing_file
    os.system(stringer)

    return

def create_mapping_from_hg19_var_id_to_hg38_pos(liftover_output_file):
	f = open(liftover_output_file)
	mapping = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		var_id = data[3]
		new_chrom = data[0].split('hr')[1]
		new_pos = data[1]
		mapping[var_id] = (new_chrom, new_pos)
	f.close()
	return mapping

##################
# Command line args
#####################
input_bim = sys.argv[1]
output_root = sys.argv[2]
liftover_directory = sys.argv[3]



# Make temp bed file
temp_hg19_bed_file = output_root + '_tmp_hg19_bed_file.txt'
make_tmp_hg19_bed_file_from_bim(input_bim, temp_hg19_bed_file)


# Run liftover
liftover_output_file = output_root + '_tmp_hg38_bed_file.txt'
liftover_missing_file = output_root + '_tmp_hg38_missing_file.txt'
run_liftover(temp_hg19_bed_file, liftover_output_file, liftover_missing_file, liftover_directory)


# Create mapping from hg19_variant_id to hg38_pos
hg19_var_id_to_hg38_pos = create_mapping_from_hg19_var_id_to_hg38_pos(liftover_output_file)

# Create new bim file
hg38_bim = output_root + '.bim'
f = open(input_bim)
t = open(hg38_bim,'w')
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	old_var_id = data[1]
	old_chrom_num = data[0]
	old_pos = data[3]
	third_col = data[2]
	a1 = data[4]
	a2 = data[5]
	# Quick error check
	if a1 != old_var_id.split(':')[3] or a2 != old_var_id.split(':')[2]:
		print('assumption erorro')
		pdb.set_trace()
	if old_var_id not in hg19_var_id_to_hg38_pos:
		# These did not map
		# These will be filtered out in the next step
		new_var_id = 'chr' + old_chrom_num + '_' + old_pos + '_missing_' + a2 + '_' + a1 + '_b38'
		t.write(old_chrom_num + '\t' + new_var_id + '\t' + third_col + '\t' + old_pos + '\t' + a1 + '\t' + a2 + '\n')
	else:
		new_chrom_num, new_pos = hg19_var_id_to_hg38_pos[old_var_id]
		new_var_id = 'chr' + new_chrom_num + '_' + new_pos + '_' + a2 + '_' + a1 + '_b38'
		t.write(new_chrom_num + '\t' + new_var_id + '\t' + third_col + '\t' + new_pos + '\t' + a1 + '\t' + a2 + '\n')
f.close()
t.close()






