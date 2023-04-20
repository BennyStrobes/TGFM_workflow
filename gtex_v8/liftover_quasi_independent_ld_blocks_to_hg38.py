import numpy as np 
import os
import sys
import pdb



def create_hg19_bed_file_from_chrom_specific_input_files(quasi_independent_dir, hg19_bed_file):
	t = open(hg19_bed_file,'w')

	for chrom_num in range(1,23):
		chrom_specific_input_file = quasi_independent_dir + 'EUR/' + 'fourier_ls-chr' + str(chrom_num) + '.bed'
		f = open(chrom_specific_input_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split()
			if head_count == 0:
				head_count = head_count + 1
				continue
			if float(data[2]) <= float(data[1]):
				print('assumption eroor')
				pdb.set_trace()
			if data[0] != 'chr' + str(chrom_num):
				print('assumption eroor')
				pdb.set_trace()
			t.write('\t'.join(np.asarray(data)) + '\n')
		f.close()
	t.close()

#Run liftOver with parameters specified by Chris Wilks (https://github.com/ChristopherWilks/snaptron/blob/master/scripts/liftover_snaptron_coords.pl)
def run_liftover(input_file, output_file, missing_file, liftover_directory):
    stringer = liftover_directory + 'liftOver -minMatch=1.00 -ends=2 ' + input_file + ' ' + liftover_directory + 'hg19ToHg38.over.chain.gz ' + output_file + ' ' + missing_file
    os.system(stringer)

def create_hg38_bed_file_from_liftover_output(liftover_output_file, hg38_bed_file):
	f = open(liftover_output_file)
	t = open(hg38_bed_file,'w')
	prev_chrom = 'chr0'
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		curr_chrom = data[0]
		curr_end = data[2]
		curr_start = data[1]
		if float(curr_start) >= float(curr_end):
			print('assumption eroror')
			pdb.set_trace()
		if prev_chrom != curr_chrom:
			prev_chrom = curr_chrom
			prev_end = curr_end
			t.write('\t'.join(np.asarray(data)) + '\n')
		else:
			if curr_start == prev_end:
				prev_chrom = curr_chrom
				prev_end = curr_end
				t.write('\t'.join(np.asarray(data)) + '\n')
			else:
				t.write(data[0] + '\t' + prev_end + '\t' + curr_end + '\n')
				prev_chrom = curr_chrom
				prev_end = curr_end
	f.close()
	t.close()

def create_large_ld_blocks_file(hg38_bed_file, large_hg38_bed_file, merge_size):
	f = open(hg38_bed_file)
	t = open(large_hg38_bed_file,'w')
	merge_counter = 0
	prev_chrom = 'chr0'
	totes = 0.0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		curr_chrom = data[0]
		curr_end = data[2]
		curr_start = data[1]
		if float(curr_start) >= float(curr_end):
			print('assumption eroror')
			pdb.set_trace()
		if prev_chrom != curr_chrom:
			if prev_chrom != 'chr0':
				totes = totes + merge_counter
				t.write(prev_chrom + '\t' + merge_start + '\t' + merge_end + '\n')
			prev_chrom = curr_chrom
			merge_start = curr_start
			merge_end = curr_end
			merge_counter = 1
		else:
			merge_end = curr_end
			merge_counter = merge_counter + 1
			if merge_counter >= merge_size:
				totes = totes + merge_counter
				t.write(prev_chrom + '\t' + merge_start + '\t' + merge_end + '\n')
				prev_chrom = 'chr0'
	if prev_chrom != 'chr0':
		t.write(prev_chrom + '\t' + merge_start + '\t' + merge_end + '\n')

	f.close()
	t.close()



liftover_directory = sys.argv[1]
quasi_independent_dir = sys.argv[2]
quasi_independent_ld_blocks_hg38_dir = sys.argv[3]


# First create bed file in hg19
hg19_bed_file = quasi_independent_ld_blocks_hg38_dir + 'quasi_independent_ld_blocks_hg19_bed.txt'
create_hg19_bed_file_from_chrom_specific_input_files(quasi_independent_dir, hg19_bed_file)

# Run liftover
liftover_output_file = quasi_independent_ld_blocks_hg38_dir + '_liftover_output.txt'
liftover_missing_file = quasi_independent_ld_blocks_hg38_dir + '_liftover_missing.txt'
run_liftover(hg19_bed_file, liftover_output_file, liftover_missing_file, liftover_directory)

# create ld blocks bed file hg38
hg38_bed_file = quasi_independent_ld_blocks_hg38_dir + 'quasi_independent_ld_blocks_hg38_bed.txt'
create_hg38_bed_file_from_liftover_output(liftover_output_file, hg38_bed_file)

# create large ld blocks bed file hg38
large_hg38_bed_file = quasi_independent_ld_blocks_hg38_dir + 'large_10_quasi_independent_ld_blocks_hg38_bed.txt'
create_large_ld_blocks_file(hg38_bed_file, large_hg38_bed_file, 10)


