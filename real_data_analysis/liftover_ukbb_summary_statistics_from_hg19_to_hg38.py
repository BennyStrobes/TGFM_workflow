import numpy as np 
import os
import sys
import pdb
import gzip

def association_file_to_bed_file(hg19_gwas_file, temp_hg19_bed):
	f = gzip.open(hg19_gwas_file)
	t = open(temp_hg19_bed,'w')
	head_count = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 16:
			print('assumption eroror')
			pdb.set_trace()
		chrom_num = data[1]
		pos = data[2]
		pos2 = int(pos) + 1
		t.write('chr' + chrom_num + '\t' + pos + '\t' + str(pos2) + '\n')
	f.close()
	t.close()


def create_ordered_list_of_gwas_studies_by_parsing_input_directory(ukbb_sumstats_hg19_dir):
	gwas_studies = []
	for file_name in os.listdir(ukbb_sumstats_hg19_dir):
		if file_name.endswith('.bgen.stats.gz') == False:
			continue
		if file_name.endswith('.interim.bgen.stats.gz'):
			continue
		study_name = file_name.split('v3.')[-1].split('.bgen')[0]
		if study_name.startswith('460K'):
			continue
		gwas_studies.append(study_name)
	if len(gwas_studies) != len(np.unique(gwas_studies)):
		print('assumption eroror')
		pdb.set_trace()

	gwas_studies = np.sort(np.asarray(gwas_studies))
	return gwas_studies


#Run liftOver with parameters specified by Chris Wilks (https://github.com/ChristopherWilks/snaptron/blob/master/scripts/liftover_snaptron_coords.pl)
def run_liftover(input_file, output_file, missing_file, liftover_directory):
    stringer = liftover_directory + 'liftOver -minMatch=1.00 -ends=2 ' + input_file + ' ' + liftover_directory + 'hg19ToHg38.over.chain.gz ' + output_file + ' ' + missing_file
    os.system(stringer)


def recreate_cafeh_association_file_in_hg38_format(hg19_association_file, hg38_association_file, liftover_output_file, liftover_missing_file):
	missing = {}
	f = open(liftover_missing_file)
	for line in f:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		data = line.split('\t')
		missing[data[0] + '_' + data[1]] = 1
	f.close()
	f = gzip.open(hg19_association_file)
	t = open(hg38_association_file,'w')
	g = open(liftover_output_file)
	head_count = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		chrom_num = data[1]
		pos = data[2]
		snp_key = 'chr' + chrom_num + '_' + pos
		if snp_key in missing:
			continue
		hg38_info = next(g).rstrip().split('\t')
		hg38_pos = hg38_info[1]	
		hg38_chrom = hg38_info[0].split('hr')[1]
		t.write(data[0] + '\t' + hg38_chrom + '\t' + hg38_pos + '\t' + '\t'.join(data[3:]) + '\n')		
	f.close()
	g.close()
	t.close()

def get_num_lines(file_name):
	counter = 0
	f = open(file_name)
	for line in f:
		counter = counter + 1
	f.close()
	return counter

def error_checking(liftover_output_file, hg19_association_file):
	num_lines_lift = get_num_lines(liftover_output_file)
	num_lines_association = get_num_lines(hg19_association_file) - 1
	if num_lines_association != num_lines_lift:
		print('assumption erororo')

def filter_sum_stats_based_on_maf(hg38_gwas_file, hg38_gwas_maf_filter_file, maf_thresh):
	f = open(hg38_gwas_file)
	t = open(hg38_gwas_maf_filter_file,'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		af = float(data[6])
		if af < maf_thresh or af > (1.0 - maf_thresh):
			continue
		t.write(line + '\n')
	f.close()
	t.close()


def create_ukbb_snp_list(hg38_gwas_maf_filter_file, ukbb_snp_list):
	f = open(hg38_gwas_file)
	t = open(ukbb_snp_list,'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rs_id = data[0]
		snp_id = 'chr' + data[1] + '_' + data[2] + '_' + data[4] + '_' + data[5]
		t.write(snp_id + '\t' + rs_id + '\n')
	f.close()
	t.close()


liftover_directory = sys.argv[1]
ukbb_sumstats_hg19_dir = sys.argv[2]
ukbb_sumstats_hg38_dir = sys.argv[3]


# Create ordered list of gwas studies
ordered_gwas_studies = create_ordered_list_of_gwas_studies_by_parsing_input_directory(ukbb_sumstats_hg19_dir)


# Now liftover each study to hg38
# Also print studies and new study files to output file
study_organizer_file = ukbb_sumstats_hg38_dir + 'ukbb_hg38_sumstat_files.txt'
t_outer = open(study_organizer_file,'w')
t_outer.write('study_name\tstudy_file\n')
# Loop through studies

for itera, gwas_study in enumerate(ordered_gwas_studies):
	print(gwas_study)
	# Original gwas file in hg19 coordinate space
	hg19_gwas_file = ukbb_sumstats_hg19_dir + 'bolt_337K_unrelStringentBrit_MAF0.001_v3.' + gwas_study + '.bgen.stats.gz'

	# First convert hg19 association to temporary bed format
	temp_hg19_bed = ukbb_sumstats_hg38_dir + gwas_study + '_temp_bed_hg19.txt'
	association_file_to_bed_file(hg19_gwas_file, temp_hg19_bed)

	# Run liftover
	liftover_output_file = ukbb_sumstats_hg38_dir + gwas_study + '_liftover_output.txt'
	liftover_missing_file = ukbb_sumstats_hg38_dir + gwas_study + '_liftover_missing.txt'
	run_liftover(temp_hg19_bed, liftover_output_file, liftover_missing_file, liftover_directory)

	# Recreate gwas association file in hg38 format
	hg38_gwas_file = ukbb_sumstats_hg38_dir + gwas_study + '_hg38_liftover.bgen.stats'
	recreate_cafeh_association_file_in_hg38_format(hg19_gwas_file, hg38_gwas_file, liftover_output_file, liftover_missing_file)

	# Quick error checking
	error_checking(liftover_output_file, hg38_gwas_file)


	# Print study and new lifted over file to output
	t_outer.write(gwas_study + '\t' + hg38_gwas_file + '\n')


	# Remove temporary files
	os.system('rm ' + liftover_output_file)
	os.system('rm ' + liftover_missing_file)
	os.system('rm ' + temp_hg19_bed)

	# For first trait (ie only do this once), print list of snps in UKBB
	if itera == 0:
		ukbb_snp_list = ukbb_sumstats_hg38_dir + 'ukbb_hg38_liftover_bgen_snps.txt'
		create_ukbb_snp_list(hg38_gwas_file, ukbb_snp_list)
t_outer.close()
