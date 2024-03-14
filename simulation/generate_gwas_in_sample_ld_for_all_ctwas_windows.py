import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os
import pdb
from pandas_plink import read_plink1_bin
import pyreadr





def extract_dictionary_list_of_window_rsids(bim_file, window_start, window_end):
	dicti = {}
	f = open(bim_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		rs_id = data[1]
		snp_pos = int(data[3])
		if snp_pos >= window_start and snp_pos < window_end: # in window
			# Quick error check
			if rs_id in dicti:
				print('assumption eroror')
				pdb.set_trace()
			dicti[rs_id] = 1
	f.close()
	return dicti



def load_in_ordered_array_of_rsids_from_bim_file(genotype_bim):
	f = open(genotype_bim)
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		arr.append(data[1])
	f.close()
	arr = np.asarray(arr)

	# Quick error check 
	if len(np.unique(arr)) != len(arr):
		print('assumption eroror')
		pdb.set_trace()

	return arr


def get_window_snp_info(gwas_plink_stem, window_rsids):
	# Load in ordered array of rsids
	genotype_bim = gwas_plink_stem + '.bim'
	ordered_rsids = load_in_ordered_array_of_rsids_from_bim_file(genotype_bim)
	# Initialize snp df
	snp_df = []

	# Stream genotype bim file
	f=open(genotype_bim)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# ignore snps not from window
		if data[1] not in window_rsids:
			continue
		tmp_arr = []
		tmp_arr.append(data[0])
		tmp_arr.append(data[1])
		tmp_arr.append(data[3])
		tmp_arr.append(data[4])
		tmp_arr.append(data[5])
		tmp_arr.append('1.0')
		tmp_arr = np.asarray(tmp_arr)

		snp_df.append(tmp_arr)
	f.close()
	return snp_df



def compute_ld_on_window(gwas_plink_stem, window_rsids):
	# Load in ordered array of rsids
	genotype_bim = gwas_plink_stem + '.bim'
	ordered_rsids = load_in_ordered_array_of_rsids_from_bim_file(genotype_bim)

	# Load in variant names in plink format
	window_variant_names = []
	for snp_iter, snp_rsid in enumerate(ordered_rsids):
		# Skip snps that are not hapmap3 snps
		if snp_rsid not in window_rsids:
			continue
		window_variant_names.append('variant' + str(snp_iter))
	window_variant_names = np.asarray(window_variant_names)

	# Load in genotype object
	genotype_obj = read_plink1_bin(gwas_plink_stem + '.bed', gwas_plink_stem + '.bim', gwas_plink_stem + '.fam', verbose=False)


	# Extract genotype data for this batch of snps
	window_variant_genotype = np.asarray(genotype_obj.sel(variant=window_variant_names))

	# Standardize and mean impute genotype
	nvary = window_variant_genotype.shape[1]
	for var_iter in range(nvary):
		variant_genotype = window_variant_genotype[:,var_iter]
		std_variant_genotype = np.copy(variant_genotype)
		nan_indices = np.isnan(variant_genotype)
		non_nan_mean = np.mean(variant_genotype[nan_indices==False])
		std_variant_genotype[nan_indices] = non_nan_mean
		std_variant_genotype = (std_variant_genotype - np.mean(std_variant_genotype))/np.std(std_variant_genotype)
		window_variant_genotype[:,var_iter] = std_variant_genotype


	# Compute LD
	LD_mat = np.corrcoef(np.transpose(window_variant_genotype))

	return LD_mat



##############################
# Command line argumemnts
##############################
chrom_num = sys.argv[1]
ctwas_window_list_file = sys.argv[2]
processed_genotype_data_dir = sys.argv[3]
processed_ctwas_genotype_data_dir = sys.argv[4]


# Bim file for this chromosomse
bim_file = processed_genotype_data_dir + 'simulated_gwas_data_' + chrom_num + '.bim'

# Loop through simulation windows
g = open(ctwas_window_list_file)
head_count = 0
for line in g:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue

	# Extract relevent info from this window
	window_chrom = data[0]
	window_start = int(data[1])
	window_end = int(data[2])

	# Extract dictionary list of window rsids
	window_rsids = extract_dictionary_list_of_window_rsids(bim_file, window_start, window_end)

	# outtput stem
	output_stem = processed_ctwas_genotype_data_dir + 'ctwas_window_chr' + str(chrom_num) + '.R_snp.' + str(window_start) + '_' + str(window_end)
	print(str(chrom_num) + '.R_snp.' + str(window_start) + '_' + str(window_end))

	####################################################
	# Save variant info
	####################################################
	gwas_plink_stem = processed_genotype_data_dir + 'simulated_gwas_data_' + str(chrom_num)  # Genotype files
	window_snp_df = get_window_snp_info(gwas_plink_stem, window_rsids)
	Rvar_file = output_stem + '.Rvar'
	t = open(Rvar_file,'w')
	t.write('chrom\tid\tpos\talt\tref\tvariance\n')
	for window_snp in window_snp_df:
		t.write('\t'.join(window_snp) + '\n')
	t.close()

	####################################################
	# extract window LD 
	####################################################
	window_ld = compute_ld_on_window(gwas_plink_stem, window_rsids)
	# Save window ld to output
	temporary_rds_file = output_stem + '_temp.RDS'
	perm_rds_file = output_stem + '.RDS'
	# Save to temporary RDS through pyreadr (though the output of this an R data frame)
	pyreadr.write_rds(temporary_rds_file, pd.DataFrame(window_ld))
	# Convert RDS object from R data frame to R matrix
	os.system('Rscript convert_rds_file_from_df_to_matrix.R ' + temporary_rds_file + ' ' + perm_rds_file)
	# Delete temporary file
	os.system('rm ' + temporary_rds_file)

g.close()