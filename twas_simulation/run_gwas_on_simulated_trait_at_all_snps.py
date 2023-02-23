import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os
import pdb
from pandas_plink import read_plink1_bin
import statsmodels.api as sm

# Extract dictionary list of hapmap3 rsids (ie regression snps in sldsc)
def extract_dictionary_list_of_hapmap3_rsids(hm3_rsid_file):
	# Initialize dictionary
	dicti = {}

	# Stream rsid file
	f = open(hm3_rsid_file)
	for line in f:
		line = line.rstrip()
		if line in dicti:
			print('assumption eroror')
			pdb.set_trace()
		dicti[line] = 1
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

# Run GWAS on hapmap3 rsids (ie regression snps in sldsc)
def run_gwas_on_hm3_rsids(hapmap3_rsids, trait_values_file, gwas_plink_stem, gwas_output_file, batch_size=500):
	# Load in trait vector
	trait_vector = np.loadtxt(trait_values_file)

	# Load in ordered array of rsids
	genotype_bim = gwas_plink_stem + '.bim'
	ordered_rsids = load_in_ordered_array_of_rsids_from_bim_file(genotype_bim)

	# Load in genotype object
	genotype_obj = read_plink1_bin(gwas_plink_stem + '.bed', gwas_plink_stem + '.bim', gwas_plink_stem + '.fam', verbose=False)

	# Open output file handle
	t = open(gwas_output_file,'w')
	t.write('rsid\tbeta\tbeta_se\tz\n')

	# Initialize arrays to perform batch based analysis
	batch_rsids = []
	batch_variant_names = []

	snp_counter = 0
	# Loop through snps and run gwas for hm3 snps
	for snp_iter, snp_rsid in enumerate(ordered_rsids):
		# Skip snps that are not hapmap3 snps
		if snp_rsid not in hapmap3_rsids:
			continue
		snp_counter = snp_counter + 1
		batch_rsids.append(snp_rsid)
		batch_variant_names.append('variant' + str(snp_iter))

		# Only extract genotypes for variants in bach
		if np.mod(len(batch_rsids), batch_size) == 0:
			# Extract genotype data for this batch of snps
			batch_variant_genotype = np.asarray(genotype_obj.sel(variant=batch_variant_names))

			# Loop through variants in batch
			for batch_iter, batch_rsid in enumerate(np.asarray(batch_rsids)):
				# Extract genotype data for this snp
				variant_genotype = batch_variant_genotype[:,batch_iter]
				std_variant_genotype = np.copy(variant_genotype)
				# Mean impute genotype
				nan_indices = np.isnan(variant_genotype)
				non_nan_mean = np.mean(variant_genotype[nan_indices==False])
				std_variant_genotype[nan_indices] = non_nan_mean
				# Standardize genotype
				std_variant_genotype = (std_variant_genotype - np.mean(std_variant_genotype))/np.std(std_variant_genotype)


				# Now get effect size, standard erorr and z-score for association between standardized genotype and trait
				# Fit model using statsmodels
				mod = sm.OLS(trait_vector, std_variant_genotype)
				res = mod.fit()
				# Extract results
				effect_size = res.params[0]
				effect_size_se = res.bse[0]
				effect_size_z = effect_size/effect_size_se

				# Print to output file
				t.write(batch_rsid + '\t' + str(effect_size) + '\t' + str(effect_size_se) + '\t' + str(effect_size_z) + '\n')

			# Reset batch_rsids and batch_variant_names
			batch_rsids = []
			batch_variant_names = []

	# Print stragglers
	if len(batch_rsids) > 0:
		# Extract genotype data for this batch of snps
		batch_variant_genotype = np.asarray(genotype_obj.sel(variant=batch_variant_names))

		# Loop through variants in batch
		for batch_iter, batch_rsid in enumerate(np.asarray(batch_rsids)):
			# Extract genotype data for this snp
			variant_genotype = batch_variant_genotype[:,batch_iter]
			std_variant_genotype = np.copy(variant_genotype)
			# Mean impute genotype
			nan_indices = np.isnan(variant_genotype)
			non_nan_mean = np.mean(variant_genotype[nan_indices==False])
			std_variant_genotype[nan_indices] = non_nan_mean
			# Standardize genotype
			std_variant_genotype = (std_variant_genotype - np.mean(std_variant_genotype))/np.std(std_variant_genotype)


			# Now get effect size, standard erorr and z-score for association between standardized genotype and trait
			# Fit model using statsmodels
			mod = sm.OLS(trait_vector, std_variant_genotype)
			res = mod.fit()
			# Extract results
			effect_size = res.params[0]
			effect_size_se = res.bse[0]
			effect_size_z = effect_size/effect_size_se

			# Print to output file
			t.write(batch_rsid + '\t' + str(effect_size) + '\t' + str(effect_size_se) + '\t' + str(effect_size_z) + '\n')


	t.close()

	return


def run_gwas_for_a_single_gene(genotype_obj, batch_rsids, batch_rsids2, gwas_trait_vector, gene_gwas_output_file):
	# Open output file handle
	t = open(gene_gwas_output_file,'w')
	t.write('rsid\tbeta\tbeta_se\tz\n')

	batch_variant_genotype = np.asarray(genotype_obj.sel(variant=batch_rsids))

	# Loop through variants in batch
	for batch_iter, batch_rsid in enumerate(np.asarray(batch_rsids)):
		# Extract genotype data for this snp
		variant_genotype = batch_variant_genotype[:,batch_iter]
		std_variant_genotype = np.copy(variant_genotype)
		# Mean impute genotype
		nan_indices = np.isnan(variant_genotype)
		non_nan_mean = np.mean(variant_genotype[nan_indices==False])
		std_variant_genotype[nan_indices] = non_nan_mean
		# Standardize genotype
		std_variant_genotype = (std_variant_genotype - np.mean(std_variant_genotype))/np.std(std_variant_genotype)


		# Now get effect size, standard erorr and z-score for association between standardized genotype and trait
		# Fit model using statsmodels
		mod = sm.OLS(gwas_trait_vector, std_variant_genotype)
		res = mod.fit()
		# Extract results
		effect_size = res.params[0]
		effect_size_se = res.bse[0]
		effect_size_z = effect_size/effect_size_se

		# Print to output file
		t.write(batch_rsids2[batch_iter] + '\t' + str(effect_size) + '\t' + str(effect_size_se) + '\t' + str(effect_size_z) + '\n')
	t.close()
	return


def run_gwas_on_all_rsids(trait_values_file, gwas_plink_stem, gwas_output_file, batch_size=1000):
	# Load in trait vector
	trait_vector = np.loadtxt(trait_values_file)

	# Load in ordered array of rsids
	genotype_bim = gwas_plink_stem + '.bim'
	ordered_rsids = load_in_ordered_array_of_rsids_from_bim_file(genotype_bim)

	# Load in genotype object
	genotype_obj = read_plink1_bin(gwas_plink_stem + '.bed', gwas_plink_stem + '.bim', gwas_plink_stem + '.fam', verbose=False)

	# Open output file handle
	t = open(gwas_output_file,'w')
	t.write('rsid\tbeta\tbeta_se\tz\n')

	# Initialize arrays to perform batch based analysis
	batch_rsids = []
	batch_variant_names = []

	snp_counter = 0
	# Loop through snps and run gwas for hm3 snps
	for snp_iter, snp_rsid in enumerate(ordered_rsids):
		# Skip snps that are not hapmap3 snps
		snp_counter = snp_counter + 1
		batch_rsids.append(snp_rsid)
		batch_variant_names.append('variant' + str(snp_iter))

		# Only extract genotypes for variants in bach
		if np.mod(len(batch_rsids), batch_size) == 0:
			# Extract genotype data for this batch of snps
			batch_variant_genotype = np.asarray(genotype_obj.sel(variant=batch_variant_names))

			# Loop through variants in batch
			for batch_iter, batch_rsid in enumerate(np.asarray(batch_rsids)):
				# Extract genotype data for this snp
				variant_genotype = batch_variant_genotype[:,batch_iter]
				std_variant_genotype = np.copy(variant_genotype)
				# Mean impute genotype
				nan_indices = np.isnan(variant_genotype)
				non_nan_mean = np.mean(variant_genotype[nan_indices==False])
				std_variant_genotype[nan_indices] = non_nan_mean
				# Standardize genotype
				std_variant_genotype = (std_variant_genotype - np.mean(std_variant_genotype))/np.std(std_variant_genotype)


				# Now get effect size, standard erorr and z-score for association between standardized genotype and trait
				# Fit model using statsmodels
				mod = sm.OLS(trait_vector, std_variant_genotype)
				res = mod.fit()
				# Extract results
				effect_size = res.params[0]
				effect_size_se = res.bse[0]
				effect_size_z = effect_size/effect_size_se

				# Print to output file
				t.write(batch_rsid + '\t' + str(effect_size) + '\t' + str(effect_size_se) + '\t' + str(effect_size_z) + '\n')

			# Reset batch_rsids and batch_variant_names
			batch_rsids = []
			batch_variant_names = []

	# Print stragglers
	if len(batch_rsids) > 0:
		# Extract genotype data for this batch of snps
		batch_variant_genotype = np.asarray(genotype_obj.sel(variant=batch_variant_names))

		# Loop through variants in batch
		for batch_iter, batch_rsid in enumerate(np.asarray(batch_rsids)):
			# Extract genotype data for this snp
			variant_genotype = batch_variant_genotype[:,batch_iter]
			std_variant_genotype = np.copy(variant_genotype)
			# Mean impute genotype
			nan_indices = np.isnan(variant_genotype)
			non_nan_mean = np.mean(variant_genotype[nan_indices==False])
			std_variant_genotype[nan_indices] = non_nan_mean
			# Standardize genotype
			std_variant_genotype = (std_variant_genotype - np.mean(std_variant_genotype))/np.std(std_variant_genotype)


			# Now get effect size, standard erorr and z-score for association between standardized genotype and trait
			# Fit model using statsmodels
			mod = sm.OLS(trait_vector, std_variant_genotype)
			res = mod.fit()
			# Extract results
			effect_size = res.params[0]
			effect_size_se = res.bse[0]
			effect_size_z = effect_size/effect_size_se

			# Print to output file
			t.write(batch_rsid + '\t' + str(effect_size) + '\t' + str(effect_size_se) + '\t' + str(effect_size_z) + '\n')


	t.close()

	return


def extract_gwas_trait_vector(gene_trait_value_file):
	f = open(gene_trait_value_file)
	arr = []
	for line in f:
		line = line.rstrip()
		arr.append(float(line))
	f.close()

	return np.asarray(arr)



##############################
# Command line argumemnts
##############################
simulation_number = int(sys.argv[1])
chrom_num = sys.argv[2]
simulation_name_string = sys.argv[3]
processed_genotype_data_dir = sys.argv[4]
simulated_trait_dir = sys.argv[5]
simulated_gwas_dir = sys.argv[6]
simulated_gene_expression_dir = sys.argv[7]


####################################################
# Extract genotype data
####################################################

gwas_plink_stem = processed_genotype_data_dir + 'simulated_gwas_data_' + str(chrom_num)  # Genotype files

# Load in ordered array of rsids
genotype_bim = gwas_plink_stem + '.bim'
ordered_rsids = load_in_ordered_array_of_rsids_from_bim_file(genotype_bim)
ordered_variant_names = []
for ii,rsid in enumerate(ordered_rsids):
	ordered_variant_names.append('variant' + str(ii))
ordered_variant_names = np.asarray(ordered_variant_names)

# Load in genotype object
genotype_obj = read_plink1_bin(gwas_plink_stem + '.bed', gwas_plink_stem + '.bim', gwas_plink_stem + '.fam', verbose=False)



####################################################
# Run GWAS on all snps (in each gene seperately)
####################################################
simulated_expression_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'  # This file contains a line for each gene, and we will use it to select which genes in which tissue are used
f = open(simulated_expression_summary_file)
head_count = 0
gene_counter = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	print(gene_counter)
	# This corresponds to a single gene
	# Extract relevent fields for this gene
	gene_name = data[0]
	causal_eqtl_effect_file = data[3]
	eqtl_snp_id_file = data[4]
	eqtl_snp_indices_file = data[5]
	eqtl_snp_indices = np.load(eqtl_snp_indices_file)

	# Extract trait file for this gene
	gene_trait_value_file = simulated_trait_dir + simulation_name_string + '_expression_mediated_trait_values_' + gene_name + '.txt'  # Trait vector
	gwas_trait_vector = extract_gwas_trait_vector(gene_trait_value_file)

	# Run gwas for all rsids in this gene
	gene_gwas_output_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results_' + gene_name +'.txt'
	run_gwas_for_a_single_gene(genotype_obj, ordered_variant_names[eqtl_snp_indices], ordered_rsids[eqtl_snp_indices], gwas_trait_vector, gene_gwas_output_file)
	gene_counter = gene_counter + 1


f.close()





