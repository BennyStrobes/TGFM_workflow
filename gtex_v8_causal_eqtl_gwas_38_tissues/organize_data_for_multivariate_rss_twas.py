import numpy as np 
import os
import sys
import pdb
import pickle

def extract_global_variant_arr(bim_file_names):
	variants = {}
	for bim_file_name in bim_file_names:
		g = open(bim_file_name)
		for line in g:
			line = line.rstrip()
			data = line.split('\t')
			variants[data[1]] = 1
		g.close()
	return np.sort([*variants])

def create_mapping_from_variant_name_to_global_variant_arr_pos(global_variant_arr):
	dicti = {}
	for i, val in enumerate(global_variant_arr):
		dicti[val] = i
	return dicti

def create_gene_local_to_global_mapping(variant_name_to_global_variant_arr_pos, bim_file_names):
	gene_local_to_global_mapping = []
	for bim_file_name in bim_file_names:
		arr = []
		g = open(bim_file_name)
		for line in g:
			line = line.rstrip()
			data = line.split('\t')
			variant_id = data[1]
			global_pos = variant_name_to_global_variant_arr_pos[variant_id]
			arr.append(global_pos)
		g.close()	
		gene_local_to_global_mapping.append(np.asarray(arr))
	return gene_local_to_global_mapping

def create_global_bim_data(gene_local_to_global_mapping, bim_file_names, num_global_variants, global_variant_arr):
	# String for initialization purposes
	init_string = 'NAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'

	# Initialize global bim data
	global_bim_data = np.reshape([init_string]*(6*num_global_variants), (num_global_variants, 6))

	for bim_counter, bim_file_name in enumerate(bim_file_names):
		# Get local bim data
		bim_data = np.loadtxt(bim_file_name, dtype=str, delimiter='\t')
		# Fill in global bim data with local bim data
		global_bim_data[gene_local_to_global_mapping[bim_counter],:] = bim_data

	# Quick error checking
	if sum(sum(global_bim_data == init_string)) != 0:
		print('assumption error creating global bim data')
		pdb.set_trace()
	# More error checking
	if np.array_equal(global_bim_data[:,1], global_variant_arr) != True:
		print('assumption error creating global bim data')
		pdb.set_trace()
	return global_bim_data

def create_global_geno_data(gene_local_to_global_mapping, geno_file_names, num_global_variants):
	# Initialize global geno data
	global_geno_data = np.zeros((489, num_global_variants))

	for geno_counter, geno_file_name in enumerate(geno_file_names):
		# Get local geno data
		local_geno_data = np.loadtxt(geno_file_name)

		# Fill in global geno data with local geno data
		global_geno_data[:, gene_local_to_global_mapping[geno_counter]] = local_geno_data
	
	# Simple error checking
	if np.min(np.var(global_geno_data,axis=0)) < .98 or np.max(np.var(global_geno_data,axis=0)) > 1.02:
		print('geno assumption eroror')
		pdb.set_trace()

	return global_geno_data

def create_global_susie_data(gene_local_to_global_mapping, eqtl_susie_mu_file_names, num_global_variants):
	# initialize global susie data
	global_susie_data = []

	for susie_counter, susie_file_name in enumerate(eqtl_susie_mu_file_names):
		# Initialize global susie data for this gene
		global_susie_data_mat = np.zeros((10, num_global_variants))
		# get local data
		local_susie_data = np.loadtxt(susie_file_name)

		# for this gene, Fill in global susie data with local susie data
		global_susie_data_mat[:, gene_local_to_global_mapping[susie_counter]] = local_susie_data

		# Now add to global arr (seperate element for each gene)
		global_susie_data.append(global_susie_data_mat)
	return global_susie_data

def create_global_gwas_data(gene_local_to_global_mapping, gwas_beta_file_names, num_global_variants):
	# Initialize global data
	global_gwas_data = np.zeros(num_global_variants)

	for gwas_counter, gwas_file_name in enumerate(gwas_beta_file_names):
		local_gwas_data = np.loadtxt(gwas_file_name)

		global_gwas_data[gene_local_to_global_mapping[gwas_counter]] = local_gwas_data

	return global_gwas_data

def create_global_fusion_weights_data(gene_local_to_global_mapping, gwas_beta_file_names, num_global_variants):
	# Initialize global data
	global_weight_data = []

	for gwas_counter, gwas_file_name in enumerate(gwas_beta_file_names):
		local_gwas_data = np.loadtxt(gwas_file_name)

		global_weight_arr = np.zeros(num_global_variants)
		global_weight_arr[gene_local_to_global_mapping[gwas_counter]] = local_gwas_data

		global_weight_data.append(global_weight_arr)

	return global_weight_data


def organize_single_components_data(component_name, trait_name, component_gene_file, gwas_sample_size):
	# load in data
	raw_data = np.loadtxt(component_gene_file,dtype=str, delimiter='\t')
	gene_data_header = raw_data[0,:]
	gene_data = raw_data[1:, :]

	# Extract relevent fields from data
	gene_names = gene_data[:,0]
	tissue_names = gene_data[:,1]
	geno_file_names = gene_data[:,8]
	bim_file_names = gene_data[:, 9]
	eqtl_susie_mu_file_names = gene_data[:,10]
	eqtl_susie_alpha_file_names = gene_data[:, 11]
	eqtl_susie_mu_sd_file_names = gene_data[:, 12]
	gwas_beta_file_names = gene_data[:, 13]
	gwas_beta_se_file_names = gene_data[:, 14]
	fusion_weight_file_names = gene_data[:,15]

	# create full gene names
	full_gene_names = []
	for index in range(len(gene_names)):
		full_gene_name = gene_names[index] + '_' + tissue_names[index]
		full_gene_names.append(full_gene_name)
	full_gene_names = np.asarray(full_gene_names)

	# Extract number of genes (this is really gene-tissue pairs. but who cares amirite)
	num_genes = len(gene_names)

	# extract global list of variants (aggregated across all genes)
	global_variant_arr = extract_global_variant_arr(bim_file_names)
	num_global_variants = len(global_variant_arr)

	# Create dictionary mapping from variant name to global-variant arr position
	variant_name_to_global_variant_arr_pos = create_mapping_from_variant_name_to_global_variant_arr_pos(global_variant_arr)

	# For each gene seperately create array of indexes corresponding to global variant arr
	gene_local_to_global_mapping = create_gene_local_to_global_mapping(variant_name_to_global_variant_arr_pos, bim_file_names)

	# create global bim data
	global_bim_data = create_global_bim_data(gene_local_to_global_mapping, bim_file_names, num_global_variants, global_variant_arr)

	# Create global geno data
	global_geno_data = create_global_geno_data(gene_local_to_global_mapping, geno_file_names, num_global_variants)
	# create global LD
	global_ld = np.corrcoef(np.transpose(global_geno_data))
	#global_ld = np.transpose(global_geno_data)

	# Create global susie mu data
	global_susie_mu_data = create_global_susie_data(gene_local_to_global_mapping, eqtl_susie_mu_file_names, num_global_variants)

	# Create global susie alpha data
	global_susie_alpha_data = create_global_susie_data(gene_local_to_global_mapping, eqtl_susie_alpha_file_names, num_global_variants)

	# Create global susie mu sd data
	global_susie_mu_sd_data = create_global_susie_data(gene_local_to_global_mapping, eqtl_susie_mu_sd_file_names, num_global_variants)

	# Create global gwas beta data
	global_gwas_beta_data = create_global_gwas_data(gene_local_to_global_mapping, gwas_beta_file_names, num_global_variants)

	# Create global gwas beta sd data
	global_gwas_beta_se_data = create_global_gwas_data(gene_local_to_global_mapping, gwas_beta_se_file_names, num_global_variants)	

	# Create global fusion_weights
	global_fusion_weights_data = create_global_fusion_weights_data(gene_local_to_global_mapping, fusion_weight_file_names, num_global_variants)	

	# Place all global data into dictionary
	global_dictionary = {'genes': full_gene_names, 'variants': global_variant_arr, 'bim': global_bim_data, 'reference_ld':global_ld, 'susie_mu':global_susie_mu_data, 'susie_alpha': global_susie_alpha_data, 'susie_mu_sd': global_susie_mu_sd_data, 'fusion_weights': global_fusion_weights_data, 'gwas_beta': global_gwas_beta_data, 'gwas_beta_se': global_gwas_beta_se_data, 'gwas_sample_size': gwas_sample_size}


	# REMOVE data
	for itera in range(len(geno_file_names)):
		os.system('rm ' + geno_file_names[itera])
		os.system('rm ' + bim_file_names[itera])
		os.system('rm ' + eqtl_susie_mu_file_names[itera])
		os.system('rm ' + eqtl_susie_alpha_file_names[itera])
		os.system('rm ' + eqtl_susie_mu_sd_file_names[itera])
		os.system('rm ' + gwas_beta_file_names[itera])
		os.system('rm ' + gwas_beta_se_file_names[itera])
		os.system('rm ' + fusion_weight_file_names[itera])

	return global_dictionary

chrom_num = sys.argv[1]
trait_name = sys.argv[2]
gtex_pseudotissue_file = sys.argv[3]
output_dir = sys.argv[4]
gwas_sample_size = int(sys.argv[5])
gene_version = sys.argv[6]


# each line in this file is a disease component (for this trait) on this chromosome 
# Each line also contains all file containing info on all the genes in cis to this component
input_file = output_dir + trait_name + '_' + gene_version + '_' + chrom_num + '_component_rss_multivariate_twas_data.txt'

# Create output file
output_file = output_dir + trait_name + '_' + gene_version + '_' + chrom_num + '_component_rss_multivariate_twas_data_organized.txt'
t = open(output_file,'w')
# Header to output file
t.write('chrom_num\tcomponent_name\tnumber_of_cis_genes\tcomponent_twas_data\n')

# Weird (rare bug) where susie identifies 2 components with same lead snp (SHOULD REALLY BE FIXED EARLIER ON)
used_components = {}

f = open(input_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	component_name = data[2] + '_' + data[4]
	component_gene_file = data[7]
	num_genes = int(data[8])

	if component_name in used_components:
		print('SKIPPING COMPONENT BECAUSE ALREADY SEEN')
		continue
	used_components[component_name] = 1

	if num_genes > 0:
		component_data = organize_single_components_data(component_name, trait_name, component_gene_file, gwas_sample_size)
		pkl_file = output_dir + trait_name + '_' + component_name + '_' + gene_version + '_multivariate_twas_data.pkl'
		g = open(pkl_file, "wb")
		pickle.dump(component_data, g)
		g.close()
		t.write(data[0] + '\t' + component_name + '\t' + str(num_genes) + '\t' + pkl_file + '\n')
	else:
		t.write(data[0] + '\t' + component_name + '\t' + str(num_genes) + '\t' + 'NA' + '\n')

	# gene_num=1
	# gwas_z = component_data['gwas_beta']/component_data['gwas_beta_se']
	# LD = component_data['reference_ld']
	# weights = np.sum(component_data['susie_mu'][gene_num]*component_data['susie_alpha'][gene_num],axis=0)  # Extract susie PMCES for this gene
	# twas_z = np.dot(weights, gwas_z)/np.sqrt(np.dot(np.dot(weights, LD), weights))

f.close()
t.close()