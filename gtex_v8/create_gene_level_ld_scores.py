import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import os
import pdb
import numpy as np
from pandas_plink import read_plink1_bin
import pickle
import pandas as pd
import pyreadr
import gzip
import time




def get_pseudotissue_names(gtex_pseudotissue_file):
	f = open(gtex_pseudotissue_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0])
	f.close()
	return np.asarray(arr)



def create_tissue_to_gene_model_df(pseudotissue_name, gtex_susie_gene_models_dir, gene_type, chrom_num):
	gene_model_summary_file = gtex_susie_gene_models_dir + pseudotissue_name + '/' + pseudotissue_name + '_' + gene_type + '_pos_file.txt'
	raw_table = np.loadtxt(gene_model_summary_file, dtype=str, delimiter='\t')[1:]
	rows_on_chromosome = raw_table[:,2] == chrom_num
	raw_table_chrom = raw_table[rows_on_chromosome, :]
	return raw_table_chrom



def create_mapping_from_rsid_to_snpid(rsids, snp_ids):
	mapping = {}
	for ii, rsid in enumerate(rsids):
		if rsid in mapping:
			print('assumption eroror')
			pdb.set_trace()
		mapping[rsid] = snp_ids[ii]
	return mapping

def create_mapping_from_snpid_to_rsid(rsids, snp_ids):
	mapping = {}
	for ii, rsid in enumerate(rsids):
		if snp_ids[ii] in mapping:
			print('assumption eororor')
			pdb.set_trace()
		mapping[snp_ids[ii]] = rsid
	return mapping

def create_mapping_snpid_to_reference_index(snp_ids, alt_snp_ids):
	mapping = {}
	for ii, snpid in enumerate(snp_ids):
		mapping[snpid] = (ii, 1.0)
		alt_snpid = alt_snp_ids[ii]
		if alt_snpid == snpid:
			print('assumption eroror')
			pdb.set_trace()
		mapping[alt_snpid] = (ii, -1.0)

	return mapping

def get_gene_model_snp_positions(variant_names):
	positions = []
	for variant_name in variant_names:
		positions.append(variant_name.split('_')[1])
	return np.asarray(positions).astype(int)

def project_gene_model_onto_full_window(gene_susie_mu, gene_snpid_names, window_snpid_to_reference_index, n_window_snps, sign_aware_model=True):
	# First get number of components
	n_components = gene_susie_mu.shape[0]

	# Initialize gene model on full window
	window_susie_mu = np.zeros((n_components, n_window_snps))

	# Loop through gene snpid names
	for ii, gene_snpid_name in enumerate(gene_snpid_names):
		# Get window position and sign
		if gene_snpid_name not in window_snpid_to_reference_index:
			print('skip snp ' + gene_snpid_name)
			continue
		window_position, sign = window_snpid_to_reference_index[gene_snpid_name]

		if sign_aware_model:
			window_susie_mu[:, window_position] = (gene_susie_mu[:, ii])*sign
		else:
			window_susie_mu[:, window_position] = (gene_susie_mu[:, ii])

	return window_susie_mu

def get_window_regression_snp_indices(G_window_snpids, regression_snp_id_to_regression_snp_position):
	window_regression_snp_indices = []
	global_regression_snp_positions = []

	for snpid in G_window_snpids:
		if snpid in regression_snp_id_to_regression_snp_position:
			window_regression_snp_indices.append(True)
			global_regression_snp_positions.append(regression_snp_id_to_regression_snp_position[snpid])
		else:
			window_regression_snp_indices.append(False)

	# Put into nice array format
	window_regression_snp_indices = np.asarray(window_regression_snp_indices)
	global_regression_snp_positions = np.asarray(global_regression_snp_positions)

	# Quick error checking
	if np.sum(window_regression_snp_indices) != len(global_regression_snp_positions):
		print('assumption eroror')
		pdb.set_trace()

	return window_regression_snp_indices, global_regression_snp_positions

def compute_gene_variance(susie_mu, susie_mu_sd, susie_alpha, ld):
	gene_var = 0.0

	# Component level eqtl effect sizes for this gene		
	gene_component_effect_sizes = (susie_mu)*susie_alpha

	# eQTL effect sizes for this gene
	gene_eqtl_effect_sizes = np.sum(gene_component_effect_sizes,axis=0)


	num_susie_components = susie_mu.shape[0]
	for k_index in range(num_susie_components):
		gene_var = gene_var + np.sum((np.square(susie_mu[k_index,:]) + np.square(susie_mu_sd[k_index,:]))*np.diag(ld)*susie_alpha[k_index,:])
		eqtl_component_pmces = (susie_mu[k_index,:])*(susie_alpha[k_index,:])
		gene_var = gene_var - np.dot(np.dot(eqtl_component_pmces,ld), eqtl_component_pmces)
	gene_var = gene_var + np.dot(np.dot(gene_eqtl_effect_sizes,ld), gene_eqtl_effect_sizes)
				
	return gene_var


def extract_ld_annotations_for_this_gene_region(ld, susie_mu, susie_mu_sd, susie_alpha, variant_indices, n_ref_panel_samples):
	# Number of variants
	num_var = ld.shape[0]

	# Component level eqtl effect sizes for this gene		
	gene_component_effect_sizes = susie_mu*susie_alpha

	# eQTL effect sizes for this gene
	gene_eqtl_effect_sizes = np.sum(gene_component_effect_sizes,axis=0)

	# Compute squared eqtl effect sizes for this gene
	gene_squared_eqtl_effect_sizes = np.sum((np.square(susie_mu) + np.square(susie_mu_sd))*susie_alpha,axis=0) + gene_eqtl_effect_sizes*gene_eqtl_effect_sizes - np.sum(gene_component_effect_sizes*gene_component_effect_sizes,axis=0)

	# E[beta_k*beta_j]
	cross_terms = np.dot(np.reshape(gene_eqtl_effect_sizes, (num_var,1)), np.reshape(gene_eqtl_effect_sizes, (1,num_var))) - np.dot(np.transpose(gene_component_effect_sizes), gene_component_effect_sizes)

	# Compute gene variance
	gene_variance = compute_gene_variance(susie_mu, susie_mu_sd, susie_alpha, ld)

	# Compute ld scores for diagonal piece
	diagonal_ld_scores = np.sum((np.square(ld[variant_indices,:])*gene_squared_eqtl_effect_sizes),axis=1)
	
	# Comptute ld scores for off diagonal elements
	np.fill_diagonal(cross_terms, np.zeros(num_var))
	non_diagonal_ld_scores = np.sum(np.dot(cross_terms, ld[:, variant_indices])*ld[:,variant_indices],axis=0)

	# Generate complete ld scores
	ld_scores = (diagonal_ld_scores + non_diagonal_ld_scores)/gene_variance

	# Get adjusted ld scores
	adj_ld_scores = ld_scores - ((1.0-ld_scores)/(n_ref_panel_samples-2.0))

	return ld_scores, adj_ld_scores

def extract_ld_annotations_for_this_gene_region_with_eqtl_point_estimate(ld, gene_eqtl_effect_sizes, variant_indices, n_ref_panel_samples, version='B'):
	gene_variance = np.dot(np.dot(gene_eqtl_effect_sizes, ld), gene_eqtl_effect_sizes)

	if version == 'A':
		# Number of variants
		num_var = ld.shape[0]

		# Compute squared eqtl effect sizes for this gene
		gene_squared_eqtl_effect_sizes = np.square(gene_eqtl_effect_sizes)

		# E[beta_k*beta_j]
		cross_terms = np.dot(np.reshape(gene_eqtl_effect_sizes, (num_var,1)), np.reshape(gene_eqtl_effect_sizes, (1,num_var)))


		diagonal_ld_scores = np.sum((np.square(ld[variant_indices,:])*gene_squared_eqtl_effect_sizes),axis=1)

		# Temp
		np.fill_diagonal(cross_terms, np.zeros(num_var))
		non_diagonal_ld_scores = np.sum(np.dot(cross_terms, ld[:, variant_indices])*ld[:,variant_indices],axis=0)

		# Generate complete ld scores
		ld_scores = (diagonal_ld_scores + non_diagonal_ld_scores)/gene_variance
	elif version == 'B':
		standardized_gene_eqtl_effect_sizes = gene_eqtl_effect_sizes/np.sqrt(gene_variance)
		ld_scores = np.square(np.dot(ld[variant_indices,:], standardized_gene_eqtl_effect_sizes))

	# Get adjusted ld scores
	adj_ld_scores = ld_scores - ((1.0-ld_scores)/(n_ref_panel_samples-2.0))

	return ld_scores, adj_ld_scores

def extract_gene_regression_weights_for_this_gene_region_with_eqtl_point_estimate(ld, gene_eqtl_effect_sizes, variant_indices, n_ref_panel_samples):
	gene_variance = np.dot(np.dot(gene_eqtl_effect_sizes, ld), gene_eqtl_effect_sizes)
	standardized_gene_eqtl_effect_sizes = gene_eqtl_effect_sizes/np.sqrt(gene_variance)

	standardized_gene_eqtl_effect_sizes_shaped = np.reshape(standardized_gene_eqtl_effect_sizes, (len(standardized_gene_eqtl_effect_sizes), 1))

	beta_beta_t = np.dot(standardized_gene_eqtl_effect_sizes_shaped, np.transpose(standardized_gene_eqtl_effect_sizes_shaped))

	eqtl_cov = np.dot(np.dot(ld[variant_indices,:], beta_beta_t), ld[:,variant_indices])

	# Remove diagonal elements
	np.fill_diagonal(eqtl_cov,0.0)

	# Square and sum
	np.sum(np.square(eqtl_cov), axis=0)

	weights = np.sum(np.square(eqtl_cov), axis=1)

	return weights

def extract_ld_annotations_for_this_gene_region_with_sample_correlation(pmces, sample_geno, variant_indices, n_ref_panel_samples):
	stand_sample_geno = np.copy(sample_geno)
	for kk in range(stand_sample_geno.shape[1]):
		stand_sample_geno[:, kk] = (sample_geno[:, kk] - np.mean(sample_geno[:,kk]))/np.std(sample_geno[:,kk])

	pred_expr = np.dot(stand_sample_geno, pmces)

	full = np.hstack((np.reshape(pred_expr, (len(pred_expr), 1)), stand_sample_geno[:, variant_indices]))

	tmper = np.corrcoef(np.transpose(full))

	ld_scores = np.square(tmper[0,1:])
	adj_ld_scores = ld_scores - ((1.0-ld_scores)/(n_ref_panel_samples-2.0))
	return ld_scores, adj_ld_scores



def compute_gene_level_ld_scores_for_single_gene(gene_tissue_model, ref_genotype_obj, snpid_to_reference_index, regression_snp_id_to_regression_snp_position, global_gene_level_pmces_ld_scores, global_gene_level_pmces_adj_ld_scores, global_gene_level_pmces_weights):
	# First get gene model ub and lb snp name
	gene_snpid_names = np.asarray(gene_tissue_model['variant_names'])[:,0]
	gene_snp_positions = get_gene_model_snp_positions(gene_snpid_names)
	gene_snpid_lb = gene_snpid_names[np.argmin(gene_snp_positions)]
	gene_snpid_ub = gene_snpid_names[np.argmax(gene_snp_positions)]

	# Now get ref positions of gene snp id lb and ub
	ref_pos_gene_snp_lb = snpid_to_reference_index[gene_snpid_lb][0]
	ref_pos_gene_snp_ub = snpid_to_reference_index[gene_snpid_ub][0]

	# Gene window cm lb and ub
	gene_window_cm_lb = ref_genotype_obj['cm'][ref_pos_gene_snp_lb] - 1.0
	gene_window_cm_ub = ref_genotype_obj['cm'][ref_pos_gene_snp_ub] + 1.0

	# Get indices corresponding to ref genotype for this gene window
	ref_genotype_indices_for_window = (ref_genotype_obj['cm'] >= gene_window_cm_lb) & (ref_genotype_obj['cm'] < gene_window_cm_ub)

	# Subset ref genotype data to only snps in ref_genotype_indices_for_window
	G_window_geno = ref_genotype_obj['G'][:, ref_genotype_indices_for_window]
	G_window_ld = np.corrcoef(np.transpose(G_window_geno))
	G_window_rsids = ref_genotype_obj['rsid'][ref_genotype_indices_for_window]
	G_window_snpids = ref_genotype_obj['snp_id'][ref_genotype_indices_for_window]
	G_window_alt_snpids = ref_genotype_obj['alt_snp_id'][ref_genotype_indices_for_window]

	# Get number of window snps
	n_window_snps = len(G_window_alt_snpids)
	# Number of reference panel samples
	n_ref_panel_samples = G_window_geno.shape[0]

	# Get indices in window corresponding to regression snps
	window_regression_snp_indices, global_regression_snp_positions = get_window_regression_snp_indices(G_window_snpids, regression_snp_id_to_regression_snp_position)

	# Calculate number of regression snps in this window
	n_regression_snps_in_window = np.sum(window_regression_snp_indices)

	# If no regression snps, can return here
	if n_regression_snps_in_window == 0:
		print('SKIP WINDOW')
		return global_gene_level_pmces_ld_scores, global_gene_level_pmces_adj_ld_scores, global_gene_level_pmces_weights

	# Create mapping from snp id to reference snp position in this window
	window_snpid_to_reference_index = create_mapping_snpid_to_reference_index(G_window_snpids, G_window_alt_snpids)

	# Project gene models onto full window
	window_susie_mu = project_gene_model_onto_full_window(np.asarray(gene_tissue_model['susie_mu']), gene_snpid_names, window_snpid_to_reference_index, n_window_snps, sign_aware_model=True)
	window_susie_alpha = project_gene_model_onto_full_window(np.asarray(gene_tissue_model['susie_alpha']), gene_snpid_names, window_snpid_to_reference_index, n_window_snps, sign_aware_model=False)
	local_susie_mu_sd = np.sqrt(np.asarray(gene_tissue_model['susie_mu2']) - np.square(np.asarray(gene_tissue_model['susie_mu'])))
	window_susie_mu_sd = project_gene_model_onto_full_window(local_susie_mu_sd, gene_snpid_names, window_snpid_to_reference_index, n_window_snps, sign_aware_model=False)

	# Extract LD scores using susie distribution of eqtl effect sizes
	#window_gene_ld_scores, window_adj_ld_scores = extract_ld_annotations_for_this_gene_region(G_window_ld, window_susie_mu, window_susie_mu_sd, window_susie_alpha, window_regression_snp_indices, n_ref_panel_samples)
	# Extract LD scores using point estimate eqtl effect sizes
	gene_pmces = np.sum(window_susie_mu*window_susie_alpha,axis=0)
	window_gene_ld_scores_pe, window_adj_ld_scores_pe = extract_ld_annotations_for_this_gene_region_with_eqtl_point_estimate(G_window_ld, gene_pmces, window_regression_snp_indices, n_ref_panel_samples, version='B')
	# Extract LD scores by correlating predicted expression with observed genotype (Identical to the above)
	#window_gene_ld_scores_corr, window_adj_ld_scores_corr = extract_ld_annotations_for_this_gene_region_with_sample_correlation(np.sum(window_susie_mu*window_susie_alpha,axis=0), G_window_geno,  window_regression_snp_indices, n_ref_panel_samples)

	# Extract gene regression weights
	#window_gene_regression_weights_pe = extract_gene_regression_weights_for_this_gene_region_with_eqtl_point_estimate(G_window_ld, gene_pmces, window_regression_snp_indices, n_ref_panel_samples)

	# Add updates to global array
	#global_gene_level_ld_scores[global_regression_snp_positions] = global_gene_level_ld_scores[global_regression_snp_positions] + window_gene_ld_scores
	#global_gene_level_adj_ld_scores[global_regression_snp_positions] = global_gene_level_adj_ld_scores[global_regression_snp_positions] + window_adj_ld_scores
	global_gene_level_pmces_ld_scores[global_regression_snp_positions] = global_gene_level_pmces_ld_scores[global_regression_snp_positions] + window_gene_ld_scores_pe
	global_gene_level_pmces_adj_ld_scores[global_regression_snp_positions] = global_gene_level_pmces_adj_ld_scores[global_regression_snp_positions] + window_adj_ld_scores_pe
	#global_gene_level_pmces_weights[global_regression_snp_positions] = global_gene_level_pmces_weights[global_regression_snp_positions] + window_gene_regression_weights_pe

	return global_gene_level_pmces_ld_scores, global_gene_level_pmces_adj_ld_scores, global_gene_level_pmces_weights


def load_in_regression_snp_ids(variant_level_ld_score_file, rsid_to_snpid, snpid_to_reference_index):
	rsids = []
	snpids = []
	snp_to_regression_snp_position = {}
	f = gzip.open(variant_level_ld_score_file)
	head_count = 0
	snp_counter = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[1]
		# This filter occurs because of weird filter at the beginnong on removing repeat snps (only allow this to happen once per genome)
		snp_id = rsid_to_snpid[rsid]
		snp_id_info = snp_id.split('_')
		snp_id_alt = snp_id_info[0] + '_' + snp_id_info[1] + '_' + snp_id_info[3] + '_' + snp_id_info[2] + '_' + snp_id_info[4]
		# Quick error check
		if snp_id not in snpid_to_reference_index:
			print('assumption eroror')
			pdb.set_trace()
	
		# Add to arrays
		rsids.append(rsid)
		snpids.append(snp_id)
		if snp_id in snp_to_regression_snp_position or snp_id_alt in snp_to_regression_snp_position:
			print('assumption eroror')
			pdb.set_trace()
		snp_to_regression_snp_position[snp_id] = snp_counter
		snp_to_regression_snp_position[snp_id_alt] = snp_counter
		snp_counter = snp_counter + 1
	f.close()

	return np.asarray(rsids), np.asarray(snpids), snp_to_regression_snp_position


def get_non_repeat_columns_from_plink_obj(G_obj, variant_level_ld_score_file):
	f = gzip.open(variant_level_ld_score_file)
	head_count = 0
	hm3_rs_ids = {}
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[1]
		if rsid in hm3_rs_ids:
			print('assumption erorro')
			pdb.set_trace()
		hm3_rs_ids[rsid] = 1
	f.close()

	G_obj_chrom = np.asarray(G_obj.chrom)
	G_obj_pos = np.asarray(G_obj.pos)
	G_obj_a0 = np.asarray(G_obj.a0)
	G_obj_a1 = np.asarray(G_obj.a1)
	G_obj_rs = np.asarray(G_obj.snp)

	n_var = len(G_obj_pos)
	used = {}
	valid_columns = []

	# Go through hm3 snps first
	for ii in range(n_var):
		snp_id1 = 'chr' + G_obj_chrom[ii] + '_' + str(G_obj_pos[ii]) + '_' + G_obj_a0[ii] + '_' + G_obj_a1[ii]
		snp_id2 = 'chr' + G_obj_chrom[ii] + '_' + str(G_obj_pos[ii]) + '_' + G_obj_a1[ii] + '_' + G_obj_a0[ii]
		rsid = G_obj_rs[ii]
		if rsid not in hm3_rs_ids:
			continue
		if snp_id1 in used or snp_id2 in used:
			print('assumption error')
			pdb.set_trace()
		used[snp_id1] = 1
		used[snp_id2] = 1
		valid_columns.append(ii)
	for ii in range(n_var):
		snp_id1 = 'chr' + G_obj_chrom[ii] + '_' + str(G_obj_pos[ii]) + '_' + G_obj_a0[ii] + '_' + G_obj_a1[ii]
		snp_id2 = 'chr' + G_obj_chrom[ii] + '_' + str(G_obj_pos[ii]) + '_' + G_obj_a1[ii] + '_' + G_obj_a0[ii]
		rsid = G_obj_rs[ii]
		if rsid in hm3_rs_ids:
			continue
		if snp_id1 in used or snp_id2 in used:
			continue
		used[snp_id1] = 1
		used[snp_id2] = 1
		valid_columns.append(ii)
	valid_columns = np.sort(np.asarray(valid_columns))
	return valid_columns


variant_level_ld_score_file = sys.argv[1]
genotype_stem = sys.argv[2]
pseudotissue_name = sys.argv[3]
gtex_susie_gene_models_dir = sys.argv[4]
chrom_num = sys.argv[5]
gene_type = sys.argv[6]
gene_level_sldsc_output_root = sys.argv[7]

print('CHROM' + str(chrom_num))


# Create dictionary from pseudotissue to array of gene models for that tissue
tissue_to_gene_model_df = create_tissue_to_gene_model_df(pseudotissue_name, gtex_susie_gene_models_dir, gene_type, chrom_num)

# Number of genes on this chromosome in this tissue
n_genes = tissue_to_gene_model_df.shape[0]

# Load in Reference Genotype data
G_obj = read_plink1_bin(genotype_stem + '.bed', genotype_stem + '.bim', genotype_stem + '.fam', verbose=False)
# For some reason, this file has some repeats. First get columns that don't correspond to repeats
valid_variant_columns = get_non_repeat_columns_from_plink_obj(G_obj, variant_level_ld_score_file)
G_obj = G_obj[:, valid_variant_columns]

G_obj_geno = G_obj.values # Numpy 2d array of dimension num samples X num snps
G_obj_chrom = np.asarray(G_obj.chrom)
G_obj_pos = np.asarray(G_obj.pos)
# For our purposes, a0 is the effect allele
# For case of plink package, a0 is the first column in the plink bim file
G_obj_a0 = np.asarray(G_obj.a0)
G_obj_a1 = np.asarray(G_obj.a1)
# RSids
G_obj_rsids = np.asarray(G_obj.snp)
# Centimorgan distances
G_obj_cm = np.asarray(G_obj.cm)
# Snp ids
G_obj_snp_ids = 'chr' + G_obj_chrom + '_' + (G_obj_pos.astype(str)) + '_' + G_obj_a0 + '_' + G_obj_a1 + '_b38'
G_obj_alt_snp_ids = 'chr' + G_obj_chrom + '_' + (G_obj_pos.astype(str)) + '_' + G_obj_a1 + '_' + G_obj_a0 + '_b38'

# Put geno into organized dictionary
genotype_obj = {'G': G_obj_geno, 'rsid': G_obj_rsids, 'snp_id': G_obj_snp_ids, 'alt_snp_id': G_obj_alt_snp_ids, 'position': G_obj_pos, 'cm': G_obj_cm}

# Quick error checks on genotype data
if len(np.unique(G_obj_rsids)) != len(G_obj_rsids):
	print('assumption eroror')
	pdb.set_trace()
if len(np.unique(G_obj_snp_ids)) != len(G_obj_snp_ids):
	print('assumption error')
	pdb.set_trace()


# Create mapping from rsid to snpid and snpid to rsid
rsid_to_snpid = create_mapping_from_rsid_to_snpid(genotype_obj['rsid'], genotype_obj['snp_id'])
snpid_to_rsid = create_mapping_from_snpid_to_rsid(genotype_obj['rsid'], genotype_obj['snp_id'])

# Create mapping from snp id to reference index
snpid_to_reference_index = create_mapping_snpid_to_reference_index(genotype_obj['snp_id'], genotype_obj['alt_snp_id'])

# Get regression rsids
regression_rsids, regression_snpids, regression_snp_id_to_regression_snp_position = load_in_regression_snp_ids(variant_level_ld_score_file, rsid_to_snpid, snpid_to_reference_index)
n_regression_snps = len(regression_rsids)

# Initialize vector to keep track of gene level ld scores for regression snps on this chromosome
#global_gene_level_ld_scores = np.zeros(n_regression_snps)
#global_gene_level_adj_ld_scores = np.zeros(n_regression_snps)
global_gene_level_pmces_ld_scores = np.zeros(n_regression_snps)
global_gene_level_pmces_adj_ld_scores = np.zeros(n_regression_snps)
global_gene_level_pmces_weights = np.zeros(n_regression_snps)


start_time = time.time()
# Loop through genes
print(n_genes)
for g_counter in range(n_genes):
	print(g_counter)
	gene_info_vec = tissue_to_gene_model_df[g_counter,:]
	rdat_gene_model_file = gene_info_vec[0]
	# Load gene-tissue model from gene-tissue weight file
	gene_tissue_model = pyreadr.read_r(rdat_gene_model_file)

	#global_gene_level_ld_scores, global_gene_level_adj_ld_scores, global_gene_level_pmces_ld_scores, global_gene_level_pmces_adj_ld_scores, global_gene_level_pmces_weights = compute_gene_level_ld_scores_for_single_gene(gene_tissue_model, genotype_obj, snpid_to_reference_index, regression_snp_id_to_regression_snp_position, global_gene_level_ld_scores, global_gene_level_adj_ld_scores, global_gene_level_pmces_ld_scores, global_gene_level_pmces_adj_ld_scores, global_gene_level_pmces_weights)
	global_gene_level_pmces_ld_scores, global_gene_level_pmces_adj_ld_scores, global_gene_level_pmces_weights = compute_gene_level_ld_scores_for_single_gene(gene_tissue_model, genotype_obj, snpid_to_reference_index, regression_snp_id_to_regression_snp_position, global_gene_level_pmces_ld_scores, global_gene_level_pmces_adj_ld_scores, global_gene_level_pmces_weights)

	end_time = time.time()
	#print(end_time-start_time)
	start_time = end_time

# Save results to output file
#output_file1 = gene_level_sldsc_output_root + pseudotissue_name + '_' + gene_type + '_gene_ld_scores'
#np.savetxt(output_file1, global_gene_level_ld_scores, fmt="%s")

#output_file2 = gene_level_sldsc_output_root + pseudotissue_name + '_' + gene_type + '_gene_adj_ld_scores'
#np.savetxt(output_file2, global_gene_level_adj_ld_scores, fmt="%s")

output_file3 = gene_level_sldsc_output_root + pseudotissue_name + '_' + gene_type + '_pmces_gene_ld_scores'
np.savetxt(output_file3, global_gene_level_pmces_ld_scores, fmt="%s")

output_file4 = gene_level_sldsc_output_root + pseudotissue_name + '_' + gene_type + '_pmces_gene_adj_ld_scores'
np.savetxt(output_file4, global_gene_level_pmces_adj_ld_scores, fmt="%s")

#output_file5 = gene_level_sldsc_output_root + pseudotissue_name + '_' + gene_type + '_pmces_gene_weights'
#np.savetxt(output_file5, global_gene_level_pmces_weights, fmt="%s")





