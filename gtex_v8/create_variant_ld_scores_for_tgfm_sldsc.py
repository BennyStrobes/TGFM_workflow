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

def fill_in_variant_level_ld_scores_using_sliding_window_approach(genotype_obj, regression_snp_id_to_regression_snp_position, window_size):
	# Initialize dictionaries to keep track of data
	regression_snpids_to_nearby_ld_vec = {}
	regression_snp_id_to_nearby_snp_ids = {}

	# Organize cm window info
	min_cm = np.min(genotype_obj['cm'])
	max_cm = np.max(genotype_obj['cm'])
	start_cm = min_cm - (2.0*window_size)

	# Initialize vector of reference snp positions
	ref_snp_positions = np.arange(len(genotype_obj['cm']))

	# Loop through windows
	while start_cm < max_cm:
		print(start_cm)
		# Define window positions
		mid_left_cm = start_cm + window_size
		mid_right_cm = start_cm + 2.0*window_size
		end_cm = start_cm + 3.0*window_size
	
		# Get indices corresponding to ref genotype for this gene window
		ref_genotype_indices_for_window = (genotype_obj['cm'] >= start_cm) & (genotype_obj['cm'] < end_cm)

		# Subset ref genotype data to only snps in ref_genotype_indices_for_window
		G_window_geno = genotype_obj['G'][:, ref_genotype_indices_for_window]
		n_ref_samp = G_window_geno.shape[0]
		G_window_ld = np.corrcoef(np.transpose(G_window_geno))

		G_window_snps = genotype_obj['snp_id'][ref_genotype_indices_for_window]
		G_window_alt_snps = genotype_obj['alt_snp_id'][ref_genotype_indices_for_window]
		G_window_cms = genotype_obj['cm'][ref_genotype_indices_for_window]
		G_window_positions = ref_snp_positions[ref_genotype_indices_for_window]

		# Now loop through snps in window and find regression snps in middle of the window
		n_window_snps = len(G_window_snps)
		for snp_iter in range(n_window_snps):
			# In middle of window test
			snp_cm = G_window_cms[snp_iter]
			if snp_cm >= mid_left_cm and snp_cm < mid_right_cm:
				# Now check if its a regression snp
				snp_id = G_window_snps[snp_iter]
				alt_snp_id = G_window_alt_snps[snp_iter]
				if snp_id in regression_snp_id_to_regression_snp_position:
					snp_name = snp_id
				elif alt_snp_id in regression_snp_id_to_regression_snp_position:
					snp_name = alt_snp_id
				else:
					continue
				nearby_snps_in_window = np.abs(G_window_cms - snp_cm) <= window_size
				variant_ld = G_window_ld[snp_iter, nearby_snps_in_window]
				variant_ld_scores = np.square(variant_ld)
				variant_adj_ld_scores = variant_ld_scores - ((1.0-variant_ld_scores)/(n_ref_samp-2.0))

				# Quick error check
				if snp_name in regression_snpids_to_nearby_ld_vec:
					print('assumption erororo')
					pdb.set_trace()

				# Add to dictionary to keep track
				regression_snpids_to_nearby_ld_vec[snp_name] = variant_adj_ld_scores
				regression_snp_id_to_nearby_snp_ids[snp_name] = G_window_positions[nearby_snps_in_window]
		# Go to next window
		start_cm = start_cm + window_size

	return regression_snpids_to_nearby_ld_vec, regression_snp_id_to_nearby_snp_ids

def create_mapping_from_snp_id_to_annotation_vector(variant_annotation_file, rsid_to_snpid):
	f = open(variant_annotation_file)
	snp_to_anno = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			anno_names = np.asarray(data[4:])
			continue
		anno_vec = np.asarray(data[4:]).astype(float)
		rsid = data[2]
		if rsid not in rsid_to_snpid:
			continue
		snp_id = rsid_to_snpid[rsid]
		snp_to_anno[snp_id] = anno_vec
	return snp_to_anno, anno_names

def save_annotation_matrix(snp_id_to_annotation_vector, snp_ids, alt_snp_ids, annotation_output_file):
	n_snps = len(snp_ids)
	anno_mat = []
	for snp_iter in range(n_snps):
		snp_id = snp_ids[snp_iter]
		alt_snp_id = alt_snp_ids[snp_iter]
		if snp_id in snp_id_to_annotation_vector and alt_snp_id in snp_id_to_annotation_vector:
			print('assumption eroror')
			pdb.set_trace()
		elif snp_id in snp_id_to_annotation_vector:
			anno_vec = snp_id_to_annotation_vector[snp_id]
		elif alt_snp_id in snp_id_to_annotation_vector:
			anno_vec = snp_id_to_annotation_vector[snp_id]
		else:
			print('assumption eroror')
			pdb.set_trace()
		anno_mat.append(anno_vec)

	anno_mat = np.asarray(anno_mat)
	np.save(annotation_output_file, anno_mat)


#####################
# Command line args
#####################
hapmap3_rsid_file = sys.argv[1]
ldsc_baseline_ld_annotation_dir = sys.argv[2]
ref_1kg_genotype_dir = sys.argv[3]
chrom_num = sys.argv[4]
preprocessed_tgfm_sldsc_data_dir = sys.argv[5]



# Load in genotype data for this chromosome
genotype_stem = ref_1kg_genotype_dir + '1000G.EUR.hg38.' + str(chrom_num)

# Existing variant ld score file
variant_level_ld_score_file = preprocessed_tgfm_sldsc_data_dir + 'baselineLD_no_qtl.' + chrom_num + '.l2.ldscore.gz'
#variant_annotation_file = ldsc_baseline_ld_annotation_dir + 'baselineLD.' + chrom_num + '.annot.gz'
variant_annotation_file = preprocessed_tgfm_sldsc_data_dir + 'baselineLD_no_qtl.' + chrom_num + '.annot'

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

# Create mapping from snp id to annotation vector
snp_id_to_annotation_vector, annotation_names = create_mapping_from_snp_id_to_annotation_vector(variant_annotation_file, rsid_to_snpid)

# Save annotation matrix
annotation_output_file = preprocessed_tgfm_sldsc_data_dir + 'genomic_annotation_' + 'baselineLD_no_qtl.' + chrom_num + '.annot.npy'
save_annotation_matrix(snp_id_to_annotation_vector, genotype_obj['snp_id'], genotype_obj['alt_snp_id'], annotation_output_file)

# In order to save memory
del snp_id_to_annotation_vector


# Initialize dictionary to keep track of variants in ld with each regression snp
window_size = 1 # IN CM
regression_snpids_to_nearby_ld_vec, regression_snp_id_to_nearby_snp_ids = fill_in_variant_level_ld_scores_using_sliding_window_approach(genotype_obj, regression_snp_id_to_regression_snp_position, window_size)

# Save dictionaries as pickle
ld_dictionary_pkl_file = preprocessed_tgfm_sldsc_data_dir + 'regression_np_id_to_variant_ld_score_' + chrom_num + '.pkl'
f = open(ld_dictionary_pkl_file, 'wb') 
pickle.dump(regression_snpids_to_nearby_ld_vec, f)
f.close()

snp_index_dictionary_pkl_file = preprocessed_tgfm_sldsc_data_dir + 'regression_snp_id_to_snp_index_' + chrom_num + '.pkl'
f = open(snp_index_dictionary_pkl_file, 'wb') 
pickle.dump(regression_snp_id_to_nearby_snp_ids, f)
f.close()

# Save regression snp ids to output
regression_snp_id_output_file = preprocessed_tgfm_sldsc_data_dir + 'regression_snp_ids_' + chrom_num + '.txt'
t = open(regression_snp_id_output_file,'w')
t.write('snp_id\trs_id\n')
n_reg_snps = len(regression_snpids)
for snp_iter in range(n_reg_snps):
	rsid = regression_rsids[snp_iter]
	snp_id = regression_snpids[snp_iter]
	snp_info = snp_id.split('_')
	alt_snp_id = snp_info[0] + '_' + snp_info[1] + '_' + snp_info[3] + '_' + snp_info[2] + '_' + snp_info[4]
	if snp_id in regression_snpids_to_nearby_ld_vec:
		t.write(snp_id + '\t' + rsid + '\n')
	elif alt_snp_id in regression_snpids_to_nearby_ld_vec:
		t.write(alt_snp_id + '\t' + rsid + '\n')
	else:
		print('assumption erororo')
		pdb.set_trace()
t.close()


#file_obj = open(ld_dictionary_pkl_file, 'rb')
#dicti_t = pickle.load(file_obj)
#file_obj.close()





