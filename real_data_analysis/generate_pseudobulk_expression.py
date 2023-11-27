import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import pdb
import scanpy as sc
import h5py
from sklearn.decomposition import PCA
from anndata import AnnData
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.linear_model import LinearRegression
from sklearn.cluster import KMeans
import pandas as pd
from scipy.sparse import csr_matrix
import scipy.io
import gzip
import rnaseqnorm





def extract_gene_name_and_type_from_info(gene_info):
	ensamble_id = 'NA'
	gene_id = 'NA'
	gene_type = 'NA'
	gene_status = 'NA'
	for stringer in gene_info:
		if stringer.startswith('ID='):
			ensamble_id = stringer.split('ID=')[1]
		if stringer.startswith('gene_type='):
			gene_type = stringer.split('gene_type=')[1]
		if stringer.startswith('gene_status='):
			gene_status = stringer.split('gene_status=')[1]
		if stringer.startswith('gene_name='):
			gene_id = stringer.split('gene_name=')[1]
	if ensamble_id == 'NA' or gene_id == 'NA' or gene_type == 'NA' or gene_status == 'NA':
		print('assumption erroror')
		pdb.set_trace()
	return ensamble_id, gene_id, gene_type, gene_status

def extract_protein_coding_known_autosomal_genes(gene_struct, gene_annotation_file):
	# Make dictionary list of genes that are protein-coding, known, and autosomal
	valid_genes = {}
	f = open(gene_annotation_file)
	head_count = 0
	types = {}
	status = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip header
		if line.startswith('#'):
			continue
		if len(data) != 9:
			print('assumption error')
			pdb.set_trace()
		sequence_class = data[2]
		if sequence_class != 'gene':
			continue
		if data[0] == 'chrX' or data[0] == 'chrY' or data[0] == 'chrM':
			continue

		gene_info = data[8].split(';')

		ensamble_id, gene_id, gene_type, gene_status = extract_gene_name_and_type_from_info(gene_info)

		if gene_status != 'KNOWN':
			continue
		if gene_type != 'protein_coding':
			continue

		cat_gene_name = gene_id + '_' + ensamble_id.split('.')[0]
		# Error check
		if cat_gene_name in valid_genes:
			print('repeated gene')
			pdb.set_trace()
		valid_genes[cat_gene_name] = 1

	f.close()
	# Make binary vectory corresponding to ordered list of whehter our genes are protein-coding, autosomal and known
	binary_vector = []
	ensamble_ids = gene_struct[gene_struct.columns[0]]
	gene_ids = gene_struct.index
	num_genes = len(gene_ids)
	for gene_num in range(num_genes):
		ensamble_id = ensamble_ids[gene_num]
		gene_id = gene_ids[gene_num]
		if gene_id + '_' + ensamble_id not in valid_genes:
			binary_vector.append(False)
		else:
			binary_vector.append(True)
	return np.asarray(binary_vector)

def get_ordered_individuals_that_we_have_genotype_and_expression_data_for(adata, genotyped_individuals_file, min_cells_per_indi):
	arr = []
	rna_seqed_indis = np.unique(adata.obs['ind_cov'])
	dicti = {}
	num_cells = []
	for rna_seqed_indi in rna_seqed_indis:
		dicti[rna_seqed_indi] = 1
	head_count = 0
	f = open(genotyped_individuals_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		indi = data[1]
		if indi not in dicti:
			print('assumption eororor')
			pdb.set_trace()
		num_cells.append(sum(adata.obs['ind_cov'] == indi))
		arr.append(indi)
	f.close()
	arr = np.asarray(arr)
	num_cells = np.asarray(num_cells)
	indices = num_cells > min_cells_per_indi
	return arr[indices]


def generate_pseudobulk_expression(cell_indi_arr, num_genes, sc_expr_mat_subset, ordered_individuals, min_cells_per_indi):
	reordered_individuals = []
	for individual_index, individual in enumerate(ordered_individuals):
		# Boolean array of cells assigned to this individual
		cells_to_individual_boolean = cell_indi_arr == individual
		if np.sum(cells_to_individual_boolean) >= min_cells_per_indi:
			reordered_individuals.append(individual)
	reordered_individuals = np.asarray(reordered_individuals)

	num_individuals = len(reordered_individuals)
	#num_genes = adata_subset.shape[1]
	pseudobulk_expr = np.zeros((num_individuals, num_genes))

	num_cells_per_individual = []

	for individual_index, individual in enumerate(reordered_individuals):
		# Boolean array of cells assigned to this individual
		cells_to_individual_boolean = cell_indi_arr == individual

		pseudobulk_expr[individual_index, :] = np.asarray(np.sum(sc_expr_mat_subset[cells_to_individual_boolean,:], axis=0))

		num_cells_per_individual.append(np.sum(cells_to_individual_boolean))

	return pseudobulk_expr, num_cells_per_individual, reordered_individuals

def normalize_expression(raw_pseudobulk_expression, sample_level_normalization, gene_level_normalization):
	# Initialize output normalized expression matrix
	normalized_expression = np.zeros(raw_pseudobulk_expression.shape)

	##################################
	# Perform sample level normalization
	##################################
	if sample_level_normalization == 'qn':
		df = pd.DataFrame(np.transpose(raw_pseudobulk_expression))
		temp_out = rnaseqnorm.normalize_quantiles(df)
		raw_pseudobulk_expression = np.transpose(np.asarray(temp_out))
	elif sample_level_normalization == 'edger_cpm':
		df = pd.DataFrame(np.transpose(raw_pseudobulk_expression))
		temp_out = rnaseqnorm.edgeR_cpm(df, log=True)
		raw_pseudobulk_expression = np.transpose(np.asarray(temp_out))

	##################################
	# Perform gene level normalization
	##################################
	if gene_level_normalization == 'zscore':
		for gene_num in range(normalized_expression.shape[1]):
			temp_expr = (raw_pseudobulk_expression[:, gene_num] - np.mean(raw_pseudobulk_expression[:, gene_num]))/np.std(raw_pseudobulk_expression[:, gene_num])
			temp_expr = temp_expr - np.mean(temp_expr)
			if np.sum(np.isnan(temp_expr)) > 0:
				pdb.set_trace()
			normalized_expression[:, gene_num] = temp_expr
	elif gene_level_normalization == 'ign':
		# Code from GTEx v8
		# Project each gene onto a gaussian
		df = pd.DataFrame(np.transpose(raw_pseudobulk_expression))
		norm_df = rnaseqnorm.inverse_normal_transform(df)
		normalized_expression = np.transpose(np.asarray(norm_df))
	else:
		print(gene_level_normalization + ' gene level normalization method currently not implemented')
		pdb.set_trace()

	return normalized_expression

# Generate expression PC loadings and variance explained of those expression PCs
def generate_pca_scores_and_variance_explained(X, num_pcs):
	# Run PCA (via SVD)
	#uuu, sss, vh = np.linalg.svd(np.transpose(X), full_matrices=False)
	#svd_loadings = np.transpose(vh)[:,:num_pcs]
	#ve = (np.square(sss)/np.sum(np.square(sss)))[:num_pcs]

	# Faster in sklearn
	_pca = PCA(n_components=num_pcs, svd_solver='arpack')
	svd_loadings = _pca.fit_transform(X)
	ve = _pca.explained_variance_ratio_

	return svd_loadings, ve


def filter_genes_in_pseudobulk_file(pseudobulk_expression, ordered_gene_names, fraction_samples_expressed_threshold):
	# Quick error checking
	if pseudobulk_expression.shape[1] != ordered_gene_names.shape[0]:
		print('assumption error')
		pdb.set_trace()
	num_genes = pseudobulk_expression.shape[1]
	new_gene_indices = []
	for gene_num in range(num_genes):
		# Compute fraction of samples for this gene that are expressed
		fraction_expressed = float(sum(pseudobulk_expression[:,gene_num] > 0.0))/len(pseudobulk_expression[:,gene_num])
		if fraction_expressed < fraction_samples_expressed_threshold:
			continue
		new_gene_indices.append(gene_num)
	new_gene_indices = np.asarray(new_gene_indices)
	return pseudobulk_expression[:, new_gene_indices], ordered_gene_names[new_gene_indices]


def extract_individual_covariates(ordered_individuals, adata, individual_to_genotype_pcs):
	# Create  header for covariates
	covariate_header = 'Age\tSex\tSLE_status'
	for pc_num in range(1,4):
		covariate_header = covariate_header + '\tGenotype_pc_' + str(pc_num)

	# Initialize dictionary to keep track of individual to covariate mapping
	dicti = {}

	# Loop through individauls
	for individual in ordered_individuals:
		# Extract list of cells corresponding to this individual
		individual_indices = adata.obs['ind_cov'] == individual
		
		# Get Age of this individual
		age_arr = adata.obs['Age'][individual_indices]
		individual_age = age_arr[0]

		# Get sex of this individaul
		sex_arr = adata.obs['Sex'][individual_indices]
		if len(np.unique(sex_arr)) != 1:
			print('assumption erororo')
			pdb.set_trace()
		individual_sex = sex_arr[0]
		individual_sex_binary = '0.0'
		if individual_sex == 'Female':
			individual_sex_binary = '1.0'

		# Get SLE status of this individual
		sle_arr = adata.obs['SLE_status'][individual_indices]
		if len(np.unique(sle_arr)) != 1:
			print('assumption eroror')
			pdb.set_trace()
		individual_sle = sle_arr[0]
		individual_sle_binary = '0.0'
		if individual_sle == 'SLE':
			individual_sle_binary = '1.0'

		# Genotype PCs
		individual_genotype_pcs = individual_to_genotype_pcs[individual]

		# Construct genotype arr
		individual_covariate_arr = np.hstack(([individual_age, individual_sex_binary, individual_sle_binary], individual_genotype_pcs))
		
		# Add to dictionary 
		if individual in dicti:
			print('asssumption eororor')
			pdb.set_trace()
		dicti[individual] = individual_covariate_arr
	return dicti, covariate_header

def create_mapping_from_individual_to_genotype_pcs(genotype_pc_file):
	dicti = {}
	head_count = 0
	f = open(genotype_pc_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		indi_id = data[0]
		if head_count == 0:
			head_count = head_count + 1
			continue
		if indi_id in dicti:
			print('assumption erroror')
			pdb.set_trace()
		dicti[indi_id] = data[1:]
	f.close()
	return dicti

def generate_full_covariate_mat(cluster_specific_individuals, individual_to_covariate_vector, covariate_header, expr_pcs):
	# Error chekcing
	if len(cluster_specific_individuals) != expr_pcs.shape[0]:
		print('assumption eroror')

	full_covariate_header = 'Individual_id\t' + covariate_header

	for pc_num in range(1,11):
		full_covariate_header = full_covariate_header + '\texpression_pc_' + str(pc_num)


	full_covariate_data = []

	for i, individual in enumerate(cluster_specific_individuals):

		indi_cov = individual_to_covariate_vector[individual]
		indi_expr_pcs = expr_pcs[i,:]

		full_covariate_data.append(np.hstack(([individual], indi_cov, indi_expr_pcs)))

	full_covariate_data = np.asarray(full_covariate_data)

	return np.vstack((np.asarray(full_covariate_header.split('\t')), full_covariate_data))

def make_num_cells_per_individual_output_file(num_cells_per_individual, cluster_specific_individuals, num_cells_output_file):
	t = open(num_cells_output_file, 'w')
	t.write('Individual_id\tNumber_of_cells\n')
	for i, individual in enumerate(cluster_specific_individuals):
		t.write(individual + '\t' + str(num_cells_per_individual[i]) + '\n')
	t.close()

#######################
# Command line args
#######################
processed_sc_expression_dir = sys.argv[1]
input_h5py_file = sys.argv[2]
processed_genotype_dir = sys.argv[3]
pseudobulk_expression_dir = sys.argv[4]
gene_annotation_file = sys.argv[5]


# Genotype PC file
genotype_pc_file = processed_genotype_dir + 'inds_v2_sample_filter_ind_snp_pcs.eigenvec'
individual_to_genotype_pcs = create_mapping_from_individual_to_genotype_pcs(genotype_pc_file)

# Create ouptut file summarizing all created pseudobulk data sets
psuedobulk_data_set_summary_file = pseudobulk_expression_dir + 'pseudobulk_data_set_summary.txt'
t = open(psuedobulk_data_set_summary_file,'w')
t.write('data_set_name\tcluster_method_name\tcluster_name\tpseudobulk_expression_file\tcovariate_file\tnum_donors\tnum_genes\tnum_cells_per_individual_file\n')

# File containing ordered list of individuals that we have genotype and expresssion data for
genotyped_individuals_file = processed_sc_expression_dir + 'individual_info_european_rna_and_dna.txt'
# NOW GONNA USE filtered_sample_info_file for this

# Load in (all genes version of ANn data)
adata = sc.read_h5ad(input_h5py_file)

# Quick hacky fix to cluster_assignments mat
# Use Perez et al assignments fro third bin
cluster_assignments = np.asarray(adata.obs['cg_cov'])


# Extract protein coding gene indices
protein_coding_gene_indices = extract_protein_coding_known_autosomal_genes(adata.raw.var, gene_annotation_file)

# Subset adata matrix to protein-coding autosomal
#adata.raw = adata.raw[:, protein_coding_gene_indices]

# Extract ordered list of individuals that we have genotype data for
min_cells_per_indi = 2500
ordered_individuals = get_ordered_individuals_that_we_have_genotype_and_expression_data_for(adata, genotyped_individuals_file, min_cells_per_indi)



# Extract individual covariates
individual_to_covariate_vector, covariate_header = extract_individual_covariates(ordered_individuals, adata, individual_to_genotype_pcs)

# Ordered gene names
ordered_gene_names = adata.raw.var['gene_ids'][protein_coding_gene_indices]

# Np matrix expression mat
#sc_expr_mat = adata.X.toarray()
sc_expr_mat = (adata.raw.X[:,protein_coding_gene_indices]).toarray()



# Name of cluster method
cluster_method_name = 'cell_types'
# Vector of unique cluster assignments
unique_cluster_assignments = np.sort(np.unique(cluster_assignments))
# Number of unique cluster assignments
num_unique_cluster_assignments = len(unique_cluster_assignments)

# Loop through unique cluster assignments
for cluster_assignment in unique_cluster_assignments:
	print(cluster_assignment)
	# Get boolean vector of cells representing whether cell is assigned to the current cluster
	cell_assigned_to_cluster_boolean = cluster_assignments == cluster_assignment
		
	# Generate pseudobulk expression
	min_cells_per_indi = 5
	pseudobulk_expression, num_cells_per_individual, cluster_specific_individuals = generate_pseudobulk_expression(adata.obs['ind_cov'][cell_assigned_to_cluster_boolean], sc_expr_mat.shape[1], sc_expr_mat[cell_assigned_to_cluster_boolean,:], ordered_individuals, min_cells_per_indi)

	# Filter genes in pseudobulk file
	fraction_samples_expressed_threshold = .8
	pseudobulk_expression, pseudobulk_gene_names = filter_genes_in_pseudobulk_file(pseudobulk_expression, ordered_gene_names, fraction_samples_expressed_threshold)


	# Normalize gene expression
	# Options for sample level normalization are currently 'none'
	sample_level_normalization = 'edger_cpm'
	# Options for gene level normalization are 'zscore' and 'ign'
	gene_level_normalization = 'zscore'
	normalized_expression = normalize_expression(pseudobulk_expression, sample_level_normalization, gene_level_normalization)

	# Run PCA on pseudobulk data
	num_pcs = 10
	expr_pcs, expr_pcs_pve = generate_pca_scores_and_variance_explained(normalized_expression, num_pcs)

	# Create full covariates
	full_covariate_mat = generate_full_covariate_mat(cluster_specific_individuals, individual_to_covariate_vector, covariate_header, expr_pcs)

	# Create expr mat with headers and column names
	temp_gene_names = np.hstack((['gene_id'], np.asarray(pseudobulk_gene_names)))
	normalized_expression_annotated = np.transpose(np.vstack((temp_gene_names, np.hstack((np.transpose(np.asmatrix(cluster_specific_individuals)),normalized_expression.astype(str))))))

	#name of data set
	data_set_name = cluster_method_name + '_' + cluster_assignment

	# Save expression mat to output file
	expression_output_file = pseudobulk_expression_dir + data_set_name + '_pseudobulk_expression.txt'
	np.savetxt(expression_output_file, normalized_expression_annotated, fmt="%s", delimiter='\t')


	# Save covariate mat to output file
	covariate_output_file = pseudobulk_expression_dir + data_set_name + '_pseudobulk_covariates.txt'
	np.savetxt(covariate_output_file, full_covariate_mat, fmt="%s", delimiter='\t')

	# Make number of cells per individual output file
	num_cells_output_file = pseudobulk_expression_dir + data_set_name + '_num_cells_per_individual.txt'
	make_num_cells_per_individual_output_file(num_cells_per_individual, cluster_specific_individuals, num_cells_output_file)

	t.write(data_set_name + '\t' + cluster_method_name + '\t' + cluster_assignment + '\t' + expression_output_file + '\t' + covariate_output_file + '\t')
	t.write(str(len(cluster_specific_individuals)) + '\t' + str(pseudobulk_gene_names.shape[0]) + '\t' + num_cells_output_file + '\n')

t.close()



################
# Filter pseudobulk data set summary file to cell types with reasonable number of samples
psuedobulk_data_set_summary_file = pseudobulk_expression_dir + 'pseudobulk_data_set_summary.txt'
psuedobulk_data_set_summary_filtered_file = pseudobulk_expression_dir + 'pseudobulk_data_set_summary_filtered.txt'
head_count = 0
f = open(psuedobulk_data_set_summary_file)
t = open(psuedobulk_data_set_summary_filtered_file,'w')
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	if float(data[5]) < 100:
		continue
	t.write(line + '\n')
f.close()
t.close()




