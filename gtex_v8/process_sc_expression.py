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



######################
# Command line args
######################
input_h5py_file = sys.argv[1]
input_h5py_pseudobulk_file = sys.argv[2]
processed_sc_expression_dir = sys.argv[3]

######################
# Filtering parameters
#######################
np.random.seed(0)
sc.settings.verbosity = 3 
num_hvg = 2000
regress_out_batch = 'True'

########################
# Get QC metrics
adata = sc.read_h5ad(input_h5py_file)

adata_pb = sc.read_h5ad(input_h5py_pseudobulk_file)

# Extract qc metrics for each cell
qc_metrics = sc.pp.calculate_qc_metrics(adata.raw)[0]

######################
# Load in ScanPy data
#######################
adata_scran = sc.read_h5ad(input_h5py_scran_file)


#######################################
# Add cell qc metrics to scran anndata object
#######################################
for column_name in qc_metrics.columns:
	adata_scran.obs[column_name] = qc_metrics[column_name]
adata_scran.var['mt'] = adata_scran.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata_scran, qc_vars=['mt'],inplace=True)

#######################################
# Filter to European ancestry individuals
#######################################
adata_scran = adata_scran[adata_scran.obs['pop_cov'] == 'European', :]

#######################################
# Gene filtering
#######################################
sc.pp.filter_cells(adata_scran, min_genes=200)
adata_scran = adata_scran[adata_scran.obs.pct_counts_mt < 5, :]
# Get genes expressed in at least XX% of cells
min_fraction_of_cells = .01
sc.pp.filter_genes(adata_scran, min_cells=(adata_scran.X.shape[0])*min_fraction_of_cells)

#######################
# Save adata obect to output
#######################
un_filtered_h5_output_file = processed_sc_expression_dir + '_scran_normalized_all_genes.h5ad'
adata_scran.write(un_filtered_h5_output_file)


# Get highly variable genes
sc.pp.highly_variable_genes(adata_scran, n_top_genes=(num_hvg))
# Filter to those highly variable genes
adata_scran = adata_scran[:, adata_scran.var.highly_variable]

#######################################
# Regress out batch
#######################################
if regress_out_batch == 'True':
	# Seemingly somewhat hacky fix recommended by theis lab, necessary to run combat below
	adata_scran = adata_scran.copy()
	# Run combat
	sc.pp.combat(adata_scran, key='batch_cov')

#######################################
# Save un-scaled data to raw
#######################################
adata_scran.raw = adata_scran

#######################################
# Standardize	
#######################################
sc.pp.scale(adata_scran, max_value=10)


#######################################	
# Run PCA
#######################################
print('pca')
sc.tl.pca(adata_scran, svd_solver='arpack')

##################
# Run UMAP on full data
##################
sc.pp.neighbors(adata_scran)
sc.tl.umap(adata_scran)

#######################
# output file stem
#######################	
output_file_stem = processed_sc_expression_dir + 'scran_normalization_hvg_' + str(num_hvg) + '_regress_batch_' + regress_out_batch
	
#######################
# Save Expression PCs
#######################
pc_output_file = output_file_stem + '_sc_expression_pcs.txt'
np.savetxt(pc_output_file, adata_scran.obsm['X_pca'], fmt="%s", delimiter='\t', comments='')

#######################
# Save Expression PCs variance
#######################
pc_pve_output_file = output_file_stem + '_sc_expression_pcs_percent_variance_explained.txt'
np.savetxt(pc_pve_output_file, adata_scran.uns['pca']['variance_ratio'], fmt="%s", delimiter='\t', comments='')


#######################
# Save UMAP loadings
#######################
umap_output_file = output_file_stem + '_sc_expression_umaps.txt'
np.savetxt(umap_output_file, adata_scran.obsm['X_umap'], fmt="%s", delimiter='\t', comments='')

#######################
# Save Covariate Info
#######################
covariate_output_file = output_file_stem + '_cell_covariates.txt'
np.savetxt(covariate_output_file, adata_scran.obs, fmt="%s", delimiter='\t', header='\t'.join(adata_scran.obs.columns), comments='')

#######################
# Save adata obect to output
#######################
temp2_h5_output_file = output_file_stem + '_2.h5ad'
adata_scran.write(temp2_h5_output_file)