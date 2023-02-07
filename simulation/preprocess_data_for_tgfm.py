import sys
import numpy as np 
import os
import sys
import pdb
import pickle
import time
import gzip









def create_dictionary_mapping_from_rsid_to_genomic_annotation_vector(annotation_file):
	dicti = {}
	rs_id_vec = []
	variant_position_vec = []
	f = open(annotation_file)
	head_count=0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[2]
		annotation = np.asarray(data[4:]).astype(float)
		variant_position_vec.append(int(data[1]))
		rs_id_vec.append(rsid)
		# Quick error check
		if rsid in dicti:
			print('assumption eroror')
			pdb.set_trace()

		dicti[rsid] = annotation

	f.close()

	return dicti, np.asarray(variant_position_vec), np.asarray(rs_id_vec)




def extract_gene_tissue_pairs_and_associated_gene_models_in_window(window_start, window_end, gene_summary_file, simulated_learned_gene_models_dir, simulation_name_string, eqtl_sample_size, window_indices):
	# Initialize output vectors
	gene_tissue_pairs = []
	weight_vectors = []
	gene_tss_arr = []

	# Loop through genes (note: not gene tissue pairs)
	f = open(gene_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent fields from line
		ensamble_id = data[0]
		gene_tss = int(data[2])
		cis_snp_indices_file = data[5]
		cis_snp_indices = np.load(cis_snp_indices_file)

		# Quick error check
		if len(cis_snp_indices) != len(window_indices):
			print('assumption eroror')
			pdb.set_trace()

		if gene_tss >= (window_start + 100000) and gene_tss <= (window_end - 100000):
			# Gene is in cis with respect to window

			# Fitted gene file
			fitted_gene_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_gene_model_pmces.npy'
			gene_model_mat = np.load(fitted_gene_file)

			# Relevent info
			n_tiss = gene_model_mat.shape[0]
			n_cis_snps = gene_model_mat.shape[1]

			# QUick error check
			if np.sum(cis_snp_indices) != n_cis_snps:
				print('assumption eroror')
				pdb.set_trace()

			# Loop through tissues
			for tiss_iter in range(n_tiss):
				# Skip gene, tissue pairs with no gene models
				if np.array_equal(gene_model_mat[tiss_iter,:], np.zeros(n_cis_snps)):
					continue

				# Convert gene-tissue eqtl effect sizes (a vector of length number of cis snps) to window-level eqtl effect sizes
				window_level_eqtl_effect_sizes = np.zeros(np.sum(window_indices))
				window_level_eqtl_effect_sizes[cis_snp_indices[window_indices]] = gene_model_mat[tiss_iter,:]

				# Add info to arrays
				weight_vectors.append(window_level_eqtl_effect_sizes)
				gene_tss_arr.append(gene_tss)
				gene_tissue_pairs.append(ensamble_id + '_' + 'tissue' + str(tiss_iter))

	f.close()

	return np.asarray(gene_tissue_pairs), weight_vectors, np.asarray(gene_tss_arr)


def create_anno_matrix_for_set_of_rsids(rsid_to_genomic_annotation, window_rsids):
	anno_mat_arr = []
	for window_rsid in window_rsids:
		anno_mat_arr.append(rsid_to_genomic_annotation[window_rsid])
	return np.asarray(anno_mat_arr)


def compute_log_expected_probability(full_anno_mat, sldsc_tau_mean, threshold):
	expected_per_ele_h2 = np.dot(full_anno_mat, sldsc_tau_mean)
	expected_per_ele_h2[expected_per_ele_h2 <= threshold] = threshold

	prob = expected_per_ele_h2/np.sum(expected_per_ele_h2)

	return np.log(prob)

def compute_expected_log_probability(full_anno_mat, sldsc_tau_mean, sldsc_tau_cov, threshold, n_samples=10000):
	# Sample a whole bunch of taus
	sampled_taus = np.random.multivariate_normal(mean=sldsc_tau_mean, cov=sldsc_tau_cov, size=n_samples)

	sampled_expected_per_ele_h2 = np.dot(full_anno_mat, np.transpose(sampled_taus))

	sampled_expected_per_ele_h2[sampled_expected_per_ele_h2 <= threshold] = threshold

	prob = np.copy(sampled_expected_per_ele_h2)
	for sample_iter in range(n_samples):
		prob[:, sample_iter] = prob[:, sample_iter]/np.sum(prob[:, sample_iter])

	log_prob = np.log(prob)

	return np.mean(log_prob,axis=1)

def compute_log_expected_probability_variant_v_gene_only(full_anno_mat, sldsc_tau_mean, threshold):
	expected_per_ele_h2 = np.dot(full_anno_mat, sldsc_tau_mean)
	variant_elements = full_anno_mat[:,0] == 1
	gene_elements = full_anno_mat[:,0] != 1
	# Quick error check
	if np.sum(variant_elements) + np.sum(gene_elements) != full_anno_mat.shape[0]:
		print('assumption eroror')
		pdb.set_trace()
	window_snp_h2 = np.sum(expected_per_ele_h2[variant_elements])
	window_gene_h2 = np.sum(expected_per_ele_h2[gene_elements])
	if window_snp_h2 < threshold:
		window_snp_h2 = threshold
	if window_gene_h2 < threshold:
		window_gene_h2 = threshold

	avg_per_ele_h2 = np.ones(len(expected_per_ele_h2))
	avg_per_ele_h2[variant_elements] = window_snp_h2/np.sum(variant_elements)
	avg_per_ele_h2[gene_elements] = window_gene_h2/np.sum(gene_elements)

	prob = avg_per_ele_h2/np.sum(avg_per_ele_h2)

	return np.log(prob)



def compute_various_versions_of_log_prior_probabilities(window_rsids, window_snp_anno_mat, gene_tissue_pairs, sldsc_tau_mean, sldsc_tau_cov, sparse_sldsc_tau, threshold=1e-30):
	# Merge together window_snp_anno_mat with tissue_anno_mat
	n_snps = len(window_rsids)
	n_genes = len(gene_tissue_pairs)
	tissue_anno_mat = np.zeros((n_genes, 10))
	for gene_iter, gene_tissue_pair in enumerate(gene_tissue_pairs):
		tissue_index = int(gene_tissue_pair.split('_')[1].split('issue')[1])
		tissue_anno_mat[gene_iter, tissue_index] = 1

	full_anno_mat_top = np.hstack((window_snp_anno_mat, np.zeros((n_snps,10))))
	full_anno_mat_bottom = np.hstack((np.zeros((n_genes, window_snp_anno_mat.shape[1])), tissue_anno_mat))
	full_anno_mat = np.vstack((full_anno_mat_top, full_anno_mat_bottom))

	# Quick error check
	if full_anno_mat.shape[1] != len(sldsc_tau_mean):
		print('assumption error')
		pdb.set_trace()

	variant_v_gene_only_ln_pi = compute_log_expected_probability_variant_v_gene_only(full_anno_mat, sldsc_tau_mean, threshold)	
	point_estimate_ln_pi = compute_log_expected_probability(full_anno_mat, sldsc_tau_mean, threshold)
	sparse_estimate_ln_pi = compute_log_expected_probability(full_anno_mat, sparse_sldsc_tau, threshold)
	distribution_estimate_ln_pi = compute_expected_log_probability(full_anno_mat, sldsc_tau_mean, sldsc_tau_cov, threshold)

	return point_estimate_ln_pi, sparse_estimate_ln_pi, distribution_estimate_ln_pi, variant_v_gene_only_ln_pi

def load_in_window_gwas_betas_and_ses(window_gwas_summary_file, window_rsids):
	f = open(window_gwas_summary_file)
	beta_vec =[]
	beta_se_vec = []
	rsid_vec = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[0]
		beta = float(data[1])
		beta_se = float(data[2])
		beta_vec.append(beta)
		beta_se_vec.append(beta_se)
		rsid_vec.append(rsid)

	f.close()

	if np.array_equal(np.asarray(rsid_vec), window_rsids) == False:
		print('assumption eroror')
		pdb.set_trace()

	return np.asarray(beta_vec), np.asarray(beta_se_vec)

# Convert gwas summary statistics to *STANDARDIZED* effect sizes
# Following SuSiE code found in these two places:
########1. https://github.com/stephenslab/susieR/blob/master/R/susie_rss.R  (LINES 277-279)
########2. https://github.com/stephenslab/susieR/blob/master/R/susie_ss.R (LINES 148-156 AND 203-205)
def convert_to_standardized_summary_statistics(gwas_beta_raw, gwas_beta_se_raw, gwas_sample_size, R, sigma2=1.0):
	gwas_z_raw = gwas_beta_raw/gwas_beta_se_raw

	XtX = (gwas_sample_size-1)*R
	Xty = np.sqrt(gwas_sample_size-1)*gwas_z_raw
	var_y = 1

	dXtX = np.diag(XtX)
	csd = np.sqrt(dXtX/(gwas_sample_size-1))
	csd[csd == 0] = 1

	XtX = (np.transpose((1/csd) * XtX) / csd)
	Xty = Xty / csd

	dXtX2 = np.diag(XtX)

	beta_scaled = (1/dXtX2)*Xty
	beta_se_scaled = np.sqrt(sigma2/dXtX2)

	return beta_scaled, beta_se_scaled, XtX


def save_ln_pi_output_file(ln_pi, output_file, window_rsids, gene_tissue_pairs):
	t2 = open(output_file,'w')
	t2.write('element_name\tln_pi\n')

	genetic_element_names = np.hstack((window_rsids, gene_tissue_pairs))

	# Quick error check
	if len(genetic_element_names) != len(ln_pi):
		print('assumption error')
		pdb.set_trace()

	for itera, genetic_element_name in enumerate(genetic_element_names):
		t2.write(genetic_element_name + '\t' + str(ln_pi[itera]) + '\n')

	t2.close()

	return


#####################
# Command line args
#####################
simulation_number = int(sys.argv[1])
chrom_num = sys.argv[2]
simulation_name_string = sys.argv[3]
n_gwas_individuals = int(sys.argv[4])
eqtl_sample_size = int(sys.argv[5])
simulation_window_list_file = sys.argv[6]
annotation_file = sys.argv[7]
simulated_gwas_dir = sys.argv[8]
simulated_gene_expression_dir = sys.argv[9]
simulated_learned_gene_models_dir = sys.argv[10]
simulated_sldsc_results_dir = sys.argv[11]
simulated_tgfm_input_data_dir = sys.argv[12]

# Load in SLDSC heritability model
sldsc_tau_mean_file = simulated_sldsc_results_dir + simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_sldsc_results_mean_jacknifed_taus.txt'
sldsc_tau_cov_file = simulated_sldsc_results_dir + simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_sldsc_results_covariance_jacknifed_taus.txt'
sldsc_tau_mean = np.loadtxt(sldsc_tau_mean_file)
sldsc_tau_cov = np.loadtxt(sldsc_tau_cov_file)
# Load in SLDSC sparse heritability model
sldsc_sparse_tau_file = simulated_sldsc_results_dir + simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_sldsc_results_organized_0.5_sparse_ard_eqtl_coefficients_mv_update_res.txt'
sparse_sldsc_tau_raw = np.loadtxt(sldsc_sparse_tau_file, dtype=str, delimiter='\t')
sparse_sldsc_tau = sparse_sldsc_tau_raw[2:,1].astype(float)

# Create dictionary mapping from rsid to genomic annotation vector
rsid_to_genomic_annotation, variant_position_vec, rsids = create_dictionary_mapping_from_rsid_to_genomic_annotation_vector(annotation_file)


# Open outputful summarizing TGFM input (one line for each window)
tgfm_input_data_summary_file = simulated_tgfm_input_data_dir + simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_tgfm_input_data_summary.txt'
t = open(tgfm_input_data_summary_file,'w')
# Write header
t.write('window_name\tLD_npy_file\tTGFM_input_pkl\tlog_prior_probability_file_stem\n')




# Now loop through windows
f = open(simulation_window_list_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	# Skip header
	if head_count == 0:
		head_count = head_count + 1
		continue

	# Extract relevent info for this window
	window_name = data[0]
	window_start = int(data[1])
	window_middle_start = int(data[2])
	window_middle_end = int(data[3])
	window_end = int(data[4])

	# Extract indices of variants in this window
	window_indices = (variant_position_vec >= window_start) & (variant_position_vec < window_end)
	window_rsids = rsids[window_indices]
	window_variant_position_vec = variant_position_vec[window_indices]

	# Create annotation matrix for window rsids
	window_anno_mat = create_anno_matrix_for_set_of_rsids(rsid_to_genomic_annotation, window_rsids)

	# Extract gene-tissue pairs and fitted models in this window
	gene_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'
	gene_tissue_pairs, gene_tissue_pair_weight_vectors, gene_tissue_pairs_tss = extract_gene_tissue_pairs_and_associated_gene_models_in_window(window_start, window_end, gene_summary_file, simulated_learned_gene_models_dir, simulation_name_string, eqtl_sample_size, window_indices)

	# Extract LD
	ld_mat_file = simulated_tgfm_input_data_dir + simulation_name_string + '_' + window_name + '_in_sample_ld.npy'
	ld_mat = np.load(ld_mat_file)

	# Get middle variant indices and middle gene indices
	middle_variant_indices = np.where((window_variant_position_vec >= window_middle_start) & (window_variant_position_vec < window_middle_end))[0]
	middle_gene_indices = np.where((gene_tissue_pairs_tss >= window_middle_start) & (gene_tissue_pairs_tss < window_middle_end))[0]

	# Load in GWAS betas and standard errors
	window_gwas_summary_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results_window_' + window_name + '.txt'
	gwas_beta, gwas_beta_se = load_in_window_gwas_betas_and_ses(window_gwas_summary_file, window_rsids)
	# Standardize gwas beta and se
	beta_scaled, beta_se_scaled, XtX = convert_to_standardized_summary_statistics(gwas_beta, gwas_beta_se, n_gwas_individuals, ld_mat)

	# Compute various ln(pi) and save to output # and save those results to output files
	ln_pi_output_stem = simulated_tgfm_input_data_dir + simulation_name_string + '_' + window_name + '_eqtl_ss_' + str(eqtl_sample_size) + '_ln_pi'
	# Uniform prior
	n_window_elements = len(window_rsids) + len(gene_tissue_pairs)
	uniform_pi = np.ones(n_window_elements)*(1.0/n_window_elements)
	uniform_ln_pi = np.log(uniform_pi)
	save_ln_pi_output_file(uniform_ln_pi, ln_pi_output_stem + '_uniform.txt', window_rsids, gene_tissue_pairs)

	thresholds = [1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-20, 1e-30]
	for threshold in thresholds:
		point_estimate_ln_pi, sparse_estimate_ln_pi, distribution_estimate_ln_pi, variant_v_gene_only_ln_pi = compute_various_versions_of_log_prior_probabilities(window_rsids, window_anno_mat, gene_tissue_pairs, sldsc_tau_mean, sldsc_tau_cov, sparse_sldsc_tau, threshold=threshold)
		save_ln_pi_output_file(point_estimate_ln_pi, ln_pi_output_stem + '_point_estimate_' + str(threshold) + '.txt', window_rsids, gene_tissue_pairs)
		save_ln_pi_output_file(sparse_estimate_ln_pi, ln_pi_output_stem + '_sparse_estimate_' + str(threshold) + '.txt', window_rsids, gene_tissue_pairs)
		save_ln_pi_output_file(distribution_estimate_ln_pi, ln_pi_output_stem + '_distribution_estimate_' + str(threshold) + '.txt', window_rsids, gene_tissue_pairs)
		save_ln_pi_output_file(variant_v_gene_only_ln_pi, ln_pi_output_stem + '_variant_v_gene_only_' + str(threshold) + '.txt', window_rsids, gene_tissue_pairs)


	# Organize TGFM data into nice data structure
	tgfm_data = {}
	tgfm_data['genes'] = gene_tissue_pairs
	tgfm_data['variants'] = window_rsids
	tgfm_data['gwas_beta'] = beta_scaled
	tgfm_data['gwas_beta_se'] = beta_se_scaled
	tgfm_data['gwas_sample_size'] = n_gwas_individuals
	tgfm_data['gene_eqtl_pmces'] = np.asarray(gene_tissue_pair_weight_vectors)
	tgfm_data['middle_gene_indices'] = middle_gene_indices
	tgfm_data['middle_variant_indices'] = middle_variant_indices

	# Save TGFM output data to pickle
	window_pickle_output_file = simulated_tgfm_input_data_dir + simulation_name_string + '_' + window_name + '_eqtl_ss_' + str(eqtl_sample_size) + '_tgfm_input_data.pkl'
	g = open(window_pickle_output_file, "wb")
	pickle.dump(tgfm_data, g)
	g.close()


	# Write to output summary file
	t.write(window_name + '\t' + ld_mat_file + '\t' + window_pickle_output_file + '\t' + ln_pi_output_stem + '\n')

t.close()
f.close()



