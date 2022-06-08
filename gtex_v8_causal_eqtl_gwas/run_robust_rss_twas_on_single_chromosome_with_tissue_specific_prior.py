import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import pickle
import tgfm_causal_twas
import tgfm_robust_causal_twas
import tgfm_fine_mapping
import tgfm_rss_fine_mapping
import rpy2
import rpy2.robjects.numpy2ri as numpy2ri
import rpy2.robjects as ro
ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
susieR_pkg = importr('susieR')


def extract_gwas_component_susie_pmces_in_same_order_as_twas_data(gwas_component_susie_pmces_file, twas_bim):
	num_twas_var = twas_bim.shape[0]
	# First create mapping from twas variants to position
	mapping = {}
	for var_iter in range(num_twas_var):
		var_name = twas_bim[var_iter, 1]
		mapping[var_name] = (var_iter, twas_bim[var_iter,4], twas_bim[var_iter,5])
		var_info = var_name.split('_')
		alt_var_name = var_info[0] + '_' + var_info[1] + '_' + var_info[3] + '_' + var_info[2] + '_' + var_info[4]
		mapping[alt_var_name] = (var_iter, twas_bim[var_iter,4], twas_bim[var_iter,5])

	# Load in gwas pmces data
	gwas_component_pmces_data = np.loadtxt(gwas_component_susie_pmces_file, dtype=str,delimiter='\t')
	gwas_component_pmces_data = gwas_component_pmces_data[1:,:]

	
	# Initialize output array
	gwas_pmces = np.zeros(num_twas_var)
	gwas_susie_mu = np.zeros(num_twas_var)
	gwas_susie_mu_sd = np.zeros(num_twas_var)
	gwas_susie_alpha = np.zeros(num_twas_var)


	used_positions = {}
	for var_iter in range(gwas_component_pmces_data.shape[0]):
		variant_id = gwas_component_pmces_data[var_iter,0] + '_b38'
		if variant_id in mapping:
			variant_info = variant_id.split('_')
			rev = 1.0
			if variant_info[2] != mapping[variant_id][1]:
				rev = -1.0
			gwas_pmces[mapping[variant_id][0]] = float(gwas_component_pmces_data[var_iter,1])*rev
			gwas_susie_mu[mapping[variant_id][0]] = float(gwas_component_pmces_data[var_iter,2])*rev
			gwas_susie_mu_sd[mapping[variant_id][0]] = float(gwas_component_pmces_data[var_iter,4])
			gwas_susie_alpha[mapping[variant_id][0]] = float(gwas_component_pmces_data[var_iter,3])

			used_positions[mapping[variant_id][0]] = 1
	if len(used_positions) != num_twas_var:
		print('assumption eroror')
		pdb.set_trace()
	return gwas_pmces, gwas_susie_mu, gwas_susie_mu_sd, gwas_susie_alpha


def get_eqtl_component_level_pmces(susie_mu, susie_alpha):
	susie_pmces = []
	for component_num in range(len(susie_mu)):
		susie_pmces.append(susie_mu[component_num]*susie_alpha[component_num])
	return susie_pmces

def extract_tissue_names(gtex_pseudotissue_file):
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

def run_variant_level_susie_on_this_region(betas, betas_se, sample_size, ref_ld, original_gwas_susie_alpha):
	# res <- susie_rss(bhat=as.numeric(beta_mat[trait_num,]), shat=as.numeric(std_err_mat[trait_num,]), R=LD, n=sample_sizes[trait_num])
	try:
		m = len(betas)
		susie_variant_obj = susieR_pkg.susie_rss(bhat=betas.reshape((m,1)), shat=betas_se.reshape((m,1)), R=ref_ld, n=sample_size)

		susie_alpha = susie_variant_obj.rx2('alpha')
		susie_mu = susie_variant_obj.rx2('mu')
		susie_mu2 = susie_variant_obj.rx2('mu2')
		susie_mu_sd = np.sqrt(susie_mu2 - np.square(susie_mu))

		component_num = -1
		correlations = []
		for temp_component in range(susie_alpha.shape[0]):
			corry = np.corrcoef(susie_alpha[temp_component, :], original_gwas_susie_alpha)[0,1]
			correlations.append(corry)
		correlations = np.asarray(correlations)

		component_num = np.nanargmax(correlations)

		if correlations[component_num] < .5:
			print('low correlation error')
			discover_component_bool = False
		else:
			discover_component_bool = True
		pass_bool = True
	except ValueError:
		pass_bool = False
		susie_alpha = 'na'
		susie_mu = 'na'
		susie_mu_sd = 'na'
		component_num = -1
		discover_component_bool = False
		pass_bool = False

	return susie_alpha, susie_mu, susie_mu_sd, component_num, discover_component_bool, pass_bool

def run_variant_level_susie_on_this_region_debug(betas, betas_se, sample_size, ref_ld, component_name, trait_name, twas_bim_df):
	'''
	variant_id_to_pos = {}
	nrow = twas_bim_df.shape[0]
	for row_num in range(nrow):
		variant_id_raw = twas_bim_df[row_num,1]
		variant_id_info = variant_id_raw.split('_')
		variant_id1 = variant_id_info[0] + '_' + variant_id_info[1] + '_' + variant_id_info[2] + '_' + variant_id_info[3]
		variant_id2 = variant_id_info[0] + '_' + variant_id_info[1] + '_' + variant_id_info[3] + '_' + variant_id_info[2]
		if variant_id1 in variant_id_to_pos or variant_id2 in variant_id_to_pos:
			print('assumption eroror')
			pdb.set_trace()
		variant_id_to_pos[variant_id1] = row_num
		variant_id_to_pos[variant_id2] = row_num
	'''


	orig_susie_data_dir= '/n/groups/price/ben/causal_eqtl_gwas/gtex_v8_causal_eqtl_gwas/ukbb_preprocessed_for_genome_wide_susie/'
	window_name = component_name.split('_')[-1]
	beta_file = orig_susie_data_dir + window_name + '_beta.txt'
	beta_se_file = orig_susie_data_dir + window_name + '_beta_std_err.txt'
	ref_1kg_genotype_file = orig_susie_data_dir + window_name + '_ref_1kg_genotype.txt'
	sample_size_file = orig_susie_data_dir + window_name + '_study_sample_sizes.txt'
	variant_name_file = orig_susie_data_dir + window_name + '_variant_ids.txt'
	study_file = orig_susie_data_dir + window_name + '_studies.txt'

	studies_orig = np.loadtxt(study_file, dtype=str)
	study_index = np.where(studies_orig == trait_name)[0][0]

	beta_orig = np.loadtxt(beta_file)[study_index,:]
	beta_se_orig = np.loadtxt(beta_se_file)[study_index,:]
	ref_geno_orig = np.loadtxt(ref_1kg_genotype_file)
	sample_size_orig = np.loadtxt(sample_size_file)[study_index]
	variant_name_orig = np.loadtxt(variant_name_file, dtype=str)
	ld_orig = np.corrcoef(np.transpose(ref_geno_orig))

	# Variant re-ordering
	variant_indices = []
	variant_id_to_pos = {}
	for i,variant_name in enumerate(variant_name_orig):
		variant_id_to_pos[variant_name] = i
	nrow = twas_bim_df.shape[0]
	for row_num in range(nrow):
		variant_id_raw = twas_bim_df[row_num,1]
		variant_id_info = variant_id_raw.split('_')
		variant_id1 = variant_id_info[0] + '_' + variant_id_info[1] + '_' + variant_id_info[2] + '_' + variant_id_info[3]
		variant_id2 = variant_id_info[0] + '_' + variant_id_info[1] + '_' + variant_id_info[3] + '_' + variant_id_info[2]
		if variant_id1 in variant_id_to_pos and variant_id2 in variant_id_to_pos:
			print('assumption eroror')
			pb.set_trace()
		elif variant_id1 in variant_id_to_pos:
			variant_indices.append(variant_id_to_pos[variant_id1])
		elif variant_id2 in variant_id_to_pos:
			variant_indices.append(variant_id_to_pos[variant_id2])
	variant_indices = np.asarray(variant_indices)


	m = len(beta_orig)
	susie_variant_obj_orig = susieR_pkg.susie_rss(bhat=beta_orig.reshape((m,1)), shat=beta_se_orig.reshape((m,1)), R=ld_orig, n=sample_size)

	# NOTE: THIS PERFECTLY RECAPTURES ORIGINAL (gwas_susie_alpha)
	susie_alpha_orig = susie_variant_obj_orig.rx2('alpha')
	susie_alpha_orig_reordered = susie_alpha_orig[:,variant_indices]

	m = len(variant_indices)
	ld_orig_data_small = np.corrcoef(np.transpose(ref_geno_orig[:, variant_indices]))
	susie_variant_obj_orig_data_small = susieR_pkg.susie_rss(bhat=beta_orig[variant_indices].reshape((m,1)), shat=beta_se_orig[variant_indices].reshape((m,1)), R=ld_orig_data_small, n=sample_size)

	pdb.set_trace()


	# res <- susie_rss(bhat=as.numeric(beta_mat[trait_num,]), shat=as.numeric(std_err_mat[trait_num,]), R=LD, n=sample_sizes[trait_num])
	m = len(betas)
	susie_variant_obj = susieR_pkg.susie_rss(bhat=betas.reshape((m,1)), shat=betas_se.reshape((m,1)), R=ref_ld, n=sample_size)

	pdb.set_trace()


def correlate_predicted_trait_effect_sizes_with_predictec_trait_eqtl_effect_sizes(gwas_susie_mu, gwas_susie_alpha, eqtl_susie_mu, eqtl_susie_alpha, twas_alpha, twas_alpha_sd):
	gwas_susie_pmces = np.sum(gwas_new_susie_mu*gwas_new_susie_alpha,axis=0)
	num_genes = len(twas_alpha)
	eqtl_susie_pmces = np.zeros(len(gwas_susie_pmces))
	for gene_num in range(num_genes):
		eqtl_susie_pmces = eqtl_susie_pmces + np.sum(eqtl_susie_mu[gene_num]*eqtl_susie_alpha[gene_num], axis=0)*twas_alpha[gene_num]

	pdb.set_trace()

def get_tissues(gene_names):
	tissues = []
	for gene_name in gene_names:
		gene_info = gene_name.split('_')
		tissue = '_'.join(gene_info[1:])
		tissues.append(tissue)
	return np.asarray(tissues)


chrom_num = sys.argv[1]
trait_name = sys.argv[2]
gtex_pseudotissue_file = sys.argv[3]
component_data_file = sys.argv[4]
ukbb_genome_wide_susie_organized_results_dir = sys.argv[5]
output_dir = sys.argv[6]
gene_version = sys.argv[7]


ordered_tissue_names = extract_tissue_names(gtex_pseudotissue_file)
tissue_to_position_mapping = {}
for i, val in enumerate(ordered_tissue_names):
	tissue_to_position_mapping[val] = i

# Extract tissue specific prior (vector of length number of tissues)
tissue_specific_prior_precision_file = output_dir + trait_name + '_' + gene_version + '_count_genes_once_null_init_robust_tissue_specific_prior_precision_temp.txt'
tissue_specific_prior_data = np.loadtxt(tissue_specific_prior_precision_file,dtype=str,delimiter='\t')
tissue_specific_prior_expected_variance = 1.0/tissue_specific_prior_data[1:,1].astype(float)

# Extract pleiotropic variance
pleiotropic_variance_file = output_dir + trait_name + '_' + gene_version + '_count_genes_once_null_init_robust_pleiotropic_prior_precision_temp.txt'
pleiotropic_variance = 1.0/float(np.loadtxt(pleiotropic_variance_file,dtype=str)[1,0])



ordered_tissue_names = extract_tissue_names(gtex_pseudotissue_file)
# open output file
output_file = output_dir + trait_name + '_' + chrom_num + '_summary_tgfm_robust_results_' + gene_version + '_tissue_specific_prior.txt'

t = open(output_file,'w')
# write header
t.write('component_name\tnum_genes\ttgfm_twas_pickle\ttgfm_rss_regression_fine_mapping_table_file\ttgfm_rss_regression_tissue_fine_mapping_table_file\tpickled_data_file\n')

# Loop through trait components
f = open(component_data_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	# Skip header
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Extract relevent fields from line corresponding to trait component
	component_name = data[1]
	num_genes = int(data[2])
	pickled_data_file = data[3]


	print(component_name + '\t' + str(num_genes))

	if pickled_data_file == 'NA':
		t.write(component_name + '\t' + str(num_genes) + '\t' 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\n')
	else:
		# Load in pickled data for this component
		f = open(pickled_data_file, "rb")
		twas_data = pickle.load(f)
		f.close()

		# Load in gwas component susie pmces
		gwas_component_susie_pmces_file = ukbb_genome_wide_susie_organized_results_dir + trait_name + '_' + component_name + '_component_posterior_mean_causal_effect_sizes.txt'
		gwas_component_susie_pmces, gwas_susie_mu, gwas_susie_mu_sd, gwas_susie_alpha = extract_gwas_component_susie_pmces_in_same_order_as_twas_data(gwas_component_susie_pmces_file, twas_data['bim'])

		# Run variant-level susie on this region
		gwas_new_susie_alpha, gwas_new_susie_mu, gwas_new_susie_mu_sd, gwas_new_susie_component_num, gwas_new_susie_discover_component_bool, pass_bool = run_variant_level_susie_on_this_region(twas_data['gwas_beta'], twas_data['gwas_beta_se'], twas_data['gwas_sample_size'], twas_data['reference_ld'], gwas_susie_alpha)

		if pass_bool == False:
			t.write(component_name + '\t' + str(num_genes) + '\t' 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\n')
			continue

		# Convert gwas summary statistics to *STANDARDIZED* effect sizes
		# Following SuSiE code found in these two places:
		########1. https://github.com/stephenslab/susieR/blob/master/R/susie_rss.R  (LINES 277-279)
		########2. https://github.com/stephenslab/susieR/blob/master/R/susie_ss.R (LINES 148-156 AND 203-205)
		beta_scaled, beta_se_scaled, XtX = tgfm_causal_twas.convert_to_standardized_summary_statistics(twas_data['gwas_beta'], twas_data['gwas_beta_se'], twas_data['gwas_sample_size'], twas_data['reference_ld'])
		twas_data['gwas_beta'] = beta_scaled
		twas_data['gwas_beta_se'] = beta_se_scaled

		pred = np.sum(gwas_new_susie_alpha*gwas_new_susie_mu,axis=0)
		resid =twas_data['gwas_beta']- np.dot(XtX, pred)*(1.0/np.diag(XtX))

		# Extract ordered tissues for this component
		component_tissues = []
		for gene_full_name in twas_data['genes']:
			component_tissues.append('_'.join(gene_full_name.split('_')[1:]))
		component_tissues = np.asarray(component_tissues)


		# Extract ordered gene prior variances for this component
		ordered_gene_prior_variances = []
		for component_tissue in component_tissues:
			prior_var = tissue_specific_prior_expected_variance[tissue_to_position_mapping[component_tissue]]
			ordered_gene_prior_variances.append(prior_var)
		ordered_gene_prior_variances = np.asarray(ordered_gene_prior_variances)

		# RUN robust TGFM CAUSAL TWAS
		twas_obj = tgfm_robust_causal_twas.TGFM_CAUSAL_TWAS(estimate_prior_variance=False, prior_variance_vector=ordered_gene_prior_variances, pleiotropic_prior_variance=pleiotropic_variance, convergence_thresh=1e-6)
		twas_obj.fit(twas_data_obj=twas_data)
		# np.corrcoef(np.dot(XtX, twas_obj.beta_mu)*(1.0/np.diag(XtX)), np.dot(XtX, pred)*(1.0/np.diag(XtX)))
		

		# Prepare data for TGFM fine mapping
		#tgfm_fine_mapping_data = {'gwas_component_pmces': gwas_component_susie_pmces, 'gwas_susie_mu': gwas_susie_mu, 'gwas_susie_mu_sd': gwas_susie_mu_sd, 'gwas_susie_alpha': gwas_susie_alpha, 'eqtl_pmces': get_eqtl_component_level_pmces(twas_data['susie_mu'], twas_data['susie_alpha']), 'eqtl_mu':twas_data['susie_mu'], 'eqtl_mu_sd': twas_data['susie_mu_sd'], 'eqtl_alpha':twas_data['susie_alpha'], 'twas_alpha': twas_obj.alpha_mu, 'twas_alpha_sd': np.sqrt(twas_obj.alpha_var), 'genes': twas_data['genes'], 'variants': twas_data['variants'], 'ordered_tissue_names': ordered_tissue_names}

		# Run TGFM fine mapping
		#tgfm_fm_obj = tgfm_fine_mapping.TGFM_FM(likelihood_version='probabilistic_full', mean_component_boolean=False)
		#tgfm_fm_obj.fit(fm_data_obj=tgfm_fine_mapping_data)

		# Run RSS TGFM fine mapping
		tgfm_rss_fm_data = {'gwas_susie_mu': gwas_new_susie_alpha, 'gwas_susie_alpha': gwas_new_susie_alpha, 'gwas_susie_mu_sd':gwas_new_susie_mu_sd, 'gwas_susie_component': gwas_new_susie_component_num, 'twas_alpha': twas_obj.alpha_mu, 'twas_alpha_sd': np.sqrt(twas_obj.alpha_var), 'ordered_tissue_names': ordered_tissue_names}

		tgfm_rss_fm_obj_regr = tgfm_rss_fine_mapping.TGFM_RSS_FM(residual_version='regress')
		tgfm_rss_fm_obj_regr.fit(twas_data_obj=twas_data, tgfm_data_obj=tgfm_rss_fm_data)

		#tgfm_rss_fm_obj_pred = tgfm_rss_fine_mapping.TGFM_RSS_FM(residual_version='predict')
		#tgfm_rss_fm_obj_pred.fit(twas_data_obj=twas_data, tgfm_data_obj=tgfm_rss_fm_data)


		#tgfm_rss_fm_data = {'gwas_susie_mu': gwas_new_susie_mu, 'gwas_susie_alpha': gwas_new_susie_alpha, 'gwas_susie_mu_sd':gwas_new_susie_mu_sd, 'gwas_susie_component': gwas_new_susie_component_num, 'gwas_susie_XtX': XtX, 'twas_alpha': twas_obj.nominal_twas_rss_alpha_mu, 'twas_alpha_sd': np.sqrt(twas_obj.nominal_twas_rss_alpha_var), 'ordered_tissue_names': ordered_tissue_names}
		#tgfm_rss_fm_obj_nom_regr = tgfm_rss_fine_mapping.TGFM_RSS_FM(residual_version='regress')
		#tgfm_rss_fm_obj_nom_regr.fit(twas_data_obj=twas_data, tgfm_data_obj=tgfm_rss_fm_data)




		# Save TWAS results to output file
		tgfm_twas_pkl_file = output_dir + trait_name + '_' + component_name + '_' + gene_version + '_tgfm_robust_twas_tissue_specific_prior_results.pkl'
		g = open(tgfm_twas_pkl_file, "wb")
		pickle.dump(twas_obj, g)
		g.close()	


		# Save Fine-mapping results to output files
		#tgfm_fm_results_table_file = output_dir + trait_name + '_' + component_name + '_' + gene_version + '_tgfm_robust_fine_mapping_table.txt'
		#tgfm_fm_obj.posterior_prob_df.to_csv(tgfm_fm_results_table_file, sep='\t', index=False)
		#tgfm_tissue_fm_results_table_file = output_dir + trait_name + '_' + component_name + '_' + gene_version + '_tgfm_robust_tissue_fine_mapping_table.txt'
		#tgfm_fm_obj.tissue_posterior_prob_df.to_csv(tgfm_tissue_fm_results_table_file, sep='\t', index=False)
		# TGFM RSS Regression
		tgfm_rss_regression_fm_results_table_file = output_dir + trait_name + '_' + component_name + '_' + gene_version + '_tgfm_robust_rss_regression_fine_mapping_tissue_specific_prior_table.txt'
		tgfm_rss_fm_obj_regr.posterior_prob_df.to_csv(tgfm_rss_regression_fm_results_table_file, sep='\t', index=False)
		tgfm_rss_regression_tissue_fm_results_table_file = output_dir + trait_name + '_' + component_name + '_' + gene_version + '_tgfm_robust_rss_regression_tissue_fine_mapping_tissue_specific_prior_table.txt'
		tgfm_rss_fm_obj_regr.tissue_posterior_prob_df.to_csv(tgfm_rss_regression_tissue_fm_results_table_file, sep='\t', index=False)
		# TGFM RSS prediction
		#tgfm_rss_prediction_fm_results_table_file = output_dir + trait_name + '_' + component_name + '_' + gene_version + '_tgfm_robust_rss_prediction_fine_mapping_table.txt'
		#tgfm_rss_fm_obj_pred.posterior_prob_df.to_csv(tgfm_rss_prediction_fm_results_table_file, sep='\t', index=False)
		#tgfm_rss_prediction_tissue_fm_results_table_file = output_dir + trait_name + '_' + component_name + '_' + gene_version + '_tgfm_robust_rss_prediction_tissue_fine_mapping_table.txt'
		#tgfm_rss_fm_obj_pred.tissue_posterior_prob_df.to_csv(tgfm_rss_prediction_tissue_fm_results_table_file, sep='\t', index=False)
		# TGFM RSS NOMINAL Regression
		#tgfm_rss_regression_nom_fm_results_table_file = output_dir + trait_name + '_' + component_name + '_' + gene_version + '_tgfm_robust_rss_regression_nom_fine_mapping_table.txt'
		#tgfm_rss_fm_obj_nom_regr.posterior_prob_df.to_csv(tgfm_rss_regression_nom_fm_results_table_file, sep='\t', index=False)
		#tgfm_rss_regression_nom_tissue_fm_results_table_file = output_dir + trait_name + '_' + component_name + '_' + gene_version + '_tgfm_robust_rss_regression_nom_tissue_fine_mapping_table.txt'
		#tgfm_rss_fm_obj_nom_regr.tissue_posterior_prob_df.to_csv(tgfm_rss_regression_nom_tissue_fm_results_table_file, sep='\t', index=False)


		# Write output files to component level output
		t.write(component_name + '\t' + str(num_genes) + '\t' + tgfm_twas_pkl_file + '\t' + tgfm_rss_regression_fm_results_table_file + '\t' + tgfm_rss_regression_tissue_fm_results_table_file + '\t' + pickled_data_file + '\n')
f.close()
t.close()