import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import pickle
import tgfm
import rpy2
import rpy2.robjects.numpy2ri as numpy2ri
import rpy2.robjects as ro
import scipy.stats
ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
susieR_pkg = importr('susieR')

def compute_pips(alpha_mat):
	LL = alpha_mat.shape[0]
	n_elements = alpha_mat.shape[1]
	anti_pips = np.ones(n_elements)

	for component_iter in range(LL):
		anti_pips = anti_pips*(1.0 - alpha_mat[component_iter,:])
	pips = 1.0 - anti_pips
	return pips


def susie_inference_shell(tgfm_data, ld_mat, init_method, est_resid_var_bool):
	# More hack: need to redo twas z
	variant_z = tgfm_data['gwas_beta']/tgfm_data['gwas_beta_se']

	susie_variant_only = susieR_pkg.susie_rss(z=variant_z.reshape((len(variant_z),1)), R=ld_mat, n=tgfm_data['gwas_sample_size'], L=15, estimate_residual_variance=est_resid_var_bool)

	susie_variant_only_obj = {'alpha':susie_variant_only.rx2('alpha'), 'mu':susie_variant_only.rx2('mu'),'mu2':susie_variant_only.rx2('mu2')}

	susie_variant_only_obj['pips'] = compute_pips(susie_variant_only.rx2('alpha'))

	return susie_variant_only_obj


def extract_gene_variant_ld(standardized_eqtl_effects, variant_ld):
	expression_covariance = np.dot(np.dot(standardized_eqtl_effects, variant_ld), np.transpose(standardized_eqtl_effects))
	np.fill_diagonal(expression_covariance, 1.0)
	dd = np.diag(1.0/np.sqrt(np.diag(expression_covariance)))
	ge_ld = np.dot(np.dot(dd, expression_covariance),dd)
	gene_variant_ld = np.dot(standardized_eqtl_effects,variant_ld) # Ngenes X n_variants
	return gene_variant_ld


def extract_variants_tagged_by_fine_mapped_genes(gene_variant_ld, tgfm_alpha_pips, gene_pip_thresh, abs_corr_thresh):
	n_var = gene_variant_ld.shape[1]
	n_gene = gene_variant_ld.shape[0]

	tagged_variants = np.asarray(['False']*n_var)

	for gene_iter in range(n_gene):
		if tgfm_alpha_pips[gene_iter] < gene_pip_thresh:
			continue
		tagged_variants[np.abs(gene_variant_ld[gene_iter,:]) > abs_corr_thresh] = 'True'
	return tagged_variants




######################
# Command line args
######################
tgfm_input_summary_file = sys.argv[1]
tgfm_output_stem = sys.argv[2]
init_method = sys.argv[3]
est_resid_var_str = sys.argv[4]
ln_pi_method_name = sys.argv[5]



window_pvalue_thresh = 1e-5



if est_resid_var_str == 'False':
	est_resid_var_bool = False
elif est_resid_var_str == 'True':
	est_resid_var_bool = True 
else:
	print('assumption eroror')
	pdb.set_trace()



# Open PIP file handle
pip_output_file = tgfm_output_stem + '_tgfm_pip_summary.txt'
t_pip = open(pip_output_file,'w')
t_pip.write('window_name\tinclusion_elements\tinclusion_probabilities\n')




# Now loop through windows
# In each window run TGFM independently
# Loop through trait components
tgfm_input_data = np.loadtxt(tgfm_input_summary_file,dtype=str,delimiter='\t')
tgfm_input_data = tgfm_input_data[1:,:]



# Get n_windows on this run
n_windows = tgfm_input_data.shape[0]

for window_iter in range(n_windows):
	data = tgfm_input_data[window_iter, :]

	##############################
	# Extract relevent fields
	###############################
	window_name = data[0]
	print(window_name)



	ld_file = data[1]
	tgfm_input_pkl = data[2]
	ln_pi_file_stem = data[3]

	##############################
	# Load in Data
	###############################
	# Load in tgfm input data
	g = open(tgfm_input_pkl, "rb")
	tgfm_data = pickle.load(g)
	g.close()

	# Extract gwas p
	gwas_z = tgfm_data['gwas_beta']/tgfm_data['gwas_beta_se']
	gwas_p = scipy.stats.norm.sf(abs(gwas_z))*2.0

	# Ignore windows with no pvalues less than some threshold
	if np.min(gwas_p) > window_pvalue_thresh:
		print('skipped because of window pvalue threshold')
		t_pip.write(window_name + '\tNA\tNA\n')
		continue
	# Skip windows with no genes
	if len(tgfm_data['genes']) == 0:
		print('skipped because of no genes')
		t_pip.write(window_name + '\tNA\tNA\n')
		continue

	# Load in LD
	ld_mat = np.load(ld_file)
	# Add ld to tgfm_data obj
	tgfm_data['reference_ld'] = ld_mat

	# Extract full ld between genes, variants, and gene-variants
	gene_variant_ld = extract_gene_variant_ld(tgfm_data['gene_eqtl_pmces'], tgfm_data['reference_ld'])

	##############################
	# Run TGFM
	###############################
	susie_variant_only_obj = susie_inference_shell(tgfm_data, ld_mat, init_method, est_resid_var_bool)
	# Add other relevent info
	susie_variant_only_obj['middle_variant_indices'] = tgfm_data['middle_variant_indices']


	# Write to credible set output
	# More of just a place holder here
	t_pip.write(window_name + '\tNA\tNA\n')
	t_pip.flush()

	# Load in TGFM results file
	tmp = tgfm_output_stem.split('susie_variant_only')[0] + 'susie_sampler_pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped_' + window_name + '_results.pkl'
	g = open(tmp, "rb")
	tmp_res = pickle.load(g)
	g.close()


	print(np.corrcoef(tmp_res['expected_beta_pips'], susie_variant_only_obj['pips'])[0,1])


	# Extract variants tagged by fine-mapped genes
	valid_indices_mat = []
	valid_indices_types = []
	abs_corr_threshs = [.25, .5, .75, .9]
	gene_pip_threshs = [.05, .1, .25, .5]
	for abs_corr_thresh in abs_corr_threshs:
		for gene_pip_thresh in gene_pip_threshs:
			tagged_variants = extract_variants_tagged_by_fine_mapped_genes(gene_variant_ld, tmp_res['expected_alpha_pips'], gene_pip_thresh, abs_corr_thresh)
			valid_indices = tagged_variants=='False'
			valid_indices_mat.append(valid_indices)
			valid_indices_types.append(str(abs_corr_thresh) + '_' + str(gene_pip_thresh))
	valid_indices_mat = np.asarray(valid_indices_mat)
	valid_indices_types = np.asarray(valid_indices_types)
	susie_variant_only_obj['valid_indices_mat'] = valid_indices_mat
	susie_variant_only_obj['valid_indices_types'] = valid_indices_types



	print(np.corrcoef(tmp_res['expected_beta_pips'][valid_indices_mat[0,:]], susie_variant_only_obj['pips'][valid_indices_mat[0,:]])[0,1])


	# Write pickle file
	window_tgfm_output_file = tgfm_output_stem + '_' + window_name + '_results.pkl'
	g = open(window_tgfm_output_file, "wb")
	pickle.dump(susie_variant_only_obj, g)
	g.close()

t_pip.close()



