import numpy as np 
import os
import sys
import pdb
import math
import scipy.stats


def log_sum(lnbfs, weights=None):
	# Taken from sumstats package
    """Sum of a sequence of bayes factors in logarithmic space
    # Basically this np.log(np.sum(np.exp(lnbfs)))
    Parameters
    ----------
    lnbfs
        sequence of (natural) log bayes factors
    weights
        sequence of weights for `lnbfs`
    
    Returns
    -------
    float
        the logarithm of the sum of the bayes factors
    """
    
    lnbfs = tuple(lnbfs)
    if not weights:
        weights = (1,) * len(lnbfs)
    max_lnbf = max(lnbfs)
    try:
        return (
            max_lnbf + math.log(
                math.fsum(
                   math.exp(lnbf - max_lnbf) * weight
                   for lnbf, weight in zip(lnbfs, weights)
                )
            )
        )
    except ValueError:
        if len(lnbfs) == 2:
            return min(lnbfs)
        else:
            raise RuntimeError(
                'The sum of absolute bayes factors may have been rounded to '
                'zero due to a pathalogically high maximum'
            )


def log_diff(x,y):
	my_max = np.max((x,y))
	res = my_max + np.log(np.exp(x-my_max) - np.exp(y-my_max))
	return res


def get_coloc_object(study_name, snp_names, sample_size, beta_vec, var_beta_vec):
	dicti = {}
	dicti['study_name'] = study_name
	dicti['snp_names'] = snp_names
	dicti['sample_size'] = sample_size
	dicti['beta'] = beta_vec
	dicti['var_beta'] = var_beta_vec
	return dicti

def get_approx_log_bf_estimates(coloc_object, sd_prior=0.15):
	z = coloc_object['beta']/np.sqrt(coloc_object['var_beta'])
	r = np.square(sd_prior)/(np.square(sd_prior) + coloc_object['var_beta'])
	lABF = 0.5*(np.log(1.0 - r) + (r*np.square(z)))
	return lABF

def get_snp_pph4_from_log_bf(lbf_1, lbf_2):
	internal_sum_lbf = lbf_1 + lbf_2
	denom_log_abf = log_sum(internal_sum_lbf)
	snp_pph4 = np.exp(internal_sum_lbf - denom_log_abf)
	return snp_pph4

def get_log_bayes_sums(l1, l2):
	lsum = l1 + l2

	lb_sum_h1 = log_sum(l1)
	lb_sum_h2 = log_sum(l2)
	lb_sum_h3 =  log_diff(log_sum(l1) + log_sum(l2), log_sum(lsum))
	lb_sum_h4 = log_sum(lsum)
	return lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4

def run_coloc_with_precomputed_log_bayes_sums(lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4, p0=1.0, p1=1e-4, p2=1e-4, p12=1e-5):
	lb_sum_h0 = 0.0

	joint_sum_h0 = lb_sum_h0 + 0.0
	joint_sum_h1 = lb_sum_h1 + np.log(p1) - np.log(p0)
	joint_sum_h2 = lb_sum_h2 + np.log(p2) - np.log(p0)
	joint_sum_h3 = lb_sum_h3 + np.log(p1) + np.log(p2) - (2.0*np.log(p0))
	joint_sum_h4 = lb_sum_h4 + np.log(p12) - np.log(p0)

	all_abf = np.asarray([joint_sum_h0, joint_sum_h1, joint_sum_h2, joint_sum_h3, joint_sum_h4])
	denom_log_abf = log_sum(all_abf)
	pps = np.exp(all_abf - denom_log_abf)
	return pps

def write_mat_to_output_file(mat, row_names, col_names, output_file):
	t = open(output_file,'w')
	t.write('\t'.join(col_names) + '\n')
	for row_index in range(len(row_names)):
		t.write(row_names[row_index] + '\t' + '\t'.join(mat[row_index,:]) + '\n')
	t.close()


def run_coloc_for_single_gene(gene_name, n_file, beta_file, std_err_file, tissue_names, tissue_name_to_position, coloc_output_dir, trait_name):
	# Coloc thresholds
	coloc_thresholds = np.asarray([.5, .7, .9, .95, .99])

	# Load in coloc data
	beta_df = np.loadtxt(beta_file, dtype=str,delimiter='\t')
	std_err_df = np.loadtxt(std_err_file, dtype=str,delimiter='\t')
	n_df = np.loadtxt(n_file, dtype=str,delimiter='\t')

	# Get studies observed for this gene
	observed_studies = beta_df[1:, 0]
	num_studies = len(observed_studies)
	eqtl_studies = observed_studies[:-1]

	# Get snp names
	snp_names = beta_df[0,1:]

	# Get object for trait for coloc
	trait_coloc_object = get_coloc_object(observed_studies[-1], snp_names, float(n_df[-1,1]), beta_df[-1,1:].astype(float), np.square(std_err_df[-1, 1:].astype(float)))


	# Make output matrix to pph values across tissues
	pph_mat = np.zeros((len(tissue_names), 5))
	pph_mat[:,0] = 1.0

	# Keep track of log bayes output matrix
	log_bayes_mat = np.reshape(['NULL']*(len(tissue_names)*5), (len(tissue_names), 5)).astype('U16')
	log_bayes_mat[:,-1] = str(len(snp_names))

	# Keep track of snp pph4s
	snp_pph4_mat = np.zeros((len(tissue_names), len(snp_names)))


	# Initilize output
	predicted_effects_list = []
	coloc_at_threshold_arr = []
	for coloc_threshold in coloc_thresholds:
		coloc_thresh_mat = np.zeros((len(snp_names), len(tissue_names)))
		predicted_effects_list.append(coloc_thresh_mat)
		coloc_at_threshold_arr.append(False)


	# Loop through eqtl studies
	for eqtl_study_num, eqtl_study_name in enumerate(eqtl_studies):
		# Global tissue position of this eqtl study
		global_tissue_position = tissue_name_to_position[eqtl_study_name]

		# Get object for eqtl study for coloc
		eqtl_coloc_object = get_coloc_object(observed_studies[eqtl_study_num], snp_names, float(n_df[(eqtl_study_num+1),1]), beta_df[(eqtl_study_num+1),1:].astype(float), np.square(std_err_df[(eqtl_study_num+1), 1:].astype(float)))

		# Estimate log bayes factors for each study
		trait_approx_lbf = get_approx_log_bf_estimates(trait_coloc_object)
		eqtl_approx_lbf = get_approx_log_bf_estimates(eqtl_coloc_object)
		
		# Esimate SNP PPH4 from log bayes factors
		snp_pph4 = get_snp_pph4_from_log_bf(eqtl_approx_lbf, trait_approx_lbf)
		snp_pph4[snp_pph4 < .01] = 0.0
		snp_pph4_mat[global_tissue_position,:] = snp_pph4

		# Get log bayes sums (basically summary stats for the gene
		lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4 = get_log_bayes_sums(eqtl_approx_lbf, trait_approx_lbf)

		# Keep track log sum bayes for this gene in each tissue
		log_bayes_mat[global_tissue_position,:-1] = np.asarray([lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4]).astype(str)

		# Keep track of Coloc probabilities
		pph_vec = run_coloc_with_precomputed_log_bayes_sums(lb_sum_h1, lb_sum_h2, lb_sum_h3, lb_sum_h4)
		pph_mat[global_tissue_position,:] = pph_vec

		# Compute predicted effects at each threshold
		for threshold_iter, coloc_threshold in enumerate(coloc_thresholds):
			if pph_vec[4] > coloc_threshold:
				coloc_at_threshold_arr[threshold_iter] = True
				predicted_effects_list[threshold_iter][:, global_tissue_position] = snp_pph4*pph_vec[4]*trait_coloc_object['beta']

	# Save snp pph4 mat to output file
	snp_pph4_output_file = coloc_output_dir + trait_name + '_' + gene_name + '_snp_pph4.txt'
	write_mat_to_output_file(snp_pph4_mat.astype(str), tissue_names, snp_names, snp_pph4_output_file)


	# Save pph mat
	pph_mat_output_file = coloc_output_dir + trait_name + '_' + gene_name + '_coloc_posterior_probabilities.txt'
	write_mat_to_output_file(pph_mat.astype(str), tissue_names, np.asarray(["PPH0", "PPH1", "PPH2", "PPH3", "PPH4"]), pph_mat_output_file)

	# Save predicted effects mat
	for threshold_iter, coloc_threshold in enumerate(coloc_thresholds):
		if coloc_at_threshold_arr[threshold_iter]:
			predicted_effect_size_file = coloc_output_dir + trait_name + '_' + gene_name + '_coloc_' + str(coloc_threshold) + '_predicted_effect_sizes.txt'
			write_mat_to_output_file(predicted_effects_list[threshold_iter].astype(str), snp_names, tissue_names, predicted_effect_size_file)

	# Save log bayes sums
	lb_mat_output_file = coloc_output_dir + trait_name + '_' + gene_name + '_coloc_log_bayes_sums.txt'
	write_mat_to_output_file(log_bayes_mat, tissue_names, np.asarray(["LB1", "LB2", "LB3", "LB4", "num_snps"]), lb_mat_output_file)


gene_file = sys.argv[1]
trait_name = sys.argv[2]
gtex_tissue_file = sys.argv[3]
coloc_output_dir = sys.argv[4]
job_number = float(sys.argv[5])
total_jobs = float(sys.argv[6])



# Load in tissue names
tissue_df = np.loadtxt(gtex_tissue_file,dtype=str,delimiter='\t')
tissue_names = tissue_df[1:,0]

# Create mapping from tissue name to position
tissue_name_to_position = {}
for index, tissue_name in enumerate(tissue_names):
	tissue_name_to_position[tissue_name] = index


# Load in gene data frame
gene_df = np.loadtxt(gene_file, dtype=str,delimiter='\t')
gene_df_header = gene_df[0,:]
gene_df = gene_df[1:,:]
# Get total number of genes
num_genes = gene_df.shape[0]


# For parallelization purposes, determine which genes to test in this thread
tasks_per_job = np.floor(num_genes/total_jobs) + 1
start_task = int(np.floor(job_number*tasks_per_job + 1) - 1)
end_task = int(np.floor((job_number + 1)*tasks_per_job) -1)
if end_task > (num_genes - 1):
	end_task = num_genes -1




# Loop through all genes and run colocalization for that gene
for gene_num in range(start_task, (end_task+1)):
	print(gene_num)
	gene_name = gene_df[gene_num, 0]
	n_file = gene_df[gene_num, 4]
	beta_file = gene_df[gene_num, 5]
	std_err_file = gene_df[gene_num, 6]

	# Run coloc for a single gene
	run_coloc_for_single_gene(gene_name, n_file, beta_file, std_err_file, tissue_names, tissue_name_to_position, coloc_output_dir, trait_name)


