import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import pandas as pd
import os
import pdb
from sklearn.linear_model import LinearRegression
from sparse_sldsc_from_multivariate_summary_statistics import SPARSE_SLDSC
from sparse_sldsc_from_multivariate_summary_statistics_with_some_fixed_effects import SPARSE_SLDSC_SOME_FIXED
from sparse_sldsc_from_multivariate_summary_statistics_fixed_genotype import SPARSE_SLDSC_FIXED_TERM
from sparse_sldsc_from_multivariate_summary_statistics_ard_with_some_fixed_effects import SPARSE_SLDSC_ARD_SOME_FIXED
from sparse_sldsc_from_multivariate_summary_statistics_ard_with_some_fixed_effects_mv_updates import SPARSE_SLDSC_ARD_SOME_FIXED_MV_UPDATES
from sparse_sldsc_from_multivariate_summary_statistics_ard_with_some_fixed_effects_mv_updates_eqtl_intercept import SPARSE_SLDSC_ARD_SOME_FIXED_MV_UPDATES_EQTL_INTERCEPT


def get_anno_names(example_anno_file):
	f = open(example_anno_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		anno_names = np.asarray(data[3:])
		break
	f.close()
	return anno_names




def extract_tau_and_tau_se_from_log_file(sldsc_log_file, anno_names):
	f = open(sldsc_log_file)
	for line in f:
		line = line.rstrip()
		if line.startswith('Coefficients:'):
			coef = np.asarray(line.split()[1:]).astype(float)
			if len(coef) < len(anno_names):
				#line = f.next()
				line = next(f)
				data = np.asarray(line.split()).astype(float)
				coef = np.hstack([coef, data])
		if line.startswith('Coefficient SE:'):
			coef_se = np.asarray(line.split()[2:]).astype(float)
			if len(coef_se) < len(anno_names):
				#line = f.next()
				line = next(f)
				data = np.asarray(line.split()).astype(float)
				coef_se = np.hstack([coef_se, data])
	return coef, coef_se

def compute_jacknifed_covariance_matrix(jacknifed_taus):
	jacknife_mean = np.mean(jacknifed_taus,axis=0)
	diff = jacknifed_taus - jacknife_mean
	num_jacknife_samples = jacknifed_taus.shape[0]

	jacknifed_cov = np.dot(np.transpose(diff),diff)*(num_jacknife_samples-1.0)/num_jacknife_samples

	#from sklearn.covariance import GraphicalLassoCV, ledoit_wolf
	#X = diff*(num_jacknife_samples-1.0)
	#model = GraphicalLassoCV(max_iter=3000)
	#model.fit(X)

	return jacknifed_cov, jacknife_mean

def extract_non_eqtl_annotations(anno_names):
	non_eqtl_indices = []
	eqtl_indices = []
	for ii, anno_name in enumerate(anno_names):
		if anno_name.endswith('L2') or anno_name.startswith('intercept'):
			non_eqtl_indices.append(ii)
		else:
			eqtl_indices.append(ii)
	return np.asarray(non_eqtl_indices), np.asarray(eqtl_indices)

def load_in_annotation_sdev_data(annotation_sd_file):
	aa = np.loadtxt(annotation_sd_file,dtype=str,delimiter='\t')
	return aa[1:,0], aa[1:,1].astype(float)

def load_in_mvec(anno_stem, suffix):
	for chrom_num in range(1,23):
		filer = anno_stem + '.' + str(chrom_num) + '.l2' + suffix
		tmp = np.loadtxt(filer)
		if chrom_num == 1:
			counter = np.zeros(len(tmp))
		counter = counter + tmp
	return counter

def print_organized_summary_file(output_file_name, anno_names, tau, tau_se):
	tau_z = tau/tau_se
	t = open(output_file_name,'w')
	t.write('Annotation\ttau\ttau_se\ttau_z\n')
	for ii in range(len(anno_names)):
		t.write(anno_names[ii] + '\t' + str(tau[ii]) + '\t' + str(tau_se[ii]) + '\t' + str(tau_z[ii]) + '\n')
	t.close()

def extract_expression_mediated_h2(jacknifed_taus, eqtl_start_index, m_vec):
	jacknifed_med_h2 = jacknifed_taus*m_vec

	jacknifed_geno_h2 = np.sum(jacknifed_med_h2[:,:eqtl_start_index],axis=1)
	jacknifed_expr_h2 = np.sum(jacknifed_med_h2[:,eqtl_start_index:],axis=1)

	jacknifed_expr_med_h2 = jacknifed_expr_h2/(jacknifed_expr_h2 + jacknifed_geno_h2)

	mean_expr_med_h2 = np.mean(jacknifed_expr_med_h2)

	diff = jacknifed_expr_med_h2 - mean_expr_med_h2
	num_jacknife_samples = jacknifed_expr_med_h2.shape[0]

	jacknifed_se = np.sqrt(np.dot(np.transpose(diff),diff)*(num_jacknife_samples-1.0)/num_jacknife_samples)
	
	return mean_expr_med_h2, jacknifed_se


def jacknife_mean_and_se(veccy):
	mean_veccy = np.mean(veccy)

	diff = veccy - mean_veccy
	num_jacknife_samples = veccy.shape[0]

	jacknifed_se = np.sqrt(np.dot(np.transpose(diff),diff)*(num_jacknife_samples-1.0)/num_jacknife_samples)


	return mean_veccy, jacknifed_se


def print_average_per_snp_and_per_gene_tissue_h2_for_sparse_model(taus, eqtl_start_index, m_vec, output_file, anno_names):
	# Compute per snp h2
	per_snp_h2 = np.sum((m_vec[:eqtl_start_index])*(taus[:eqtl_start_index]))/(m_vec[0])

	tissue_names = anno_names[eqtl_start_index:]
	n_gene_tiss = m_vec[eqtl_start_index:]

	per_gene_tissue_h2 = taus[eqtl_start_index:]

	# Print to output
	t = open(output_file,'w')
	t.write('genetic_element_name\tper_genetic_element_h2\tper_genetic_element_h2_se\n')
	t.write('variant\t' + str(per_snp_h2) + '\t' + str('NA') + '\n')
	for tissue_iter, tissue_name in enumerate(tissue_names):
		t.write(tissue_name + '\t' + str(per_gene_tissue_h2[tissue_iter]) + '\t' + str('NA') + '\n')
	t.close()
	return




def print_average_per_snp_and_per_gene_tissue_h2(jacknifed_taus, eqtl_start_index, m_vec, output_file, anno_names):
	# Compute total heritability going through snps and through gene expression
	jacknifed_med_h2 = jacknifed_taus*m_vec
	jacknifed_geno_h2 = np.sum(jacknifed_med_h2[:,:eqtl_start_index],axis=1)
	jacknifed_tiss_expr_h2 = jacknifed_med_h2[:,eqtl_start_index:]

	# Get number of snps
	n_snps = m_vec[0]
	n_gene_tiss = m_vec[eqtl_start_index:]
	tissue_names = anno_names[eqtl_start_index:]

	# Per snp and per-gene h2 across jacknifed samples
	per_snp_h2 = jacknifed_geno_h2/n_snps
	per_gene_tissue_h2 = jacknifed_tiss_expr_h2/n_gene_tiss

	# Get mean and se from jacknife samples
	avg_per_snp_h2, avg_per_snp_h2_se = jacknife_mean_and_se(per_snp_h2)

	# Print to output
	t = open(output_file,'w')
	t.write('genetic_element_name\tper_genetic_element_h2\tper_genetic_element_h2_se\n')
	t.write('variant\t' + str(avg_per_snp_h2) + '\t' + str(avg_per_snp_h2_se) + '\n')
	for tissue_iter, tissue_name in enumerate(tissue_names):
		avg_per_gene_h2, avg_per_gene_h2_se = jacknife_mean_and_se(per_gene_tissue_h2[:, tissue_iter])
		t.write(tissue_name + '\t' + str(avg_per_gene_h2) + '\t' + str(avg_per_gene_h2_se) + '\n')
	t.close()
	return


def print_average_per_snp_and_per_gene_h2(jacknifed_taus, eqtl_start_index, m_vec, output_file):
	# Compute total heritability going through snps and through gene expression
	jacknifed_med_h2 = jacknifed_taus*m_vec
	jacknifed_geno_h2 = np.sum(jacknifed_med_h2[:,:eqtl_start_index],axis=1)
	jacknifed_expr_h2 = np.sum(jacknifed_med_h2[:,eqtl_start_index:],axis=1)

	# Get number of snps
	n_snps = m_vec[0]
	n_genes = np.sum(m_vec[eqtl_start_index:])

	# Per snp and per-gene h2 across jacknifed samples
	per_snp_h2 = jacknifed_geno_h2/n_snps
	per_gene_h2 = jacknifed_expr_h2/n_genes

	# Get mean and se from jacknife samples
	avg_per_snp_h2, avg_per_snp_h2_se = jacknife_mean_and_se(per_snp_h2)
	avg_per_gene_h2, avg_per_gene_h2_se = jacknife_mean_and_se(per_gene_h2)

	# Print to output
	t = open(output_file,'w')
	t.write('genetic_element_name\tper_genetic_element_h2\tper_genetic_element_h2_se\n')
	t.write('variant\t' + str(avg_per_snp_h2) + '\t' + str(avg_per_snp_h2_se) + '\n')
	t.write('gene\t' + str(avg_per_gene_h2) + '\t' + str(avg_per_gene_h2_se) + '\n')
	t.close()
	return


def print_organized_h2_mediated(output_file_name, h2_med, h2_med_se):
	t = open(output_file_name,'w')
	t.write('h2_med\th2_med_se\n')
	t.write(str(h2_med) + '\t' + str(h2_med_se) + '\n')
	t.close()

def extract_and_print_per_element_h2_for_nonnegative_bootstrapped_models(nonnegative_bootstrapped_coefficient_file, m_vec, nonnegative_bootstrapped_per_element_h2_file):
	# Extract data
	global_taus = []
	bootstrapped_taus = []
	anno_names = []
	f = open(nonnegative_bootstrapped_coefficient_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		anno_names.append(data[0].split('_0')[0])
		global_taus.append(float(data[1]))
		bootstrapped_taus.append(np.asarray(data[4].split(';')).astype(float))
	f.close()
	anno_names = np.asarray(anno_names)
	global_taus = np.asarray(global_taus)
	bootstrapped_taus = np.asarray(bootstrapped_taus)

	# Extract eqtl start index
	for ii, ele_name in enumerate(anno_names):
		if ele_name.endswith('L2') == False:
			eqtl_start_index = ii
			break

	# Extract per element h2
	# For variants
	global_per_snp_h2 = np.sum(global_taus[:eqtl_start_index]*m_vec[:eqtl_start_index])/m_vec[0]
	bootstrapped_per_snp_h2 = np.dot(np.transpose(bootstrapped_taus[:eqtl_start_index,:]), m_vec[:eqtl_start_index])/m_vec[0]
	# For genes
	global_per_gene_h2 = np.copy(global_taus[eqtl_start_index:])
	bootstrapped_per_gene_h2 = np.copy(bootstrapped_taus[eqtl_start_index:,:])
	tissue_names = anno_names[eqtl_start_index:]

	# Print to output
	t = open(nonnegative_bootstrapped_per_element_h2_file,'w')
	# Header
	t.write('element_name\tglobal_per_element_h2\tbootstrapped_per_element_h2\n')
	# For variants
	t.write('variant\t' + str(global_per_snp_h2) + '\t' + ';'.join(bootstrapped_per_snp_h2.astype(str)) + '\n')
	for tissue_iter, tissue_name in enumerate(tissue_names):
		t.write(tissue_name + '\t' + str(global_per_gene_h2[tissue_iter]) + '\t' + ';'.join(bootstrapped_per_gene_h2[tissue_iter,:].astype(str)) + '\n')
	t.close()
	return

sldsc_output_root = sys.argv[1]
anno_stem = sys.argv[2]

# Extract relevent info from log file
sldsc_log_file = sldsc_output_root + '.log'
example_anno_file = anno_stem + '.21.l2.ldscore'
anno_names = get_anno_names(example_anno_file)
tau, tau_se = extract_tau_and_tau_se_from_log_file(sldsc_log_file, anno_names)
tau_z = tau/tau_se
# Print to output
print_organized_summary_file(sldsc_output_root + 'organized_res.txt', anno_names, tau, tau_se)


# Total counts
m_vec = load_in_mvec(anno_stem, '.M')
m_5_50_vec = load_in_mvec(anno_stem, '.M_5_50')


# Get organized per element h2 for non-negative bootstrapped model
nonnegative_bootstrapped_coefficient_file = sldsc_output_root + 'nonnegative_eqtl_bootstrapped_sldsc_coefficients.txt'
nonnegative_bootstrapped_per_element_h2_file = sldsc_output_root + 'nonnegative_eqtl_bootstrapped_sldsc_per_element_h2.txt'
extract_and_print_per_element_h2_for_nonnegative_bootstrapped_models(nonnegative_bootstrapped_coefficient_file, m_vec, nonnegative_bootstrapped_per_element_h2_file)
nonnegative_bootstrapped_per_element_h2_5_50_file = sldsc_output_root + 'nonnegative_eqtl_bootstrapped_sldsc_per_element_h2_5_50.txt'
extract_and_print_per_element_h2_for_nonnegative_bootstrapped_models(nonnegative_bootstrapped_coefficient_file, m_5_50_vec, nonnegative_bootstrapped_per_element_h2_5_50_file)

# Print partioned h2 to output
print_organized_summary_file(sldsc_output_root + 'organized_mediated_h2.txt', anno_names, tau*m_vec, tau_se*m_vec)
print_organized_summary_file(sldsc_output_root + 'organized_mediated_h2_5_50.txt', anno_names, tau*m_5_50_vec, tau_se*m_5_50_vec)


# Load in annotation names and annotation sdevs
#annotation_sd_file = anno_stem + '_annotation_sdev.txt'
#annotation_names, annotation_sdev = load_in_annotation_sdev_data(annotation_sd_file)


# Load in jacknifed data
jacknifed_taus_file = sldsc_output_root + '.part_delete'
jacknifed_taus = np.loadtxt(jacknifed_taus_file)

# Load in Jacknifed intercept data
jacknifed_intercept_file = sldsc_output_root + '.intercept_delete'
jacknifed_intercept = np.loadtxt(jacknifed_intercept_file)

# Compute standardized jacknifed mean tau and covarianc
#annotation_sdev = annotation_sdev*300000.0
#annotation_sdev = (annotation_sdev*0.0 + 1.0)*300000.0
#standardized_jacknifed_taus = jacknifed_taus*annotation_sdev

######################
# Add intercept
#standardized_jacknifed_taus = np.hstack((jacknifed_intercept.reshape(len(jacknifed_intercept),1), standardized_jacknifed_taus))
#annotation_sdev = np.hstack([np.asarray([1.0]), annotation_sdev])
#annotation_names = np.hstack([np.asarray(['intercept']), annotation_names])
anno_names = np.hstack([np.asarray(['intercept']), anno_names])
######################



#jacknifed_tau_covariance, jacknifed_tau_mean = compute_jacknifed_covariance_matrix(standardized_jacknifed_taus)

#jacknifed_tau_z = jacknifed_tau_mean/np.sqrt(np.diag(jacknifed_tau_covariance))


# Extract annotations to be fixed (ie the non-eqtl annotations)
non_eqtl_annotations, eqtl_annotations = extract_non_eqtl_annotations(anno_names)
eqtl_start_index = np.min(eqtl_annotations)

# Jacknife expression mediated h2
h2_med, h2_med_se = extract_expression_mediated_h2(jacknifed_taus, eqtl_start_index-1, m_vec)
h2_5_50_med, h2_5_50_med_se = extract_expression_mediated_h2(jacknifed_taus, eqtl_start_index-1, m_5_50_vec)
print_organized_h2_mediated(sldsc_output_root + 'h2_med.txt', h2_med, h2_med_se)
print_organized_h2_mediated(sldsc_output_root + 'h2_5_50_med.txt', h2_5_50_med, h2_5_50_med_se)


# Compute average per-snp and per-gene h2
print_average_per_snp_and_per_gene_h2(jacknifed_taus, eqtl_start_index-1, m_vec, sldsc_output_root + 'avg_per_snp_and_gene_h2.txt')
print_average_per_snp_and_per_gene_h2(jacknifed_taus, eqtl_start_index-1, m_5_50_vec, sldsc_output_root + 'avg_per_snp_and_gene_h2_5_50.txt')

# Compute average per-snp and per-gene-tissue h2
print_average_per_snp_and_per_gene_tissue_h2(jacknifed_taus, eqtl_start_index-1, m_vec, sldsc_output_root + 'avg_per_snp_and_gene_tissue_h2.txt', anno_names[1:])
print_average_per_snp_and_per_gene_tissue_h2(jacknifed_taus, eqtl_start_index-1, m_5_50_vec, sldsc_output_root + 'avg_per_snp_and_gene_tissue_h2_5_50.txt', anno_names[1:])


# TEMP HACK
#jacknifed_tau_mean2 = np.loadtxt(sldsc_output_root + 'coef_estimate.txt')
#jacknifed_tau_covariance2 = np.loadtxt(sldsc_output_root + 'coef_cov_estimate.txt')
#annotation_sdev_full = np.hstack((annotation_sdev, [1]))
#jacknifed_tau_mean2 = jacknifed_tau_mean2*annotation_sdev_full
#jacknifed_tau_covariance2 = np.dot(np.dot(np.diag(annotation_sdev_full), jacknifed_tau_covariance2), np.diag(annotation_sdev_full))


# Susie sparse updates
#sparse_sldsc_obj = SPARSE_SLDSC_SOME_FIXED(max_iter=20000, L=10)
#sparse_sldsc_obj.fit(tau=jacknifed_tau_mean, tau_cov=jacknifed_tau_covariance, fixed_coefficients=non_eqtl_annotations)

#################
# Multivariate updates (MV) UPDATES ACROSS range of Regularization parameters
#for reg_param in [1e-2, 1e-1, 5e-1, 1.0, 2.0, 5.0, 10.0]:
'''
for reg_param in [5e-1]:
	print(reg_param)

	non_eqtl_annotations_include_intercept = np.hstack([non_eqtl_annotations])
	sparse_sldsc_obj = SPARSE_SLDSC_ARD_SOME_FIXED_MV_UPDATES_EQTL_INTERCEPT(max_iter=10000, L=10, nonneg=False, nonneg_int=eqtl_start_index, regularization_param=reg_param)
	sparse_sldsc_obj.fit(tau=jacknifed_tau_mean, tau_cov=jacknifed_tau_covariance, fixed_coefficients=non_eqtl_annotations_include_intercept)
	model_beta_mu = sparse_sldsc_obj.beta_mu
	model_beta_var = np.diag(sparse_sldsc_obj.beta_cov)
	
	eqtl_coef = (sparse_sldsc_obj.alpha_mu+model_beta_mu[eqtl_annotations])/annotation_sdev[eqtl_annotations]
	print(np.sort(eqtl_coef))
	print(anno_names[eqtl_annotations][np.argsort(eqtl_coef)])
	pdb.set_trace()
	#print_organized_summary_file(sldsc_output_root + 'organized_' + str(reg_param) + '_sparse_ard_eqtl_coefficients_mv_update_res.txt', anno_names, model_beta_mu/annotation_sdev, np.sqrt(model_beta_var)/annotation_sdev)


	non_eqtl_annotations_include_intercept = np.hstack([non_eqtl_annotations])
	sparse_sldsc_obj = SPARSE_SLDSC_ARD_SOME_FIXED_MV_UPDATES(max_iter=20000, L=10, nonneg=False, nonneg_int=eqtl_start_index, regularization_param=reg_param)
	sparse_sldsc_obj.fit(tau=jacknifed_tau_mean, tau_cov=jacknifed_tau_covariance, fixed_coefficients=non_eqtl_annotations_include_intercept)
	model_beta_mu = sparse_sldsc_obj.beta_mu
	model_beta_var = np.diag(sparse_sldsc_obj.beta_cov)
	print_organized_summary_file(sldsc_output_root + 'organized_' + str(reg_param) + '_sparse_ard_eqtl_coefficients_mv_update_res.txt', anno_names, model_beta_mu/annotation_sdev, np.sqrt(model_beta_var)/annotation_sdev)
	# Compute average per-snp and per-gene-tissue h2 for sparse model
	print_average_per_snp_and_per_gene_tissue_h2_for_sparse_model((model_beta_mu/annotation_sdev)[1:], eqtl_start_index-1, m_vec, sldsc_output_root + str(reg_param) + '_sparse_ard_eqtl_coefficients_mv_update_avg_per_snp_and_gene_tissue_h2.txt', anno_names[1:])

	sparse_sldsc_obj = SPARSE_SLDSC_ARD_SOME_FIXED_MV_UPDATES(max_iter=20000, L=10, nonneg=False, nonneg_int=eqtl_start_index, regularization_param=reg_param)
	sparse_sldsc_obj.fit(tau=jacknifed_tau_mean, tau_cov=jacknifed_tau_covariance, fixed_coefficients=np.asarray([]))
	model_beta_mu = sparse_sldsc_obj.beta_mu
	model_beta_var = np.diag(sparse_sldsc_obj.beta_cov)
	print_organized_summary_file(sldsc_output_root + 'organized_' + str(reg_param) + '_sparse_ard_all_coefficients_mv_update_res.txt', anno_names, model_beta_mu/annotation_sdev, np.sqrt(model_beta_var)/annotation_sdev)
	# Compute average per-snp and per-gene-tissue h2 for sparse model
	print_average_per_snp_and_per_gene_tissue_h2_for_sparse_model((model_beta_mu/annotation_sdev)[1:], eqtl_start_index-1, m_vec, sldsc_output_root + str(reg_param) + '_sparse_ard_all_coefficients_mv_update_avg_per_snp_and_gene_tissue_h2.txt', anno_names[1:])

'''



'''
# Univariate UPDATES ACROSS range of Regularization parameters
for reg_param in [1e-20, 1e-2, 1e-1, 1.0, 2.0, 5.0, 10.0, 50.0]:
	print(reg_param)
	non_eqtl_annotations_include_intercept = np.hstack([non_eqtl_annotations,[131]])
	sparse_sldsc_obj = SPARSE_SLDSC_ARD_SOME_FIXED(max_iter=30000, L=10, nonneg=False, nonneg_int=eqtl_start_index, regularization_param=reg_param)
	sparse_sldsc_obj.fit(tau=jacknifed_tau_mean2, tau_cov=jacknifed_tau_covariance2, fixed_coefficients=non_eqtl_annotations_include_intercept)
	model_beta_mu = np.hstack((sparse_sldsc_obj.fixed_beta_mu[:-1], sparse_sldsc_obj.beta_mu))
	model_beta_var = np.hstack((sparse_sldsc_obj.fixed_beta_var[:-1], sparse_sldsc_obj.beta_var))
	print_organized_summary_file(sldsc_output_root + 'organized_' + str(reg_param) + '_sparse_ard_eqtl_coefficients_res.txt', anno_names, model_beta_mu/annotation_sdev, np.sqrt(model_beta_var)/annotation_sdev)


	sparse_sldsc_obj = SPARSE_SLDSC_ARD_SOME_FIXED(max_iter=30000, L=10, nonneg=False, nonneg_int=eqtl_start_index, regularization_param=reg_param)
	sparse_sldsc_obj.fit(tau=jacknifed_tau_mean2, tau_cov=jacknifed_tau_covariance2, fixed_coefficients=np.asarray([131]))
	print_organized_summary_file(sldsc_output_root + 'organized_' + str(reg_param) + '_sparse_ard_all_coefficients_res.txt', anno_names, sparse_sldsc_obj.beta_mu/annotation_sdev, np.sqrt(sparse_sldsc_obj.beta_var)/annotation_sdev)
'''


'''
sparse_sldsc_obj = SPARSE_SLDSC_ARD_SOME_FIXED(max_iter=40000, L=10, nonneg=False, nonneg_int=eqtl_start_index)
sparse_sldsc_obj.fit(tau=jacknifed_tau_mean, tau_cov=jacknifed_tau_covariance, fixed_coefficients=non_eqtl_annotations)
model_beta_mu = np.hstack((sparse_sldsc_obj.fixed_beta_mu, sparse_sldsc_obj.beta_mu))
model_beta_var = np.hstack((sparse_sldsc_obj.fixed_beta_var, sparse_sldsc_obj.beta_var))
print(np.sort((model_beta_mu/annotation_sdev)[eqtl_annotations]))
print(anno_names[eqtl_annotations][np.argsort((model_beta_mu/annotation_sdev)[eqtl_annotations])])
pdb.set_trace()
print_organized_summary_file(sldsc_output_root + 'organized_sparse_ard_no_geno_regularization_res.txt', anno_names, model_beta_mu/annotation_sdev, np.sqrt(model_beta_var)/annotation_sdev)
'''
'''
sparse_sldsc_obj = SPARSE_SLDSC_ARD_SOME_FIXED(max_iter=20000, L=10, nonneg=True, nonneg_int=-1.0)
sparse_sldsc_obj.fit(tau=jacknifed_tau_mean, tau_cov=jacknifed_tau_covariance, fixed_coefficients=non_eqtl_annotations)
model_beta_mu = np.hstack((sparse_sldsc_obj.fixed_beta_mu, sparse_sldsc_obj.beta_mu))
model_beta_var = np.hstack((sparse_sldsc_obj.fixed_beta_var, sparse_sldsc_obj.beta_var))
print_organized_summary_file(sldsc_output_root + 'organized_sparse_nonneg_ard_no_geno_regularization_res.txt', anno_names, model_beta_mu/annotation_sdev, np.sqrt(model_beta_var)/annotation_sdev)



sparse_sldsc_obj = SPARSE_SLDSC_ARD_SOME_FIXED(max_iter=20000, L=10, nonneg=True, nonneg_int=eqtl_start_index)
sparse_sldsc_obj.fit(tau=jacknifed_tau_mean, tau_cov=jacknifed_tau_covariance, fixed_coefficients=np.asarray([]))
print_organized_summary_file(sldsc_output_root + 'organized_sparse_nonneg_ard_res.txt', anno_names, sparse_sldsc_obj.beta_mu/annotation_sdev, np.sqrt(sparse_sldsc_obj.beta_var)/annotation_sdev)

sparse_sldsc_obj = SPARSE_SLDSC_ARD_SOME_FIXED(max_iter=20000, L=10, nonneg=False)
sparse_sldsc_obj.fit(tau=jacknifed_tau_mean, tau_cov=jacknifed_tau_covariance, fixed_coefficients=np.asarray([]))
print_organized_summary_file(sldsc_output_root + 'organized_sparse_ard_res.txt', anno_names, sparse_sldsc_obj.beta_mu/annotation_sdev, np.sqrt(sparse_sldsc_obj.beta_var)/annotation_sdev)
'''
