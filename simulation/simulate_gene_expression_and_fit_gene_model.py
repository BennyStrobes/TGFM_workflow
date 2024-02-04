import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import pickle
from pandas_plink import read_plink1_bin
import rpy2
import rpy2.robjects.numpy2ri as numpy2ri
import rpy2.robjects as ro
import statsmodels.api as sm
ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
susieR_pkg = importr('susieR')

def generate_tissue_covariance_structure_across_causal_effects(ge_h2):
	ge_per_snp_h2 = ge_h2/5
	cov_mat = np.zeros((10,10)) + 0.737*ge_per_snp_h2
	cov_mat[:3,:3] = cov_mat[:3,:3]*0.0 + 0.789*ge_per_snp_h2
	cov_mat[3:6,3:6] = cov_mat[3:6,3:6]*0.0 + 0.789*ge_per_snp_h2
	cov_mat[6:,6:] = cov_mat[6:,6:]*0.0 + 0.789*ge_per_snp_h2
	np.fill_diagonal(cov_mat, ge_per_snp_h2)
	return cov_mat


def simulate_causal_eqtl_effect_sizes_across_tissues(n_cis_snps, ge_h2, n_tiss=10):
	ge_per_snp_h2 = ge_h2/5
	# Initialize matrix of causal eqtl effect sizes across tissues
	causal_effect_sizes = np.zeros((n_cis_snps, n_tiss))

	# Generate cross tissue covariance structure matrix
	tissue_covariance_mat = generate_tissue_covariance_structure_across_causal_effects(ge_h2)

	# First simulate 3 variants with shared effects
	shared_variant_indices = np.random.choice(np.arange(n_cis_snps), size=3, replace=False, p=None)
	causal_effect_sizes[shared_variant_indices, :] =  np.random.multivariate_normal(np.zeros(10), tissue_covariance_mat,size=3)
	
	# Now simulate 2 variants with tissue-specific effects
	remaining_variant_indices = np.delete(np.arange(n_cis_snps), shared_variant_indices)
	for tiss_iter in range(n_tiss):
		tissue_specific_indices = np.random.choice(remaining_variant_indices, size=2, replace=False, p=None)
		causal_effect_sizes[tissue_specific_indices, tiss_iter] = np.random.normal(loc=0.0, scale=np.sqrt(ge_per_snp_h2),size=2)
	return causal_effect_sizes

def simulate_causal_eqtl_effect_sizes_across_tissues_shell(n_cis_snps, fraction_genes_cis_h2, ge_h2):
	min_h2 = 0
	while min_h2 < .01:
		causal_eqtl_effects = simulate_causal_eqtl_effect_sizes_across_tissues(n_cis_snps, ge_h2)
		min_h2 = np.min(np.sum(np.square(causal_eqtl_effects),axis=0))
	# Make (1.0 - fraction_genes_cis_h2)% of genes not cis heritable
	gene_cis_h2_boolean = np.random.binomial(n=1, p=fraction_genes_cis_h2, size=causal_eqtl_effects.shape[1])
	for tiss_iter, boolean_value in enumerate(gene_cis_h2_boolean):
		if boolean_value == 0:
			causal_eqtl_effects[:, tiss_iter] = (causal_eqtl_effects[:, tiss_iter])*0.0
	return causal_eqtl_effects

def simulate_causal_eqtl_effect_sizes(cis_window, simulated_gene_position_file, simulated_gene_expression_dir,simulation_name_string, ref_eqtl_sample_size, processed_genotype_data_dir, chrom_num, simulated_causal_eqtl_effect_summary_file, fraction_genes_cis_h2, ge_h2):
	# Load in genotype data across chromosome for single eQTL data set (note: it doesn't matter what eqtl data set we are using because we are just using snp positions here and all eqtl data sets have the same snps)
	genotype_stem = processed_genotype_data_dir + 'simulated_eqtl_' + str(ref_eqtl_sample_size) + '_data_' + chrom_num
	G_obj = read_plink1_bin(genotype_stem + '.bed', genotype_stem + '.bim', genotype_stem + '.fam', verbose=False)
	G_obj_chrom = np.asarray(G_obj.chrom)
	G_obj_pos = np.asarray(G_obj.pos)
	# For our purposes, a0 is the effect allele
	# For case of plink package, a0 is the first column in the plink bim file
	G_obj_a0 = np.asarray(G_obj.a0)
	G_obj_a1 = np.asarray(G_obj.a1)
	# RSids
	G_obj_rsids = np.asarray(G_obj.snp)
	# Snp ids
	G_obj_snp_ids = 'chr' + G_obj_chrom + '_' + (G_obj_pos.astype(str)) + '_' + G_obj_a0 + '_' + G_obj_a1

	t = open(simulated_causal_eqtl_effect_summary_file,'w')
	t.write('gene_id\tchr\ttss\tcausal_eqtl_effect_file\tcis_snp_id_file\tcis_snp_index_file\ttotal_n_snps\n')

	# Loop through genes (run analysis for each gene seperately)
	f = open(simulated_gene_position_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Here we define a single gene
		ensamble_id = data[1]
		tss = int(data[2])
		# Get snp names in cis window fo this gene
		cis_window_start = tss - cis_window
		cis_window_end = tss + cis_window
		cis_snp_indices = (G_obj_pos >= cis_window_start) & (G_obj_pos < cis_window_end)
		# Ignore genes with no or very few cis snps
		if np.sum(cis_snp_indices) < 10:
			continue
		cis_rsids = G_obj_rsids[cis_snp_indices]
		cis_snp_ids = G_obj_snp_ids[cis_snp_indices]
		n_cis_snps = len(cis_rsids)

		# Simulate causal eqtl effects across tissues
		causal_eqtl_effects = simulate_causal_eqtl_effect_sizes_across_tissues_shell(n_cis_snps, fraction_genes_cis_h2, ge_h2)

		# Save results to output
		# Causal eqtl effects
		gene_causal_effect_file = simulated_gene_expression_dir + simulation_name_string + '_' + ensamble_id + '_causal_eqtl_effects.npy'
		np.save(gene_causal_effect_file, causal_eqtl_effects)
		# SNP ids
		gene_cis_snpid_file = simulated_gene_expression_dir + simulation_name_string + '_' + ensamble_id + '_cis_snpids.npy'
		np.save(gene_cis_snpid_file, np.hstack((cis_snp_ids.reshape(n_cis_snps,1), cis_rsids.reshape(n_cis_snps,1))))
		# SNP indices
		gene_cis_snp_indices_file = simulated_gene_expression_dir + simulation_name_string + '_' + ensamble_id + '_cis_snp_indices.npy'
		np.save(gene_cis_snp_indices_file, np.where(cis_snp_indices==True)[0])

		# Write to output file
		t.write(ensamble_id + '\t' + chrom_num + '\t' + str(tss) + '\t' + gene_causal_effect_file + '\t' + gene_cis_snpid_file + '\t' + gene_cis_snp_indices_file + '\t' + str(len(cis_snp_indices)) + '\n')

	f.close()
	t.close()

def mean_impute_and_standardize_genotype(G_obj_geno):
	# Fill in missing values
	G_obj_geno_stand = np.copy(G_obj_geno)
	ncol = G_obj_geno_stand.shape[1]
	n_missing = []
	for col_iter in range(ncol):
		nan_indices = np.isnan(G_obj_geno[:,col_iter])
		non_nan_mean = np.mean(G_obj_geno[nan_indices==False, col_iter])
		G_obj_geno_stand[nan_indices, col_iter] = non_nan_mean
		n_missing.append(np.sum(nan_indices))
	n_missing = np.asarray(n_missing)

	# Quick error check
	if np.sum(np.std(G_obj_geno_stand,axis=0) == 0) > 0:
		print('no variance genotype assumption error')
		pdb.set_trace()

	# Standardize genotype
	G_obj_geno_stand = (G_obj_geno_stand -np.mean(G_obj_geno_stand,axis=0))/np.std(G_obj_geno_stand,axis=0)

	return G_obj_geno_stand


def simulate_gene_expression(G, sim_beta):
	#gene_heritability = np.sum(np.square(sim_beta))
	genetic_pred_ge = np.dot(G, sim_beta)
	gene_heritability = np.var(genetic_pred_ge)
	if gene_heritability > 1:
		gene_heritability = .99
	ge = np.random.normal(loc=genetic_pred_ge, scale=np.sqrt(1.0-gene_heritability))
	return ge


def compute_marginal_regression_coefficients(sim_stand_expr, gene_geno):
	'''
	n_snps = gene_geno.shape[1]
	marginal_effects = np.zeros(n_snps)
	marginal_effects_se = np.zeros(n_snps)
	pvalue_arr = np.zeros(n_snps)
	for snp_iter in range(n_snps):
		# Now get effect size, standard erorr and z-score for association between standardized genotype and trait
		# Fit model using statsmodels
		mod = sm.OLS(sim_stand_expr, gene_geno[:, snp_iter])
		res = mod.fit()
		# Extract results
		effect_size = res.params[0]
		effect_size_se = res.bse[0]
		pvalue = res.pvalues[0]
		effect_size_z = effect_size/effect_size_se

		#marginal_effects.append(effect_size)
		#marginal_effects_se.append(effect_size_se)
		#pvalue_arr.append(pvalue)
		marginal_effects[snp_iter] = effect_size
		marginal_effects_se[snp_iter] = effect_size_se
		pvalue_arr[snp_iter] = pvalue
	'''

	marginal_effects2 = np.dot(np.transpose(gene_geno), sim_stand_expr)/len(sim_stand_expr)
	marginal_effects_se2 = np.sqrt(np.var(np.transpose((gene_geno*marginal_effects2)) - sim_stand_expr, ddof=1,axis=1)/len(sim_stand_expr))


	return marginal_effects2, marginal_effects_se2

def extract_h2_from_greml_hsq_file(hsq_file):
	f = open(hsq_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if line.startswith('V(G)/Vp'):
			hsq = float(data[1])
			hsq_se = float(data[2])
		if line.startswith('Pval'):
			pval = float(data[1])
	f.close()
	return hsq, hsq_se, pval

def run_greml_no_covariate_h2_analysis_and_FUSION(chrom_num, gene_tss, expr_vec, genotype_stem, cis_window, expr_ind_ids, tmp_output_stem, cis_rsids, genotype_mat):
	gcta_path="/n/groups/price/tiffany/subpheno/fusion_twas-master/gcta_nr_robust"
	# Create gene pheno file
	gene_pheno_file = tmp_output_stem + 'gene_pheno'
	n_samples = len(expr_vec)
	#gene_pheno_data = np.hstack((np.zeros((n_samples,1)).astype(int).astype(str), expr_ind_ids.reshape(n_samples,1), expr_vec.astype(str).reshape(n_samples,1)))
	gene_pheno_data = np.hstack((expr_ind_ids.reshape(n_samples,1), expr_ind_ids.reshape(n_samples,1), expr_vec.astype(str).reshape(n_samples,1)))
	np.savetxt(gene_pheno_file, gene_pheno_data, fmt="%s", delimiter='\t')

	# Run PLINK to get plink file specifically consisting of cis snps
	start_pos = gene_tss - int(cis_window)
	end_pos = gene_tss + int(cis_window)
	plink_window_stem = tmp_output_stem + 'window_plink'
	command_string = 'plink --bfile ' + genotype_stem + ' --keep-allele-order --pheno ' + gene_pheno_file + ' --threads 1 --make-bed --out ' + plink_window_stem + ' --keep ' + gene_pheno_file + ' --chr ' + chrom_num + ' --from-bp ' + str(start_pos) + ' --to-bp ' + str(end_pos) +' --allow-no-sex'
	os.system(command_string)

	# MAKE GRM WITH PLINK
	command_string = 'plink --allow-no-sex --bfile ' + plink_window_stem + ' --make-grm-bin --threads 1 --out ' + plink_window_stem
	os.system(command_string)

	# estimate heritability with GREML
	greml_h2_res_file = tmp_output_stem + 'h2_res' 
	#arg = paste( gcta_path ," --grm ",temp_tissue_specific_stem," --pheno ",raw.pheno.file," --qcovar ",covariate_file," --out ",temp_tissue_specific_stem," --reml --reml-no-constrain --reml-lrt 1",sep='')
	command_string = gcta_path + ' --grm ' + plink_window_stem + ' --pheno ' + gene_pheno_file + ' --out ' + greml_h2_res_file + ' --reml --reml-no-constrain --reml-lrt 1'
	os.system(command_string)
	# Now extract heritabilities
	hsq, hsq_se, hsq_p = extract_h2_from_greml_hsq_file(greml_h2_res_file + '.hsq')
	causal_effects = np.zeros(len(cis_rsids))
	if hsq > 0 and hsq_p < .05:
		command_string = 'plink --allow-no-sex --bfile ' + plink_window_stem + ' --keep-allele-order --threads 1 --lasso ' + str(hsq) + ' --out ' + plink_window_stem + 'lasso'
		os.system(command_string)
		f = open(plink_window_stem + 'lasso.lasso')
		mapping = {}
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			mapping[data[1]] = -float(data[3])
		f.close()
		for ii,snp_id in enumerate(cis_rsids):
			if snp_id in mapping:
				causal_effects[ii] = mapping[snp_id]

	# Clear temporary files from directory
	os.system('rm ' + tmp_output_stem + '*')
	return hsq, hsq_se, hsq_p, causal_effects


def bootstrap_marginal_regression_coefficients(sim_stand_expr, gene_geno, n_bootstraps=500):
	bootstrapped_marginals = []
	n_samples = len(sim_stand_expr)
	eqtl_sample_indices = np.arange(n_samples)
	for bootstrap_iter in range(n_bootstraps):
		bootstrapped_eqtl_sample_indices = np.random.choice(eqtl_sample_indices, size=n_samples, replace=True, p=None)
		#bs_coef, bs_coef_se = compute_marginal_regression_coefficients(sim_stand_expr[bootstrapped_eqtl_sample_indices], gene_geno[bootstrapped_eqtl_sample_indices,:])
		bs_coef2 = np.dot(sim_stand_expr[bootstrapped_eqtl_sample_indices], gene_geno[bootstrapped_eqtl_sample_indices,:])/n_samples
		bootstrapped_marginals.append(bs_coef2)
	bootstrapped_marginals = np.asarray(bootstrapped_marginals)
	return bootstrapped_marginals


def extract_indices_of_pcs_that_explain_specified_fraction_of_variance(svd_lambda, fraction):
	if np.sum(svd_lambda < 0.0) != 0:
		print('negative lambda assumption error')
		pdb.set_trace()
	pve = svd_lambda/np.sum(svd_lambda)
	cum_sum = []
	cum = 0.0
	for pve_val in pve:
		cum = cum + pve_val
		cum_sum.append(cum)
	cum_sum = np.asarray(cum_sum)

	valid_pcs = cum_sum <= fraction
	if np.sum(valid_pcs) < 2:
		print('assumption eororor')
		pdb.set_trace()

	return valid_pcs

def reduce_dimensionality_of_genotype_data_with_pca(window_variant_genotype):
	# Reduce dimensionality of window variant genotype with svd
	n_ref = window_variant_genotype.shape[0]
	uu, ss, vh = np.linalg.svd(window_variant_genotype, full_matrices=False)
	svd_lambda = ss*ss/(n_ref-1)
	svd_Q = np.transpose(vh)
	# Prune to PCs that explain 99% of variance
	pc_indices = extract_indices_of_pcs_that_explain_specified_fraction_of_variance(svd_lambda, .99)
	pruned_svd_lambda = svd_lambda[pc_indices]
	pruned_svd_Q = svd_Q[:, pc_indices]

	n_pcs = len(pruned_svd_lambda)
	max_num_pcs = int(np.floor(n_ref*.75))
	if n_pcs >= max_num_pcs:
		pc_filter = np.asarray([False]*n_pcs)
		pc_filter[:max_num_pcs] = True
		pruned_svd_lambda = pruned_svd_lambda[pc_filter]
		pruned_svd_Q = pruned_svd_Q[:, pc_filter]
		n_pcs = len(pruned_svd_lambda)

	return pruned_svd_lambda, pruned_svd_Q


def unbiased_r_squared(marginal_effects, bootstrapped_marginal_effects, pruned_svd_ld_inv_mat, pruned_svd_lambda, pruned_svd_Q, eqtl_sample_size):
	n_pcs = len(pruned_svd_lambda)
	n_bootstraps = bootstrapped_marginal_effects.shape[0]

	alpha_mat = np.dot(pruned_svd_ld_inv_mat, np.transpose(bootstrapped_marginal_effects))

	bootstrapped_gene_variances = np.sum(alpha_mat*np.transpose(bootstrapped_marginal_effects),axis=0)

	bootstrapped_standardized_r_squared = np.transpose(np.square(bootstrapped_marginal_effects))/bootstrapped_gene_variances
	expected_bootstrapped_standardized_r_squared = np.mean(bootstrapped_standardized_r_squared,axis=1)

	#expected_pmces_standardized_r_squared = np.square(np.mean(bootstrapped_marginal_effects,axis=0))/np.dot(np.mean(alpha_mat,axis=1), np.mean(bootstrapped_marginal_effects,axis=0))
	expected_pmces_standardized_r_squared = np.square(marginal_effects)/np.dot(np.dot(pruned_svd_ld_inv_mat, np.transpose(marginal_effects)), marginal_effects)

	noise_terms = expected_bootstrapped_standardized_r_squared - expected_pmces_standardized_r_squared

	unbiased_standardized_r_squared = expected_pmces_standardized_r_squared - noise_terms

	return unbiased_standardized_r_squared


def simulate_gene_expression_and_fit_gene_model_shell(simulated_causal_eqtl_effect_summary_file, eqtl_sample_size, simulation_name_string, processed_genotype_data_dir, simulated_learned_gene_models_dir, chrom_num):
	# Load in genotype data across chromosome for eQTL data set
	genotype_stem = processed_genotype_data_dir + 'simulated_eqtl_' + str(eqtl_sample_size) + '_data_' + chrom_num
	G_obj = read_plink1_bin(genotype_stem + '.bed', genotype_stem + '.bim', genotype_stem + '.fam', verbose=False)
	G_obj_geno = G_obj.values # Numpy 2d array of dimension num samples X num snps
	G_obj_chrom = np.asarray(G_obj.chrom)
	G_obj_pos = np.asarray(G_obj.pos)
	# For our purposes, a0 is the effect allele
	# For case of plink package, a0 is the first column in the plink bim file
	G_obj_a0 = np.asarray(G_obj.a0)
	G_obj_a1 = np.asarray(G_obj.a1)
	# RSids
	G_obj_rsids = np.asarray(G_obj.snp)
	G_obj_sample_names = np.asarray(G_obj.sample)
	# Snp ids
	G_obj_snp_ids = 'chr' + G_obj_chrom + '_' + (G_obj_pos.astype(str)) + '_' + G_obj_a0 + '_' + G_obj_a1

	# Mean impute and standardize genotype
	G_obj_geno_stand = mean_impute_and_standardize_genotype(G_obj_geno)
	# Now loop through genes
	f = open(simulated_causal_eqtl_effect_summary_file)
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		counter = counter + 1
		print(counter)
		# AT A SINGLE GENE
		# Extract relevent information from the line
		ensamble_id = data[0]
		gene_tss = int(data[2])
		gene_causal_eqtl_effect_file = data[3]
		gene_cis_snp_indices_file = data[5]
		total_n_genome_snps = int(data[6])

		# Extract simulated causal effect sizes for this gene (cis_snpsXn_tissues)
		sim_causal_eqtl_effect_sizes = np.load(gene_causal_eqtl_effect_file)

		# Extract indices of cis_snps
		cis_snp_indices_raw = np.load(gene_cis_snp_indices_file)
		cis_snp_indices = np.asarray([False]*total_n_genome_snps)
		cis_snp_indices[cis_snp_indices_raw] = True


		n_cis_snps = np.sum(cis_snp_indices)
		cis_rsids = G_obj_rsids[cis_snp_indices]

		# Quick error check
		if np.sum(cis_snp_indices) < 10:
			print('assumption eroror')
			pdb.set_trace()

		# Extract standardized matrix of cis snps around the gene
		gene_geno = G_obj_geno_stand[:, cis_snp_indices]

		'''
		# Reduce dimensionality of genotype data
		pruned_svd_lambda, pruned_svd_Q = reduce_dimensionality_of_genotype_data_with_pca(gene_geno)
		pruned_svd_ld_inv_mat = np.dot(np.dot(pruned_svd_Q, np.diag(1.0/pruned_svd_lambda)), np.transpose(pruned_svd_Q))
		'''
		# Now loop through tissues (run simulation and gene-modeling seperately in each tissue)
		n_tiss = sim_causal_eqtl_effect_sizes.shape[1]
		pmces_cross_tissues = []
		marginal_beta_cross_tissues = []
		marginal_beta_var_cross_tissues = []
		#marginal_effect_sizes_cross_tissues = []
		#marginal_effect_se_cross_tissues = []
		#fusion_pmces_cross_tissues = []
		#greml_h2_cross_tissues = []
		#greml_h2_se_cross_tissues = []
		#unbiased_r_squared_cross_tissues = []
		for tissue_iter in range(n_tiss):
			# Simulate gene expression in this gene, tissue pair
			sim_expr = simulate_gene_expression(gene_geno, sim_causal_eqtl_effect_sizes[:, tissue_iter])
			# Standardize simulated gene expression
			sim_stand_expr = (sim_expr - np.mean(sim_expr))/np.std(sim_expr)
			# Run eQTL variant fine-mapping with SuSiE
			susie_fitted = susieR_pkg.susie(gene_geno, sim_stand_expr,L=10)
			# Test whether there exist any identified susie components for this gene
			if type(susie_fitted.rx2('sets').rx2('cs_index')) == rpy2.rinterface_lib.sexp.NULLType:
				# no susie components identified
				pmces = np.zeros(n_cis_snps)
			else:
				susie_components = np.asarray(susie_fitted.rx2('sets').rx2('cs_index')) - 1
				pmces = np.sum(susie_fitted.rx2('alpha')*susie_fitted.rx2('mu'),axis=0)
				# Also need to save individual data objects
				# alpha
				alpha_model_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_tissue_' + str(tissue_iter) + '_gene_model_susie_alpha.npy'
				np.save(alpha_model_output_file, susie_fitted.rx2('alpha'))
				# mu
				mu_model_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_tissue_' + str(tissue_iter) + '_gene_model_susie_mu.npy'
				np.save(mu_model_output_file, susie_fitted.rx2('mu'))
				# mu_var
				mu_var = susie_fitted.rx2('mu2') - np.square(susie_fitted.rx2('mu'))
				mu_var_model_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_tissue_' + str(tissue_iter) + '_gene_model_susie_mu_var.npy'
				np.save(mu_var_model_output_file, mu_var)
			'''
			# Run greml and FUSION
			tmp_output_stem = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size)  + '_tmp_h2_'
			hsq, hsq_se, hsq_p, fusion_causal_effects = run_greml_no_covariate_h2_analysis_and_FUSION(str(chrom_num), gene_tss, sim_stand_expr, genotype_stem, 100000, G_obj_sample_names, tmp_output_stem, cis_rsids, gene_geno)
			try:
				hsq, hsq_se, hsq_p, fusion_causal_effects = run_greml_no_covariate_h2_analysis_and_FUSION(str(chrom_num), gene_tss, sim_stand_expr, genotype_stem, 100000, G_obj_sample_names, tmp_output_stem, cis_rsids, gene_geno)
			except:
				fusion_causal_effects = np.zeros(len(cis_rsids))
				hsq = 'nan'
				hsq_se = 'nan'
				hsq_p = 'nan'
			greml_h2_cross_tissues.append(str(hsq))
			greml_h2_se_cross_tissues.append(str(hsq_se))
			fusion_pmces_cross_tissues.append(fusion_causal_effects)
			'''
			marginal_effects, marginal_effects_se = compute_marginal_regression_coefficients(sim_stand_expr, gene_geno)
			'''
			boostrapped_marginal_effects = bootstrap_marginal_regression_coefficients(sim_stand_expr, gene_geno, n_bootstraps=1000)
			marginal_bootstrapped_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_tissue_' + str(tissue_iter) + '_marginal_effects_bootstrapped.npy'
			np.save(marginal_bootstrapped_output_file, boostrapped_marginal_effects)

			# Extract unbiased r_squared
			unbiased_r_squared_res = unbiased_r_squared(marginal_effects, boostrapped_marginal_effects, pruned_svd_ld_inv_mat, pruned_svd_lambda, pruned_svd_Q, eqtl_sample_size)
			unbiased_r_squared_cross_tissues.append(unbiased_r_squared_res)

			marginal_effect_sizes_cross_tissues.append(marginal_effects)
			marginal_effect_se_cross_tissues.append(marginal_effects_se)
			'''
			pmces_cross_tissues.append(pmces)
			marginal_beta_cross_tissues.append(marginal_effects)
			marginal_beta_var_cross_tissues.append(np.square(marginal_effects_se))

		# Convert pmces to numpy array
		pmces_cross_tissues = np.asarray(pmces_cross_tissues)
		marginal_beta_cross_tissues = np.asarray(marginal_beta_cross_tissues)
		marginal_beta_var_cross_tissues = np.asarray(marginal_beta_var_cross_tissues)
		#fusion_pmces_cross_tissues = np.asarray(fusion_pmces_cross_tissues)
		#greml_h2_cross_tissues = np.asarray(greml_h2_cross_tissues)
		#greml_h2_se_cross_tissues = np.asarray(greml_h2_se_cross_tissues)


		# Save gene-model PMCES across tissues to output
		gene_model_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_gene_model_pmces.npy'
		np.save(gene_model_output_file, pmces_cross_tissues)

		# Save marginal pvalue across tissues to output
		marginal_beta_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_marginal_beta.npy'
		np.save(marginal_beta_output_file, marginal_beta_cross_tissues)
		marginal_beta_var_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_marginal_beta_var.npy'
		np.save(marginal_beta_var_output_file, marginal_beta_var_cross_tissues)

		# Save gene-model marginal effects across tissues to output
		#gene_model_marg_effects_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_marginal_effects_gene_model.npy'
		#np.save(gene_model_marg_effects_output_file, marginal_effect_sizes_cross_tissues)

		# Save gene-model marginal effects se across tissues to output
		#gene_model_marg_effects_se_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_marginal_effects_se_gene_model.npy'
		#np.save(gene_model_marg_effects_se_output_file, marginal_effect_se_cross_tissues)

		# Save gene-model marginal effects se across tissues to output
		#unbiased_r_sq_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_unbiased_r_squared.npy'
		#np.save(unbiased_r_sq_output_file, unbiased_r_squared_cross_tissues)

		# Save gene-model PMCES across tissues to output
		#gene_model_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_fusion_lasso_pmces_gene_model.npy'
		#np.save(gene_model_output_file, fusion_pmces_cross_tissues)

		# Save greml estimates to output
		#greml_estimate_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_greml_estimates.npy'
		#np.save(greml_estimate_output_file, np.vstack((greml_h2_cross_tissues,greml_h2_se_cross_tissues)))

	f.close()
	return


############################
# COMMAND LINE ARGUMENTS
############################
simulation_number = int(sys.argv[1])
chrom_num = sys.argv[2]
cis_window = int(sys.argv[3])
simulated_gene_position_file = sys.argv[4]
simulated_gene_expression_dir = sys.argv[5]
simulated_learned_gene_models_dir = sys.argv[6]
simulation_name_string = sys.argv[7]
processed_genotype_data_dir = sys.argv[8]
ge_h2_str = sys.argv[9]

# Get true gene expression h2
ge_h2 = float('.' + ge_h2_str)

# Set seed
np.random.seed(simulation_number)

# Define vector of eQTL sample sizes
eqtl_sample_sizes = np.asarray([300, 500, 1000])


# Fraction of genes cis heritable for a given tissue
fraction_genes_cis_h2 = .5


############################
# Simulate causal eQTL effect sizes across tissues
############################
# Create file to keep track of causal eqtl effect sizes across genes
simulated_causal_eqtl_effect_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'
simulate_causal_eqtl_effect_sizes(cis_window, simulated_gene_position_file, simulated_gene_expression_dir, simulation_name_string, eqtl_sample_sizes[0], processed_genotype_data_dir, chrom_num, simulated_causal_eqtl_effect_summary_file, fraction_genes_cis_h2, ge_h2)

############################
# Simulate Gene expression and fit gene models for each data-set (eqtl sample-size), tissue
############################
# Do seperately in each data set (eqtl sample size) because genotypes are different
for eqtl_sample_size in eqtl_sample_sizes:
	print(eqtl_sample_size)
	simulate_gene_expression_and_fit_gene_model_shell(simulated_causal_eqtl_effect_summary_file, eqtl_sample_size, simulation_name_string, processed_genotype_data_dir, simulated_learned_gene_models_dir, chrom_num)





