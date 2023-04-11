import sys
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import time

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

def extract_sparse_genotype_matrix_for_this_gene_from_ld(gene_ld_mat, eqtl_sample_size):
	svd_lambda, svd_Q = np.linalg.eig(gene_ld_mat)
	# Filter out complex eigenvalues
	real_valued_eigs = np.iscomplex(svd_lambda) == False
	svd_lambda = svd_lambda[real_valued_eigs].astype(float)
	if np.sum(np.iscomplex(svd_Q[:, real_valued_eigs])) != 0.0:
		print('asssumption eroror')
		pdb.set_trace()
	svd_Q = svd_Q[:, real_valued_eigs].astype(float)
	# Filter out negative eigen values
	positive_eigs = svd_lambda > 0.0
	svd_lambda = svd_lambda[positive_eigs]
	svd_Q = svd_Q[:, positive_eigs]


	# Prune to PCs that explain 99% of variance
	pc_indices = extract_indices_of_pcs_that_explain_specified_fraction_of_variance(svd_lambda, .99)
	pruned_svd_lambda = svd_lambda[pc_indices]
	pruned_svd_Q = svd_Q[:, pc_indices]

	n_pcs = len(pruned_svd_lambda)
	max_num_pcs = int(np.floor(eqtl_sample_size*.75))
	if n_pcs >= max_num_pcs:
		pc_filter = np.asarray([False]*n_pcs)
		pc_filter[:max_num_pcs] = True
		pruned_svd_lambda = pruned_svd_lambda[pc_filter]
		pruned_svd_Q = pruned_svd_Q[:, pc_filter]
		n_pcs = len(pruned_svd_lambda)

	pruned_svd_ld_inv_mat = np.dot(np.dot(pruned_svd_Q, np.diag(1.0/pruned_svd_lambda)), np.transpose(pruned_svd_Q))

	return pruned_svd_ld_inv_mat


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


def calculate_distance_between_two_vectors(vec1, vec2):
	#dist = np.sqrt(np.sum(np.square(vec1-vec2)))
	dist = np.mean(np.abs(vec1-vec2))
	return dist





class SUSIE_IMPUTE(object):
	def __init__(self, L=10, max_iter=200, convergence_thresh=1e-5, estimate_prior_variance=True, alpha_0=.1, ard_a_prior=0.0, ard_b_prior=0.0, variant_init_log_pi=None):
		# Prior on gamma distributions defining residual variance and
		self.L = L
		self.max_iter = max_iter
		self.convergence_thresh = convergence_thresh
		self.estimate_prior_variance = estimate_prior_variance
		self.ard_a_prior = ard_a_prior
		self.ard_b_prior = ard_b_prior
		self.variant_init_log_pi = variant_init_log_pi
		self.alpha_0 = alpha_0
	def fit(self, beta, beta_se, LD, noise_ratios):
		""" Fit the model.
			Args:
			twas_data_obj
		"""

		print('###############################')
		print('SuSiE IMPUTE')
		print('###############################')

		#####################
		# Initialize variables
		self.initialize_variables(beta, beta_se, LD, noise_ratios)


		print('Initialization complete')
		self.iter = 0
		# Compute residual
		self.residual = np.copy(beta)

		# BEGIN CAVI
		for itera in range(self.max_iter):
			# Update alpha (ie predicted effect of each gene on the trait)
			self.update_susie_effects(self.variant_log_pi, self.component_variances)


			if self.estimate_prior_variance:
				self.update_component_variances()

			diff = calculate_distance_between_two_vectors(self.beta_mu, self.prev_beta_mu)
			self.convergence_tracker.append(diff)
			self.prev_beta_mu = np.copy(self.beta_mu)
			self.iter = self.iter + 1
			print(diff)
			if diff <= self.convergence_thresh:
				self.converged = True
				break

		if self.converged == False:
			print('Did not converge after ' + str(self.max_iter) + ' iterations')


	def update_prior_probabilities(self):
		for l_index in range(self.L):
			alpha_l = self.alpha_0 + self.beta_phi[l_index,:]
			self.variant_log_pi[l_index,:] = scipy.special.digamma(alpha_l) - scipy.special.digamma(np.sum(alpha_l))


	def update_susie_effects(self, expected_log_variant_pi, component_variances):
		# Loop through components
		for l_index in range(self.L):
			# Include current component (as it is currently removed)
			component_non_med_pred = self.beta_mu[l_index,:]*self.beta_phi[l_index,:]
			self.residual = self.residual + np.dot(self.srs_inv, component_non_med_pred)


			#################
			# Beta effects
			#################
			variant_b_terms = self.residual*self.s_inv_2_diag
			variant_b_terms[self.noisy_predictors] = variant_b_terms[self.noisy_predictors] + np.dot(self.residual, self.s_inv_r_inv_s_inv)*self.noise_ratios[self.noisy_predictors]
			variant_a_terms = (self.precomputed_a_terms) - .5*(1.0/component_variances[l_index])

			mixture_beta_var = -1.0/(2.0*variant_a_terms)
			mixture_beta_mu = variant_b_terms*mixture_beta_var
			
			################
			# Normalization
			###############
			un_normalized_lv_beta_weights = expected_log_variant_pi[l_index,:] - (.5*np.log(component_variances[l_index])) + (.5*np.square(mixture_beta_mu)/mixture_beta_var) + (.5*np.log(mixture_beta_var))

			normalizing_term = scipy.special.logsumexp(np.hstack((un_normalized_lv_beta_weights)))

			# Save results to global model parameters
			self.beta_phi[l_index,:] = np.exp(un_normalized_lv_beta_weights-normalizing_term)
			self.beta_mu[l_index,:] = mixture_beta_mu
			self.beta_var[l_index,:] = mixture_beta_var


			# Remove current component (as it is currently removed)
			component_non_med_pred = self.beta_mu[l_index,:]*self.beta_phi[l_index,:]
			self.residual = self.residual - np.dot(self.srs_inv, component_non_med_pred)



	def update_component_variances(self):
		variance_a_term = self.ard_a_prior + (1.0/2.0)
		variance_b_term = self.ard_b_prior + np.sum((np.square(self.beta_mu) + self.beta_var)*self.beta_phi,axis=1)/np.sum(self.beta_phi,axis=1)/2.0
		self.component_variances = variance_b_term/variance_a_term
	def update_gamma_alpha(self):
		alpha_squared_expected_val = np.square(self.alpha_mu) + self.alpha_var
		# VI updates
		self.gamma_alpha_a = self.ard_a_prior + (self.G/2.0)
		self.gamma_alpha_b = self.ard_b_prior + (np.sum(alpha_squared_expected_val)/2.0)

		self.expected_gamma_alpha = np.ones(self.G)*(self.gamma_alpha_a/self.gamma_alpha_b)

	def update_gamma_beta(self):
		beta_squared_expected_val = np.square(self.beta_mu) + self.beta_var
		# VI updates
		self.gamma_beta_a = self.ard_a_prior + (self.K/2.0)
		self.gamma_beta_b = self.ard_b_prior + (np.sum(beta_squared_expected_val)/2.0)

		self.expected_gamma_beta = (self.gamma_beta_a/self.gamma_beta_b)



	def nominal_twas_rss_updates(self, twas_data_obj):
		self.nominal_twas_rss_z = np.zeros(self.G)
		self.nominal_twas_rss_alpha_mu = np.zeros(self.G)
		self.nominal_twas_rss_alpha_var = np.zeros(self.G)
		for g_index in range(self.G):
			# Remove effect of the gene corresponding to g_index from the residaul
			if self.fusion_weights == False:
				eqtl_pmces = np.sum((twas_data_obj['susie_mu'][g_index])*(twas_data_obj['susie_alpha'][g_index]),axis=0)
			else:
				eqtl_pmces = twas_data_obj['fusion_weights'][g_index]
				
			b_term = np.dot(np.multiply(twas_data_obj['gwas_beta'], self.s_inv_2_diag), eqtl_pmces)
			#a_term = self.precomputed_a_terms[g_index]
			a_term = self.precomputed_a_terms[g_index]

			alpha_var = -1.0/(2.0*a_term)
			alpha_mu = b_term*alpha_var
			self.nominal_twas_rss_z[g_index] = alpha_mu/np.sqrt(alpha_var)
			self.nominal_twas_rss_alpha_mu[g_index] = alpha_mu
			self.nominal_twas_rss_alpha_var[g_index] = alpha_var


	def initialize_variables(self, beta, beta_se, LD, noise_ratios):

		# Number of variants
		self.K = len(beta)

		# Set log pi
		if self.variant_init_log_pi is None:
			variant_pi = np.ones(self.K)/(self.K)
			self.variant_log_pi = np.log(variant_pi)
		else:
			self.variant_log_pi = np.copy(self.variant_init_log_pi)
		# Repeat rows
		self.variant_log_pi = np.transpose(np.repeat(np.reshape(self.variant_log_pi,(len(self.variant_log_pi),1)),self.L,axis=1))

		# Noise ratio data
		self.noise_ratios = np.copy(noise_ratios)
		self.noisy_predictors = self.noise_ratios!=0.0


		# Initialize variational distributions defining component variances
		self.component_variances = np.ones(self.L)*1e4

		# Initialize variational distribution defining betas (the causal effect of pleiotropic genotype on the trait)
		self.beta_mu = np.zeros((self.L, self.K))
		self.beta_var = np.ones((self.L, self.K))
		self.beta_phi = np.ones((self.L, self.K))/(self.K)

		self.prev_beta_mu = np.copy(self.beta_mu)

		# Generate S matrix
		s_squared_vec = np.square(beta_se)
		s_vec = np.sqrt(s_squared_vec)
		S_mat = np.diag(s_vec)
		S_inv_mat = np.diag(1.0/s_vec)
		S_inv_2_mat = np.diag(1.0/np.square(s_vec))

		# Compute (S^-1)R(S^-1) taking advantage of fact that S^-1 is a diagonal matrix
		D_mat = np.multiply(np.multiply(np.diag(S_inv_mat)[:, None], LD), np.diag(S_inv_mat))
		# Compute (S)R(S^-1) taking advantage of fact that S and S^-1 is a diagonal matrix
		srs_inv_mat = np.multiply(np.multiply(np.diag(S_mat)[:, None], LD), np.diag(S_inv_mat)) - np.diag(self.noise_ratios)


		# Generate data object containing statistics that are precomputed
		self.srs_inv = srs_inv_mat
		self.s_inv_2_diag = np.diag(S_inv_2_mat)

		self.D_diag = np.diag(D_mat)

		self.converged = False
		self.convergence_tracker = []


		# Invert LD matrix
		#inverse_ld = np.linalg.inv(LD + np.eye(LD.shape[0])*1e-05)
		inverse_ld = extract_sparse_genotype_matrix_for_this_gene_from_ld(LD, 100000.0)
		#inverse_ld2 = np.linalg.inv(LD + np.eye(LD.shape[0])*1e-09)
		#inverse_ld = np.linalg.pinv(LD)
		#pdb.set_trace()

		# Precompute some terms
		self.precomputed_a_terms = (-0.5*self.D_diag) - (self.noise_ratios*self.D_diag) - (0.5*np.square(self.noise_ratios)*self.D_diag*np.diag(inverse_ld))
		self.s_inv_r_inv_s_inv = np.multiply(np.multiply(np.diag(S_inv_mat)[:, None], inverse_ld), np.diag(S_inv_mat))[:, self.noisy_predictors]


