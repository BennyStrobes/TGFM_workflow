import sys
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import time


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





class TGFM_CAUSAL_TWAS(object):
	def __init__(self, L=10, max_iter=200, convergence_thresh=1e-5, estimate_prior_variance=False, estimate_pleiotropic_prior_variance=False,fusion_weights=False,robust=True, pleiotropic_prior_variance=1000, ard_a_prior=1e-16, ard_b_prior=1e-16, init_pi=None):
		# Prior on gamma distributions defining residual variance and
		self.L = L
		self.max_iter = max_iter
		self.convergence_thresh = convergence_thresh
		self.estimate_prior_variance = estimate_prior_variance
		self.estimate_pleiotropic_prior_variance = estimate_pleiotropic_prior_variance
		self.pleiotropic_prior_variance = pleiotropic_prior_variance
		self.ard_a_prior = ard_a_prior
		self.ard_b_prior = ard_b_prior
		self.fusion_weights = fusion_weights
		self.robust = robust
		self.init_pi = init_pi
	def fit(self, twas_data_obj):
		""" Fit the model.
			Args:
			twas_data_obj
		"""

		print('###############################')
		print('TGFM-SUSIE CAUSAL TWAS')
		print('###############################')

		#####################
		# Initialize variables
		self.initialize_variables(twas_data_obj)
		self.nominal_twas_rss_updates(twas_data_obj)

		print('Initialization complete')
		self.iter = 0
		# Compute residual
		self.residual = twas_data_obj['gwas_beta'] - np.dot(self.srs_inv, self.global_pmces)
		for temper in range(1):
			if self.robust:
				# Update betea (ie predicted pleiotropic gwas effect sizes)
				self.update_beta(self.expected_gamma_beta)
		# BEGIN CAVI
		for itera in range(self.max_iter):
			# Update alpha (ie predicted effect of each gene on the trait)
			self.update_susie_alpha(np.log(self.pi), self.component_variances)
			if self.robust:
				# Update betea (ie predicted pleiotropic gwas effect sizes)
				self.update_beta(self.expected_gamma_beta)

			if self.estimate_prior_variance:
				self.update_component_variances()
			if self.estimate_pleiotropic_prior_variance:
				if self.robust:
					self.update_gamma_beta()
			diff = calculate_distance_between_two_vectors(self.alpha_mu, self.prev_alpha_mu)
			self.convergence_tracker.append(diff)
			self.prev_alpha_mu = np.copy(self.alpha_mu)
			#diff = calculate_distance_between_two_vectors(self.beta_mu, self.prev_beta_mu)
			#self.prev_beta_mu = np.copy(self.beta_mu)
			#print(diff)
			self.iter = self.iter + 1
			if diff <= self.convergence_thresh:
				self.converged = True
				self.update_susie_alpha_no_pi(self.component_variances)
				break

		self.update_susie_alpha_no_pi(self.component_variances)
		if self.converged == False:
			print('Did not converge after ' + str(self.max_iter) + ' iterations')

	def update_susie_alpha(self, expected_log_pi, component_variances):
		if False:
			self.update_susie_alpha_fusion_weights(expected_log_pi, component_variances)
		else:
			self.update_susie_alpha_susie_weights(expected_log_pi, component_variances)
	def update_susie_alpha_no_pi(self, component_variances):
		if False:
			self.update_susie_alpha_fusion_weights_no_pi(expected_log_pi, component_variances)
		else:
			self.update_susie_alpha_susie_weights_no_pi(component_variances)

	def update_susie_alpha_susie_weights_no_pi(self, component_variances):
		# Include effects of all components in residual
		# residual should only contain beta-effects
		gene_trait_pred = np.dot(np.sum(self.alpha_mu*self.phi,axis=0), self.gene_eqtl_pmces)
		self.residual = self.residual + np.dot(self.srs_inv, gene_trait_pred)

		self.phi_no_weights = np.copy(self.phi)
		# Loop through components
		for l_index in range(self.L):
			# Calculate terms inolved in update
			gene_weights = np.sum(self.alpha_mu*self.phi,axis=0) - (self.alpha_mu[l_index,:]*self.phi[l_index,:])
			b_terms = np.dot(self.gene_eqtl_pmces, np.multiply(self.residual, self.s_inv_2_diag)) + np.dot(self.precomputed_gene_gene_terms, gene_weights)
			a_terms = self.precomputed_a_terms - .5*(1.0/component_variances[l_index])

			mixture_alpha_var = -1.0/(2.0*a_terms)
			mixture_alpha_mu = b_terms*mixture_alpha_var
			
			un_normalized_lv_weights = - (.5*np.log(component_variances[l_index])) + (.5*np.square(mixture_alpha_mu)/mixture_alpha_var) + (.5*np.log(mixture_alpha_var))
			un_normalized_lv_weights = un_normalized_lv_weights - scipy.special.logsumexp(un_normalized_lv_weights)
			self.phi_no_weights[l_index,:] = np.exp(un_normalized_lv_weights)

		gene_trait_pred = np.dot(np.sum(self.alpha_mu*self.phi,axis=0), self.gene_eqtl_pmces)
		self.residual = self.residual - np.dot(self.srs_inv, gene_trait_pred)


	def update_susie_alpha_susie_weights(self, expected_log_pi, component_variances):
		# Include effects of all components in residual
		# residual should only contain beta-effects
		gene_trait_pred = np.dot(np.sum(self.alpha_mu*self.phi,axis=0), self.gene_eqtl_pmces)
		self.residual = self.residual + np.dot(self.srs_inv, gene_trait_pred)
		# Loop through components
		for l_index in range(self.L):
			# Calculate terms inolved in update
			gene_weights = np.sum(self.alpha_mu*self.phi,axis=0) - (self.alpha_mu[l_index,:]*self.phi[l_index,:])
			b_terms = np.dot(self.gene_eqtl_pmces, np.multiply(self.residual, self.s_inv_2_diag)) + np.dot(self.precomputed_gene_gene_terms, gene_weights)
			a_terms = self.precomputed_a_terms - .5*(1.0/component_variances[l_index])

			mixture_alpha_var = -1.0/(2.0*a_terms)
			mixture_alpha_mu = b_terms*mixture_alpha_var
			
			un_normalized_lv_weights = expected_log_pi - (.5*np.log(component_variances[l_index])) + (.5*np.square(mixture_alpha_mu)/mixture_alpha_var) + (.5*np.log(mixture_alpha_var))
			un_normalized_lv_weights = un_normalized_lv_weights - scipy.special.logsumexp(un_normalized_lv_weights)
			self.phi[l_index,:] = np.exp(un_normalized_lv_weights)
			self.alpha_mu[l_index,:] = mixture_alpha_mu
			self.alpha_var[l_index,:] = mixture_alpha_var

		gene_trait_pred = np.dot(np.sum(self.alpha_mu*self.phi,axis=0), self.gene_eqtl_pmces)
		self.residual = self.residual - np.dot(self.srs_inv, gene_trait_pred)


	def update_susie_alpha_fusion_weights(self, expected_log_pi, component_variances):
		# Include deffect of the first component (as it is currently removed)
		component_gene_trait_pred = np.dot(self.alpha_mu[0,:]*self.phi[0,:], self.gene_eqtl_pmces)
		self.residual = self.residual + np.dot(self.srs_inv, component_gene_trait_pred)

		# Loop through components
		for l_index in range(self.L):
			# Calculate terms inolved in update
			b_terms = np.dot(self.gene_eqtl_pmces, np.multiply(self.residual, self.s_inv_2_diag))
			a_terms = self.precomputed_a_terms - .5*(1.0/component_variances[l_index])

			mixture_alpha_var = -1.0/(2.0*a_terms)
			mixture_alpha_mu = b_terms*mixture_alpha_var
			
			un_normalized_lv_weights = expected_log_pi - (.5*np.log(component_variances[l_index])) + (.5*np.square(mixture_alpha_mu)/mixture_alpha_var) + (.5*np.log(mixture_alpha_var))
			un_normalized_lv_weights = un_normalized_lv_weights - scipy.special.logsumexp(un_normalized_lv_weights)
			self.phi[l_index,:] = np.exp(un_normalized_lv_weights)
			self.alpha_mu[l_index,:] = mixture_alpha_mu
			self.alpha_var[l_index,:] = mixture_alpha_var

			# Update resid for next round (after this resid includes effects of all components)
			if l_index == (self.L - 1):
				gene_trait_pred = -np.dot(self.alpha_mu[l_index,:]*self.phi[l_index,:], self.gene_eqtl_pmces)
			else:
				gene_trait_pred = np.dot(self.alpha_mu[(l_index+1),:]*self.phi[(l_index+1),:], self.gene_eqtl_pmces) - np.dot(self.alpha_mu[l_index,:]*self.phi[l_index,:], self.gene_eqtl_pmces)
			self.residual = self.residual + np.dot(self.srs_inv, gene_trait_pred)


	def update_alpha(self, expected_gamma_alpha):
		# Include effect of the gene corresponding to 0 from the residaul
		gene_trait_pred = self.gene_eqtl_pmces[0]*self.alpha_mu[0]
		self.residual = self.residual + np.dot(self.srs_inv, gene_trait_pred)

		for g_index in range(self.G):
			# Calculate terms involved in update	
			b_term = np.dot(np.multiply(self.residual, self.s_inv_2_diag), self.gene_eqtl_pmces[g_index])
			a_term = self.precomputed_a_terms[g_index] - .5*expected_gamma_alpha[g_index]

			# VI Updates
			self.alpha_var[g_index] = -1.0/(2.0*a_term)
			self.alpha_mu[g_index] = b_term*self.alpha_var[g_index]

			# Update resid for next round (after this resid includes effects of all genes)
			if g_index == (self.G - 1):
				gene_trait_pred = -self.gene_eqtl_pmces[g_index]*self.alpha_mu[g_index]
			else:
				gene_trait_pred = self.gene_eqtl_pmces[(g_index+1)]*self.alpha_mu[(g_index+1)] - self.gene_eqtl_pmces[g_index]*self.alpha_mu[g_index]
			self.residual = self.residual + np.dot(self.srs_inv, gene_trait_pred)


	def update_beta(self, expected_gamma_beta):
		precomputed_a_terms = (-.5*self.D_diag) - (.5*expected_gamma_beta)

		# Remove the pleiotropic effect of the variant corresponding to 0 from the residaul
		self.residual = self.residual + self.srs_inv[:,0]*self.beta_mu[0]

		# Update each snp in parallel
		for k_index in range(self.K):

			# Remove the pleiotropic effect of the variant corresponding to k_index from the residaul
			#self.residual = self.residual + self.srs_inv[:,k_index]*self.beta_mu[k_index]

			# Calculate terms involved in update
			b_term = self.residual[k_index]*self.s_inv_2_diag[k_index]
			#a_term = (-.5*self.D_diag[k_index]) - (.5*expected_gamma_beta)

			# VI Updates
			self.beta_var[k_index] = -1.0/(2.0*precomputed_a_terms[k_index])
			self.beta_mu[k_index] = b_term*self.beta_var[k_index]

			# Update resid for next round (after this resid includes effects of all genes)
			if k_index == (self.K - 1):
				self.residual = self.residual - self.srs_inv[:,k_index]*self.beta_mu[k_index]
			else:
				self.residual = self.residual + self.srs_inv[:,(k_index+1)]*self.beta_mu[(k_index+1)]- self.srs_inv[:,k_index]*self.beta_mu[k_index]

	def update_component_variances(self):
		#for l_index in range(self.L):
			#self.component_variances[l_index] = np.sum((np.square(self.alpha_mu[l_index,:]) + self.alpha_var[l_index,:])*self.phi[l_index,:])
		self.component_variances = np.sum((np.square(self.alpha_mu) + self.alpha_var)*self.phi,axis=1)

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


	def initialize_variables(self, twas_data_obj):
		# Number of genes
		self.G = len(twas_data_obj['genes'])
		# Number of variants
		self.K = len(twas_data_obj['variants'])
		# Gene names
		self.genes = twas_data_obj['genes']

		if self.init_pi is None:
			self.pi = np.ones(self.G)/self.G
		else:
			self.pi = np.copy(self.init_pi)

		self.expected_gamma_beta = 1.0/self.pleiotropic_prior_variance

		# Compute nominal twas z-scores
		self.nominal_twas_z = np.zeros(self.G)
		gwas_z = twas_data_obj['gwas_beta']/twas_data_obj['gwas_beta_se']
		if self.fusion_weights == False:
			for g_index in range(self.G):
				weights = np.sum(twas_data_obj['susie_mu'][g_index]*twas_data_obj['susie_alpha'][g_index],axis=0)  # Extract susie PMCES for this gene
				twas_z = np.dot(weights, gwas_z)/np.sqrt(np.dot(np.dot(weights, twas_data_obj['reference_ld']), weights))
				self.nominal_twas_z[g_index] = twas_z
		else:
			for g_index in range(self.G):
				weights = twas_data_obj['fusion_weights'][g_index]  # Extract susie PMCES for this gene
				twas_z = np.dot(weights, gwas_z)/np.sqrt(np.dot(np.dot(weights, twas_data_obj['reference_ld']), weights))
				self.nominal_twas_z[g_index] = twas_z		


		# Initialize variational distributions defining alphas (the causal effect of genetically-predicted expression in each gene on the trait)
		# Currently using null intitialization
		self.alpha_mu = np.zeros((self.L, self.G))
		self.alpha_var = np.ones((self.L, self.G))
		self.phi = np.ones((self.L, self.G))/self.G
		self.component_variances = np.ones(self.L)*1e-4

		self.prev_alpha_mu = np.copy(self.alpha_mu)

		# Initialize variational distribution defining betas (the causal effect of pleiotropic genotype on the trait)
		self.beta_mu = np.zeros(self.K)
		self.beta_var = np.ones(self.K)

		# Generate S matrix
		s_squared_vec = np.square(twas_data_obj['gwas_beta_se']) + (np.square(twas_data_obj['gwas_beta'])/twas_data_obj['gwas_sample_size'])
		s_vec = np.sqrt(s_squared_vec)
		S_mat = np.diag(s_vec)
		S_inv_mat = np.diag(1.0/s_vec)
		S_inv_2_mat = np.diag(1.0/np.square(s_vec))

		# Compute (S^-1)R(S^-1) taking advantage of fact that S^-1 is a diagonal matrix
		D_mat = np.multiply(np.multiply(np.diag(S_inv_mat)[:, None], twas_data_obj['reference_ld']), np.diag(S_inv_mat))
		# Compute (S)R(S^-1) taking advantage of fact that S and S^-1 is a diagonal matrix
		srs_inv_mat = np.multiply(np.multiply(np.diag(S_mat)[:, None], twas_data_obj['reference_ld']), np.diag(S_inv_mat))

		# Generate data object containing statistics that are precomputed
		self.srs_inv = srs_inv_mat
		self.s_inv_2_diag = np.diag(S_inv_2_mat)

		self.D_diag = np.diag(D_mat)

		# Get eqtl PMCES
		self.gene_eqtl_pmces = []
		if self.fusion_weights == False:
			for g_index in range(self.G):
				eqtl_pmces = np.sum((twas_data_obj['susie_mu'][g_index])*(twas_data_obj['susie_alpha'][g_index]),axis=0)
				self.gene_eqtl_pmces.append(eqtl_pmces)
		else:
			for g_index in range(self.G):
				eqtl_pmces = twas_data_obj['fusion_weights'][g_index]
				self.gene_eqtl_pmces.append(eqtl_pmces)	

		self.converged = False

		# CAN PRECOMPUTE a-terms for variational updates (will not change iteration to iteration)
		self.precomputed_a_terms = np.zeros(self.G)

		self.global_pmces = np.zeros(self.K) + self.beta_mu

		self.convergence_tracker = []

		if self.fusion_weights == False:
			self.precomputed_gene_gene_terms = np.zeros((self.G, self.G))
			for g_index in range(self.G):
				num_susie_components = twas_data_obj['susie_mu'][g_index].shape[0]
				for k_index in range(num_susie_components):
					self.precomputed_a_terms[g_index] = self.precomputed_a_terms[g_index] - np.sum(.5*(np.square(twas_data_obj['susie_mu'][g_index][k_index,:]) + np.square(twas_data_obj['susie_mu_sd'][g_index][k_index,:]))*np.diag(D_mat)*twas_data_obj['susie_alpha'][g_index][k_index,:])
					eqtl_component_pmces = (twas_data_obj['susie_mu'][g_index][k_index,:])*(twas_data_obj['susie_alpha'][g_index][k_index,:])
					self.precomputed_a_terms[g_index] = self.precomputed_a_terms[g_index] + .5*np.dot(np.dot(eqtl_component_pmces,D_mat), eqtl_component_pmces)
				self.precomputed_a_terms[g_index] = self.precomputed_a_terms[g_index] - .5*np.dot(np.dot(self.gene_eqtl_pmces[g_index],D_mat), self.gene_eqtl_pmces[g_index])
				self.global_pmces = self.global_pmces + self.gene_eqtl_pmces[g_index]*np.sum(self.alpha_mu[:,g_index]*self.phi[:, g_index])
				
				# Cross terms
				self.precomputed_gene_gene_terms[g_index, g_index] = 2.0*self.precomputed_a_terms[g_index]
				for g_prime_index in range((g_index+1), self.G):
					temp_val = -np.dot(np.dot(self.gene_eqtl_pmces[g_prime_index],D_mat), self.gene_eqtl_pmces[g_index])
					self.precomputed_gene_gene_terms[g_index, g_prime_index] = temp_val
					self.precomputed_gene_gene_terms[g_prime_index, g_index] = temp_val
		else:
			self.precomputed_gene_gene_terms = np.zeros((self.G, self.G))
			for g_index in range(self.G):
				self.precomputed_a_terms[g_index] = - .5*np.dot(np.dot(self.gene_eqtl_pmces[g_index],D_mat), self.gene_eqtl_pmces[g_index])
				self.global_pmces = self.global_pmces + self.gene_eqtl_pmces[g_index]*np.sum(self.alpha_mu[:,g_index]*self.phi[:, g_index])
				# Cross terms
				self.precomputed_gene_gene_terms[g_index, g_index] = 2.0*self.precomputed_a_terms[g_index]
				for g_prime_index in range((g_index+1), self.G):
					temp_val = -np.dot(np.dot(self.gene_eqtl_pmces[g_prime_index],D_mat), self.gene_eqtl_pmces[g_index])
					self.precomputed_gene_gene_terms[g_index, g_prime_index] = temp_val
					self.precomputed_gene_gene_terms[g_prime_index, g_index] = temp_val

		# Converg gene_eqtl_pmces to more matrix friendly format
		self.gene_eqtl_pmces = np.asarray(self.gene_eqtl_pmces)

		# Compute expression correlation matrix
		expression_covariance = np.dot(np.dot(self.gene_eqtl_pmces, twas_data_obj['reference_ld']), np.transpose(self.gene_eqtl_pmces))
		dd = np.diag(1.0/np.sqrt(np.diag(expression_covariance)))
		self.ge_ld = np.dot(np.dot(dd, expression_covariance),dd)
