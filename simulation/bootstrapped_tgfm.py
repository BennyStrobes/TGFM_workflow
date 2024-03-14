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


def fill_in_causal_effect_size_matrix(null_mat, bs_eqtls_pmces_sparse):
	null_mat[bs_eqtls_pmces_sparse[:,0].astype(int), bs_eqtls_pmces_sparse[:,1].astype(int)] = bs_eqtls_pmces_sparse[:,2]
	return null_mat




class TGFM(object):
	def __init__(self, L=10, max_iter=200, convergence_thresh=1e-5, estimate_prior_variance=False, fusion_weights=False, ard_a_prior=1e-16, ard_b_prior=1e-16, gene_init_log_pi=None, variant_init_log_pi=None, single_variance_component=True, bootstrap_prior=False):
		# Prior on gamma distributions defining residual variance and
		self.L = L
		self.max_iter = max_iter
		self.convergence_thresh = convergence_thresh
		self.estimate_prior_variance = estimate_prior_variance
		self.ard_a_prior = ard_a_prior
		self.ard_b_prior = ard_b_prior
		self.fusion_weights = fusion_weights
		self.gene_init_log_pi = gene_init_log_pi
		self.variant_init_log_pi = variant_init_log_pi
		self.single_variance_component = single_variance_component
		self.bootstrap_prior = bootstrap_prior
	def fit(self, twas_data_obj, phi_init=None, mu_init=None, mu_var_init=None):
		""" Fit the model.
			Args:
			twas_data_obj
		"""

		print('###############################')
		print('Bootstrapped TGFM')
		print('###############################')

		#####################
		# Initialize variables
		self.initialize_variables(twas_data_obj, phi_init, mu_init, mu_var_init)

		print('Initialization complete')
		self.iter = 0
		# Compute residual
		self.residual = []
		self.residual_incl_alpha = []
		for bs_iter in range(self.n_bs):
			self.residual.append(twas_data_obj['gwas_beta'])
			self.residual_incl_alpha.append(twas_data_obj['gwas_beta'])
		self.residual = np.asarray(self.residual)
		self.residual_incl_alpha = np.asarray(self.residual_incl_alpha)

		# Remove initialized non-mediated effects
		# Currently not implemented to remove initialized mediated efffects (but will return an  errror if we try to)
		for l_index in range(self.L - 1):
			component_non_med_pred = self.beta_mus[l_index]*self.beta_phis[l_index]
			self.residual = self.residual - np.dot(component_non_med_pred, self.srs_inv)
			self.residual_incl_alpha = self.residual_incl_alpha - np.dot(component_non_med_pred, self.srs_inv)


		# BEGIN CAVI
		for itera in range(self.max_iter):
			# Update alpha (ie predicted effect of each gene on the trait)
			self.update_susie_effects(self.gene_log_pi, self.variant_log_pi, self.component_variances, self.component_variances)

			if self.estimate_prior_variance:
				self.update_component_variances()

			self.iter = self.iter + 1

		# DONE WITH CAVI
		# COMPUTE PIPS for each bootstrap
		self.compute_pips()


	def compute_pips(self):
		self.alpha_pips = np.ones((self.n_bs, self.G))
		self.beta_pips = np.ones((self.n_bs, self.K))

		for component_iter in range(self.L):
			self.alpha_pips = self.alpha_pips*(1.0 - self.alpha_phis[component_iter])
			self.beta_pips = self.beta_pips*(1.0 - self.beta_phis[component_iter])

		self.alpha_pips = 1.0 - self.alpha_pips
		self.beta_pips = 1.0 - self.beta_pips

		self.expected_alpha_pips = np.mean(self.alpha_pips,axis=0)
		self.expected_beta_pips = np.mean(self.beta_pips,axis=0)




	def update_susie_effects(self, expected_log_pi, expected_log_variant_pi, alpha_component_variances, beta_component_variances):
		# Initialize array to keep track of causal eqtl effect sizes across bootstraps
		gene_eqtl_pmces = np.zeros((self.G, self.K))
		# Loop through components
		for l_index in range(self.L):
			#print(l_index)

			if l_index == 0:
				prev_l_index = self.L - 1
			else:
				prev_l_index = l_index - 1

			# Include current component (as it is currently removed)
			component_non_med_pred = self.beta_mus[l_index]*self.beta_phis[l_index] - self.beta_mus[prev_l_index]*self.beta_phis[prev_l_index]
			self.residual = self.residual + np.dot(self.gene_trait_pred[l_index] - self.gene_trait_pred[prev_l_index] + component_non_med_pred, self.srs_inv)
			self.residual_incl_alpha = self.residual_incl_alpha + np.dot(component_non_med_pred, self.srs_inv)
			for bs_iter in range(self.n_bs):
				# Load in gene eQTL PMCES for this bootstrap
				gene_eqtl_pmces = gene_eqtl_pmces*0.0
				#bs_eqtls_pmces_sparse = np.load(self.pmces_files[bs_iter])
				bs_eqtls_pmces_sparse = self.sparse_sampled_gene_eqtl_pmces[bs_iter]
				gene_eqtl_pmces = fill_in_causal_effect_size_matrix(gene_eqtl_pmces, bs_eqtls_pmces_sparse)
				#gene_eqtl_pmces = np.load(self.pmces_files[bs_iter])

				#################
				# Alpha effects
				#################
				# Calculate terms inolved in update
				#gene_weights = np.sum(self.alpha_mu*self.alpha_phi,axis=0) - (self.alpha_mu[l_index,:]*self.alpha_phi[l_index,:])
				gene_weights = self.agg_gene_trait[bs_iter] - (self.alpha_mus[l_index][bs_iter,:]*self.alpha_phis[l_index][bs_iter,:])
				b_terms = np.dot(gene_eqtl_pmces, np.multiply(self.residual_incl_alpha[bs_iter,:], self.s_inv_2_diag)) + np.dot(self.precomputed_gene_gene_terms[bs_iter], gene_weights)
				a_terms = self.precomputed_a_terms[bs_iter] - .5*(1.0/alpha_component_variances[bs_iter][l_index])

				mixture_alpha_var = -1.0/(2.0*a_terms)
				mixture_alpha_mu = b_terms*mixture_alpha_var
				# HACK TO DEAL WITH MISSING gene models (mixture alpha_mu should never be precisely zero unless there is no gene model)
				mixture_alpha_var[mixture_alpha_mu == 0.0] = 1e-30

				#################
				# Beta effects
				#################
				variant_b_terms = self.residual[bs_iter,:]*self.s_inv_2_diag
				variant_a_terms = (-.5*self.D_diag) - .5*(1.0/beta_component_variances[bs_iter][l_index])

				mixture_beta_var = -1.0/(2.0*variant_a_terms)
				mixture_beta_mu = variant_b_terms*mixture_beta_var

			
				################
				# Normalization (across beta and alpha)
				###############
				if self.bootstrap_prior:
					un_normalized_lv_alpha_weights = expected_log_pi[:,bs_iter] - (.5*np.log(alpha_component_variances[bs_iter][l_index])) + (.5*np.square(mixture_alpha_mu)/mixture_alpha_var) + (.5*np.log(mixture_alpha_var))
					un_normalized_lv_beta_weights = expected_log_variant_pi[:,bs_iter] - (.5*np.log(beta_component_variances[bs_iter][l_index])) + (.5*np.square(mixture_beta_mu)/mixture_beta_var) + (.5*np.log(mixture_beta_var))
				else:
					un_normalized_lv_alpha_weights = expected_log_pi - (.5*np.log(alpha_component_variances[bs_iter][l_index])) + (.5*np.square(mixture_alpha_mu)/mixture_alpha_var) + (.5*np.log(mixture_alpha_var))
					un_normalized_lv_beta_weights = expected_log_variant_pi - (.5*np.log(beta_component_variances[bs_iter][l_index])) + (.5*np.square(mixture_beta_mu)/mixture_beta_var) + (.5*np.log(mixture_beta_var))

				normalizing_term = scipy.special.logsumexp(np.hstack((un_normalized_lv_alpha_weights, un_normalized_lv_beta_weights)))


				# Update agg_gene_trait pt 1
				self.agg_gene_trait[bs_iter] = self.agg_gene_trait[bs_iter] - (self.alpha_mus[l_index][bs_iter,:]*self.alpha_phis[l_index][bs_iter,:])

				# Save results to global model parameters
				mixture_alpha_phi = np.exp(un_normalized_lv_alpha_weights-normalizing_term)
				self.alpha_phis[l_index][bs_iter,:] = mixture_alpha_phi
				self.alpha_mus[l_index][bs_iter,:] = mixture_alpha_mu
				self.alpha_vars[l_index][bs_iter,:] = mixture_alpha_var

				mixture_beta_phi = np.exp(un_normalized_lv_beta_weights-normalizing_term)
				self.beta_phis[l_index][bs_iter,:] = mixture_beta_phi
				self.beta_mus[l_index][bs_iter,:] = mixture_beta_mu
				self.beta_vars[l_index][bs_iter,:] = mixture_beta_var


				# Update KL Terms
				if self.bootstrap_prior:
					lbf = np.hstack((un_normalized_lv_alpha_weights-expected_log_pi[:, bs_iter], un_normalized_lv_beta_weights-expected_log_variant_pi[:, bs_iter]))
					# Save LBF
					self.alpha_lbfs[l_index][bs_iter,:] = un_normalized_lv_alpha_weights-expected_log_pi[:, bs_iter]
					self.beta_lbfs[l_index][bs_iter,:] = un_normalized_lv_beta_weights-expected_log_variant_pi[:, bs_iter]
					# Update KL terms
					maxlbf = np.max(lbf)
					ww = np.exp(lbf - maxlbf)
					ww_weighted = ww*np.exp(np.hstack((expected_log_pi[:, bs_iter], expected_log_variant_pi[:, bs_iter])))
				else:
					lbf = np.hstack((un_normalized_lv_alpha_weights-expected_log_pi, un_normalized_lv_beta_weights-expected_log_variant_pi))
					# Save LBF
					self.alpha_lbfs[l_index][bs_iter,:] = un_normalized_lv_alpha_weights-expected_log_pi
					self.beta_lbfs[l_index][bs_iter,:] = un_normalized_lv_beta_weights-expected_log_variant_pi
					# Update KL terms
					maxlbf = np.max(lbf)
					ww = np.exp(lbf - maxlbf)
					ww_weighted = ww*np.exp(np.hstack((expected_log_pi, expected_log_variant_pi)))
				kl_term1 = -(np.log(np.sum(ww_weighted)) + maxlbf) # THIS TERM IS CORRECT
				kl_term2 = np.sum((mixture_beta_mu*mixture_beta_phi)*variant_b_terms) + np.sum((mixture_alpha_mu*mixture_alpha_phi)*b_terms)
				kl_term3 = -.5*(np.sum((self.NN - 1)*mixture_beta_phi*(np.square(mixture_beta_mu) + mixture_beta_var)) + np.sum((self.NN - 1)*mixture_alpha_phi*(np.square(mixture_alpha_mu) + mixture_alpha_var)))
				self.KL_terms[l_index][bs_iter] = kl_term1 + kl_term2 + kl_term3

				# Update agg_gene_trait pt 2
				self.agg_gene_trait[bs_iter] = self.agg_gene_trait[bs_iter] + (self.alpha_mus[l_index][bs_iter,:]*self.alpha_phis[l_index][bs_iter,:])

				# Update gene trait pred
				self.gene_trait_pred[l_index][bs_iter,:] = np.dot((self.alpha_mus[l_index][bs_iter,:])*(self.alpha_phis[l_index][bs_iter,:]), gene_eqtl_pmces)


	def update_component_variances(self):
		# NOTE: COMPONENT VARIANCE IS SHARED ACROSS BETA AND ALPHA: Maybe a bad idea??
		for bs_iter in range(self.n_bs):
			for l_iter in range(self.L):
				self.component_variances[bs_iter][l_iter] = np.sum((np.square(self.alpha_mus[l_iter][bs_iter,:]) + self.alpha_vars[l_iter][bs_iter,:])*self.alpha_phis[l_iter][bs_iter,:]) + np.sum((np.square(self.beta_mus[l_iter][bs_iter,:]) + self.beta_vars[l_iter][bs_iter,:])*self.beta_phis[l_iter][bs_iter,:])


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


	def initialize_variables(self, twas_data_obj, phi_init, mu_init, mu_var_init):
		# Number of genes
		self.G = len(twas_data_obj['genes'])
		# Number of variants
		self.K = len(twas_data_obj['variants'])
		# Gene names
		self.genes = twas_data_obj['genes']
		# GWAS sample size
		self.NN = twas_data_obj['gwas_sample_size']

		if self.variant_init_log_pi is None:
			gene_pi = np.ones(self.G)/(self.G + self.K)
			variant_pi = np.ones(self.K)/(self.G + self.K)
			self.gene_log_pi = np.log(gene_pi)
			self.variant_log_pi = np.log(variant_pi)

		else:
			self.gene_log_pi = np.copy(self.gene_init_log_pi)
			self.variant_log_pi = np.copy(self.variant_init_log_pi)


		# Generate S matrix
		self.gwas_variant_z = twas_data_obj['gwas_beta']/twas_data_obj['gwas_beta_se']
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

		# Number of bootstrapped
		self.n_bs = len(twas_data_obj['sparse_sampled_gene_eqtl_pmces'])

		# Get eqtl PMCES
		self.precomputed_a_terms = []
		self.precomputed_gene_gene_terms = []

		self.component_variances = []
		self.agg_gene_trait = []
		self.nominal_twas_z = []

		bs_eqtls_pmces = np.zeros((self.G, self.K))
		for bs_iter in range(self.n_bs):
			bs_eqtls_pmces = bs_eqtls_pmces*0.0
			#bs_eqtls_pmces_sparse = np.load(twas_data_obj['bs_gene_eqtl_pmces_files'][bs_iter])
			bs_eqtls_pmces_sparse = twas_data_obj['sparse_sampled_gene_eqtl_pmces'][bs_iter]
			bs_eqtls_pmces = fill_in_causal_effect_size_matrix(bs_eqtls_pmces, bs_eqtls_pmces_sparse)
			#bs_precomputed_gene_gene_terms_old = -np.dot(np.dot(bs_eqtls_pmces,D_mat), np.transpose(bs_eqtls_pmces))

			# Only perform computation for snps with non-zero eqtl effect in at least one gene
			non_zero_snp_indices = np.sum(bs_eqtls_pmces!=0.0,axis=0) != 0.0

			bs_precomputed_gene_gene_terms = -np.dot(np.dot(bs_eqtls_pmces[:, non_zero_snp_indices],D_mat[:, non_zero_snp_indices][non_zero_snp_indices,:]), np.transpose(bs_eqtls_pmces[:, non_zero_snp_indices]))

			bs_precomputed_a_terms = .5*np.diag(bs_precomputed_gene_gene_terms)
			self.precomputed_a_terms.append(bs_precomputed_a_terms)
			self.precomputed_gene_gene_terms.append(bs_precomputed_gene_gene_terms)

			self.component_variances.append(np.ones(self.L)*1e4)
			self.agg_gene_trait.append(np.zeros(self.G))

			# Compute nominal twas zs
			'''
			temp_arr = []
			for gene_iter in range(self.G):
				weights = bs_eqtls_pmces[gene_iter,:]
				#twas_z = np.dot(weights, gwas_z)/np.sqrt(np.dot(np.dot(weights, twas_data_obj['reference_ld']), weights))
				twas_z = np.dot(weights, gwas_z) # Don't need denominator because weights are already standardized
				temp_arr.append(twas_z)
			temp_arr = np.asarray(temp_arr)
			'''
			bs_nominal_twas_z = np.dot(bs_eqtls_pmces, self.gwas_variant_z)
			self.nominal_twas_z.append(bs_nominal_twas_z)


		self.alpha_mus = []
		self.alpha_vars = []
		self.alpha_phis = []
		self.beta_mus = []
		self.beta_vars = []
		self.beta_phis = []
		self.gene_trait_pred = []
		self.KL_terms = []
		self.alpha_lbfs = []
		self.beta_lbfs = []

		if phi_init is not None and mu_init is not None and mu_var_init is not None:
			# Pre-specified initialization
			print('pre-specified')
			for l_iter in range(self.L):
				# Initialize variational distributions defining alphas (the causal effect of genetically-predicted expression in each gene on the trait)
				self.alpha_mus.append(np.repeat(np.reshape(mu_init[l_iter,:],(1, mu_init.shape[1])),self.n_bs,axis=0)[:,:self.G])
				self.alpha_vars.append(np.repeat(np.reshape(mu_var_init[l_iter,:],(1, mu_var_init.shape[1])),self.n_bs,axis=0)[:,:self.G])
				self.alpha_phis.append(np.repeat(np.reshape(phi_init[l_iter,:],(1, phi_init.shape[1])),self.n_bs,axis=0)[:,:self.G])

				# Initialize variational distribution defining betas (the causal effect of pleiotropic genotype on the trait)
				self.beta_mus.append(np.repeat(np.reshape(mu_init[l_iter,:],(1, mu_init.shape[1])),self.n_bs,axis=0)[:,self.G:])
				self.beta_vars.append(np.repeat(np.reshape(mu_var_init[l_iter,:],(1, mu_var_init.shape[1])),self.n_bs,axis=0)[:,self.G:])
				self.beta_phis.append(np.repeat(np.reshape(phi_init[l_iter,:],(1, phi_init.shape[1])),self.n_bs,axis=0)[:,self.G:])
				
				# Initialize quantity to keep track of predicted gene-trait effects
				self.gene_trait_pred.append(np.zeros((self.n_bs, self.K)))
				if np.sum(np.sum(np.repeat(np.reshape(mu_init[l_iter,:],(1, mu_init.shape[1])),self.n_bs,axis=0)[:,:self.G])) != 0.0:
					print('assumption erroror: non-zedro initialized gene effect')
					pdb.set_trace()

				# Initialize KL terms
				self.KL_terms.append(np.zeros(self.n_bs))
				# Initialize lbf
				self.alpha_lbfs.append(np.zeros((self.n_bs, self.G)))
				self.beta_lbfs.append(np.zeros((self.n_bs, self.K)))
		else:
			for l_iter in range(self.L):
				# Initialize variational distributions defining alphas (the causal effect of genetically-predicted expression in each gene on the trait)
				# Currently using null intitialization
				self.alpha_mus.append(np.zeros((self.n_bs, self.G)))
				self.alpha_vars.append(np.ones((self.n_bs, self.G)))
				self.alpha_phis.append(np.ones((self.n_bs, self.G))/(self.G + self.K))

				# Initialize variational distribution defining betas (the causal effect of pleiotropic genotype on the trait)
				self.beta_mus.append(np.zeros((self.n_bs, self.K)))
				self.beta_vars.append(np.ones((self.n_bs, self.K)))
				self.beta_phis.append(np.ones((self.n_bs, self.K))/(self.G + self.K))

				# Initialize quantity to keep track of predicted gene-trait effects
				self.gene_trait_pred.append(np.zeros((self.n_bs, self.K)))

				# Initialize KL terms
				self.KL_terms.append(np.zeros(self.n_bs))
				# Initialize lbf
				self.alpha_lbfs.append(np.zeros((self.n_bs, self.G)))
				self.beta_lbfs.append(np.zeros((self.n_bs, self.K)))


		#self.pmces_files = twas_data_obj['bs_gene_eqtl_pmces_files']
		self.sparse_sampled_gene_eqtl_pmces = twas_data_obj['sparse_sampled_gene_eqtl_pmces']




