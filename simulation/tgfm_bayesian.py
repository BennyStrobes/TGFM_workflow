import sys
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import time

def sample_eqtl_effects_from_susie_distribution(gene_susie_mu, gene_susie_alpha, gene_susie_mu_var):

	n_components = gene_susie_mu.shape[0]
	n_snps = gene_susie_mu.shape[1]

	sampled_eqtl_effects = np.zeros(n_snps)

	for component_iter in range(n_components):
		# Randomly draw snp index for component
		random_snp_index = np.random.choice(np.arange(n_snps).astype(int), replace=False, p=gene_susie_alpha[component_iter,:])
		
		effect_size_mean = gene_susie_mu[component_iter,random_snp_index]
		effect_size_var = gene_susie_mu_var[component_iter, random_snp_index]

		random_effect_size = np.random.normal(loc=effect_size_mean, scale=np.sqrt(effect_size_var))

		sampled_eqtl_effects[random_snp_index] = sampled_eqtl_effects[random_snp_index] + random_effect_size

	return sampled_eqtl_effects


def calculate_gene_variance_according_to_susie_distribution(susie_mu, susie_alpha, susie_mu_sd, ld):
	gene_var = 0.0

	# Component level eqtl effect sizes for this gene		
	gene_component_effect_sizes = (susie_mu)*susie_alpha

	# eQTL effect sizes for this gene
	gene_eqtl_effect_sizes = np.sum(gene_component_effect_sizes,axis=0)


	num_susie_components = susie_mu.shape[0]
	for k_index in range(num_susie_components):
		gene_var = gene_var + np.sum((np.square(susie_mu[k_index,:]) + np.square(susie_mu_sd[k_index,:]))*np.diag(ld)*susie_alpha[k_index,:])
		eqtl_component_pmces = (susie_mu[k_index,:])*(susie_alpha[k_index,:])
		gene_var = gene_var - np.dot(np.dot(eqtl_component_pmces,ld), eqtl_component_pmces)
	gene_var = gene_var + np.dot(np.dot(gene_eqtl_effect_sizes,ld), gene_eqtl_effect_sizes)
				
	return gene_var


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

def temp_elbo_calc(z_vec, LD, samp_size, alpha, mu, mu2, KL_terms):
	bb = alpha*mu
	b_bar = np.sum(bb,axis=0)
	postb2 = alpha*mu2
	elbo_term1 = samp_size -1
	elbo_term2 = -2.0*np.sum(np.sqrt(samp_size-1)*b_bar*z_vec)
	elbo_term3 = np.sum(b_bar*np.dot((samp_size-1.0)*LD, b_bar))
	elbo_term4 = - np.sum(np.dot(bb, (samp_size-1.0)*LD)*bb)
	elbo_term5 = np.sum(np.dot(np.diag(LD*(samp_size-1)), np.transpose(postb2)))
	elbo_term7 = (-samp_size/2.0)*np.log(2.0*np.pi)

	elbo = elbo_term7 - .5*(elbo_term1 + elbo_term2 + elbo_term3 + elbo_term4 + elbo_term5) - np.sum(KL_terms)
	return elbo




class TGFM(object):
	def __init__(self, L=10, max_iter=200, convergence_thresh=1e-5, ard_a_prior=1e-16, ard_b_prior=1e-16, gene_init_log_pi=None, variant_init_log_pi=None):
		# Prior on gamma distributions defining residual variance and
		self.L = L
		self.max_iter = max_iter
		self.convergence_thresh = convergence_thresh
		self.ard_a_prior = ard_a_prior
		self.ard_b_prior = ard_b_prior
		self.gene_init_log_pi = gene_init_log_pi
		self.variant_init_log_pi = variant_init_log_pi
	def fit(self, twas_data_obj):
		""" Fit the model.
			Args:
			twas_data_obj
		"""

		print('###############################')
		print('TGFM-BAYESIAN CAUSAL TWAS')
		print('###############################')

		#####################
		# Initialize variables
		self.initialize_variables(twas_data_obj)

		print('Initialization complete')
		self.iter = 0

		# BEGIN CAVI
		for itera in range(self.max_iter):
			# Update alpha (ie predicted effect of each gene on the trait)
			self.update_susie_effects(self.gene_log_pi, self.variant_log_pi, self.component_variances, self.component_variances)

			# Update component variances
			self.update_component_variances()

			# Update eQTL effects in each gene
			self.update_child_eqtl_effects()

			diff = calculate_distance_between_two_vectors(self.alpha_mu, self.prev_alpha_mu)
			print(diff)
			self.convergence_tracker.append(diff)
			self.prev_alpha_mu = np.copy(self.alpha_mu)

			self.iter = self.iter + 1
			if diff <= self.convergence_thresh:
				self.converged = True
				break
		if self.converged == False:
			print('Did not converge after ' + str(self.max_iter) + ' iterations')


	def update_child_eqtl_effects(self):
		# Precompute some quantities of interest
		# Expectation of alpha
		E_alpha = np.sum(self.alpha_mu*self.alpha_phi, axis=0)
		# Expectation of beta
		E_beta = np.sum(self.beta_mu*self.beta_phi, axis=0)

		# Expectation of alpha_i*alpha_j matrix
		E_alpha_alpha = np.dot(E_alpha.reshape(len(E_alpha), 1),E_alpha.reshape(1, len(E_alpha)))
		diag_term = np.zeros(len(E_alpha_alpha))
		for gwas_component in range(self.L):
			E_alpha_component = self.alpha_mu[gwas_component,:]*self.alpha_phi[gwas_component,:]
			E_alpha_alpha = E_alpha_alpha - np.dot(E_alpha_component.reshape(len(E_alpha_component), 1),E_alpha_component.reshape(1, len(E_alpha_component)))
			diag_term = diag_term + self.alpha_phi[gwas_component,:]*(np.square(self.alpha_mu[gwas_component,:]) + self.alpha_var[gwas_component,:])
		np.fill_diagonal(E_alpha_alpha, np.diag(E_alpha_alpha) + diag_term)

		# Expectation of alpha_i,beta_j matrix
		E_alpha_beta = np.dot(E_alpha.reshape(len(E_alpha), 1),E_beta.reshape(1, len(E_beta)))
		for gwas_component in range(self.L):
			E_alpha_component = self.alpha_mu[gwas_component,:]*self.alpha_phi[gwas_component,:]
			E_beta_component = self.beta_mu[gwas_component, :]*self.beta_phi[gwas_component,:]
			E_alpha_beta = E_alpha_beta - np.dot(E_alpha_component.reshape(len(E_alpha_component), 1),E_beta_component.reshape(1, len(E_beta_component)))


		# Now loop through genes and perform update in each gene
		for gene_iter in range(self.G):
			self.update_child_eqtl_effects_for_single_gene(gene_iter, E_alpha, E_alpha_alpha, E_alpha_beta)

	def update_child_eqtl_effects_for_single_gene(self, gene_iter, E_alpha, E_alpha_alpha, E_alpha_beta):
		#############################################
		# First need to precompute some quantities
		#############################################
		# Number of components for this eqtl
		n_eqtl_components =  self.child_eqtl_mus[gene_iter].shape[0]

		# Indices with respect to this eqtl
		gene_eqtl_indices = self.eqtl_indices[gene_iter,:]

		# Precompute correction terms for other elements
		# Note this output is limited to effects on eQTL indices
		tmp = np.transpose(self.gene_eqtl_pmces)*E_alpha_alpha[gene_iter,:]
		correction_term = np.dot(E_alpha_beta[gene_iter,:] + np.sum(tmp,axis=1) - tmp[:, gene_iter], self.D_mat[:, gene_eqtl_indices])

		# Get predicted gene eqtl effects from all components
		gene_eqtl_effects = np.sum(self.child_eqtl_mus[gene_iter]*self.child_eqtl_alphas[gene_iter],axis=0)

		# prior prob
		n_var = np.sum(gene_eqtl_indices)
		eqtl_variant_pi = np.ones(np.sum(n_var))/n_var

		# Loop through eqtl components and perform update
		for l_index in range(n_eqtl_components):
			# Remove effect from current component
			gene_eqtl_effects = gene_eqtl_effects - self.child_eqtl_mus[gene_iter][l_index,:]*self.child_eqtl_alphas[gene_iter][l_index,:]

			# Begin calculating updates
			correction_term_from_other_eqtl_comp = np.dot(E_alpha_alpha[gene_iter,gene_iter]*gene_eqtl_effects, self.D_mat[gene_eqtl_indices,:][:,gene_eqtl_indices])
			b_term1 = self.gwas_beta[gene_eqtl_indices]*self.s_inv_2_diag[gene_eqtl_indices]*E_alpha[gene_iter]
			b_term2 = self.parent_eqtl_mus[gene_iter][l_index,:]/self.parent_eqtl_mu_vars[gene_iter][l_index,:]
			b_terms = b_term1 + b_term2 - correction_term_from_other_eqtl_comp - correction_term
			a_terms = (-.5*self.D_diag[gene_eqtl_indices]*E_alpha_alpha[gene_iter,gene_iter]) - .5*(1.0/self.parent_eqtl_mu_vars[gene_iter][l_index,:])

			mixture_alpha_var = -1.0/(2.0*a_terms)
			mixture_alpha_mu = b_terms*mixture_alpha_var

			# Normalization
			un_normalized_lv_weights = np.log(self.parent_eqtl_alphas[gene_iter][l_index,:]) - (.5*np.log(self.parent_eqtl_mu_vars[gene_iter][l_index,:])) + (.5*np.square(mixture_alpha_mu)/mixture_alpha_var) + (.5*np.log(mixture_alpha_var))
			#un_normalized_lv_weights = np.log(eqtl_variant_pi) - (.5*np.log(self.parent_eqtl_mu_vars[gene_iter][l_index,:])) + (.5*np.square(mixture_alpha_mu)/mixture_alpha_var) + (.5*np.log(mixture_alpha_var))
			normalizing_term = scipy.special.logsumexp(un_normalized_lv_weights)

			# Update global variables
			self.child_eqtl_alphas[gene_iter][l_index,:] = np.exp(un_normalized_lv_weights-normalizing_term)
			self.child_eqtl_mus[gene_iter][l_index,:] = mixture_alpha_mu
			self.child_eqtl_mu_vars[gene_iter][l_index,:] = mixture_alpha_var

			# Re-include updated eqtl effects from current component
			gene_eqtl_effects = gene_eqtl_effects + self.child_eqtl_mus[gene_iter][l_index,:]*self.child_eqtl_alphas[gene_iter][l_index,:]

		# Get current PMCES for this gene
		gene_pmces = np.sum(self.child_eqtl_mus[gene_iter]*self.child_eqtl_alphas[gene_iter],axis=0)
		# Fill in PMCES matrix
		self.gene_eqtl_pmces[gene_iter, gene_eqtl_indices] = gene_pmces

	def update_susie_effects(self, expected_log_pi, expected_log_variant_pi, alpha_component_variances, beta_component_variances):
		#############################################
		# First need to precompute some quantities
		#############################################
		# Loop through genees
		for g_index in range(self.G):

			# Get scaled gene variance to add to precomputed a term
			scaled_gene_variance = calculate_gene_variance_according_to_susie_distribution(self.child_eqtl_mus[g_index], self.child_eqtl_alphas[g_index], np.sqrt(self.child_eqtl_mu_vars[g_index]), self.D_mat[self.eqtl_indices[g_index,:], :][:, self.eqtl_indices[g_index,:]])
			# Add precomputed a term to global counter
			self.precomputed_a_terms[g_index] = -.5*scaled_gene_variance

		# Precompute gene gene terms
		self.precomputed_gene_gene_terms = -np.dot(np.dot(self.gene_eqtl_pmces,self.D_mat), np.transpose(self.gene_eqtl_pmces))
		# Fix diagonal to account for additional variance
		np.fill_diagonal(self.precomputed_gene_gene_terms, self.precomputed_a_terms*2.0)

		# Residualize out existing effects
		gene_trait_pred = np.sum(np.dot(self.alpha_mu*self.alpha_phi, self.gene_eqtl_pmces),axis=0)
		non_med_pred = np.sum(self.beta_mu*self.beta_phi,axis=0)
		self.residual = self.gwas_beta - np.dot(self.srs_inv, gene_trait_pred + non_med_pred)
		self.residual_incl_alpha = self.gwas_beta - np.dot(self.srs_inv, non_med_pred)


		#############################################
		# Now run standard iterative updates
		#############################################
		# Loop through components
		for l_index in range(self.L):

			# Include current component (as it is currently removed)
			component_gene_trait_pred = np.dot(self.alpha_mu[l_index,:]*self.alpha_phi[l_index,:], self.gene_eqtl_pmces)
			component_non_med_pred = self.beta_mu[l_index,:]*self.beta_phi[l_index,:]
			self.residual = self.residual + np.dot(self.srs_inv, component_gene_trait_pred + component_non_med_pred)
			self.residual_incl_alpha = self.residual_incl_alpha + np.dot(self.srs_inv, component_non_med_pred)

			#################
			# Alpha effects
			#################
			# Calculate terms inolved in update
			gene_weights = np.sum(self.alpha_mu*self.alpha_phi,axis=0) - (self.alpha_mu[l_index,:]*self.alpha_phi[l_index,:])
			b_terms = np.dot(self.gene_eqtl_pmces, np.multiply(self.residual_incl_alpha, self.s_inv_2_diag)) + np.dot(self.precomputed_gene_gene_terms, gene_weights)
			a_terms = self.precomputed_a_terms - .5*(1.0/alpha_component_variances[l_index])

			mixture_alpha_var = -1.0/(2.0*a_terms)
			mixture_alpha_mu = b_terms*mixture_alpha_var

			#################
			# Beta effects
			#################
			variant_b_terms = self.residual*self.s_inv_2_diag
			variant_a_terms = (-.5*self.D_diag) - .5*(1.0/beta_component_variances[l_index])

			mixture_beta_var = -1.0/(2.0*variant_a_terms)
			mixture_beta_mu = variant_b_terms*mixture_beta_var


			
			################
			# Normalization (across beta and alpha)
			###############
			un_normalized_lv_alpha_weights = expected_log_pi - (.5*np.log(alpha_component_variances[l_index])) + (.5*np.square(mixture_alpha_mu)/mixture_alpha_var) + (.5*np.log(mixture_alpha_var))
			un_normalized_lv_beta_weights = expected_log_variant_pi - (.5*np.log(beta_component_variances[l_index])) + (.5*np.square(mixture_beta_mu)/mixture_beta_var) + (.5*np.log(mixture_beta_var))

			normalizing_term = scipy.special.logsumexp(np.hstack((un_normalized_lv_alpha_weights, un_normalized_lv_beta_weights)))

			# Save results to global model parameters
			self.alpha_phi[l_index,:] = np.exp(un_normalized_lv_alpha_weights-normalizing_term)
			self.alpha_mu[l_index,:] = mixture_alpha_mu
			self.alpha_var[l_index,:] = mixture_alpha_var

			self.beta_phi[l_index,:] = np.exp(un_normalized_lv_beta_weights-normalizing_term)
			self.beta_mu[l_index,:] = mixture_beta_mu
			self.beta_var[l_index,:] = mixture_beta_var

			# Update KL Terms
			lbf = np.hstack((un_normalized_lv_alpha_weights-expected_log_pi, un_normalized_lv_beta_weights-expected_log_variant_pi))
			#scipy.stats.norm.logpdf(, loc=0, scale=1)
			#pdb.set_trace()
			#betahat = np.hstack((b_terms, variant_b_terms))/(self.NN-1)
			#shat2 = 1.0/(self.NN-1)
			#lbf2 = scipy.stats.norm.logpdf(betahat, loc=0, scale=np.sqrt(alpha_component_variances[l_index] + shat2)) - scipy.stats.norm.logpdf(betahat, loc=0, scale=np.sqrt(shat2))
			#pdb.set_trace()


			maxlbf = np.max(lbf)
			ww = np.exp(lbf - maxlbf)
			ww_weighted = ww*np.exp(np.hstack((expected_log_pi, expected_log_variant_pi)))
			kl_term1 = -(np.log(np.sum(ww_weighted)) + maxlbf) # THIS TERM IS CORRECT
			kl_term2 = np.sum((self.beta_mu[l_index,:]*self.beta_phi[l_index,:])*variant_b_terms) + np.sum((self.alpha_mu[l_index,:]*self.alpha_phi[l_index,:])*b_terms)
			kl_term3 = -.5*(np.sum((self.NN - 1)*self.beta_phi[l_index,:]*(np.square(self.beta_mu[l_index,:]) + self.beta_var[l_index,:])) + np.sum((self.NN - 1)*self.alpha_phi[l_index,:]*(np.square(self.alpha_mu[l_index,:]) + self.alpha_var[l_index,:])))
			self.KL_terms[l_index] = kl_term1 + kl_term2 + kl_term3
			self.LBF_terms[l_index] = -kl_term1

			# Remove current component (as it is currently removed)
			component_gene_trait_pred = np.dot(self.alpha_mu[l_index,:]*self.alpha_phi[l_index,:], self.gene_eqtl_pmces)
			component_non_med_pred = self.beta_mu[l_index,:]*self.beta_phi[l_index,:]
			self.residual = self.residual - np.dot(self.srs_inv, component_gene_trait_pred + component_non_med_pred)
			self.residual_incl_alpha = self.residual_incl_alpha - np.dot(self.srs_inv, component_non_med_pred)



	def update_component_variances(self):
		# NOTE: COMPONENT VARIANCE IS SHARED ACROSS BETA AND ALPHA: Maybe a bad idea??
		self.component_variances = np.sum((np.square(self.alpha_mu) + self.alpha_var)*self.alpha_phi,axis=1) + np.sum((np.square(self.beta_mu) + self.beta_var)*self.beta_phi,axis=1)
	
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


		# Compute nominal twas z-scores
		self.nominal_twas_z = np.zeros(self.G)
		gwas_z = twas_data_obj['gwas_beta']/twas_data_obj['gwas_beta_se']
		self.gwas_z = gwas_z
		for g_index in range(self.G):
			weights = twas_data_obj['gene_eqtl_pmces'][g_index,:]
			twas_z = np.dot(weights, gwas_z)/np.sqrt(np.dot(np.dot(weights, twas_data_obj['reference_ld']), weights))
			self.nominal_twas_z[g_index] = twas_z


		# Initialize variational distributions defining alphas (the causal effect of genetically-predicted expression in each gene on the trait)
		# Currently using null intitialization
		self.alpha_mu = np.zeros((self.L, self.G))
		self.alpha_var = np.ones((self.L, self.G))
		self.alpha_phi = np.ones((self.L, self.G))/(self.G + self.K)
		
		# Initialize component variances
		self.component_variances = np.ones(self.L)*1e4

		self.prev_alpha_mu = np.copy(self.alpha_mu)

		# Initialize variational distribution defining betas (the causal effect of pleiotropic genotype on the trait)
		self.beta_mu = np.zeros((self.L, self.K))
		self.beta_var = np.ones((self.L, self.K))
		self.beta_phi = np.ones((self.L, self.K))/(self.G + self.K)

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
		self.D_mat = D_mat
		#self.LD = twas_data_obj['reference_ld']



		# Initialize and
		# Standardize parent eqtl data
		self.parent_eqtl_mus = []
		self.parent_eqtl_mu_vars = []
		self.parent_eqtl_alphas = []
		self.child_eqtl_mus = []
		self.child_eqtl_mu_vars = []
		self.child_eqtl_alphas = []
		self.eqtl_indices = np.copy(twas_data_obj['gene_susie_indices'])
		self.eqtl_components = []
		# Loop through genes
		for gene_iter in range(self.G):
			# Get valid components for the gene
			valid_components = []
			for component_iter in range(twas_data_obj['gene_susie_mu'][gene_iter].shape[0]):
				if np.array_equal(twas_data_obj['gene_susie_mu'][gene_iter][component_iter,:], np.zeros(len(twas_data_obj['gene_susie_mu'][gene_iter][component_iter,:]))):
					continue
				valid_components.append(component_iter)
			valid_components = np.asarray(valid_components)
			# Add to global variable
			self.eqtl_components.append(valid_components)

			if len(valid_components) == 0:
				print('assumption eroror')
				pdb.set_trace()

			# Get gene's unormalized pmces
			gene_unnormalized_pmces = np.sum(twas_data_obj['gene_susie_mu'][gene_iter][valid_components, :]*twas_data_obj['gene_susie_alpha'][gene_iter][valid_components, :],axis=0)
			# Indices corresponding to gene
			gene_indices = self.eqtl_indices[gene_iter]
			# Compute standard deviation of the gene based on pmces
			pmces_gene_sd = np.sqrt(np.dot(np.dot(gene_unnormalized_pmces, twas_data_obj['reference_ld'][gene_indices,:][:, gene_indices]), gene_unnormalized_pmces))
			#full_gene_sd = np.sqrt(calculate_gene_variance_according_to_susie_distribution(twas_data_obj['gene_susie_mu'][gene_iter], twas_data_obj['gene_susie_alpha'][gene_iter], np.sqrt(twas_data_obj['gene_susie_mu_var'][gene_iter]), twas_data_obj['reference_ld'][gene_indices,:][:, gene_indices]))


			# Get gene normalized pmces  # Version 1
			gene_normalized_pmces = gene_unnormalized_pmces/pmces_gene_sd

			# Standardize parent data and set to global variable)
			self.parent_eqtl_mus.append(twas_data_obj['gene_susie_mu'][gene_iter][valid_components, :]/pmces_gene_sd)
			self.parent_eqtl_mu_vars.append(twas_data_obj['gene_susie_mu_var'][gene_iter][valid_components, :]/np.square(pmces_gene_sd))
			self.parent_eqtl_alphas.append(twas_data_obj['gene_susie_alpha'][gene_iter][valid_components, :]) # No need to normalize parent proportions

			# Also initialize child eqtl variables to the same
			self.child_eqtl_mus.append(twas_data_obj['gene_susie_mu'][gene_iter][valid_components, :]/pmces_gene_sd)
			self.child_eqtl_mu_vars.append(twas_data_obj['gene_susie_mu_var'][gene_iter][valid_components, :]/np.square(pmces_gene_sd))
			self.child_eqtl_alphas.append(twas_data_obj['gene_susie_alpha'][gene_iter][valid_components, :]) # No need to normalize parent proportions		

		# Intialize eQTL PMCES
		self.gene_eqtl_pmces = np.zeros(twas_data_obj['gene_eqtl_pmces'].shape)
		# Initialize precomputed a-terms
		self.precomputed_a_terms = np.zeros(self.G)
		# Initialize precomputed gene-gene terms
		self.precomputed_gene_gene_terms = np.zeros((self.G, self.G))

		# Initialize Converged to not converged
		self.converged = False

		# Initialize KL terms
		self.KL_terms = np.zeros(self.L)
		self.LBF_terms = np.zeros(self.L)
		self.elbo = 0.0

		# Compute residual
		self.residual = twas_data_obj['gwas_beta']
		self.residual_incl_alpha = twas_data_obj['gwas_beta']
		self.gwas_beta = np.copy(twas_data_obj['gwas_beta'])

		# Compute gene eqtl pmces
		for g_index in range(self.G):
			# Get current PMCES for this gene
			gene_pmces = np.sum(self.child_eqtl_mus[g_index]*self.child_eqtl_alphas[g_index],axis=0)
			# Fill in PMCES matrix
			self.gene_eqtl_pmces[g_index, self.eqtl_indices[g_index,:]] = gene_pmces

		self.convergence_tracker = []


		# PRECOMPUTE a-terms for variational updates (will not change iteration to iteration)
		'''
		# Compute expression correlation matrix
		expression_covariance = np.dot(np.dot(self.gene_eqtl_pmces, twas_data_obj['reference_ld']), np.transpose(self.gene_eqtl_pmces))
		dd = np.diag(1.0/np.sqrt(np.diag(expression_covariance)))
		self.ge_ld = np.dot(np.dot(dd, expression_covariance),dd)


		self.precomputed_a_terms = np.zeros(self.G)

		self.global_pmces = np.zeros(self.K) 

		self.convergence_tracker = []


		self.precomputed_gene_gene_terms = -np.dot(np.dot(self.gene_eqtl_pmces,D_mat), np.transpose(self.gene_eqtl_pmces))
		self.precomputed_a_terms = .5*np.diag(self.precomputed_gene_gene_terms)
		'''


