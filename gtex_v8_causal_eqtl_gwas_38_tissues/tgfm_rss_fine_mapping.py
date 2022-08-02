import sys
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special



def compute_independent_multivariate_gaussian_log_likelihood(xx, mean_vec, residual_variance, sample_size):
	log_like = -sample_size*np.log(2*np.pi*(residual_variance))/2 - np.sum(((xx-mean_vec)**2)/(2 * (residual_variance)))
	return log_like

def compute_probabilistic_independent_multivariate_gaussian_log_likelihood(xx, mean_vec, xx_squared, mean_vec_squared, residual_variance, sample_size):
	log_like = -sample_size*np.log(2*np.pi*(residual_variance))/2 - np.sum((xx_squared - (2.0*xx*mean_vec) + mean_vec_squared)/(2 * (residual_variance)))
	return log_like


class TGFM_RSS_FM(object):
	def __init__(self, residual_version='regress'):
		self.residual_version = residual_version
	def fit(self, twas_data_obj, tgfm_data_obj):
		""" Fit the model.
			Args:
			twas_data_obj
			tgfm_data_obj
		"""
		print('###############################')
		print('TGFM RSS Fine Mapping')
		print('###############################')
		self.initialize_variables(twas_data_obj, tgfm_data_obj)
		self.update_posterior_probs()
		self.posterior_prob_df = pd.DataFrame(data={'tissue': self.tissue_names, 'gene': self.gene_names, 'eqtl_component_number': self.eqtl_component_numbers, 'posterior_probability':self.posterior_probs, 'log_likelihood':self.log_likelihoods})
		self.generate_tissue_level_posterior_probs(tgfm_data_obj['ordered_tissue_names'])

	def generate_tissue_level_posterior_probs(self, ordered_tissue_names):
		tiss_posterior_probs = np.zeros(len(ordered_tissue_names))
		for ii, tissue_name in enumerate(ordered_tissue_names):
			tissue_indices = self.tissue_names == tissue_name
			if np.sum(tissue_indices) == 0:
				tiss_posterior_probs[ii] = np.nan
			else:
				tiss_posterior_probs[ii] = np.sum(self.posterior_probs[tissue_indices])
		self.tissue_posterior_prob_df = pd.DataFrame(data={'tissue': ordered_tissue_names, 'posterior_probability':tiss_posterior_probs})

	def update_posterior_probs(self):

		log_likelihoods = []
		for nn in range(self.N):
			ll_a_term = np.dot(np.multiply(self.residual, self.s_inv_2_diag), self.pred_trait_effects[nn,:])
			ll_b_term = -.5*self.bDb_effects[nn]
			log_likelihood = ll_a_term + ll_b_term
			log_likelihoods.append(log_likelihood)
		self.posterior_probs = np.exp(log_likelihoods - scipy.special.logsumexp(log_likelihoods))
		self.log_likelihoods = np.asarray(log_likelihoods)

	def initialize_variables(self, twas_data_obj, tgfm_data_obj):
		# Generate S matrix
		s_squared_vec = np.square(twas_data_obj['gwas_beta_se']) + (np.square(twas_data_obj['gwas_beta'])/twas_data_obj['gwas_sample_size'])
		s_vec = np.sqrt(s_squared_vec)
		S_mat = np.diag(s_vec)
		S_inv_mat = np.diag(1.0/s_vec)
		S_inv_2_mat = np.diag(1.0/np.square(s_vec))

		# Compute (S^-1)R(S^-1) taking advantage of fact that S^-1 is a diagonal matrix
		D_mat = np.multiply(np.multiply(np.diag(S_inv_mat)[:, None], twas_data_obj['reference_ld']), np.diag(S_inv_mat))
		d_vec = np.diag(D_mat)
		# Compute (S)R(S^-1) taking advantage of fact that S and S^-1 is a diagonal matrix
		srs_inv_mat = np.multiply(np.multiply(np.diag(S_mat)[:, None], twas_data_obj['reference_ld']), np.diag(S_inv_mat))

		# Generate data object containing statistics that are precomputed
		self.srs_inv = srs_inv_mat
		self.s_inv_2_diag = np.diag(S_inv_2_mat)

		# Generate residual
		if self.residual_version == 'regress':
			gwas_susie_effects_to_regress = np.zeros(tgfm_data_obj['gwas_susie_mu'].shape[1])
			for temp_component in range(tgfm_data_obj['gwas_susie_mu'].shape[0]):
				if temp_component == tgfm_data_obj['gwas_susie_component']:
					continue
				gwas_susie_effects_to_regress = gwas_susie_effects_to_regress + (tgfm_data_obj['gwas_susie_mu'][temp_component,:]*tgfm_data_obj['gwas_susie_alpha'][temp_component,:])
			self.residual = twas_data_obj['gwas_beta'] - np.dot(self.srs_inv, gwas_susie_effects_to_regress)
			self.gwas_component_pmces = tgfm_data_obj['gwas_susie_mu'][tgfm_data_obj['gwas_susie_component'],:]*tgfm_data_obj['gwas_susie_alpha'][tgfm_data_obj['gwas_susie_component'],:]
		elif self.residual_version == 'predict':
			self.gwas_component_pmces = tgfm_data_obj['gwas_susie_mu'][tgfm_data_obj['gwas_susie_component'],:]*tgfm_data_obj['gwas_susie_alpha'][tgfm_data_obj['gwas_susie_component'],:]
			self.residual = np.dot(self.srs_inv, self.gwas_component_pmces)


		# Number of genes
		self.G = len(twas_data_obj['genes'])
		# First get predicted trait effects in each eQTL component
		predicted_trait_effects = []
		bDb_effects = []
		gene_name_arr = []
		tissue_name_arr = []
		eqtl_component_num_arr = []
		twas_z_scores = []

		total_twas_pred_effects = np.zeros(len(twas_data_obj['gwas_beta']))
		for gene_num in range(self.G):
			# Predicted trait effects for this gene (mediated through expression)
			gene_pred_trait_effects = (twas_data_obj['susie_mu'][gene_num]*twas_data_obj['susie_alpha'][gene_num])*tgfm_data_obj['twas_alpha'][gene_num]

			total_twas_pred_effects = total_twas_pred_effects + np.sum(gene_pred_trait_effects,axis=0)
			num_components = gene_pred_trait_effects.shape[0]

			gene_name = twas_data_obj['genes'][gene_num].split('_')[0]
			tissue_name = '_'.join(twas_data_obj['genes'][gene_num].split('_')[1:])

			alpha_g_squared = np.square(tgfm_data_obj['twas_alpha'][gene_num]) + np.square(tgfm_data_obj['twas_alpha_sd'][gene_num])

			for component_num in range(num_components):
				if np.var(gene_pred_trait_effects[component_num,:]) > 0:
					predicted_trait_effects.append(gene_pred_trait_effects[component_num,:])
					gene_name_arr.append(gene_name)
					tissue_name_arr.append(tissue_name)
					eqtl_component_num_arr.append(component_num)

					bdb_effect = alpha_g_squared*np.sum((d_vec)*(twas_data_obj['susie_alpha'][gene_num][component_num,:])*(np.square(twas_data_obj['susie_mu'][gene_num][component_num,:]) + np.square(twas_data_obj['susie_mu_sd'][gene_num][component_num,:])))
					bDb_effects.append(bdb_effect)
					twas_z_scores.append(tgfm_data_obj['twas_alpha'][gene_num]/tgfm_data_obj['twas_alpha_sd'][gene_num])

		# aa=twas_data_obj['gwas_beta'] - np.dot(self.srs_inv,total_twas_pred_effects)
		self.pred_trait_effects = np.asarray(predicted_trait_effects)
		self.bDb_effects = np.asarray(bDb_effects)
		self.gene_names = np.asarray(gene_name_arr)
		self.tissue_names = np.asarray(tissue_name_arr)
		self.eqtl_component_numbers = np.asarray(eqtl_component_num_arr)
		self.twas_z_scores = np.asarray(twas_z_scores)
		# Get number of eQTL components
		self.N = len(self.gene_names)

		# Extract correlations between predicted eqtl trait effect sizes and gwas effect sizes
		self.correlations = []
		for nn in range(self.N):
			self.correlations.append(np.corrcoef(self.residual, self.pred_trait_effects[nn,:])[0,1])
		self.correlations = np.asarray(self.correlations)

		# Initialize fine mapped probs
		self.posterior_probs = np.ones(self.N)/self.N


