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


class TGFM_FM(object):
	def __init__(self, max_iter=100, convergence_thresh=1e-5, likelihood_version='probabilistic_full', mean_component_boolean=True):
		# Prior on gamma distributions defining residual variance and
		self.max_iter = max_iter
		self.convergence_thresh = convergence_thresh
		self.likelihood_version = likelihood_version
		self.mean_component_boolean = mean_component_boolean
	def fit(self, fm_data_obj):
		""" Fit the model.
			Args:
			fm_data_obj
		"""
		print('###############################')
		print('TGFM Fine Mapping')
		print('###############################')
		self.initialize_variables(fm_data_obj)
		for em_iter in range(self.max_iter):
			self.update_residual_variance()
			self.update_posterior_probabilities()
		self.posterior_prob_df = pd.DataFrame(data={'tissue': self.tissue_names, 'gene': self.gene_names, 'eqtl_component_number': self.eqtl_component_numbers, 'posterior_probability':self.posterior_probs})
		self.generate_tissue_level_posterior_probs(fm_data_obj['ordered_tissue_names'])

	def generate_tissue_level_posterior_probs(self, ordered_tissue_names):
		tiss_posterior_probs = np.zeros(len(ordered_tissue_names))
		for ii, tissue_name in enumerate(ordered_tissue_names):
			tissue_indices = self.tissue_names == tissue_name
			if np.sum(tissue_indices) == 0:
				tiss_posterior_probs[ii] = np.nan
			else:
				tiss_posterior_probs[ii] = np.sum(self.posterior_probs[tissue_indices])
		self.tissue_posterior_prob_df = pd.DataFrame(data={'tissue': ordered_tissue_names, 'posterior_probability':tiss_posterior_probs})
	def debug_posterior_probs(self):
		likelihood_version = 'expected_value'
		log_likelihoods = []
		for nn in range(self.N):
			if likelihood_version == 'expected_value':
				log_like = compute_independent_multivariate_gaussian_log_likelihood(self.gwas_component_pmces, self.pred_trait_effects[nn,:], self.residual_variance, self.V)
			elif likelihood_version == 'probabilistic_eqtl':
				log_like = compute_probabilistic_independent_multivariate_gaussian_log_likelihood(self.gwas_component_pmces, self.pred_trait_effects[nn,:], np.square(self.gwas_component_pmces), self.pred_trait_effects_squared[nn,:], self.residual_variance, self.V)
			elif likelihood_version == 'probabilistic_full':
				log_like = compute_probabilistic_independent_multivariate_gaussian_log_likelihood(self.gwas_component_pmces, self.pred_trait_effects[nn,:], self.gwas_component_effects_squared, self.pred_trait_effects_squared[nn,:], self.residual_variance, self.V)
			log_likelihoods.append(log_like)
		log_likelihoods2 = []
		likelihood_version = 'probabilistic_full'
		for nn in range(self.N):
			if likelihood_version == 'expected_value':
				log_like = compute_independent_multivariate_gaussian_log_likelihood(self.gwas_component_pmces, self.pred_trait_effects[nn,:], self.residual_variance, self.V)
			elif likelihood_version == 'probabilistic_eqtl':
				log_like = compute_probabilistic_independent_multivariate_gaussian_log_likelihood(self.gwas_component_pmces, self.pred_trait_effects[nn,:], np.square(self.gwas_component_pmces), self.pred_trait_effects_squared[nn,:], self.residual_variance, self.V)
			elif likelihood_version == 'probabilistic_full':
				log_like = compute_probabilistic_independent_multivariate_gaussian_log_likelihood(self.gwas_component_pmces, self.pred_trait_effects[nn,:], self.gwas_component_effects_squared, self.pred_trait_effects_squared[nn,:], self.residual_variance, self.V)
			log_likelihoods2.append(log_like)
		return np.asarray(log_likelihoods), np.asarray(log_likelihoods2)

	def update_posterior_probabilities(self):
		log_likelihoods = []
		for nn in range(self.N):
			if self.likelihood_version == 'expected_value':
				log_like = compute_independent_multivariate_gaussian_log_likelihood(self.gwas_component_pmces, self.pred_trait_effects[nn,:], self.residual_variance, self.V)
			elif self.likelihood_version == 'probabilistic_eqtl':
				log_like = compute_probabilistic_independent_multivariate_gaussian_log_likelihood(self.gwas_component_pmces, self.pred_trait_effects[nn,:], np.square(self.gwas_component_pmces), self.pred_trait_effects_squared[nn,:], self.residual_variance, self.V)
			elif self.likelihood_version == 'probabilistic_full':
				log_like = compute_probabilistic_independent_multivariate_gaussian_log_likelihood(self.gwas_component_pmces, self.pred_trait_effects[nn,:], self.gwas_component_effects_squared, self.pred_trait_effects_squared[nn,:], self.residual_variance, self.V)
			log_likelihoods.append(log_like)
		if self.mean_component_boolean:
			if self.likelihood_version == 'probabilistic_full' or self.likelihood_version == 'probabilistic_eqtl':
				mean_log_like = compute_probabilistic_independent_multivariate_gaussian_log_likelihood(self.gwas_component_pmces, self.mean_pmces, self.gwas_component_effects_squared, np.square(self.mean_pmces), self.residual_variance, self.V)
			elif self.likelihood_version == 'expected_value':
				mean_log_like = compute_independent_multivariate_gaussian_log_likelihood(self.gwas_component_pmces, self.mean_pmces, self.residual_variance, self.V)
			log_likelihoods.append(mean_log_like)
			temp_post_prob = np.exp(log_likelihoods - scipy.special.logsumexp(log_likelihoods))
			self.posterior_probs = temp_post_prob[:self.N]
			self.mean_posterior_prob = temp_post_prob[-1]
		else:
			self.posterior_probs = np.exp(log_likelihoods - scipy.special.logsumexp(log_likelihoods))
	def update_residual_variance(self):
		numerator = 0
		denominator = 0
		if self.likelihood_version == 'expected_value':
			squared_diff = np.square(self.pred_trait_effects - self.gwas_component_pmces)
		elif self.likelihood_version == 'probabilistic_eqtl' or self.likelihood_version == 'probabilistic_full':
			squared_diff = self.pred_trait_effects_squared + self.gwas_component_effects_squared - 2.0*self.pred_trait_effects*self.gwas_component_pmces
		for nn in range(self.N):
			numerator = numerator + np.sum(self.posterior_probs[nn]*squared_diff[nn,:])
			denominator = denominator + self.posterior_probs[nn]*len(squared_diff[nn,:])
		if self.mean_component_boolean:
			mean_squared_diff = np.square(self.gwas_component_pmces - self.mean_pmces)
			numerator = numerator + np.sum(self.mean_posterior_prob*mean_squared_diff)
			denominator = denominator + self.mean_posterior_prob*len(mean_squared_diff)
		self.residual_variance = numerator/denominator


	def initialize_variables(self, fm_data_obj):
		# Number of genes
		self.G = len(fm_data_obj['genes'])
		# First get predicted trait effects in each eQTL component
		predicted_trait_effects = []
		predicted_trait_effects_squared = []
		gene_name_arr = []
		tissue_name_arr = []
		eqtl_component_num_arr = []
		for gene_num in range(self.G):
			# Predicted trait effects for this gene (mediated through expression)
			gene_pred_trait_effects = fm_data_obj['eqtl_pmces'][gene_num]*fm_data_obj['twas_alpha'][gene_num]
			
			gene_pred_trait_effects_squared = (np.square(fm_data_obj['twas_alpha'][gene_num]) + np.square(fm_data_obj['twas_alpha_sd'][gene_num]))*((fm_data_obj['eqtl_alpha'][gene_num])*(np.square(fm_data_obj['eqtl_mu'][gene_num]) + np.square(fm_data_obj['eqtl_mu_sd'][gene_num])))


			num_components = gene_pred_trait_effects.shape[0]

			gene_name = fm_data_obj['genes'][gene_num].split('_')[0]
			tissue_name = '_'.join(fm_data_obj['genes'][gene_num].split('_')[1:])


			for component_num in range(num_components):
				if np.var(gene_pred_trait_effects[component_num,:]) > 0:
					predicted_trait_effects.append(gene_pred_trait_effects[component_num,:])
					predicted_trait_effects_squared.append(gene_pred_trait_effects_squared[component_num,:])
					gene_name_arr.append(gene_name)
					tissue_name_arr.append(tissue_name)
					eqtl_component_num_arr.append(component_num) 

		self.pred_trait_effects = np.asarray(predicted_trait_effects)
		self.pred_trait_effects_squared = np.asarray(predicted_trait_effects_squared)
		self.gene_names = np.asarray(gene_name_arr)
		self.tissue_names = np.asarray(tissue_name_arr)
		self.eqtl_component_numbers = np.asarray(eqtl_component_num_arr)
		# Get number of eQTL components
		self.N = len(self.gene_names)

		# extract gwas componetn pmces
		self.gwas_component_pmces = fm_data_obj['gwas_component_pmces']
		# extract gwas componetn effects_squared
		self.gwas_component_effects_squared = (np.square(fm_data_obj['gwas_susie_mu']) + np.square(fm_data_obj['gwas_susie_mu_sd']))*fm_data_obj['gwas_susie_alpha']
		# get number of variants
		self.V = len(self.gwas_component_pmces)

		# Extract correlations between predicted eqtl trait effect sizes and gwas effect sizes
		self.correlations = []
		for nn in range(self.N):
			self.correlations.append(np.corrcoef(self.gwas_component_pmces, self.pred_trait_effects[nn,:])[0,1])
		self.correlations = np.asarray(self.correlations)

		# Mean effects
		self.mean_pmces = np.mean(self.gwas_component_pmces)*np.ones(len(self.gwas_component_pmces))

		# Initialize fine mapped probs
		self.posterior_probs = np.ones(self.N)/self.N
		if self.mean_component_boolean:
			self.posterior_probs = np.ones(self.N)/(self.N+1)
			self.mean_posterior_prob = 1.0/(self.N+1)

		# Initialize residual variance (currently not necessary as it is updated first)
		self.residual_variance = 1.0
