import sys
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import time




def calculate_distance_between_two_vectors(vec1, vec2):
	#dist = np.sqrt(np.sum(np.square(vec1-vec2)))
	dist = np.mean(np.abs(vec1-vec2))
	return dist





class SPARSE_SLDSC(object):
	def __init__(self, L=10, max_iter=10000, convergence_thresh=1e-20, estimate_prior_variance=True, ard_a_prior=1e-16, ard_b_prior=1e-16):
		# Prior on gamma distributions defining residual variance and
		self.L = L
		self.max_iter = max_iter
		self.convergence_thresh = convergence_thresh
		self.estimate_prior_variance = estimate_prior_variance
		self.ard_a_prior = ard_a_prior
		self.ard_b_prior = ard_b_prior
	def fit(self, tau, tau_cov):
		""" Fit the model.
			Args:
			twas_data_obj
		"""

		print('###############################')
		print('SPARSE SLDSC')
		print('###############################')

		#####################
		# Initialize variables
		self.initialize_variables(tau, tau_cov)


		print('Initialization complete')
		self.iter = 0
		# Compute residual
		self.residual = np.copy(tau)


		# BEGIN CAVI
		for itera in range(self.max_iter):
			# Update beta (ie predicted causal effect)
			self.update_susie_effects(np.log(self.pi), self.component_variances)


			if self.estimate_prior_variance:
				self.update_component_variances()

			# Convergence stuff
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



	def update_susie_effects(self, expected_log_pi, component_variances):
		# Loop through components
		for l_index in range(self.L):

			# Include current component (as it is currently removed)
			component_pred = self.beta_mu[l_index,:]*self.beta_phi[l_index,:]
			self.residual = self.residual + component_pred


			#################
			# Beta effects
			#################
			beta_b_terms = np.dot(self.residual,self.tau_cov_inv)
			beta_a_terms = (-.5*np.diag(self.tau_cov_inv)) - .5*(1.0/component_variances[l_index])

			mixture_beta_var = -1.0/(2.0*beta_a_terms)
			mixture_beta_mu = beta_b_terms*mixture_beta_var

			
			################
			# Normalization (across beta and alpha)
			###############
			un_normalized_lv_beta_weights = expected_log_pi - (.5*np.log(component_variances[l_index])) + (.5*np.square(mixture_beta_mu)/mixture_beta_var) + (.5*np.log(mixture_beta_var))

			normalizing_term = scipy.special.logsumexp(np.hstack((un_normalized_lv_beta_weights)))

			# Save results to global model parameters
			self.beta_phi[l_index,:] = np.exp(un_normalized_lv_beta_weights-normalizing_term)

			self.beta_mu[l_index,:] = mixture_beta_mu
			self.beta_var[l_index,:] = mixture_beta_var


			# Remove current component (as it is currently removed)
			component_pred = self.beta_mu[l_index,:]*self.beta_phi[l_index,:]
			self.residual = self.residual - component_pred



	def update_component_variances(self):
		self.component_variances = np.sum((np.square(self.beta_mu) + self.beta_var)*self.beta_phi,axis=1)/np.sum(self.beta_phi,axis=1)



	def initialize_variables(self, tau, tau_cov):
		# Number of predictors
		self.K = len(tau)

		# Initialize prior on pi
		self.pi = np.ones(self.K)/self.K

		# Initialize component variances
		self.component_variances = np.ones(self.L)

		# Initialize variational distribution defining betas (the causal effects)
		self.beta_mu = np.zeros((self.L, self.K))
		self.beta_var = np.ones((self.L, self.K))
		self.beta_phi = np.ones((self.L, self.K))/(self.K)

		# Initialize parameter vector to keep of past rounds variables (for convergence purposes)
		self.prev_beta_mu = np.copy(self.beta_mu)

		# Keep track of convergence
		self.converged = False
		self.convergence_tracker = []

		# Compute invrse tau_cov
		self.tau_cov_inv = np.linalg.inv(tau_cov)
