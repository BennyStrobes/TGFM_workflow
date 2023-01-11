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





class SPARSE_SLDSC_ARD_SOME_FIXED(object):
	def __init__(self, L=10, max_iter=10000, convergence_thresh=1e-20, estimate_prior_variance=True, ard_a_prior=1e-16, ard_b_prior=1e-16, nonneg=True, nonneg_int=93):
		# Prior on gamma distributions defining residual variance and
		self.L = L
		self.max_iter = max_iter
		self.convergence_thresh = convergence_thresh
		self.estimate_prior_variance = estimate_prior_variance
		self.ard_a_prior = ard_a_prior
		self.ard_b_prior = ard_b_prior
		self.nonneg = nonneg
		self.nonneg_int = nonneg_int
	def fit(self, tau, tau_cov, fixed_coefficients):
		""" Fit the model.
			Args:
			twas_data_obj
		"""

		print('###############################')
		print('SPARSE SLDSC')
		print('###############################')

		#####################
		# Initialize variables
		self.initialize_variables(tau, tau_cov, fixed_coefficients)


		print('Initialization complete')
		self.iter = 0
		# Compute residual
		self.residual = np.copy(tau)
		

		# BEGIN CAVI
		for itera in range(self.max_iter):
			# Update beta (ie predicted causal effect)
			self.update_susie_effects(self.component_variances)


			if self.estimate_prior_variance and itera > 500:
				#if self.estimate_prior_variance:
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



	def update_susie_effects(self, component_variances):

		# Include current effect (as it is currently removed)
		for fixed_effect_index in np.random.permutation(self.fixed_coefficients):

			self.residual[fixed_effect_index] = self.residual[fixed_effect_index] + self.fixed_beta_mu[fixed_effect_index]

			# Compute variational updates
			b_term = np.dot(self.residual,self.tau_cov_inv[fixed_effect_index,:])
			a_term = (-.5*self.tau_cov_inv[fixed_effect_index,fixed_effect_index]) - .5*(1.0/self.fixed_effect_variance)
			self.fixed_beta_var[fixed_effect_index] = -1.0/(2.0*a_term)
			self.fixed_beta_mu[fixed_effect_index] = b_term*self.fixed_beta_var[fixed_effect_index]

			# Remove updated current effect
			self.residual[fixed_effect_index] = self.residual[fixed_effect_index] - self.fixed_beta_mu[fixed_effect_index]

		#for random_effect_iter, random_effect_index in enumerate(self.random_coefficients):
		lb_index = np.min(self.random_coefficients)
		for random_effect_index in np.random.permutation(self.random_coefficients):
			random_effect_iter = random_effect_index - lb_index
			self.residual[random_effect_index] = self.residual[random_effect_index] + self.beta_mu[random_effect_iter]

			b_term = np.dot(self.residual,self.tau_cov_inv[random_effect_index,:])
			a_term = (-.5*self.tau_cov_inv[random_effect_index,random_effect_index]) - .5*(1.0/self.component_variances[random_effect_iter])
			self.beta_var[random_effect_iter] = -1.0/(2.0*a_term)
			self.beta_mu[random_effect_iter] = b_term*self.beta_var[random_effect_iter]
			self.residual[random_effect_index] = self.residual[random_effect_index] - self.beta_mu[random_effect_iter]

	def update_component_variances(self):
		for kk in range(len(self.component_variances)):
			self.component_variances[kk] = np.square(self.beta_mu[kk]) + self.beta_var[kk]
			#if kk >= self.nonneg_int and self.iter > 200 and self.nonneg:
			if kk >= self.nonneg_int and self.nonneg:
				if self.beta_mu[kk] < 0.0:
					self.component_variances[kk] = self.component_variances[kk]/(1000.0)



	def initialize_variables(self, tau, tau_cov, fixed_coefficients):
		self.fixed_coefficients = fixed_coefficients

		random_coefficients =[]
		for ii in range(len(tau)):
			if ii not in self.fixed_coefficients:
				random_coefficients.append(ii)
		self.random_coefficients = np.asarray(random_coefficients)

		# Number of predictors
		self.K = len(tau) - len(fixed_coefficients)

		# Initialize component variances
		self.component_variances = np.ones(self.K)

		# Fixed effect variance
		self.fixed_effect_variance = 1.0

		# Initialize variational distribution defining betas (the causal effects)
		self.beta_mu = np.zeros((self.K))
		#self.beta_mu = tau[self.random_coefficients]
		self.beta_var = np.ones((self.K))

		# Intitialize variational distribution defining fixed genotypes
		if len(self.fixed_coefficients) == 0:
			self.fixed_beta_mu = np.zeros(len(self.fixed_coefficients))
		else:
			#self.fixed_beta_mu = tau[self.fixed_coefficients]
			self.fixed_beta_mu = np.zeros(len(self.fixed_coefficients))
		self.fixed_beta_var = np.ones(len(self.fixed_coefficients))


		# Initialize parameter vector to keep of past rounds variables (for convergence purposes)
		self.prev_beta_mu = np.copy(self.beta_mu)

		# Keep track of convergence
		self.converged = False
		self.convergence_tracker = []

		# Compute invrse tau_cov
		self.tau_cov_inv = np.linalg.inv(tau_cov)
