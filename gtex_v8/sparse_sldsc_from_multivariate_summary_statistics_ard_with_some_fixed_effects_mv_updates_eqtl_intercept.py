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





class SPARSE_SLDSC_ARD_SOME_FIXED_MV_UPDATES_EQTL_INTERCEPT(object):
	def __init__(self, L=10, regularization_param=1e-20, max_iter=10000, convergence_thresh=1e-20, learn_fixed_variance=False, estimate_prior_variance=True, ard_a_prior=1e-16, ard_b_prior=1e-16, nonneg=True, nonneg_int=93):
		# Prior on gamma distributions defining residual variance and
		self.L = L
		self.max_iter = max_iter
		self.convergence_thresh = convergence_thresh
		self.estimate_prior_variance = estimate_prior_variance
		self.ard_a_prior = ard_a_prior
		self.ard_b_prior = ard_b_prior
		self.nonneg = nonneg
		self.nonneg_int = nonneg_int
		self.regularization_param = regularization_param
		self.learn_fixed_variance = learn_fixed_variance
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
		self.residual = np.copy(tau) - self.beta_mu - self.alpha_mu
		self.tau = np.copy(tau)

		# BEGIN CAVI
		for itera in range(self.max_iter):
			# Update beta (ie predicted causal effect)
			self.update_susie_effects_mv(self.component_variances)


			#if self.estimate_prior_variance and itera > 10000:
			if self.estimate_prior_variance:
				#if self.estimate_prior_variance:
				self.update_component_variances()

			# Convergence stuff
			diff = calculate_distance_between_two_vectors(self.beta_mu, self.prev_beta_mu)
			self.convergence_tracker.append(diff)
			self.prev_beta_mu = np.copy(self.beta_mu)
			self.iter = self.iter + 1
			print(diff)
			#if diff <= self.convergence_thresh:
				#self.converged = True
				#break
		print(diff)
		if self.converged == False:
			print('Did not converge after ' + str(self.max_iter) + ' iterations')


	def update_susie_effects_mv(self, component_variances):
		# Include effects of alpha
		eqtl_indicator = np.zeros(len(self.residual))
		eqtl_indicator[self.random_coefficients] = eqtl_indicator[self.random_coefficients] + 1.0


		n_anno = self.tau_cov_inv.shape[0]
		big_tau_inv_cov = np.zeros((n_anno+1,n_anno+1))
		big_tau_inv_cov[:n_anno,:n_anno] = self.tau_cov_inv
		anno_eqtl_intercept_precision = np.dot(eqtl_indicator, self.tau_cov_inv)
		big_tau_inv_cov[-1,-1] = np.dot(anno_eqtl_intercept_precision, eqtl_indicator)
		big_tau_inv_cov[-1,:-1] = anno_eqtl_intercept_precision
		big_tau_inv_cov[:-1,-1] = anno_eqtl_intercept_precision

		# Get coefficient prior variances
		coef_prior_variances = np.zeros(len(self.beta_mu))
		if len(self.fixed_coefficients) > 0:
			coef_prior_variances[self.fixed_coefficients] = self.fixed_effect_variance
		coef_prior_variances[self.random_coefficients] = self.component_variances
		coef_prior_precisions = 1.0/coef_prior_variances
		coef_prior_precisions = np.hstack((coef_prior_precisions, np.asarray([1.0/self.shared_eqtl_variance])))

		# Learn covariance
		coef_cov = np.linalg.inv(np.diag(coef_prior_precisions) + big_tau_inv_cov)

		tmp_term1 = np.dot(self.tau_cov_inv, self.tau)
		tmp_term2 = np.dot(np.dot(eqtl_indicator, self.tau_cov_inv), self.tau)
		tmp_term = np.hstack((tmp_term1, np.asarray([tmp_term2])))
		coef_mu = np.dot(coef_cov, tmp_term)

		self.alpha_mu = coef_mu[-1]
		self.alpha_var = coef_cov[-1,-1]
		self.beta_mu = coef_mu[:-1]
		self.beta_cov = coef_cov[:-1,:-1]



	def update_susie_effects_non_mv(self, component_variances):
		# Include effects of alpha
		eqtl_indicator = np.zeros(len(self.residual))
		eqtl_indicator[self.random_coefficients] = eqtl_indicator[self.random_coefficients] + 1.0
		self.residual = self.residual + eqtl_indicator*self.alpha_mu

		self.alpha_var = 1.0/(1.0/self.shared_eqtl_variance + np.dot(np.dot(eqtl_indicator,self.tau_cov_inv), eqtl_indicator))
		self.alpha_mu = self.alpha_var*np.dot(np.dot(eqtl_indicator, self.tau_cov_inv), self.residual)

		# Remove effects alpha
		self.residual = self.residual - eqtl_indicator*self.alpha_mu
		# Include previously removed effects of beta
		self.residual = self.residual + self.beta_mu
	

		# First get coefficient prior variances
		coef_prior_variances = np.zeros(len(self.beta_mu))
		if len(self.fixed_coefficients) > 0:
			coef_prior_variances[self.fixed_coefficients] = self.fixed_effect_variance
		coef_prior_variances[self.random_coefficients] = self.component_variances
		coef_prior_precisions = 1.0/coef_prior_variances


		#self.beta_cov = np.linalg.inv(np.diag(coef_prior_precisions) + self.tau_cov_inv)
		self.beta_cov = np.linalg.solve(np.diag(coef_prior_precisions) + self.tau_cov_inv, np.eye(self.tau_cov_inv.shape[0]))
		self.beta_mu = np.dot(np.dot(self.beta_cov, self.tau_cov_inv), self.residual)


		# Remove effects of beta
		self.residual = self.residual - self.beta_mu


		'''
		# Include current effect (as it is currently removed)
		for fixed_effect_index, fixed_effect_element in enumerate(self.fixed_coefficients):

			self.residual[fixed_effect_element] = self.residual[fixed_effect_element] + self.fixed_beta_mu[fixed_effect_index]

			# Compute variational updates
			b_term = np.dot(self.residual,self.tau_cov_inv[fixed_effect_element,:])
			a_term = (-.5*self.tau_cov_inv[fixed_effect_element,fixed_effect_element]) - .5*(1.0/self.fixed_effect_variance)
			self.fixed_beta_var[fixed_effect_index] = -1.0/(2.0*a_term)
			self.fixed_beta_mu[fixed_effect_index] = b_term*self.fixed_beta_var[fixed_effect_index]

			# Remove updated current effect
			self.residual[fixed_effect_element] = self.residual[fixed_effect_element] - self.fixed_beta_mu[fixed_effect_index]

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
		'''

	def update_component_variances(self):
		#expected_a_terms_shared = ((np.square(self.alpha_mu) + self.alpha_var)/2.0) + 1e-30
		#expected_b_terms_shared = (1.0/2.0) + self.regularization_param
		#self.shared_eqtl_variance = expected_a_terms_shared/expected_b_terms_shared

		if len(self.fixed_coefficients) > 0 and self.learn_fixed_variance:
			fixed_a_term = (np.sum((np.square(self.beta_mu) + np.diag(self.beta_cov))[self.fixed_coefficients])/2.0) + 1e-16
			fixed_b_term = (len(self.fixed_coefficients)/2.0) + 1e-16
			self.fixed_effect_variance = fixed_a_term/fixed_b_term


		expected_a_terms = np.square(self.beta_mu) + np.diag(self.beta_cov)
		expected_a_terms_random_coef = expected_a_terms[self.random_coefficients]
		for kk in range(len(self.component_variances)):
			component_variance_a_term = ((expected_a_terms_random_coef[kk])/2.0) + 1e-30
			component_variance_b_term = (1.0/2.0) + self.regularization_param

			#self.component_variances[kk] = (np.square(self.beta_mu[kk]) + self.beta_var[kk] + 1e-16)/(1.0 + 1e0)
			self.component_variances[kk] = component_variance_a_term/component_variance_b_term
			#if kk >= self.nonneg_int and self.iter > 200 and self.nonneg:
			if kk >= self.nonneg_int and self.nonneg:
				if self.beta_mu[kk] < 0.0:
					self.component_variances[kk] = self.component_variances[kk]/(100000.0)


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
		self.component_variances = np.ones(self.K)*1e10

		# Fixed effect variance
		self.fixed_effect_variance = 1.0
		#self.fixed_effect_variance = np.max(np.square(tau[self.fixed_coefficients])+ (np.diag(tau_cov)[self.fixed_coefficients]))*10

		# Initialize variational distribution defining betas (the causal effects)
		#self.beta_mu = np.zeros((self.K))
		self.beta_mu = tau
		#self.beta_var = np.ones((self.K))
		self.beta_cov = tau_cov

		#self.update_component_variances()
		#pdb.set_trace()

		# Initialize parameter vector to keep of past rounds variables (for convergence purposes)
		self.prev_beta_mu = np.copy(self.beta_mu)

		# Keep track of convergence
		self.converged = False
		self.convergence_tracker = []

		# Compute invrse tau_cov
		#self.tau_cov_inv = np.linalg.inv(tau_cov)
		self.tau_cov_inv = np.linalg.solve(tau_cov, np.eye(tau_cov.shape[0]))
		###########
		###########
		self.alpha_mu = 0.0
		self.alpha_var = 1.0
		self.shared_eqtl_variance = 1.0
		###########
		###########


