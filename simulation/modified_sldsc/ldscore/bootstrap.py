from __future__ import division
import numpy as np
from scipy.optimize import lsq_linear
from scipy.optimize import nnls
import pdb
np.seterr(divide='raise', invalid='raise')


def _check_shape(x, y):
	'''Check that x and y have the correct shapes (for regression jackknives).'''
	if len(x.shape) != 2 or len(y.shape) != 2:
		raise ValueError('x and y must be 2D arrays.')
	if x.shape[0] != y.shape[0]:
		raise ValueError(
			'Number of datapoints in x != number of datapoints in y.')
	if y.shape[1] != 1:
		raise ValueError('y must have shape (n_snp, 1)')
	n, p = x.shape
	if p > n:
		raise ValueError('More dimensions than datapoints.')

	return (n, p)


class partiallyNonNegLstsqBStrap(object):

	def __init__(self, x, y, n_blocks, nonnegative_coefficients, n_bootstraps):
		self.N, self.p = _check_shape(x, y)
		self.n_blocks = n_blocks
		self.separators = self.get_separators(self.N, self.n_blocks)
		self.n_bootstraps = n_bootstraps

		###########################
		# Get bounds
		###########################
		ub = []
		lb = []

		for nonnegative_coefficient in nonnegative_coefficients:
			if nonnegative_coefficient:
				ub.append(np.inf)
				lb.append(0.0)
			else:
				ub.append(np.inf)
				lb.append(-np.inf)
		ub = np.asarray(ub)
		lb = np.asarray(lb)


		###########################
		# Run global regression
		###########################
		res = lsq_linear(x, y[:,0], bounds=(lb, ub), lsmr_tol='auto', verbose=1)
		# Save to object
		self.global_coef = res['x']
		self.global_convergence = res['success']
		if res['success'] != True:
			print('assumption error: partially nonnegative least squares did not converge')


		###########################
		# Split x and y into chunks (blocks)
		###########################
		x_chunks = []
		y_chunks = []
		counter = 0
		for block_iter in range(n_blocks):
			seperator_start = self.separators[block_iter]  # Inclusive
			seperator_end = self.separators[(block_iter+1)]  # Inclusive
			x_chunk = x[seperator_start:seperator_end,:]
			y_chunk = y[seperator_start:seperator_end, :]
			x_chunks.append(x_chunk)
			y_chunks.append(y_chunk)
			if x_chunk.shape[0] <= 1:
				print('assumption eroror')
				pdb.set_trace()
			counter = counter + x_chunk.shape[0]
		# Quick error check
		if counter != x.shape[0]:
			print('assumption eroror')
			pdb.set_trace()
		if len(x_chunks) != self.n_blocks:
			print('assumption eroror')
			pdb.set_trace()

		###########################
		# Run regression for each bootstrapped sample
		###########################
		self.bootstrapped_coefs = []
		self.bootstrapped_convergences = []
		for bs_iter in range(self.n_bootstraps):
			print(bs_iter)
			# Get bootstrap sample indices
			bs_indices = np.random.choice(np.arange(self.n_blocks), size=self.n_blocks, replace=True)

			# Extract bootstrapped data
			bs_x = []
			bs_y = []
			for bs_index in bs_indices:
				bs_x.append(x_chunks[bs_index])
				bs_y.append(y_chunks[bs_index])
			bs_x = np.vstack(bs_x)
			bs_y = np.vstack(bs_y)

			# Run the regression
			bs_res = lsq_linear(bs_x, bs_y[:,0], bounds=(lb, ub), lsmr_tol='auto', verbose=1)

			# Save data
			self.bootstrapped_coefs.append(bs_res['x'])
			self.bootstrapped_convergences.append(bs_res['success'])
			# Report error if model does not converge
			if bs_res['success'] == False:
				print('assumption error: partially nonnegative least squares did not converge')
		# Put data in nice compact data frame
		self.bootstrapped_coefs = np.asarray(self.bootstrapped_coefs)
		self.bootstrapped_convergences = np.asarray(self.bootstrapped_convergences)


	@classmethod
	def get_separators(cls, N, n_blocks):
		'''Define evenly-spaced block boundaries.'''
		return np.floor(np.linspace(0, N, n_blocks + 1)).astype(int)
