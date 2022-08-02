import numpy as np 
import os
import sys
import pdb
import scipy.special




def generate_variance_grid(grid_start, grid_end, multiplicative_grid_scaling_factor):
	# Create fixed variance grid
	variance_grid_arr = []
	grid_value =grid_start
	grid_converged = False
	while grid_converged == False:
		variance_grid_arr.append(grid_value)
		new_grid_value = grid_value*multiplicative_grid_scaling_factor
		if new_grid_value > grid_end:
			grid_converged = True 
		else:
			grid_value = new_grid_value
	variance_grid_arr = np.asarray(variance_grid_arr)
	return variance_grid_arr

def simulate_data(num_samples, num_features, variance_grid, dirichlet_prior=.5, residual_variance=1.0):
	# Simulate features
	X = np.random.normal(size=((num_samples, num_features)))
	for feature_num in range(num_features):
		X[:, feature_num] = (X[:, feature_num] - np.mean(X[:, feature_num]))/np.std(X[:, feature_num])

	# Draw global deltas
	delta = np.random.dirichlet(alpha=np.ones(len(variance_grid))*dirichlet_prior)
	
	# Draw betas and latent variables
	betas = []
	z = []
	for feature_num in range(num_features):
		category = np.random.choice(np.arange(len(variance_grid)), p=delta)
		variance = variance_grid[category]
		beta = np.random.normal(loc=0.0, scale=np.sqrt(variance))
		betas.append(beta)
		z_vec = np.zeros(len(variance_grid))
		z_vec[category] = 1.0
		z.append(z_vec)
	betas = np.asarray(betas)
	z = np.asarray(z)

	# Draw Y (output)
	pred_y = np.dot(X,betas)
	Y = np.random.normal(loc=pred_y, scale=np.sqrt(residual_variance))

	return Y, X, betas, z, delta

def update_beta_and_Z(Y, X, beta_mu, beta_var, Z, variance_grid, expected_ln_pi, residual_precision):
	resid = Y - np.dot(X, beta_mu)
	P = X.shape[1]
	beta_squared = np.zeros(P)

	for pp in range(P):
		resid = resid + (X[:, pp]*beta_mu[pp])

		#a_terms = -.5*np.sum(np.square(X[:,pp])) - (.5/variance_grid)
		#b_term = np.sum(resid*X[:,pp])
		a_terms = -.5*residual_precision*np.sum(np.square(X[:,pp])) - (.5/variance_grid)
		b_term = np.sum(residual_precision*resid*X[:,pp])


		mixture_alpha_var = -1.0/(2.0*a_terms)
		mixture_alpha_mu = b_term*mixture_alpha_var


		un_normalized_lv_weights = expected_ln_pi - (.5*np.log(variance_grid)) + (.5*np.square(mixture_alpha_mu)/mixture_alpha_var) + (.5*np.log(mixture_alpha_var))
		un_normalized_lv_weights = un_normalized_lv_weights - scipy.special.logsumexp(un_normalized_lv_weights)
		Z[pp,:] = np.exp(un_normalized_lv_weights)


		beta_mu[pp] = np.sum(mixture_alpha_mu*Z[pp,:])
		beta_var[pp] = np.sum(mixture_alpha_var*Z[pp,:]) + np.sum(np.square(mixture_alpha_mu)*Z[pp,:]) - np.sum(np.square(mixture_alpha_mu*Z[pp,:]))

		beta_squared[pp] = np.sum(Z[pp,:]*(np.square(mixture_alpha_mu) + mixture_alpha_var))


		resid = resid - (X[:, pp]*beta_mu[pp])

	return beta_mu, beta_var, Z, beta_squared

def update_residual_precision_v3(Y, G, beta_mu, beta_squared, hyper_param=1e-16):

	#beta_squared = np.square(beta_mu) + beta_var
	beta_pred = np.dot(G, beta_mu)
	G_squared = np.square(G)

	resid = np.square(Y) + np.dot(G_squared, beta_squared)
	resid = resid - (2.0*Y*beta_pred)
	# Now add terms with interactions between factors
	resid = resid + (beta_pred*beta_pred - np.dot(G_squared, np.square(beta_mu)))

	new_alpha = hyper_param + (len(Y)/2.0)
	new_beta = hyper_param + (np.sum(resid)/2.0)

	precision = new_alpha/new_beta
	return precision

def inference(Y, X, variance_grid, max_iters=200):
	# Initialize data
	N = X.shape[0]
	P = X.shape[1]
	M = len(variance_grid)
	beta_mu = np.zeros(P)
	beta_var = np.ones(P)
	Z = np.ones((P, M))/M
	delta_alpha = np.ones(M)
	residual_precision = 1.0

	for itera in range(max_iters):
		expected_ln_pi = scipy.special.psi(delta_alpha) - scipy.special.psi(np.sum(delta_alpha))
		beta_mu, beta_var, Z, beta_squared = update_beta_and_Z(Y, X, beta_mu, beta_var, Z, variance_grid, expected_ln_pi, residual_precision)

		delta_alpha = np.sum(Z,axis=0)
		
		residual_precision = update_residual_precision_v3(Y, X, beta_mu, beta_squared, hyper_param=0.0)


		print(delta_alpha/np.sum(delta_alpha))
		print(1.0/residual_precision)



	return Z, delta_alpha, beta_mu, beta_var


def inference_v2(Y, X, variance_grid, max_iters=200):
	# Initialize data
	N = X.shape[0]
	P = X.shape[1]
	M = len(variance_grid)
	beta_mu = np.zeros(P)
	beta_var = np.ones(P)
	Z = np.ones((P, M))/M
	delta_alpha = np.ones(M)
	residual_precision = 1.0

	for itera in range(max_iters):
		expected_ln_pi = scipy.special.psi(delta_alpha) - scipy.special.psi(np.sum(delta_alpha))
		beta_mu, beta_var, Z, beta_squared = update_beta_and_Z(Y, X, beta_mu, beta_var, Z, variance_grid, expected_ln_pi, residual_precision)

		delta_alpha = np.sum(Z,axis=0)
		
		residual_precision = update_residual_precision_v3(Y, X, beta_mu, np.square(beta_mu) + beta_var, hyper_param=0.0)


		print(delta_alpha/np.sum(delta_alpha))
		print(residual_precision)



	return Z, delta_alpha, beta_mu, beta_var


def multivariate_inference(Y, X, variance_grid, max_iters=200):
	# Initialize data
	N = X.shape[0]
	P = X.shape[1]
	M = len(variance_grid)
	beta_mu = np.zeros(P)
	beta_var = np.ones(P)
	Z = np.ones((P, M))/M
	delta_alpha = np.ones(M)

	for itera in range(max_iters):
		expected_ln_pi = scipy.special.psi(delta_alpha) - scipy.special.psi(np.sum(delta_alpha))
		beta_mu, beta_var, Z = multivariate_update_beta_and_Z(Y, X, beta_mu, beta_var, Z, variance_grid, expected_ln_pi)
	pdb.set_trace()

# Create variance grid
grid_start=1e-3
grid_end=100
multiplicative_grid_scaling_factor=5.0
variance_grid = generate_variance_grid(grid_start, grid_end, multiplicative_grid_scaling_factor)

print(len(variance_grid))
# Simulate data
num_samples = 10000
num_features = 3000
residual_variance = 0.2
np.random.seed(1)
Y, X, betas_sim, z_sim, delta_sim = simulate_data(num_samples, num_features, variance_grid, residual_variance=residual_variance)


print(delta_sim)

#Z, delta_alpha, beta_mu, beta_var = multivariate_inference(Y, X, variance_grid)
Z, delta_alpha, beta_mu, beta_var = inference(Y, X, variance_grid, max_iters=3000)

#Z, delta_alpha, beta_mu, beta_var = inference_v2(Y, X, variance_grid, max_iters=3000)

pdb.set_trace()