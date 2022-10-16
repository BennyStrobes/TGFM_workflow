import numpy as np 
import os
import sys
import pdb




def update_alpha(alpha_mu, alpha_var, residual, gene_eqtl_pmces, srs_inv, num_genes, s_inv_2_diag, precomputed_a_terms, expected_gamma_alpha):
	# Hack to deal with case where window has no genes
	if len(alpha_mu) == 0:
		return alpha_mu, alpha_var, residual

	# Include effect of the gene corresponding to 0 from the residaul
	gene_trait_pred = gene_eqtl_pmces[0]*alpha_mu[0]
	residual = residual + np.dot(srs_inv, gene_trait_pred)

	for g_index in range(num_genes):
			
		# Calculate terms involved in update	
		b_term = np.dot(np.multiply(residual, s_inv_2_diag), gene_eqtl_pmces[g_index])
		a_term = precomputed_a_terms[g_index] - .5*expected_gamma_alpha[g_index]

		# VI Updates
		alpha_var[g_index] = -1.0/(2.0*a_term)
		alpha_mu[g_index] = b_term*alpha_var[g_index]

		# Update resid for next round (after this resid includes effects of all genes)
		if g_index == (num_genes - 1):
			gene_trait_pred = -gene_eqtl_pmces[g_index]*alpha_mu[g_index]
		else:
			gene_trait_pred = gene_eqtl_pmces[(g_index+1)]*alpha_mu[(g_index+1)] - gene_eqtl_pmces[g_index]*alpha_mu[g_index]
		residual = residual + np.dot(srs_inv, gene_trait_pred)

	return alpha_mu, alpha_var, residual


def update_beta(beta_mu, beta_var, residual, srs_inv, num_variants, D_diag, s_inv_2_diag, expected_gamma_beta):
	precomputed_a_terms = (-.5*D_diag) - (.5*expected_gamma_beta)

	# Remove the pleiotropic effect of the variant corresponding to 0 from the residaul
	residual = residual + srs_inv[:,0]*beta_mu[0]

	# Update each snp in parallel
	for k_index in range(num_variants):

		# Calculate terms involved in update
		b_term = residual[k_index]*s_inv_2_diag[k_index]

		# VI Updates
		beta_var[k_index] = -1.0/(2.0*precomputed_a_terms[k_index])
		beta_mu[k_index] = b_term*beta_var[k_index]

		# Update resid for next round (after this resid includes effects of all genes)
		if k_index == (num_variants - 1):
			residual = residual - srs_inv[:,k_index]*beta_mu[k_index]
		else:
			residual = residual + srs_inv[:,(k_index+1)]*beta_mu[(k_index+1)]- srs_inv[:,k_index]*beta_mu[k_index]

	return beta_mu, beta_var, residual

