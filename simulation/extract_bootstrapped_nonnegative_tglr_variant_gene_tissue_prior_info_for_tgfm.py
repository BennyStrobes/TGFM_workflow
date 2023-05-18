import numpy as np 
import os
import sys
import pdb
import pickle






tgfm_input_summary_file = sys.argv[1]
eqtl_sample_size = sys.argv[2]
simulation_name_string = sys.argv[3]
simulated_sldsc_results_dir = sys.argv[4]
output_file = sys.argv[5]


# GET M-vec
tglr_tau_file = simulated_sldsc_results_dir + simulation_name_string + '_eqtl_ss_' + eqtl_sample_size + '_susie_pmces_sldsc_results_organized_res.txt'
tglr_h2_file = simulated_sldsc_results_dir + simulation_name_string + '_eqtl_ss_' + eqtl_sample_size + '_susie_pmces_sldsc_results_organized_mediated_h2.txt'
tglr_tau_data = np.loadtxt(tglr_tau_file, dtype=str,delimiter='\t')
tglr_h2_data = np.loadtxt(tglr_h2_file, dtype=str, delimiter='\t')

# Extract taus and heritabilities
anno_names = tglr_tau_data[1:,0]
taus = tglr_tau_data[1:,1].astype(float)
h2s = tglr_h2_data[1:,1].astype(float)
# Extract m-vec
m_vec = h2s/taus
n_snps = m_vec[0]



# Get non-negative bootstrapped taus
bootstrapped_nonnegative_tglr_file = simulated_sldsc_results_dir + simulation_name_string + '_eqtl_ss_' + eqtl_sample_size + '_susie_pmces_sldsc_results_nonnegative_eqtl_bootstrapped_sldsc_coefficients.txt'
bs_taus = []
anno_names = []
f = open(bootstrapped_nonnegative_tglr_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	anno_names.append(data[0])
	bs_taus.append(np.asarray(data[4].split(';')).astype(float))
f.close()
anno_names = np.asarray(anno_names)
bs_taus = np.asarray(bs_taus)

# Eqtl start index
eqtl_start_index = np.where(anno_names=='tissue0_0')[0][0]

bs_h2s = np.transpose(bs_taus)*m_vec

bs_per_snp_h2 = np.sum(bs_h2s[:,:eqtl_start_index],axis=1)/n_snps
bs_per_gene_tissue_h2 = bs_taus[eqtl_start_index:,:]

t = open(output_file,'w')
t.write('element_name\tprior\texp_E_ln_prior\tprior_distribution\n')
t.write('variant\t' + str(np.mean(bs_per_snp_h2)) + '\t' + str(np.exp(np.mean(np.log(bs_per_snp_h2)))) + '\t' + ';'.join(bs_per_snp_h2.astype(str)) + '\n')
for tiss_iter in range(10):
	t.write('tissue' + str(tiss_iter) + '\t' + str(np.mean(bs_per_gene_tissue_h2[tiss_iter,:])) + '\t' + str(np.exp(np.mean(np.log(bs_per_gene_tissue_h2[tiss_iter,:])))) + '\t' + ';'.join(bs_per_gene_tissue_h2[tiss_iter,:].astype(str)) + '\n')
t.close()



