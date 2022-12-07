import numpy as np 
import os
import sys
import pdb



def get_number_of_snps_per_annotation(ldsc_annotation_dir, suffix):
	for chrom_num in range(1,23):
		m_file = ldsc_annotation_dir + str(chrom_num) + suffix
		tmp = np.loadtxt(m_file)
		if chrom_num == 1:
			m_vec = np.copy(tmp)
		else:
			m_vec = m_vec + tmp
	return m_vec


def get_sldsc_tau(ldsc_res_file):
	aa = np.loadtxt(ldsc_res_file, dtype=str, delimiter='\t')
	return aa[1:,-3].astype(float)

ldsc_annotation_dir = sys.argv[1]
sldsc_results_root = sys.argv[2]



num_snps_m = get_number_of_snps_per_annotation(ldsc_annotation_dir, ".l2.M")
num_snps_m_5_50 = get_number_of_snps_per_annotation(ldsc_annotation_dir, ".l2.M_5_50")

ldsc_res_file = sldsc_results_root + '.results'
sldsc_tau = get_sldsc_tau(ldsc_res_file)


pdb.set_trace()

