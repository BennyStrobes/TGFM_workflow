import numpy as np 
import os
import sys
import pdb
import scipy.stats
import pickle

def get_scaler():
	dir_name = '/n/scratch3/users/b/bes710/causal_eqtl_gwas/simulation/simulated_ld_scores/'
	m_5_50_file = dir_name + 'simulation_9_chrom21_cis_window_100000_baseline.21.l2.M_5_50'
	m_file = dir_name + 'simulation_9_chrom21_cis_window_100000_baseline.21.l2.M'

	m_5_50_vec = np.loadtxt(m_5_50_file)
	m_vec = np.loadtxt(m_file)

	scaler = m_vec/m_5_50_vec
	return np.hstack((scaler, np.ones(10)))


def create_file_containing_total_h2_across_simulation_runs(simulated_sldsc_results_dir, global_simulation_name_string, total_heritability, eqtl_sample_sizes, simulation_runs, organized_total_h2_output_file):
	# Open output file and print header
	t = open(organized_total_h2_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\testimated_h2\tsimulated_h2\n')


	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_distr' + '_sldsc_results_organized_mediated_h2.txt'
			# Load in per anno med h2
			h2_data_raw = np.loadtxt(h2_file, dtype=str, delimiter='\t')

			# Extract total estimated heritability
			total_h2_est = np.sum(h2_data_raw[1:,1].astype(float))
			arr.append(total_h2_est)

			# print to output file
			t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(total_h2_est) + '\t' + str(total_heritability) + '\n')
		print(np.mean(arr))
	# Close file handle
	t.close()

	return

def create_file_containing_fraction_expression_mediated_h2_across_simulation_runs(simulated_sldsc_results_dir, global_simulation_name_string, fraction_expression_mediated_heritability, eqtl_sample_sizes, simulation_runs, organized_fraction_h2_output_file):
	# Open output file and print header
	t = open(organized_fraction_h2_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\testimated_fraction_expression_med_h2\tsimulated_expression_mediated_h2\n')	

	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			frac_h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_distr' + '_sldsc_results_h2_med.txt'
			frac_h2_raw = np.loadtxt(frac_h2_file,dtype=str)
			frac_h2 = frac_h2_raw[1,0]

			arr.append(float(frac_h2))

			# print to output file
			t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(frac_h2) + '\t' + str(fraction_expression_mediated_heritability) + '\n')
		print(np.mean(arr))
	# Close file handle
	t.close()
	return

def create_file_containing_mediated_nonnegative_tissue_taus_in_causal_and_non_causal_tissues(simulated_sldsc_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_nonnegative_taus_per_tissue_output_file, eqtl_type, anno_type):
	# Open output file and print header
	t = open(mediated_nonnegative_taus_per_tissue_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\ttissue_number\testimated_tau\tsimulated_tau\n')	

	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			if anno_type == 'genotype_intercept':
				frac_h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + eqtl_type + '_sldsc_results_' + anno_type + '_nonnegative_eqtl_bootstrapped_sldsc_coefficients.txt'
			else:
				frac_h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + eqtl_type + '_sldsc_results_nonnegative_eqtl_bootstrapped_sldsc_coefficients.txt'
			frac_h2_raw = np.loadtxt(frac_h2_file,dtype=str)
			taus = frac_h2_raw[-10:,1]

			for tissue_iter in range(10):

				# print to output file
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(tissue_iter) + '\t' + str(taus[tissue_iter]) + '\t' + '.0001' + '\n')
	# Close file handle
	t.close()
	return



def create_file_containing_mediated_h2_sparse_binary_est_in_causal_and_non_causal_tissues(simulated_sldsc_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_binary_variable_by_tissue_output_file):
	# Open output file and print header
	t = open(mediated_binary_variable_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\ttissue_number\tcoef_est\tpvalue\tcausal_status\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_distr' '_sldsc_results_organized_0.5_sparse_ard_eqtl_coefficients_mv_update_res.txt'
			# Load in per anno med h2
			h2_data_raw = np.loadtxt(h2_file, dtype=str, delimiter='\t')
			med_h2_coef = h2_data_raw[-10:,1].astype(float)
			pvalue = (med_h2_coef < 1e-13)*1.0
			for tissue_iter in range(10):
				if tissue_iter == 0 or tissue_iter == 3:
					causal_status = 'causal'
				else:
					causal_status = 'null'
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(tissue_iter) + '\t' + str(med_h2_coef[tissue_iter]) + '\t' + str(pvalue[tissue_iter]) + '\t' + causal_status + '\n')
		#print(np.mean(arr))
	t.close()

def create_file_containing_iterative_bootstrapped_sampler_prior_mediated_h2_pvalue_in_causal_and_non_causal_tissues_across_thresholds(simulated_tgfm_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_iterative_sampler_pvalue_by_tissue_output_file, thresholds, model_name):
	# Open output file and print header
	t = open(mediated_iterative_sampler_pvalue_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\ttissue_number\tz_score\tpvalue\tcausal_status\tthreshold\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			h2_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + model_name + '_bootstrapped.txt'
			# Load in per anno med h2
			h2_data_raw = np.loadtxt(h2_file, dtype=str, delimiter='\t')
			expr_data = h2_data_raw[-10:,:]
			for tissue_iter in range(10):
				if tissue_iter == 0 or tissue_iter == 3:
					causal_status = 'causal'
				else:
					causal_status = 'null'
				for threshold in thresholds:
					bootstrapped_coefficients = np.asarray(expr_data[tissue_iter,-1].split(';')).astype(float)
					emperical_pvalue = np.sum(bootstrapped_coefficients < threshold)/len(bootstrapped_coefficients)
					t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(tissue_iter) + '\t' + str(0.0) + '\t' + str(emperical_pvalue) + '\t' + causal_status + '\t' + str(threshold) + '\n')
		#print(np.mean(arr))
	t.close()		


def create_file_containing_iterative_bootstrapped_sampler_prior_mediated_h2_pvalue_in_causal_and_non_causal_tissues(simulated_tgfm_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_iterative_sampler_pvalue_by_tissue_output_file):
	# Open output file and print header
	t = open(mediated_iterative_sampler_pvalue_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\ttissue_number\tz_score\tpvalue\tcausal_status\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			h2_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_pmces_uniform_iterative_variant_gene_prior_bootstrapped.txt'
			# Load in per anno med h2
			h2_data_raw = np.loadtxt(h2_file, dtype=str, delimiter='\t')
			expr_data = h2_data_raw[-10:,:]
			for tissue_iter in range(10):
				if tissue_iter == 0 or tissue_iter == 3:
					causal_status = 'causal'
				else:
					causal_status = 'null'
				bootstrapped_coefficients = np.asarray(expr_data[tissue_iter,-1].split(';')).astype(float)
				emperical_pvalue = np.sum(bootstrapped_coefficients < 1e-9)/len(bootstrapped_coefficients)

				t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(tissue_iter) + '\t' + str(0.0) + '\t' + str(emperical_pvalue) + '\t' + causal_status + '\n')
		#print(np.mean(arr))
	t.close()		

def create_file_containing_nonnegative_bootstrapped_mediated_h2_pvalue_in_causal_and_non_causal_tissues(simulated_sldsc_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_nonnegative_pvalue_by_tissue_output_file, eqtl_type, anno_type):
	# Open output file and print header
	t = open(mediated_nonnegative_pvalue_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\ttissue_number\tz_score\tpvalue\tcausal_status\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			if anno_type == 'full_anno':
				h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + eqtl_type + '_sldsc_results_nonnegative_eqtl_bootstrapped_sldsc_coefficients.txt'
			elif anno_type == 'genotype_intercept':
				h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + eqtl_type + '_sldsc_results_genotype_intercept_nonnegative_eqtl_bootstrapped_sldsc_coefficients.txt'
			# Load in per anno med h2
			h2_data_raw = np.loadtxt(h2_file, dtype=str, delimiter='\t')
			expr_data = h2_data_raw[-10:,:]
			for tissue_iter in range(10):
				if tissue_iter == 0 or tissue_iter == 3:
					causal_status = 'causal'
				else:
					causal_status = 'null'
				bootstrapped_coefficients = np.asarray(expr_data[tissue_iter,-1].split(';')).astype(float)
				emperical_pvalue = np.sum(bootstrapped_coefficients < 1e-10)/len(bootstrapped_coefficients)

				t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(tissue_iter) + '\t' + str(0.0) + '\t' + str(emperical_pvalue) + '\t' + causal_status + '\n')
		#print(np.mean(arr))
	t.close()

def create_file_containing_nonnegative_bootstrapped_mediated_h2_pvalue_in_causal_and_non_causal_tissues_across_thresholds(simulated_sldsc_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_nonnegative_pvalue_by_tissue_output_file, eqtl_type, anno_type, thresholds):
	# Open output file and print header
	t = open(mediated_nonnegative_pvalue_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\ttissue_number\tz_score\tpvalue\tcausal_status\tthreshold\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			if anno_type == 'full_anno':
				h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + eqtl_type + '_sldsc_results_nonnegative_eqtl_bootstrapped_sldsc_coefficients.txt'
			elif anno_type == 'genotype_intercept':
				h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + eqtl_type + '_sldsc_results_genotype_intercept_nonnegative_eqtl_bootstrapped_sldsc_coefficients.txt'
			# Load in per anno med h2
			h2_data_raw = np.loadtxt(h2_file, dtype=str, delimiter='\t')
			expr_data = h2_data_raw[-10:,:]
			for tissue_iter in range(10):
				if tissue_iter == 0 or tissue_iter == 3:
					causal_status = 'causal'
				else:
					causal_status = 'null'
				bootstrapped_coefficients = np.asarray(expr_data[tissue_iter,-1].split(';')).astype(float)
				for threshold in thresholds:
					emperical_pvalue = np.sum(bootstrapped_coefficients < threshold)/len(bootstrapped_coefficients)
					t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(tissue_iter) + '\t' + str(0.0) + '\t' + str(emperical_pvalue) + '\t' + causal_status + '\t' + str(threshold) + '\n')
		#print(np.mean(arr))
	t.close()	


def create_file_containing_mediated_h2_pvalue_in_causal_and_non_causal_tissues(simulated_sldsc_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_pvalue_by_tissue_output_file, eqtl_type, anno_type):
	# Open output file and print header
	t = open(mediated_pvalue_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\ttissue_number\tz_score\tpvalue\tcausal_status\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			if anno_type == 'full_anno':
				h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + eqtl_type + '_sldsc_results_organized_res.txt'
			elif anno_type == 'genotype_intercept':
				h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + eqtl_type + '_sldsc_results_genotype_intercept_organized_res.txt'
			# Load in per anno med h2
			h2_data_raw = np.loadtxt(h2_file, dtype=str, delimiter='\t')
			med_h2_z = h2_data_raw[-10:,3].astype(float)
			pvalue = scipy.stats.norm.sf(abs(med_h2_z))*2
			for tissue_iter in range(10):
				if tissue_iter == 0 or tissue_iter == 3:
					causal_status = 'causal'
				else:
					causal_status = 'null'
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(tissue_iter) + '\t' + str(med_h2_z[tissue_iter]) + '\t' + str(pvalue[tissue_iter]) + '\t' + causal_status + '\n')
		#print(np.mean(arr))
	t.close()

def create_file_containing_expected_fraction_of_expression_mediated_disease_components(simulated_tgfm_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, model_name, simulated_ld_scores_dir, fraction_expr_med_disease_components_output_file):
	# Open output file and print header
	t = open(fraction_expr_med_disease_components_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\testimated_fraction_expression_mediated\tsimulated_fraction_expression_mediated\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			h2_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + model_name + '_bootstrapped.txt'
			# Load in per anno med h2
			h2_data_raw = np.loadtxt(h2_file, dtype=str, delimiter='\t')
			expr_med_tau = h2_data_raw[-10:,1].astype(float)
			non_mediated_tau = float(h2_data_raw[1,1])

			# Load in M-vec
			m_vec_file = simulated_ld_scores_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_joint_baseline_variant_' + str(eqtl_sample_size) + '_gene_ld_scores.l2.M'
			m_vec = np.loadtxt(m_vec_file)
			n_var = m_vec[0]
			n_genes = m_vec[-10:]

			n_var_comp = n_var*non_mediated_tau
			n_gene_comp = np.sum(n_genes*expr_med_tau)

			fraction_mediated = n_gene_comp/(n_var_comp + n_gene_comp)

			t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(fraction_mediated) + '\t' +  '.1' + '\n')
	t.close()
	return

def create_file_containing_fraction_causal_in_causal_and_non_causal_tissues(simulated_tgfm_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, model_name, fraction_causal_by_tissue_output_file):
	# Open output file and print header
	t = open(fraction_causal_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\ttissue_number\testimated_fraction_causal\tsimulated_fraction_causal\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			h2_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + model_name + '_bootstrapped.txt'
			# Load in per anno med h2
			h2_data_raw = np.loadtxt(h2_file, dtype=str, delimiter='\t')
			med_tau = h2_data_raw[-10:,1].astype(float)

			for tissue_iter in range(10):
				if tissue_iter == 0 or tissue_iter == 3:
					sim_med_h2 = 30/2000.0
				else:
					sim_med_h2 = 0.0
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(tissue_iter) + '\t' + str(med_tau[tissue_iter]) + '\t' + str(sim_med_h2) + '\n')
		#print(np.mean(arr))
	t.close()		


def create_file_containing_mediated_h2_in_causal_and_non_causal_tissues_sparse_estimate(simulated_sldsc_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_h2_by_tissue_output_file):
	# Open output file and print header
	t = open(mediated_h2_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\ttissue_number\testimated_med_h2\tsimulated_mediated_h2\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_distr' + '_sldsc_results_organized_0.5_sparse_ard_eqtl_coefficients_mv_update_res.txt'
			# Load in per anno med h2
			h2_data_raw = np.loadtxt(h2_file, dtype=str, delimiter='\t')
			med_tau = h2_data_raw[-10:,1].astype(float)

			gene_m_file = '/n/scratch3/users/b/bes710/causal_eqtl_gwas/simulation/simulated_ld_scores/' + 'simulation_9_chrom1_cis_window_100000_gene_weighted_ld_scores_eqtlss_' + str(eqtl_sample_size) + '_M.txt'
			m_vec = np.loadtxt(gene_m_file)

			med_h2 = m_vec*med_tau

			arr.append(np.sum(med_h2))
			for tissue_iter in range(10):
				if tissue_iter == 0 or tissue_iter == 3:
					sim_med_h2 = .3*.1*.5
				else:
					sim_med_h2 = 0.0
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(tissue_iter) + '\t' + str(med_h2[tissue_iter]) + '\t' + str(sim_med_h2) + '\n')
		#print(np.mean(arr))
	t.close()	


def create_file_containing_mediated_h2_in_causal_and_non_causal_tissues(simulated_sldsc_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_h2_by_tissue_output_file):
	# Open output file and print header
	t = open(mediated_h2_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\ttissue_number\testimated_med_h2\tsimulated_mediated_h2\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_distr' + '_sldsc_results_organized_mediated_h2.txt'
			# Load in per anno med h2
			h2_data_raw = np.loadtxt(h2_file, dtype=str, delimiter='\t')
			med_h2 = h2_data_raw[-10:,1].astype(float)
			arr.append(np.sum(med_h2))
			for tissue_iter in range(10):
				if tissue_iter == 0 or tissue_iter == 3:
					sim_med_h2 = .3*.1*.5
				else:
					sim_med_h2 = 0.0
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(tissue_iter) + '\t' + str(med_h2[tissue_iter]) + '\t' + str(sim_med_h2) + '\n')
		#print(np.mean(arr))
	t.close()
	'''
	eqtl_sample_size = 1000
	aa = np.loadtxt(mediated_h2_by_tissue_output_file,dtype=str,delimiter='\t')
	bb = aa[np.where(aa[:,0]==str(eqtl_sample_size))[0],:]
	for tissue_iter in range(10):
		cc = bb[np.where(bb[:,2] == str(tissue_iter))[0],:]
		print(np.mean(cc[:,3].astype(float)))
	'''

def create_file_containing_mediated_h2_power_to_detect_causal_tissues(mediated_pvalue_by_tissue_output_file, power_output_file, eqtl_sample_sizes):
	raw_pvalue_file = np.loadtxt(mediated_pvalue_by_tissue_output_file,dtype=str, delimiter='\t')

	t = open(power_output_file,'w')
	t.write('eqtl_sample_size\tpower\tpower_lb\tpower_ub\n')

	for eqtl_sample_size in eqtl_sample_sizes:
		raw_pvalue_subset = raw_pvalue_file[raw_pvalue_file[:,0] == str(eqtl_sample_size), :]
		raw_pvalue_subset = raw_pvalue_subset[raw_pvalue_subset[:,-1] == 'causal',:]
		causal_pvalues = raw_pvalue_subset[:,4].astype(float)
		power = np.sum(causal_pvalues < .05)/len(causal_pvalues)
		se = np.sqrt((power)*(1.0-power))/np.sqrt(len(causal_pvalues))
		t.write(str(eqtl_sample_size) + '\t' + str(power) + '\t' + str(power-(se*1.96)) + '\t' + str(power+(se*1.96)) + '\n')
	t.close()
	return


def create_file_containing_mediated_h2_power_to_detect_causal_tissues_across_thresholds(mediated_pvalue_by_tissue_output_file, power_output_file, eqtl_sample_sizes, thresholds):
	raw_pvalue_file = np.loadtxt(mediated_pvalue_by_tissue_output_file,dtype=str, delimiter='\t')

	t = open(power_output_file,'w')
	t.write('eqtl_sample_size\tthreshold\tpower\tpower_lb\tpower_ub\n')

	for eqtl_sample_size in eqtl_sample_sizes:
		for threshold in thresholds:
			raw_pvalue_subset = raw_pvalue_file[raw_pvalue_file[:,0] == str(eqtl_sample_size), :]
			raw_pvalue_subset = raw_pvalue_subset[raw_pvalue_subset[:,-2] == 'causal',:]
			raw_pvalue_subset = raw_pvalue_subset[raw_pvalue_subset[:,-1].astype(float)==threshold, :]
			causal_pvalues = raw_pvalue_subset[:,4].astype(float)
			power = np.sum(causal_pvalues < .05)/len(causal_pvalues)
			se = np.sqrt((power)*(1.0-power))/np.sqrt(len(causal_pvalues))
			t.write(str(eqtl_sample_size) + '\t' + str(threshold) + '\t' + str(power) + '\t' + str(power-(se*1.96)) + '\t' + str(power+(se*1.96)) + '\n')
	t.close()
	return


def create_file_containing_mediated_h2_type_1_error_across_thresholds(mediated_pvalue_by_tissue_output_file, type_1_error_output_file, eqtl_sample_sizes, thresholds):
	raw_pvalue_file = np.loadtxt(mediated_pvalue_by_tissue_output_file,dtype=str, delimiter='\t')

	t = open(type_1_error_output_file,'w')
	t.write('eqtl_sample_size\tthreshold\ttype_1_error\ttype_1_error_lb\ttype_1_error_ub\n')

	for eqtl_sample_size in eqtl_sample_sizes:
		for threshold in thresholds:
			raw_pvalue_subset = raw_pvalue_file[raw_pvalue_file[:,0] == str(eqtl_sample_size), :]
			raw_pvalue_subset = raw_pvalue_subset[raw_pvalue_subset[:,-2] == 'null',:]
			raw_pvalue_subset = raw_pvalue_subset[raw_pvalue_subset[:,-1].astype(float)==threshold, :]
			null_pvalues = raw_pvalue_subset[:,4].astype(float)
			type_1_error = np.sum(null_pvalues < .05)/len(null_pvalues)
			se = np.sqrt((type_1_error)*(1.0-type_1_error))/np.sqrt(len(null_pvalues))

			t.write(str(eqtl_sample_size) + '\t' + str(threshold) + '\t' + str(type_1_error) + '\t' + str(type_1_error-(se*1.96)) + '\t' + str(type_1_error+(se*1.96)) + '\n')
	t.close()
	return


def create_file_containing_mediated_h2_type_1_error(mediated_pvalue_by_tissue_output_file, type_1_error_output_file, eqtl_sample_sizes):
	raw_pvalue_file = np.loadtxt(mediated_pvalue_by_tissue_output_file,dtype=str, delimiter='\t')

	t = open(type_1_error_output_file,'w')
	t.write('eqtl_sample_size\ttype_1_error\ttype_1_error_lb\ttype_1_error_ub\n')

	for eqtl_sample_size in eqtl_sample_sizes:
		raw_pvalue_subset = raw_pvalue_file[raw_pvalue_file[:,0] == str(eqtl_sample_size), :]
		raw_pvalue_subset = raw_pvalue_subset[raw_pvalue_subset[:,-1] == 'null',:]
		null_pvalues = raw_pvalue_subset[:,4].astype(float)
		type_1_error = np.sum(null_pvalues < .05)/len(null_pvalues)
		se = np.sqrt((type_1_error)*(1.0-type_1_error))/np.sqrt(len(null_pvalues))

		t.write(str(eqtl_sample_size) + '\t' + str(type_1_error) + '\t' + str(type_1_error-(se*1.96)) + '\t' + str(type_1_error+(se*1.96)) + '\n')
	t.close()
	return

def create_file_containing_avg_total_h2_across_simulation_runs(organized_total_h2_output_file, organized_avg_total_h2_output_file, eqtl_sample_sizes):
	raw_file = np.loadtxt(organized_total_h2_output_file, dtype=str,delimiter='\t')

	t = open(organized_avg_total_h2_output_file,'w')
	t.write('eqtl_sample_size\ttotal_h2\ttotal_h2_lb\ttotal_h2_ub\n')

	for eqtl_sample_size in eqtl_sample_sizes:
		raw_subset = raw_file[raw_file[:,0] == str(eqtl_sample_size), :]

		est = raw_subset[:,2].astype(float)

		meaner = np.mean(est)
		se = np.std(est)/np.sqrt(len(est))

		t.write(str(eqtl_sample_size) + '\t' + str(meaner) + '\t' + str(meaner - (1.96*se)) + '\t' + str(meaner + (1.96*se)) + '\n')

	t.close()
	return

def create_file_containing_avg_fraction_h2_across_simulation_runs(organized_fraction_h2_output_file, organized_avg_fraction_h2_output_file, eqtl_sample_sizes):
	raw_file = np.loadtxt(organized_fraction_h2_output_file, dtype=str,delimiter='\t')

	t = open(organized_avg_fraction_h2_output_file,'w')
	t.write('eqtl_sample_size\tfrac_h2\tfrac_h2_lb\tfrac_h2_ub\n')

	for eqtl_sample_size in eqtl_sample_sizes:
		raw_subset = raw_file[raw_file[:,0] == str(eqtl_sample_size), :]

		est = raw_subset[:,2].astype(float)

		meaner = np.mean(est)
		se = np.std(est)/np.sqrt(len(est))

		t.write(str(eqtl_sample_size) + '\t' + str(meaner) + '\t' + str(meaner - (1.96*se)) + '\t' + str(meaner + (1.96*se)) + '\n')

	t.close()
	return

def create_file_containing_avg_fraction_causal_by_tissue_across_simulation_runs(mediated_h2_by_tissue_output_file, organized_avg_mediated_h2_by_tissue_output_file, eqtl_sample_sizes):
	raw_file = np.loadtxt(mediated_h2_by_tissue_output_file, dtype=str,delimiter='\t')

	t = open(organized_avg_mediated_h2_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\ttissue_number\tcausal_status\tfraction_causal\tfraction_causal_lb\tfraction_causal_ub\n')

	for eqtl_sample_size in eqtl_sample_sizes:
		raw_subset = raw_file[raw_file[:,0] == str(eqtl_sample_size), :]
		for tiss_iter in range(10):
			if tiss_iter == 0 or tiss_iter == 3:
				causal_status = 'causal'
			else:
				causal_status = 'null'
			raw_subset2 = raw_subset[raw_subset[:,2] == str(tiss_iter),:]

			est = raw_subset2[:,3].astype(float)

			meaner = np.mean(est)
			se = np.std(est)/np.sqrt(len(est))

			t.write(str(eqtl_sample_size) + '\t' + str(tiss_iter) + '\t' + causal_status + '\t' + str(meaner) + '\t' + str(meaner - (1.96*se)) + '\t' + str(meaner + (1.96*se)) + '\n')
	t.close()
	return
def create_file_containing_avg_med_taus_by_tissue_across_simulation_runs(mediated_h2_by_tissue_output_file, organized_avg_mediated_h2_by_tissue_output_file, eqtl_sample_sizes):
	raw_file = np.loadtxt(mediated_h2_by_tissue_output_file, dtype=str,delimiter='\t')

	t = open(organized_avg_mediated_h2_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\ttissue_number\tcausal_status\ttau\ttau_lb\ttau_ub\n')

	for eqtl_sample_size in eqtl_sample_sizes:
		raw_subset = raw_file[raw_file[:,0] == str(eqtl_sample_size), :]
		for tiss_iter in range(10):
			if tiss_iter == 0 or tiss_iter == 3:
				causal_status = 'causal'
			else:
				causal_status = 'null'
			raw_subset2 = raw_subset[raw_subset[:,2] == str(tiss_iter),:]

			est = raw_subset2[:,3].astype(float)

			meaner = np.mean(est)
			se = np.std(est)/np.sqrt(len(est))

			t.write(str(eqtl_sample_size) + '\t' + str(tiss_iter) + '\t' + causal_status + '\t' + str(meaner) + '\t' + str(meaner - (1.96*se)) + '\t' + str(meaner + (1.96*se)) + '\n')
	t.close()
	return


def create_file_containing_avg_med_h2_by_tissue_across_simulation_runs(mediated_h2_by_tissue_output_file, organized_avg_mediated_h2_by_tissue_output_file, eqtl_sample_sizes):
	raw_file = np.loadtxt(mediated_h2_by_tissue_output_file, dtype=str,delimiter='\t')

	t = open(organized_avg_mediated_h2_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\ttissue_number\tcausal_status\tmed_h2\tmed_h2_lb\tmed_h2_ub\n')

	for eqtl_sample_size in eqtl_sample_sizes:
		raw_subset = raw_file[raw_file[:,0] == str(eqtl_sample_size), :]
		for tiss_iter in range(10):
			if tiss_iter == 0 or tiss_iter == 3:
				causal_status = 'causal'
			else:
				causal_status = 'null'
			raw_subset2 = raw_subset[raw_subset[:,2] == str(tiss_iter),:]

			est = raw_subset2[:,3].astype(float)

			meaner = np.mean(est)
			se = np.std(est)/np.sqrt(len(est))

			t.write(str(eqtl_sample_size) + '\t' + str(tiss_iter) + '\t' + causal_status + '\t' + str(meaner) + '\t' + str(meaner - (1.96*se)) + '\t' + str(meaner + (1.96*se)) + '\n')
	t.close()
	return

def extract_successfully_completed_simulation_runs(simulated_tgfm_results_dir, global_simulation_name_string):
	valid_sims = []
	for sim_iter in range(300):
		valid_sim1 = False 
		valid_sim2 = False 
		for file_name in os.listdir(simulated_tgfm_results_dir):
			if file_name.startswith('simulation_' + str(sim_iter) + '_' + global_simulation_name_string + '_eqtl_ss_100_ln_pi_distribution_estimate_1e-10') and file_name.endswith('_results.pkl'):
				valid_sim1 = True
			if file_name.startswith('simulation_' + str(sim_iter) + '_' + global_simulation_name_string + '_eqtl_ss_1000_ln_pi_shared_variant_distribution_estimate_1e-30') and file_name.endswith('_results.pkl'):
				valid_sim2 = True
		if valid_sim1 and valid_sim2:
			valid_sims.append(sim_iter)
	return np.asarray(valid_sims)

def extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file):
	causal_genes = {}
	causal_variants = {}
	causal_genetic_elements = {}

	# Load in causal non-mediated variants
	f = open(causal_variant_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		rsid = data[0]
		causal_effect_size = float(data[1])
		if causal_effect_size != 0:
			# Quick error check
			if rsid in causal_genetic_elements:
				print('assumption eroror')
				pdb.set_trace()
			causal_variants[rsid] = 1
			causal_genetic_elements[rsid] = 1
	f.close()
	
	# Load in causal genes
	f = open(causal_gene_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		ensamble_id = data[0]
		causal_effect_sizes = np.asarray(data[1:]).astype(float)
		for tiss_iter, causal_effect_size in enumerate(causal_effect_sizes):
			if causal_effect_size != 0.0:
				gene_tissue_name = ensamble_id + '_tissue' + str(tiss_iter)
				# Quick error check
				if gene_tissue_name in causal_genetic_elements:
					print('assumption eroror')
					pdb.set_trace()

				causal_genes[gene_tissue_name] = 1
				causal_genetic_elements[gene_tissue_name] = 1
	f.close()

	return causal_genetic_elements, causal_variants, causal_genes

def check_if_at_least_one_causal_genetic_element_is_in_cs(cs_genetic_elements, causal_genetic_elements):
	booleaner = 0.0
	for cs_genetic_element in cs_genetic_elements:
		if cs_genetic_element in causal_genetic_elements:
			booleaner = 1.0
	return booleaner


def extract_all_variants_and_their_positions(bim_file):
	all_variants = []
	var_pos = []
	f = open(bim_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		rsid = data[1]
		pos = int(data[3])
		all_variants.append(rsid)
		var_pos.append(pos)

	f.close()
	return np.asarray(all_variants), np.asarray(var_pos)

def extract_all_genes_and_their_positions(simulated_gene_position_file):
	f = open(simulated_gene_position_file)
	all_genes = []
	all_pos = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		all_genes.append(data[1])
		all_pos.append(int(data[2]))
	f.close()
	return np.asarray(all_genes), np.asarray(all_pos)

# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
def create_file_containing_tgfm_cs_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, initialization_versions):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	# Open output file
	t = open(cs_power_per_component_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\tinitialization_version\tsimulation_number\twindow_name\tgenetic_element_class\tgenetic_element_name\tcausal_variant_in_cs\n')
	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		for eqtl_sample_size in eqtl_sample_sizes:
			# Create mapping from ln_pi_method to dictionary of elements in 95% cs for this simulation
			discovered_dicti = {}
			for ln_pi_method in ln_pi_methods:
				for initialization_version in initialization_versions:
					discovered_pi_dicti = {}
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_distr' + '_ln_pi_' + ln_pi_method + '_init_' + initialization_version + '_tgfm_component_cs_summary.txt'
					f = open(cs_file)
					head_count = 0
					for line in f:
						line = line.rstrip()
						data = line.split('\t')
						if head_count == 0:
							head_count = head_count + 1
							continue
						cs_elements = data[3].split(';')
						for cs_element in cs_elements:
							discovered_pi_dicti[cs_element] = 1
					f.close()
					discovered_dicti[ln_pi_method + '_' + initialization_version] = discovered_pi_dicti

			# TGFM window file
			# File containing which windows we ran TGFM on
			tgfm_window_file = simulated_tgfm_input_data_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_tgfm_input_data_summary.txt'
			head_count = 0
			f = open(tgfm_window_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				# Extract relevent fields
				window_name = data[0]
				# Need to extract middle_genes and middle_variants in the middle of this window
				# Extract middle start and middle end
				window_start = int(window_name.split('_')[1])
				window_middle_start = window_start + 1000000
				window_middle_end = window_start + 2000000
				window_end = window_start + 3000000
				if window_end != int(window_name.split('_')[2]):
					print('assumption erororo')
					pdb.set_trace()

				variant_indices = (all_variants_positions >= window_middle_start) & (all_variants_positions < window_middle_end)
				middle_variants = all_variants[variant_indices]
				gene_indices = (all_genes_positions >= window_middle_start) & (all_genes_positions < window_middle_end)
				middle_genes = all_genes[gene_indices]


				for gene_name_stem in middle_genes:
					for tissue_iter in range(10):
						gene_name = gene_name_stem + '_tissue' + str(tissue_iter)
						if gene_name not in causal_genetic_elements:
							continue
						# THis is a causal gene
						for ln_pi_method in ln_pi_methods:
							for initialization_version in initialization_versions:
								booler = 0.0
								if gene_name in discovered_dicti[ln_pi_method + '_' + initialization_version]:
									booler = 1.0
								t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + initialization_version + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gene_name + '\t' + str(booler) + '\n')
				for variant_name in middle_variants:
					if variant_name not in causal_genetic_elements:
						continue
					# This is a causal variant
					for ln_pi_method in ln_pi_methods:
						for initialization_version in initialization_versions:
							booler = 0.0
							if variant_name in discovered_dicti[ln_pi_method + '_' + initialization_version]:
								booler = 1.0
							t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + initialization_version + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'variant\t' + variant_name + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return

def create_mapping_from_gene_to_tss(simulated_gene_position_file):
	f = open(simulated_gene_position_file)
	dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		dicti[data[1]] = float(data[2])
	f.close()
	return dicti

def extract_ld_file_and_tgfm_pickle_corresponding_to_a_gene_based_on_its_tss(gene_tss, tgfm_input_data_file):
	hitbool = False
	head_count = 0
	f = open(tgfm_input_data_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1 
			continue
		window_start = int(data[0].split('_')[1])
		window_middle_start = window_start + 1000000
		window_middle_end = window_start + 2000000
		if gene_tss >= window_middle_start and gene_tss < window_middle_end:
			if hitbool:
				print('assumption errorr')
				pdb.set_trace()
			hitbool = True
			ld_file = data[1]
			pkl_file = data[2]
	f.close()

	if hitbool == False:
		ld_file = 'null'
		pkl_file = 'null'

	return ld_file, pkl_file

def create_causal_eqtl_window_effects(gene_causal_effects, gene_snp_names, window_variant_names):
	window_effects = np.zeros(len(window_variant_names))
	mapping = {}
	for i,val in enumerate(window_variant_names):
		mapping[val] = i

	for ii,gene_snp_name in enumerate(gene_snp_names):
		window_effects[mapping[gene_snp_name]] = gene_causal_effects[ii]
	return window_effects

def extract_best_tagged_causal_genes(causal_genes_arr, simulation_number, global_simulation_name_string, eqtl_sample_size, simulated_gene_expression_dir, simulated_gene_position_file, simulated_learned_gene_models_dir, simulated_tgfm_input_data_dir):
	# First create mapping from gene to tss
	ensamble_id_to_tss = create_mapping_from_gene_to_tss(simulated_gene_position_file)

	# Filename for tgfm input data
	tgfm_input_data_file = simulated_tgfm_input_data_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_tgfm_input_data_summary.txt'

	# Initialize vector to keep track of best tagged genes
	best_tagged_causal_genes = []
	# Loop through causal genes
	for causal_gene in causal_genes_arr:
		# Tss of gene
		gene_tss = ensamble_id_to_tss[causal_gene.split('_')[0]]

		# Extract ld file and tgfm input pickle corresponding to this gene
		ld_file, tgfm_input_pkl = extract_ld_file_and_tgfm_pickle_corresponding_to_a_gene_based_on_its_tss(gene_tss, tgfm_input_data_file)
		if ld_file == 'null':
			continue

		# Load in pkl file
		g = open(tgfm_input_pkl, "rb")
		tgfm_data = pickle.load(g)
		g.close()
		tagging_gene_names = tgfm_data['genes']
		window_variant_names = tgfm_data['variants']
		tagging_gene_eqtl_effects = tgfm_data['gene_eqtl_pmces']

		# Load in true causal eqtl effects for gene
		gene_ensamble_id = causal_gene.split('_')[0]
		gene_tissue_name = causal_gene.split('_')[1]
		gene_tissue_int = int(gene_tissue_name.split('issue')[1])
		gene_causal_effect_file = simulated_gene_expression_dir + 'simulation_' + str(simulation_number) +'_' + global_simulation_name_string + '_' + gene_ensamble_id + '_causal_eqtl_effects.npy'
		gene_causal_effects_x_tissue = np.load(gene_causal_effect_file)
		gene_causal_effects = gene_causal_effects_x_tissue[:, gene_tissue_int]
		gene_snp_names_file =  simulated_gene_expression_dir + 'simulation_' + str(simulation_number) +'_' + global_simulation_name_string + '_' + gene_ensamble_id + '_cis_snpids.npy'
		gene_snp_names = np.load(gene_snp_names_file, allow_pickle=True)[:,1]

		causal_eqtl_window_effects = create_causal_eqtl_window_effects(gene_causal_effects, gene_snp_names, window_variant_names)
		if np.var(causal_eqtl_window_effects) == 0:
			print('assumption erroro')
			pdb.set_trace()

		ld_mat = np.load(ld_file)

		concatenated_effects = np.vstack((causal_eqtl_window_effects, tagging_gene_eqtl_effects))
		var_cov_mat = np.dot(np.dot(concatenated_effects, ld_mat), np.transpose(concatenated_effects))
		vary = np.diag(var_cov_mat)
		correlations = np.abs(var_cov_mat[0,1:]/np.sqrt(vary[0]*vary[1:]))

		'''
		# Standardize eqtl effects
		causal_effect_var = np.dot(np.dot(causal_eqtl_window_effects, ld_mat), causal_eqtl_window_effects)
		causal_eqtl_window_effects = causal_eqtl_window_effects/np.sqrt(causal_effect_var)

		correlations = []
		for tagging_gene_iter, tagging_gene_name in enumerate(tagging_gene_names):
			tagging_eqtl_effect = tagging_gene_eqtl_effects[tagging_gene_iter,:]
			effect_var = np.dot(np.dot(tagging_eqtl_effect, ld_mat), tagging_eqtl_effect)
			tagging_eqtl_effect = tagging_eqtl_effect/np.sqrt(effect_var)

			correlation = np.dot(np.dot(tagging_eqtl_effect, ld_mat), causal_eqtl_window_effects)
			correlations.append(np.abs(correlation))
		correlations = np.asarray(correlations)
		'''

		best_tagging_gene = tagging_gene_names[np.argmax(correlations)]


		best_tagged_causal_genes.append(best_tagging_gene)


	return np.asarray(best_tagged_causal_genes)


def create_file_containing_tgfm_cs_calibration_given_best_tagged_genetic_element_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, simulated_gene_expression_dir, simulated_gene_position_file, simulated_learned_gene_models_dir, simulated_tgfm_input_data_dir, cs_coverage_per_high_pip_snp_output_file):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		print(simulation_number)
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
			best_tagged_causal_genes = extract_best_tagged_causal_genes(np.asarray([*causal_genes]), simulation_number, global_simulation_name_string, eqtl_sample_size, simulated_gene_expression_dir, simulated_gene_position_file, simulated_learned_gene_models_dir, simulated_tgfm_input_data_dir)
			causal_genetic_elements = {}
			for causal_variant in [*causal_variants]:
				causal_genetic_elements[causal_variant] =1
			for causal_gene in best_tagged_causal_genes:
				causal_genetic_elements[causal_gene] =1
			for ln_pi_method in ln_pi_methods:

				# Credible set file for this run
				cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_ln_pi_' + ln_pi_method + '_tgfm_component_cs_summary.txt'
				# Loop through cs in cs file
				head_count = 0
				f = open(cs_file)
				for line in f:
					line = line.rstrip()
					data = line.split('\t')
					if head_count == 0:
						head_count = head_count + 1
						continue
					# Extract relevent fields
					component_window_name = data[0]
					component_num = data[1]
					all_cs_genetic_elements = data[3].split(';')
					cs_probs = np.asarray(data[4].split(';')).astype(float)
					cs_genetic_elements = []
					for element_iter, cs_prob in enumerate(cs_probs):
						if cs_prob >= pip_threshold:
							cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
					cs_genetic_elements = np.asarray(cs_genetic_elements)
					if len(cs_genetic_elements) == 0:
						continue

					# Quick error check
					if np.sum(cs_probs) < .95:
						print('assumption error')
						pdb.set_trace()
					genetic_element_name = cs_genetic_elements[0]
					if genetic_element_name.startswith('ENSG'):
						class_name = 'gene'
					elif genetic_element_name.startswith('rs'):
						class_name = 'variant'
					else:
						print('assumptino eroror')
						pdb.set_trace()

					# Check if cs contains at least one causal genetic element
					causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs(cs_genetic_elements, causal_genetic_elements)
					t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\n')
				
					# DEBUGGING 
					if class_name == 'gene' and causal_genetic_element_in_cs_boolean == 0.0:
						pdb.set_trace()

				f.close()
		t.flush()
	t.close()
	return

def debugging(pkl_results_file, all_cs_genetic_elements, cs_prob, causal_genes, causal_variants, component_window_name, pkl_input_file):
		ensamble_id_to_tss = create_mapping_from_gene_to_tss(simulated_gene_position_file)
		# Load in pkl file
		g = open(pkl_results_file, "rb")
		tgfm_res = pickle.load(g)
		g.close()
		g = open(pkl_input_file, "rb")
		tgfm_input = pickle.load(g)
		g.close()
		window_start = int(component_window_name.split('_')[1])
		window_middle_start = window_start + 1000000
		window_middle_end = window_start + 2000000
		window_end = int(component_window_name.split('_')[2])
		genes_in_window = []
		for gene_name in [*causal_genes]:
			ensamble_id = gene_name.split('_')[0]
			tss = ensamble_id_to_tss[ensamble_id]
			if tss >= window_start and tss < window_end:
				genes_in_window.append(gene_name)
		variants_in_window = []
		for variant in tgfm_input['variants']:
			if variant in causal_variants:
				variants_in_window.append(variant)
		pdb.set_trace()
def extract_causal_genes_in_window(simulated_gene_position_file, component_window_name, causal_genes):
	ensamble_id_to_tss = create_mapping_from_gene_to_tss(simulated_gene_position_file)
	window_start = int(component_window_name.split('_')[1])
	window_middle_start = window_start + 1000000
	window_middle_end = window_start + 2000000
	window_end = int(component_window_name.split('_')[2])
	genes_in_window = []
	for gene_name in [*causal_genes]:
		ensamble_id = gene_name.split('_')[0]
		tss = ensamble_id_to_tss[ensamble_id]
		if tss >= window_start and tss < window_end:
			genes_in_window.append(gene_name)
	return np.asarray(genes_in_window)

def get_detected_genes_in_window_from_pkl_file(pkl_results_file):
	g = open(pkl_results_file, "rb")
	tgfm_res = pickle.load(g)
	g.close()
	detected_window_genes = tgfm_res['genes']
	genes_dicti = {}
	for gene_name in detected_window_genes:
		genes_dicti[gene_name] = 1
	return genes_dicti





# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
def create_file_containing_tgfm_cs_calibration_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, cs_coverage_per_component_output_file, initialization_versions):
	# Open output file handle
	t = open(cs_coverage_per_component_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\tinitialization_version\tsimulation_number\twindow_name\tcomponent_number\tcausal_genetic_element_in_cs\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
			for ln_pi_method in ln_pi_methods:
				for initialization_version in initialization_versions:

					# Credible set file for this run
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_distr' + '_ln_pi_' + ln_pi_method + '_init_' + initialization_version + '_tgfm_component_cs_summary.txt'
					# Loop through cs in cs file
					head_count = 0
					f = open(cs_file)
					for line in f:
						line = line.rstrip()
						data = line.split('\t')
						if head_count == 0:
							head_count = head_count + 1
							continue
						# Extract relevent fields
						component_window_name = data[0]
						component_num = data[1]
						cs_genetic_elements = data[3].split(';')
						cs_probs = np.asarray(data[4].split(';')).astype(float)
						# Quick error check
						if np.sum(cs_probs) < .95:
							print('assumption error')
							pdb.set_trace()
						# Check if cs contains at least one causal genetic element
						causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs(cs_genetic_elements, causal_genetic_elements)
						t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + initialization_version + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\n')
					f.close()
	t.close()
	return


# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
def create_file_containing_averaged_tgfm_cs_calibration(cs_coverage_per_component_output_file, cs_coverage_output_file, eqtl_sample_sizes, ln_pi_methods):
	per_component_coverage = np.loadtxt(cs_coverage_per_component_output_file, dtype=str)[1:,:]
	t = open(cs_coverage_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\tcoverage\tcoverage_lb\tcoverage_ub\n')

	for eqtl_sample_size in eqtl_sample_sizes:
		for ln_pi_method in ln_pi_methods:
			subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size)) & (per_component_coverage[:,1] == ln_pi_method)
			components_covered = (per_component_coverage[subset_indices,:][:,-1]).astype(float)
			prop = np.sum(components_covered)/len(components_covered)
			prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
			prop_lb = prop - (1.96*prop_se)
			prop_ub = prop + (1.96*prop_se)

			t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\n')

	t.close()


	return 


def extract_dictionary_lists_of_heritable_and_non_heritable_genes(simulation_gene_summary_file):
	# Initialze output dictionaries
	her_genes = {}
	non_her_genes = {}
	used_genes = {}
	gene_names_arr = []
	# Loop through genes
	f = open(simulation_gene_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent fields
		gene_id = data[0]
		gene_names_arr.append(gene_id)
		causal_eqtl_effect_npy_file = data[3]
		# load in causal eqtl effects to extract which tissues (if any) are cis heritable
		causal_eqtl_effect_mat = np.load(causal_eqtl_effect_npy_file)
		causal_tissues = np.var(causal_eqtl_effect_mat,axis=0) != 0
		# Loop through tissues
		for tissue_iter in range(len(causal_tissues)):
			full_gene_name = gene_id + '_tissue' + str(tissue_iter)
			if full_gene_name in used_genes:
				print('assumption erororo')
				pdb.set_trace()
			used_genes[full_gene_name] = 1
			if causal_tissues[tissue_iter]:
				her_genes[full_gene_name] = 1
			else:
				non_her_genes[full_gene_name] = 1
	f.close()
	return her_genes, non_her_genes, np.asarray(gene_names_arr)

def compute_number_of_heritable_and_non_heritable_genes_with_a_model(ordered_gene_names, eqtl_sample_size, simulation_run, global_simulation_name_string, simulated_learned_gene_models_dir, heritable_genes_dicti, not_heritable_genes_dicti):
	n_her_genes = 0.0
	n_non_her_genes = 0.0
	for gene_name in ordered_gene_names:
		learned_gene_pmces_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_sample_size) + '_gene_model_pmces.npy'
		learned_gene_pmces = np.load(learned_gene_pmces_file)
		detected_genes = np.var(learned_gene_pmces,axis=1) != 0.0
		for tiss_iter in range(len(detected_genes)):
			full_gene_name = gene_name + '_tissue' + str(tiss_iter)
			if detected_genes[tiss_iter]:
				if full_gene_name in heritable_genes_dicti:
					n_her_genes = n_her_genes + 1.0
				elif full_gene_name in not_heritable_genes_dicti:
					n_non_her_genes = n_non_her_genes + 1.0
				else:
					print('assumption erororo')
	return n_her_genes, n_non_her_genes

def create_file_containing_number_of_detected_genes(global_simulation_name_string, simulation_runs,eqtl_sample_sizes, simulated_gene_expression_dir, simulated_learned_gene_models_dir, organized_number_detected_genes_output_file):
	# Open output file handle and print header
	t = open(organized_number_detected_genes_output_file,'w')
	t.write('simulation_number\teQTL_sample_size\tn_detected_heritable_genes\tn_heritable_genes\tn_detected_non_heritable_genes\tn_non_heritable_genes\n')
	# Loop through simulation runs
	for simulation_run in simulation_runs:
		print(simulation_run)
		# Within a simulation_run, the true causal heritable genes are the same across eqtl sample sizes
		# Extract dictionary list of heritable_genes and not heritable_genes
		simulation_gene_summary_file = simulated_gene_expression_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_causal_eqtl_effect_summary.txt'
		heritable_genes_dicti, not_heritable_genes_dicti, ordered_gene_names = extract_dictionary_lists_of_heritable_and_non_heritable_genes(simulation_gene_summary_file)

		# Loop through eqtl sample sizes
		for eqtl_sample_size in eqtl_sample_sizes:
			# At this sample size, compute the number of heritable genes with a model and the number of non-heritable genes with a model
			n_heritable_genes_modeled, n_non_heritable_genes_modeled = compute_number_of_heritable_and_non_heritable_genes_with_a_model(ordered_gene_names, eqtl_sample_size, simulation_run, global_simulation_name_string, simulated_learned_gene_models_dir, heritable_genes_dicti, not_heritable_genes_dicti)

			t.write(str(simulation_run) + '\t' + str(eqtl_sample_size) + '\t' + str(n_heritable_genes_modeled) + '\t' + str(len(heritable_genes_dicti)) + '\t' + str(n_non_heritable_genes_modeled) + '\t' + str(len(not_heritable_genes_dicti)) + '\n')
		t.flush()
	t.close()
	return

def create_file_containing_coloc_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_coloc_results_dir, pip_threshold, coloc_cs_coverage_per_high_pip_snp_output_file):
	# Open output file handle
	t = open(coloc_cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tsimulation_number\twindow_name\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
			coloc_results_file = simulated_coloc_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_coloc_results.txt'
			f = open(coloc_results_file)
			used_genes = {}
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				gene_id = data[0]
				region = 'NA'
				pip = float(data[1])
				if gene_id in used_genes:
					print('assumption eroror')
				used_genes[gene_id] = (pip, region)
			f.close()

			for gene_id in [*used_genes]:
				if used_genes[gene_id][0] < pip_threshold:
					continue
				gene_pip = used_genes[gene_id][0]
				region_name = used_genes[gene_id][1]
				if gene_id in causal_genes:
					causal_genetic_element_in_cs_boolean = 1
				else:
					causal_genetic_element_in_cs_boolean = 0
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_number) + '\t' + region_name + '\t' + gene_id + '\t' + 'gene' + '\t' + str(causal_genetic_element_in_cs_boolean) + '\n')
	t.close()
	return	

def create_file_containing_focus_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, pip_threshold, focus_cs_coverage_per_high_pip_snp_output_file):
	# Open output file handle
	t = open(focus_cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tsimulation_number\twindow_name\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
			focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_focus_res.focus.tsv'
			f = open(focus_results_file)
			used_genes = {}
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				gene_id = data[0]
				if gene_id == 'NULL.MODEL':
					continue
				region = data[-1]
				pip = float(data[-3])
				if gene_id in used_genes:
					if used_genes[gene_id][0] < pip:
						used_genes[gene_id] = (pip, region)
				used_genes[gene_id] = (pip, region)
			f.close()

			for gene_id in [*used_genes]:
				if used_genes[gene_id][0] < pip_threshold:
					continue
				gene_pip = used_genes[gene_id][0]
				region_name = used_genes[gene_id][1]
				if gene_id in causal_genes:
					causal_genetic_element_in_cs_boolean = 1
				else:
					causal_genetic_element_in_cs_boolean = 0
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_number) + '\t' + region_name + '\t' + gene_id + '\t' + 'gene' + '\t' + str(causal_genetic_element_in_cs_boolean) + '\n')
	t.close()
	return	

def create_file_containing_tgfm_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
			for ln_pi_method in ln_pi_methods:
				for twas_method in twas_methods:
					if ln_pi_method.endswith('bootstrapped') and twas_method == 'susie_pmces':
						continue
					if ln_pi_method.startswith('tglr_bootstrapped') and twas_method == 'susie_pmces':
						continue
					if ln_pi_method == 'pmces_uniform_iterative_variant_gene_prior_pip_level_pmces' and twas_method == 'susie_pmces':
						continue					
					# Credible set file for this run
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
					f = open(cs_file)
					head_count = 0
					for line in f:
						line = line.rstrip()
						data = line.split('\t')
						if head_count == 0:
							head_count = head_count + 1
							continue
						component_window_name = data[0]
						component_num = 'NA'
						if len(data) < 3:
							continue
						if data[2] == 'NA':
							continue
						all_cs_genetic_elements = data[1].split(';')
						cs_probs = np.asarray(data[2].split(';')).astype(float)
						cs_genetic_elements = []
						for element_iter, cs_prob in enumerate(cs_probs):
							if cs_prob >= pip_threshold:
								cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
						cs_genetic_elements = np.asarray(cs_genetic_elements)
						if len(cs_genetic_elements) == 0:
							continue


						# Extract n_causal_components
						full_window_name = str(simulation_number) + '_' + component_window_name
						#n_causal_elements = full_window_name_to_n_causal_elements[full_window_name]

						for genetic_element_name in cs_genetic_elements:
							if genetic_element_name.startswith('ENSG'):
								class_name = 'gene'
							elif genetic_element_name.startswith('rs'):
								class_name = 'variant'
							else:
								print('assumptino eroror')
								pdb.set_trace()

							# Check if cs contains at least one causal genetic element
							causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)

							t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\n')
					f.close()



	t.close()
	return

def create_file_containing_averaged_focus_high_pip_calibration(focus_cs_coverage_per_high_pip_snp_output_file, focus_cs_high_pip_coverage_output_file, eqtl_sample_sizes):
	per_component_coverage = np.loadtxt(focus_cs_coverage_per_high_pip_snp_output_file, dtype=str)[1:,:]
	t = open(focus_cs_high_pip_coverage_output_file,'w')
	t.write('eQTL_sample_size\tgenetic_element_class\tn_detected_elements\tcoverage\tcoverage_lb\tcoverage_ub\n')

	for eqtl_sample_size in eqtl_sample_sizes:

		subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size))
		components_covered = (per_component_coverage[subset_indices,:][:,-1]).astype(float)
		prop = np.sum(components_covered)/len(components_covered)
		prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
		prop_lb = prop - (1.96*prop_se)
		prop_ub = prop + (1.96*prop_se)
		n_elements = np.sum(subset_indices)
		t.write(str(eqtl_sample_size) + '\t' + 'gene' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\n')

	t.close()
	return


def create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_component_output_file, cs_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods):
	per_component_coverage = np.loadtxt(cs_coverage_per_component_output_file, dtype=str)[1:,:]
	t = open(cs_coverage_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tgenetic_element_class\tn_detected_elements\tcoverage\tcoverage_lb\tcoverage_ub\n')

	for eqtl_sample_size in eqtl_sample_sizes:
		for ln_pi_method in ln_pi_methods:
			for twas_method in twas_methods:
				subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method)
				components_covered = (per_component_coverage[subset_indices,:][:,-1]).astype(float)
				prop = np.sum(components_covered)/len(components_covered)
				prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
				prop_lb = prop - (1.96*prop_se)
				prop_ub = prop + (1.96*prop_se)
				n_elements = np.sum(subset_indices)
				t.write(str(eqtl_sample_size) + '\t' + ln_pi_method +'\t' + twas_method  + '\t' + 'all' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\n')

				subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method) & (per_component_coverage[:, 7] == 'variant')
				components_covered = (per_component_coverage[subset_indices,:][:,-1]).astype(float)
				prop = np.sum(components_covered)/len(components_covered)
				prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
				prop_lb = prop - (1.96*prop_se)
				prop_ub = prop + (1.96*prop_se)
				n_elements = np.sum(subset_indices)
				t.write(str(eqtl_sample_size) + '\t' + ln_pi_method +'\t' + twas_method  + '\t' + 'variant' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\n')

				subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method) & (per_component_coverage[:, 7] == 'gene')
				components_covered = (per_component_coverage[subset_indices,:][:,-1]).astype(float)
				prop = np.sum(components_covered)/len(components_covered)
				prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
				prop_lb = prop - (1.96*prop_se)
				prop_ub = prop + (1.96*prop_se)
				n_elements = np.sum(subset_indices)			
				t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + 'gene' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\n')
	t.close()
	return

def create_file_containing_focus_high_pip_snp_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	# Open output file
	t = open(focus_cs_power_per_component_output_file,'w')
	t.write('eQTL_sample_size\tsimulation_number\twindow_name\tgenetic_element_class\tgenetic_element_name\tcausal_variant_in_cs\n')
	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		for eqtl_sample_size in eqtl_sample_sizes:
			# Create mapping from ln_pi_method to dictionary of elements in 95% cs for this simulation
			discovered_dicti = {}
			
			focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_focus_res.focus.tsv'
			f = open(focus_results_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				gene_id = data[0]
				if gene_id == 'NULL.MODEL':
					continue
				region = data[-1]
				pip = float(data[-3])
				if pip >= pip_threshold:
					discovered_dicti[gene_id] = 1
			f.close()			

			# TGFM window file
			# File containing which windows we ran TGFM on
			head_count = 0
			f = open(global_window_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				# Extract relevent fields
				window_name = data[0]
				# Need to extract middle_genes and middle_variants in the middle of this window
				# Extract middle start and middle end
				window_start = int(window_name.split('_')[1])
				window_middle_start = window_start + 1000000
				window_middle_end = window_start + 2000000
				window_end = window_start + 3000000
				if window_end != int(window_name.split('_')[2]):
					print('assumption erororo')
					pdb.set_trace()

				variant_indices = (all_variants_positions >= window_middle_start) & (all_variants_positions < window_middle_end)
				middle_variants = all_variants[variant_indices]
				gene_indices = (all_genes_positions >= window_middle_start) & (all_genes_positions < window_middle_end)
				middle_genes = all_genes[gene_indices]


				for gene_name_stem in middle_genes:
					for tissue_iter in range(10):
						gene_name = gene_name_stem + '_tissue' + str(tissue_iter)
						if gene_name not in causal_genetic_elements:
							continue
						# THis is a causal gene
						booler = 0.0
						if gene_name in discovered_dicti:
							booler = 1.0
						t.write(str(eqtl_sample_size) + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gene_name + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return


def create_file_containing_coloc_high_pip_snp_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	# Open output file
	t = open(focus_cs_power_per_component_output_file,'w')
	t.write('eQTL_sample_size\tsimulation_number\twindow_name\tgenetic_element_class\tgenetic_element_name\tcausal_variant_in_cs\n')
	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		for eqtl_sample_size in eqtl_sample_sizes:
			# Create mapping from ln_pi_method to dictionary of elements in 95% cs for this simulation
			discovered_dicti = {}
			
			focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_coloc_results.txt'
			f = open(focus_results_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				gene_id = data[0]
				region = 'NA'
				pip = float(data[1])
				if pip >= pip_threshold:
					discovered_dicti[gene_id] = 1
			f.close()			

			# TGFM window file
			# File containing which windows we ran TGFM on
			head_count = 0
			f = open(global_window_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				# Extract relevent fields
				window_name = data[0]
				# Need to extract middle_genes and middle_variants in the middle of this window
				# Extract middle start and middle end
				window_start = int(window_name.split('_')[1])
				window_middle_start = window_start + 1000000
				window_middle_end = window_start + 2000000
				window_end = window_start + 3000000
				if window_end != int(window_name.split('_')[2]):
					print('assumption erororo')
					pdb.set_trace()

				variant_indices = (all_variants_positions >= window_middle_start) & (all_variants_positions < window_middle_end)
				middle_variants = all_variants[variant_indices]
				gene_indices = (all_genes_positions >= window_middle_start) & (all_genes_positions < window_middle_end)
				middle_genes = all_genes[gene_indices]


				for gene_name_stem in middle_genes:
					for tissue_iter in range(10):
						gene_name = gene_name_stem + '_tissue' + str(tissue_iter)
						if gene_name not in causal_genetic_elements:
							continue
						# THis is a causal gene
						booler = 0.0
						if gene_name in discovered_dicti:
							booler = 1.0
						t.write(str(eqtl_sample_size) + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gene_name + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return


def create_file_containing_tgfm_high_pip_snp_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	# Open output file
	t = open(cs_power_per_component_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\twas_method\tsimulation_number\twindow_name\tgenetic_element_class\tgenetic_element_name\tcausal_variant_in_cs\n')
	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		for eqtl_sample_size in eqtl_sample_sizes:
			# Create mapping from ln_pi_method to dictionary of elements in 95% cs for this simulation
			discovered_dicti = {}
			for ln_pi_method in ln_pi_methods:
				for twas_method in twas_methods:
					if ln_pi_method.endswith('bootstrapped') and twas_method == 'susie_pmces':
						continue
					if ln_pi_method.startswith('tglr_bootstrapped') and twas_method == 'susie_pmces':
						continue
					if ln_pi_method == 'pmces_uniform_iterative_variant_gene_prior_pip_level_pmces' and twas_method == 'susie_pmces':
						continue
					discovered_pi_dicti = {}
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
					f = open(cs_file)
					head_count = 0
					for line in f:
						line = line.rstrip()
						data = line.split('\t')
						if head_count == 0:
							head_count = head_count + 1
							continue
						component_window_name = data[0]
						component_num = 'NA'
						if len(data) < 3:
							continue
						if data[2] == 'NA':
							continue
						all_cs_genetic_elements = data[1].split(';')
						cs_probs = np.asarray(data[2].split(';')).astype(float)
						cs_genetic_elements = []
						for element_iter, cs_prob in enumerate(cs_probs):
							if cs_prob >= pip_threshold:
								cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
						cs_genetic_elements = np.asarray(cs_genetic_elements)
						if len(cs_genetic_elements) == 0:
							continue
						for cs_element in cs_genetic_elements:
							discovered_pi_dicti[cs_element] = 1
					f.close()
					discovered_dicti[ln_pi_method + '_' + twas_method] = discovered_pi_dicti

			# TGFM window file
			# File containing which windows we ran TGFM on
			head_count = 0
			f = open(global_window_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				# Extract relevent fields
				window_name = data[0]
				# Need to extract middle_genes and middle_variants in the middle of this window
				# Extract middle start and middle end
				window_start = int(window_name.split('_')[1])
				window_middle_start = window_start + 1000000
				window_middle_end = window_start + 2000000
				window_end = window_start + 3000000
				if window_end != int(window_name.split('_')[2]):
					print('assumption erororo')
					pdb.set_trace()

				variant_indices = (all_variants_positions >= window_middle_start) & (all_variants_positions < window_middle_end)
				middle_variants = all_variants[variant_indices]
				gene_indices = (all_genes_positions >= window_middle_start) & (all_genes_positions < window_middle_end)
				middle_genes = all_genes[gene_indices]


				for gene_name_stem in middle_genes:
					for tissue_iter in range(10):
						gene_name = gene_name_stem + '_tissue' + str(tissue_iter)
						if gene_name not in causal_genetic_elements:
							continue
						# THis is a causal gene
						for ln_pi_method in ln_pi_methods:
							for twas_method in twas_methods:
								if ln_pi_method.endswith('bootstrapped') and twas_method == 'susie_pmces':
									continue
								if ln_pi_method.startswith('tglr_bootstrapped') and twas_method == 'susie_pmces':
									continue
								if ln_pi_method == 'pmces_uniform_iterative_variant_gene_prior_pip_level_pmces' and twas_method == 'susie_pmces':
									continue
								booler = 0.0
								if gene_name in discovered_dicti[ln_pi_method + '_' + twas_method]:
									booler = 1.0
								t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gene_name + '\t' + str(booler) + '\n')
				for variant_name in middle_variants:
					if variant_name not in causal_genetic_elements:
						continue
					# This is a causal variant
					for ln_pi_method in ln_pi_methods:
						for twas_method in twas_methods:
							if ln_pi_method.endswith('bootstrapped') and twas_method == 'susie_pmces':
								continue
							if ln_pi_method.startswith('tglr_bootstrapped') and twas_method == 'susie_pmces':
								continue
							if ln_pi_method == 'pmces_uniform_iterative_variant_gene_prior_pip_level_pmces' and twas_method == 'susie_pmces':
								continue
							booler = 0.0
							if variant_name in discovered_dicti[ln_pi_method + '_' + twas_method]:
								booler = 1.0
							t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'variant\t' + variant_name + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return

def create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods):
	per_component_power = np.loadtxt(cs_power_per_component_output_file, dtype=str)[1:,:]
	t = open(cs_power_output_file,'w')
	t.write(('eQTL_sample_size\tln_pi_method\ttwas_methods\tgenetic_element_class\tpower\tpower_lb\tpower_ub\n'))
	genetic_element_classes = np.asarray(['gene', 'variant'])
	for eqtl_sample_size in eqtl_sample_sizes:
		for ln_pi_method in ln_pi_methods:
			for twas_method in twas_methods:
				for element_class in genetic_element_classes:

					subset_indices = (per_component_power[:,0] == str(eqtl_sample_size)) & (per_component_power[:,1] == ln_pi_method) & (per_component_power[:,2] == twas_method) & ( per_component_power[:,5] == element_class)
					components_covered = (per_component_power[subset_indices,:][:,-1]).astype(float)
					prop = np.sum(components_covered)/len(components_covered)
					prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
					prop_lb = prop - (1.96*prop_se)
					prop_ub = prop + (1.96*prop_se)

					t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + element_class + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\n')

	t.close()
	return

def create_file_containing_averaged_focus_cs_power(focus_cs_power_per_component_output_file, focus_cs_power_output_file, eqtl_sample_sizes):
	per_component_power = np.loadtxt(focus_cs_power_per_component_output_file, dtype=str)[1:,:]
	t = open(focus_cs_power_output_file,'w')
	t.write(('eQTL_sample_size\tgenetic_element_class\tpower\tpower_lb\tpower_ub\n'))
	genetic_element_classes = np.asarray(['gene', 'variant'])
	for eqtl_sample_size in eqtl_sample_sizes:

		subset_indices = (per_component_power[:,0] == str(eqtl_sample_size))
		components_covered = (per_component_power[subset_indices,:][:,-1]).astype(float)
		prop = np.sum(components_covered)/len(components_covered)
		prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
		prop_lb = prop - (1.96*prop_se)
		prop_ub = prop + (1.96*prop_se)

		t.write(str(eqtl_sample_size) + '\t' + 'gene' + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\n')

	t.close()
	return

def create_file_containing_expected_fraction_of_tgfm_elements_mediated_by_gene_expression(global_window_file, eqtl_sample_sizes, simulation_runs, twas_method, ln_pi_method, global_simulation_name_string, simulated_tgfm_results_dir, expected_fraction_of_elements_mediated_by_gene_expression_summary_file):
	# Get window names
	window_names = np.loadtxt(global_window_file,dtype=str,delimiter='\t')[1:,0]

	# Keep track of number of high pip elements across windows
	n_high_pip_elements = {}
	for eqtl_sample_size in eqtl_sample_sizes:
		n_high_pip_elements[eqtl_sample_size] = (0.0, 0.0)  # gene, variant

	for simulation_run in simulation_runs:
		print(simulation_run)
		for eqtl_sample_size in eqtl_sample_sizes:
			for window_name in window_names:
				tgfm_results_pkl = simulated_tgfm_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_' + window_name + '_results.pkl'

				# Skip window if run not successs
				if os.path.isfile(tgfm_results_pkl) == False:
					continue

				# Load in pkl file
				g = open(tgfm_results_pkl, "rb")
				tgfm_res = pickle.load(g)
				g.close()

				# Extract pips in middle
				middle_gene_pips = tgfm_res['expected_alpha_pips'][tgfm_res['middle_gene_indices']]
				middle_variant_pips = tgfm_res['expected_beta_pips'][tgfm_res['middle_variant_indices']]

				# Update global counter
				old_tuple = n_high_pip_elements[eqtl_sample_size]
				new_tuple = (old_tuple[0] + np.sum(middle_gene_pips), old_tuple[1] + np.sum(middle_variant_pips))
				n_high_pip_elements[eqtl_sample_size] = new_tuple


	# Print to output file
	t = open(expected_fraction_of_elements_mediated_by_gene_expression_summary_file,'w')
	# Header 
	t.write('eQTL_sample_size\tnum_genes\tnum_snps\texpression_mediated_fraction\texpression_mediated_fraction_lb\texpression_mediated_fraction_ub\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		num_genes = n_high_pip_elements[eqtl_sample_size][0]
		num_variants = n_high_pip_elements[eqtl_sample_size][1]
		prop = num_genes/(num_genes + num_variants)

		prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(num_genes + num_variants)
		prop_lb = prop - (1.96*prop_se)
		prop_ub = prop + (1.96*prop_se)

		t.write(str(eqtl_sample_size) + '\t' + str(num_genes) + '\t' + str(num_variants) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\n')
	t.close()
	return	

def create_file_containing_fraction_of_high_pip_tgfm_elements_mediated_by_gene_expression(global_window_file, eqtl_sample_sizes, simulation_runs, pip_thresholds, twas_method, ln_pi_method, global_simulation_name_string, simulated_tgfm_results_dir, fraction_high_pip_elements_mediated_by_gene_expression_summary_file):
	# Get window names
	window_names = np.loadtxt(global_window_file,dtype=str,delimiter='\t')[1:,0]

	# Keep track of number of high pip elements across windows
	n_high_pip_elements = {}
	for eqtl_sample_size in eqtl_sample_sizes:
		n_high_pip_elements[eqtl_sample_size] = {}
		for pip_threshold in pip_thresholds:
			n_high_pip_elements[eqtl_sample_size][pip_threshold] = (0.0, 0.0)

	for simulation_run in simulation_runs:
		print(simulation_run)
		for eqtl_sample_size in eqtl_sample_sizes:
			#tgfm_results_pkl = simulated_tgfm_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_' + window_name + '_results.pkl'
			tgfm_results_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method  + '_tgfm_pip_summary.txt'
			f = open(tgfm_results_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				if len(data) != 3:
					continue
				if data[1] == 'NA':
					continue
				ele_names = np.asarray(data[1].split(';'))
				ele_pips = np.asarray(data[2].split(';')).astype(float)

				for pip_threshold in pip_thresholds:
					indices = ele_pips >= pip_threshold
					high_pip_eles = ele_names[indices]
					for high_pip_ele in high_pip_eles:
						if high_pip_ele.startswith('ENSG'):
							prev_tuple = n_high_pip_elements[eqtl_sample_size][pip_threshold]
							new_tuple = (prev_tuple[0] + 1, prev_tuple[1])
							n_high_pip_elements[eqtl_sample_size][pip_threshold] = new_tuple
						else:
							prev_tuple = n_high_pip_elements[eqtl_sample_size][pip_threshold]
							new_tuple = (prev_tuple[0], prev_tuple[1] + 1)
							n_high_pip_elements[eqtl_sample_size][pip_threshold] = new_tuple
			f.close()

	# Print to output file
	t = open(fraction_high_pip_elements_mediated_by_gene_expression_summary_file,'w')
	# Header 
	t.write('eQTL_sample_size\tPIP_threshold\tnum_genes\tnum_snps\texpression_mediated_fraction\texpression_mediated_fraction_lb\texpression_mediated_fraction_ub\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		for pip_threshold in pip_thresholds:
			num_genes = n_high_pip_elements[eqtl_sample_size][pip_threshold][0]
			num_variants = n_high_pip_elements[eqtl_sample_size][pip_threshold][1]
			prop = num_genes/(num_genes + num_variants)

			prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(num_genes + num_variants)
			prop_lb = prop - (1.96*prop_se)
			prop_ub = prop + (1.96*prop_se)

			t.write(str(eqtl_sample_size) + '\t' + str(pip_threshold) + '\t' + str(num_genes) + '\t' + str(num_variants) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\n')
	t.close()
	return


#######################
# Command line args
#######################
chrom_num = sys.argv[1]
cis_window = sys.argv[2]
n_gwas_individuals = sys.argv[3]
global_simulation_name_string = sys.argv[4]
total_heritability = sys.argv[5]
fraction_expression_mediated_heritability = sys.argv[6]
simulated_sldsc_results_dir = sys.argv[7]
simulated_organized_results_dir = sys.argv[8]
simulated_tgfm_results_dir = sys.argv[9]
simulated_trait_dir = sys.argv[10]
simulated_gene_expression_dir = sys.argv[11]
simulated_learned_gene_models_dir = sys.argv[12]
simulated_tgfm_input_data_dir = sys.argv[13]
simulated_gene_position_file = sys.argv[14]
processed_genotype_data_dir = sys.argv[15]
simulated_ld_scores_dir = sys.argv[16]
simulated_focus_results_dir = sys.argv[17]
simulated_coloc_results_dir = sys.argv[18]


#############################################################
# Genome-wide analysis
#############################################################
#Bim file
bim_file = processed_genotype_data_dir + 'simulated_gwas_data_' + chrom_num + '.bim'

# Used eQTL sample sizes
eqtl_sample_sizes = np.asarray([300,500,1000])

'''
# Simulation runs
# Currently hacky because had some failed simulations
simulation_runs = np.arange(1,21)
thresholds = [1e-1, 1e-2, 1e-3, 5e-4, 2.5e-4, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-14, 1e-16, 1e-18]

model_names = ['susie_pmces_uniform_iterative_variant_gene_prior_pip_level']

for model_name in model_names:

	##########################
	# Power and Type 1 Error
	##########################
	# Create file showing p-value in causal tissues and non-causal tissues
	mediated_iterative_sampler_pvalue_by_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + model_name + '_mediated_h2_pvalue_by_tissue_est_across_thresholds.txt'
	create_file_containing_iterative_bootstrapped_sampler_prior_mediated_h2_pvalue_in_causal_and_non_causal_tissues_across_thresholds(simulated_tgfm_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_iterative_sampler_pvalue_by_tissue_output_file, thresholds, model_name)
	# Create file showing power to detect causal tissues
	power_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + model_name + '_mediated_h2_power_across_thresholds.txt'
	create_file_containing_mediated_h2_power_to_detect_causal_tissues_across_thresholds(mediated_iterative_sampler_pvalue_by_tissue_output_file, power_output_file, eqtl_sample_sizes, thresholds)
	# Create file showing type 1 error for null tissues
	type_1_error_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + model_name + '_h2_type_1_error_across_thresholds.txt'
	create_file_containing_mediated_h2_type_1_error_across_thresholds(mediated_iterative_sampler_pvalue_by_tissue_output_file, type_1_error_output_file, eqtl_sample_sizes, thresholds)


	##########################
	# Bias in fraction of causal components going through gene expression
	##########################
	# Create file showing mediated h2 in causal tissues and non-causal tissues
	fraction_expr_med_disease_components_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + model_name + '_expected_fraction_expression_mediated_disease_components.txt'
	create_file_containing_expected_fraction_of_expression_mediated_disease_components(simulated_tgfm_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, model_name, simulated_ld_scores_dir, fraction_expr_med_disease_components_output_file)
	# Average estimates across simulations
	organized_avg_fraction_expr_med_disease_components_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + model_name + '_avg_expected_fraction_expression_mediated_disease_components.txt'
	create_file_containing_avg_fraction_h2_across_simulation_runs(fraction_expr_med_disease_components_output_file, organized_avg_fraction_expr_med_disease_components_output_file, eqtl_sample_sizes)
	print(organized_avg_fraction_expr_med_disease_components_output_file)


	##########################
	# Bias in fraction of causal genes per tissue
	##########################
	# Create file showing mediated h2 in causal tissues and non-causal tissues
	fraction_causal_by_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + model_name + '_fraction_causal_by_tissue.txt'
	create_file_containing_fraction_causal_in_causal_and_non_causal_tissues(simulated_tgfm_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, model_name, fraction_causal_by_tissue_output_file)
	# Average estimates across simulations
	organized_avg_fraction_causal_by_tissue_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + model_name + '_avg_fraction_causal_by_tissue.txt'
	#create_file_containing_avg_med_h2_by_tissue_across_simulation_runs(fraction_causal_by_tissue_output_file, organized_avg_fraction_causal_by_tissue_output_file, eqtl_sample_sizes)
	create_file_containing_avg_fraction_causal_by_tissue_across_simulation_runs(fraction_causal_by_tissue_output_file, organized_avg_fraction_causal_by_tissue_output_file, eqtl_sample_sizes)
'''

#############################################################
# Fraction of TGFM elements mediated by gene expression
#############################################################
# Used eQTL sample sizes
eqtl_sample_sizes = np.asarray([300,500,1000])


# Simulation runs
# Currently hacky because had some failed simulations
simulation_runs = np.arange(1,21)

# Thresholds to consider
pip_thresholds=np.asarray([.1, .3, .5, .7, .9])

# Model specification
twas_method = 'susie_sampler'
ln_pi_method = 'pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'

# Windows
global_window_file = processed_genotype_data_dir + 'chromosome_' + str(chrom_num) + '_windows_3_mb.txt'

# Create file containing fraction of high pip elememnts mediated by gene expresssion
fraction_high_pip_elements_mediated_by_gene_expression_summary_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + twas_method + '_' + ln_pi_method + '_fraction_high_pip_elements_mediated_by_gene_expresssion.txt'
#create_file_containing_fraction_of_high_pip_tgfm_elements_mediated_by_gene_expression(global_window_file, eqtl_sample_sizes, simulation_runs, pip_thresholds, twas_method, ln_pi_method, global_simulation_name_string, simulated_tgfm_results_dir, fraction_high_pip_elements_mediated_by_gene_expression_summary_file)

# Create file containing expected fraction of elememnts mediated by gene expresssion
expected_fraction_of_elements_mediated_by_gene_expression_summary_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + twas_method + '_' + ln_pi_method + '_expected_fraction_elements_mediated_by_gene_expresssion.txt'
create_file_containing_expected_fraction_of_tgfm_elements_mediated_by_gene_expression(global_window_file, eqtl_sample_sizes, simulation_runs, twas_method, ln_pi_method, global_simulation_name_string, simulated_tgfm_results_dir, expected_fraction_of_elements_mediated_by_gene_expression_summary_file)



'''
#############################################################
# Fine-mapping evaluation metrics
#############################################################
# Used eQTL sample sizes
eqtl_sample_sizes = np.asarray([300,500,1000])


# Simulation runs
# Currently hacky because had some failed simulations
simulation_runs = np.arange(1,21)
#simulation_runs = np.delete(simulation_runs, [1, 6, 17, 18])


# ln_pi methods used
#ln_pi_methods = np.asarray(['uniform', 'tglr_variant_gene', 'tglr_sparse_variant_gene_tissue', 'tglr_bootstrapped_nonnegative_sampler', 'iterative_variant_gene_tissue', 'iterative_variant_gene_tissue_bootstrapped', 'pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped', 'sampler_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])
#ln_pi_methods = np.asarray(['uniform', 'tglr_variant_gene', 'tglr_bootstrapped_nonnegative_pmces', 'tglr_bootstrapped_nonnegative_sampler', 'pmces_uniform_iterative_variant_gene_prior_pip_level_pmces', 'pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped', 'sampler_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])
ln_pi_methods = np.asarray(['uniform', 'pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])


# twas method
twas_methods = np.asarray(['susie_pmces', 'susie_sampler'])


#Bim file
bim_file = processed_genotype_data_dir + 'simulated_gwas_data_' + chrom_num + '.bim'

# Windows
global_window_file = processed_genotype_data_dir + 'chromosome_' + str(chrom_num) + '_windows_3_mb.txt'

##################################
# Coverage/Calibration to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]

# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods)


##################################
# Power to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]

for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_tgfm_high_pip_snp_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods,twas_methods)



##################################
# FOCUS Coverage/Calibration to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, simulation_num, window_num, boolean_causal_variant_in_cs)
	focus_cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_focus_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_focus_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, pip_threshold, focus_cs_coverage_per_high_pip_snp_output_file)

	focus_cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_focus_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_focus_high_pip_calibration(focus_cs_coverage_per_high_pip_snp_output_file, focus_cs_high_pip_coverage_output_file, eqtl_sample_sizes)

##################################
# FOCUS Power to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]

for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	focus_cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_focus_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_focus_high_pip_snp_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file)

	focus_cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_focus_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_focus_cs_power(focus_cs_power_per_component_output_file, focus_cs_power_output_file, eqtl_sample_sizes)



##################################
# coloc Coverage/Calibration to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, simulation_num, window_num, boolean_causal_variant_in_cs)
	coloc_cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_coloc_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_coloc_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_coloc_results_dir, pip_threshold, coloc_cs_coverage_per_high_pip_snp_output_file)

	coloc_cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_coloc_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_focus_high_pip_calibration(coloc_cs_coverage_per_high_pip_snp_output_file, coloc_cs_high_pip_coverage_output_file, eqtl_sample_sizes)

##################################
# coloc Power to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]

for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	coloc_cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_coloc_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_coloc_high_pip_snp_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_coloc_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, coloc_cs_power_per_component_output_file, pip_threshold, global_window_file)

	coloc_cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_coloc_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_focus_cs_power(coloc_cs_power_per_component_output_file, coloc_cs_power_output_file, eqtl_sample_sizes)


'''






















































































##################################
# Fine-mapping evaluation metrics
##################################
# (i) Coverage/Calibration: the proportion of credible sets that include at least one true causal genetic element across simulation replicates;
# (ii) Power: the number of true causal variants identified (i.e., covered by a credible set)
# (iii) Resolution: the size of credible sets and the number of fine-mapped variants with high confidence (e.g., PIP >95%);
# These should be done for both genes and variants

'''
##################################
# Coverage/Calibration to detect snps with PIP > threshold
##################################
pip_thresholds = [.1, .3, .5, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, initialization_versions, cs_coverage_per_high_pip_snp_output_file)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, initialization_versions, twas_methods)

##################################
# Power to detect snps with PIP > threshold
##################################
pip_thresholds = [.1, .3, .5, .9, .95, .99]

for pip_threshold in pip_thresholds:
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_tgfm_high_pip_snp_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, initialization_versions, twas_methods)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods, initialization_versions,twas_methods)
'''
'''
##################################
# Coverage/Calibration to detect snps with PIP > threshold in windows where causal genes are detected
##################################
pip_thresholds = [.1, .3, .5, .9, .95, .99]
pip_thresholds = [.9]

# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_calibration_per_component_where_causal_gene_is_detected.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp_where_causal_gene_is_detected(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, initialization_versions, cs_coverage_per_high_pip_snp_output_file, simulated_gene_position_file)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_calibration_where_causal_gene_is_detected.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, initialization_versions, twas_methods)
'''



'''

# Simulation runs
# Currently hacky because had some failed simulations
simulation_runs = np.arange(1,100)
simulation_runs = np.delete(simulation_runs, [17,88])

##############################
# Calculate number of of cis-heritable genes
##############################
organized_number_detected_genes_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_fraction_of_detected_heritable_genes.txt'
#create_file_containing_number_of_detected_genes(global_simulation_name_string, simulation_runs, eqtl_sample_sizes, simulated_gene_expression_dir, simulated_learned_gene_models_dir, organized_number_detected_genes_output_file)

##############################
# Estimate bias in TGFM-SLDSC
##############################
# Create file comparing total heritabilities across simulation runs
organized_total_h2_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_total_est_h2.txt'
create_file_containing_total_h2_across_simulation_runs(simulated_sldsc_results_dir, global_simulation_name_string, total_heritability, eqtl_sample_sizes, simulation_runs, organized_total_h2_output_file)
# Average estimates across simulations
organized_avg_total_h2_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_avg_total_est_h2.txt'
create_file_containing_avg_total_h2_across_simulation_runs(organized_total_h2_output_file, organized_avg_total_h2_output_file, eqtl_sample_sizes)


# Create file comparing fraction mediated heritabilities across simulation runs
organized_fraction_h2_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_fraction_expression_med_h2_est.txt'
create_file_containing_fraction_expression_mediated_h2_across_simulation_runs(simulated_sldsc_results_dir, global_simulation_name_string, fraction_expression_mediated_heritability, eqtl_sample_sizes, simulation_runs, organized_fraction_h2_output_file)
# Average estimates across simulations
organized_avg_fraction_h2_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_avg_fraction_expression_med_h2.txt'
create_file_containing_avg_fraction_h2_across_simulation_runs(organized_fraction_h2_output_file, organized_avg_fraction_h2_output_file, eqtl_sample_sizes)

# Create file showing mediated h2 in causal tissues and non-causal tissues
mediated_h2_by_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_mediated_h2_by_tissue_est.txt'
create_file_containing_mediated_h2_in_causal_and_non_causal_tissues(simulated_sldsc_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_h2_by_tissue_output_file)
# Average estimates across simulations
organized_avg_mediated_h2_by_tissue_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_avg_mediated_h2_by_tissue.txt'
create_file_containing_avg_med_h2_by_tissue_across_simulation_runs(mediated_h2_by_tissue_output_file, organized_avg_mediated_h2_by_tissue_output_file, eqtl_sample_sizes)


# Create file showing mediated h2 in causal tissues and non-causal tissues based on sparse model
mediated_h2_by_tissue_sparse_est_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_mediated_h2_by_tissue_sparse_est.txt'
create_file_containing_mediated_h2_in_causal_and_non_causal_tissues_sparse_estimate(simulated_sldsc_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_h2_by_tissue_sparse_est_output_file)
# Average estimates across simulations
organized_avg_mediated_h2_by_tissue_sparse_est_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_avg_mediated_h2_by_tissue_sparse_est.txt'
create_file_containing_avg_med_h2_by_tissue_across_simulation_runs(mediated_h2_by_tissue_sparse_est_output_file, organized_avg_mediated_h2_by_tissue_sparse_est_output_file, eqtl_sample_sizes)


##############################
# Power and type 1 error in sparse TGFM-SLDSC
##############################
# Create file showing binary status in causal tissues and non-causal tissues
mediated_binary_variable_by_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_mediated_h2_sparse_binary_est_by_tissue.txt'
create_file_containing_mediated_h2_sparse_binary_est_in_causal_and_non_causal_tissues(simulated_sldsc_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_binary_variable_by_tissue_output_file)
# Create file showing power to detect causal tissues
power_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_sparse_mediated_h2_power.txt'
create_file_containing_mediated_h2_power_to_detect_causal_tissues(mediated_binary_variable_by_tissue_output_file, power_output_file, eqtl_sample_sizes)
# Create file showing type 1 error for null tissues
type_1_error_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_sparse_mediated_h2_type_1_error.txt'
create_file_containing_mediated_h2_type_1_error(mediated_binary_variable_by_tissue_output_file, type_1_error_output_file, eqtl_sample_sizes)


##############################
# Power and type 1 error in TGFM-SLDSC
##############################
eqtl_types = ['susie_pmces', 'susie_distr']
anno_types = ['genotype_intercept', 'full_anno']

for eqtl_type in eqtl_types:
	for anno_type in anno_types:
		# Create file showing p-value in causal tissues and non-causal tissues
		mediated_pvalue_by_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + eqtl_type + '_' + anno_type + '_mediated_h2_pvalue_by_tissue_est.txt'
		create_file_containing_mediated_h2_pvalue_in_causal_and_non_causal_tissues(simulated_sldsc_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_pvalue_by_tissue_output_file, eqtl_type, anno_type)
		# Create file showing power to detect causal tissues
		power_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + eqtl_type + '_' + anno_type + '_mediated_h2_power.txt'
		create_file_containing_mediated_h2_power_to_detect_causal_tissues(mediated_pvalue_by_tissue_output_file, power_output_file, eqtl_sample_sizes)
		# Create file showing type 1 error for null tissues
		type_1_error_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + eqtl_type + '_' + anno_type + '_mediated_h2_type_1_error.txt'
		create_file_containing_mediated_h2_type_1_error(mediated_pvalue_by_tissue_output_file, type_1_error_output_file, eqtl_sample_sizes)
'''

'''
##############################
# Power and type 1 error in TGFM-SLDSC non-negative bootstrapped
##############################
eqtl_types = ['susie_pmces', 'susie_distr']
anno_types = ['genotype_intercept', 'full_anno']
thresholds = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12]

for eqtl_type in eqtl_types:
	for anno_type in anno_types:
		# Create file showing p-value in causal tissues and non-causal tissues
		mediated_nonnegative_pvalue_by_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + eqtl_type + '_' + anno_type + '_nonnegative_bootstrapped_mediated_h2_pvalue_by_tissue_est_across_thresholds.txt'
		create_file_containing_nonnegative_bootstrapped_mediated_h2_pvalue_in_causal_and_non_causal_tissues_across_thresholds(simulated_sldsc_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_nonnegative_pvalue_by_tissue_output_file, eqtl_type, anno_type, thresholds)
		# Create file showing power to detect causal tissues
		power_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + eqtl_type + '_' + anno_type + '_nonnegative_bootstrapped_mediated_h2_power_across_thresholds.txt'
		create_file_containing_mediated_h2_power_to_detect_causal_tissues_across_thresholds(mediated_nonnegative_pvalue_by_tissue_output_file, power_output_file, eqtl_sample_sizes, thresholds)
		# Create file showing type 1 error for null tissues
		type_1_error_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + eqtl_type + '_' + anno_type + '_nonnegative_bootstrapped_mediated_h2_type_1_error_across_thresholds.txt'
		create_file_containing_mediated_h2_type_1_error_across_thresholds(mediated_nonnegative_pvalue_by_tissue_output_file, type_1_error_output_file, eqtl_sample_sizes, thresholds)

		# Assess bias in looking at per-tissue taus
		mediated_nonnegative_taus_per_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + eqtl_type + '_' + anno_type + '_mediated_nonnegative_tau_tissue_est.txt'
		create_file_containing_mediated_nonnegative_tissue_taus_in_causal_and_non_causal_tissues(simulated_sldsc_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_nonnegative_taus_per_tissue_output_file, eqtl_type, anno_type)
		organized_avg_mediated_nonnegative_taus_by_tissue_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_' + eqtl_type + '_' + anno_type + '_avg_mediated_nonnegative_tau_tissue_est.txt'
		create_file_containing_avg_med_taus_by_tissue_across_simulation_runs(mediated_nonnegative_taus_per_tissue_output_file, organized_avg_mediated_nonnegative_taus_by_tissue_output_file, eqtl_sample_sizes)
'''

