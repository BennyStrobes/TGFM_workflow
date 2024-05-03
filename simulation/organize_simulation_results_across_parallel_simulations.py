import numpy as np 
import os
import sys
import pdb
import scipy.stats
import pickle
#from sklearn.linear_model import LinearRegression

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


def create_file_containing_average_prior_value_in_causal_and_non_causal_tissues(simulated_tgfm_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file, model_name, gene_type):
	if gene_type.startswith('max'):
		simulation_runs = simulation_runs + 100
	# Open output file and print header
	t = open(mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\ttissue_number\taverage_prior_value\tcausal_status\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			if simulation_run == 60 or simulation_run == 81 or simulation_run == 104 or simulation_run == 117:
				continue
			# Extract per annotation mediated heritability file for this run
			h2_file = simulated_tgfm_results_dir + 'simulation_new_' + str(simulation_run) + '_' + global_simulation_name_string + '_all_t_' + gene_type + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + model_name + '_bootstrapped.txt'
			pdb.set_trace()
			# Load in per anno med h2
			h2_data_raw = np.loadtxt(h2_file, dtype=str, delimiter='\t')
			expr_data = h2_data_raw[-10:,:]
			for tissue_iter in range(10):
				if tissue_iter == 0 or tissue_iter == 3:
					causal_status = 'causal'
				else:
					causal_status = 'null'
				bootstrapped_coefficients = np.asarray(expr_data[tissue_iter,-1].split(';')).astype(float)
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(tissue_iter) + '\t' + str(np.mean(bootstrapped_coefficients)) + '\t' + causal_status + '\n')
	t.close()		
	return



def create_file_containing_iterative_bootstrapped_sampler_prior_mediated_h2_gaussian_approximation_pvalue_in_causal_and_non_causal_tissues(simulated_tgfm_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file, model_name, gene_type='component_gene'):
	# Open output file and print header
	t = open(mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\ttissue_number\tz_score\tpvalue\tcausal_status\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			h2_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_nsamp_100_' + gene_type + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + model_name + '_bootstrapped.txt'
			# Load in per anno med h2
			h2_data_raw = np.loadtxt(h2_file, dtype=str, delimiter='\t')
			expr_data = h2_data_raw[-10:,:]
			for tissue_iter in range(10):
				if tissue_iter == 0 or tissue_iter == 3:
					causal_status = 'causal'
				else:
					causal_status = 'null'

				bootstrapped_coefficients = np.asarray(expr_data[tissue_iter,-1].split(';')).astype(float)
				if np.std(bootstrapped_coefficients, ddof=1) == 0.0:
					t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(tissue_iter) + '\t' + str(0.0) + '\t' + str(1.0) + '\t' + causal_status + '\n')
				else:
					zz = np.mean(bootstrapped_coefficients)/np.std(bootstrapped_coefficients, ddof=1)
					p_value = scipy.stats.norm.sf(abs(zz))  # 1-sided pvalue
					t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(tissue_iter) + '\t' + str(0.0) + '\t' + str(p_value) + '\t' + causal_status + '\n')
	t.close()		
	return


def create_file_containing_iterative_bootstrapped_sampler_prior_mediated_h2_pvalue_in_causal_and_non_causal_tissues_across_thresholds(simulated_tgfm_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_iterative_sampler_pvalue_by_tissue_output_file, thresholds, model_name):
	# Open output file and print header
	t = open(mediated_iterative_sampler_pvalue_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\ttissue_number\tz_score\tpvalue\tcausal_status\tthreshold\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			h2_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_nsamp_100_component_gene' + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + model_name + '_bootstrapped.txt'
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
			h2_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_nsamp_100_component_gene' + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + model_name + '_bootstrapped.txt'
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
			h2_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string+ '_nsamp_100_component_gene' + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + model_name + '_bootstrapped.txt'
			# Load in per anno med h2
			h2_data_raw = np.loadtxt(h2_file, dtype=str, delimiter='\t')
			med_tau = h2_data_raw[-10:,1].astype(float)

			for tissue_iter in range(10):
				if tissue_iter == 0 or tissue_iter == 3:
					sim_med_h2 = 30/2000.0
				else:
					sim_med_h2 = 0.0
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(tissue_iter) + '\t' + str(med_tau[tissue_iter]) + '\t' + str(sim_med_h2) + '\n')
			t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + 'nm_variant' + '\t' + str(h2_data_raw[1,1]) + '\t' + str(sim_med_h2) + '\n')
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


def create_file_containing_mediated_h2_power_to_detect_causal_tissues_gaussian_approximation(mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file, power_output_file, eqtl_sample_sizes):
	raw_pvalue_file = np.loadtxt(mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file,dtype=str, delimiter='\t')

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

def create_file_containing_mediated_h2_type_1_error_gaussian_approximation(mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file, type_1_error_output_file, eqtl_sample_sizes):
	raw_pvalue_file = np.loadtxt(mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file,dtype=str, delimiter='\t')

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
		causal_status = 'causal'
		raw_subset2 = raw_subset[raw_subset[:,2] == 'nm_variant',:]
		est = raw_subset2[:,3].astype(float)
		meaner = np.mean(est)
		se = np.std(est)/np.sqrt(len(est))
		t.write(str(eqtl_sample_size) + '\t' + 'nm_variant' + '\t' + causal_status + '\t' + str(meaner) + '\t' + str(meaner - (1.96*se)) + '\t' + str(meaner + (1.96*se)) + '\n')

	t.close()
	return


def create_file_containing_avg_fraction_causal_by_tissue_causal_status_across_simulation_runs(mediated_h2_by_tissue_output_file, organized_avg_mediated_h2_by_tissue_output_file, eqtl_sample_sizes):
	raw_file = np.loadtxt(mediated_h2_by_tissue_output_file, dtype=str,delimiter='\t')

	t = open(organized_avg_mediated_h2_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\tgenetic_element_type\tcausal_status\tfraction_causal\tfraction_causal_lb\tfraction_causal_ub\n')

	for eqtl_sample_size in eqtl_sample_sizes:
		raw_subset = raw_file[raw_file[:,0] == str(eqtl_sample_size), :]
		
		causal_status = 'causal'
		raw_subset2 = raw_subset[(raw_subset[:,2] == str(0)) | (raw_subset[:,2] == str(3)), :]
		est = raw_subset2[:,3].astype(float)
		meaner = np.mean(est)
		se = np.std(est)/np.sqrt(len(est))
		t.write(str(eqtl_sample_size) + '\t' + 'gene_tissue' + '\t' + causal_status + '\t' + str(meaner) + '\t' + str(meaner - (1.96*se)) + '\t' + str(meaner + (1.96*se)) + '\n')

		causal_status = 'null'
		raw_subset2 = raw_subset[(raw_subset[:,2] != str(0)) & (raw_subset[:,2] != str(3)) & (raw_subset[:,2] != 'nm_variant'), :]
		est = raw_subset2[:,3].astype(float)
		meaner = np.mean(est)
		se = np.std(est)/np.sqrt(len(est))
		t.write(str(eqtl_sample_size) + '\t' + 'gene_tissue' + '\t' + causal_status + '\t' + str(meaner) + '\t' + str(meaner - (1.96*se)) + '\t' + str(meaner + (1.96*se)) + '\n')


		causal_status = 'causal'
		raw_subset2 = raw_subset[raw_subset[:,2] == 'nm_variant',:]
		est = raw_subset2[:,3].astype(float)
		meaner = np.mean(est)
		se = np.std(est)/np.sqrt(len(est))
		t.write(str(eqtl_sample_size) + '\t' + 'nm_variant' + '\t' + causal_status + '\t' + str(meaner) + '\t' + str(meaner - (1.96*se)) + '\t' + str(meaner + (1.96*se)) + '\n')

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

def extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file, gene_level=False):
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
				if gene_level:
					gene_tissue_name = gene_tissue_name.split('_')[0]
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

def compute_number_of_heritable_and_non_heritable_genes_with_a_model(ordered_gene_names, eqtl_sample_size, simulation_run, global_simulation_name_string, simulated_learned_gene_models_dir, heritable_genes_dicti, not_heritable_genes_dicti, gene_type):
	n_her_genes = 0.0
	n_non_her_genes = 0.0
	for gene_name in ordered_gene_names:
		if gene_type == 'all_non_zero_gene':
			learned_gene_pmces_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_run) + '/simulation_new_' + str(simulation_run) + '_' + global_simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_sample_size) + '_gene_model_pmces.npy'
			learned_gene_pmces = np.load(learned_gene_pmces_file)
			detected_genes = np.var(learned_gene_pmces,axis=1) != 0.0
		elif gene_type == 'component_gene':
			valid_comp_file = learned_gene_pmces_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_run) + '/simulation_new_' + str(simulation_run) + '_' + global_simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_sample_size) + '_gene_valid_susie_comp.npy'
			detected_genes_raw = np.load(valid_comp_file)
			detected_genes = []
			for ele in detected_genes_raw:
				if ele == 'True':
					detected_genes.append(True)
				else:
					detected_genes.append(False)
			detected_genes = np.asarray(detected_genes)
		elif gene_type.startswith('max_min_ratio_'):
			max_min_ratio_file_file = learned_gene_pmces_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_run) + '/simulation_new_' + str(simulation_run) + '_' + global_simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_sample_size) + '_gene_best_to_worst_susie_pi_ratios.npy'
			max_min_ratio = np.load(max_min_ratio_file_file)
			max_min_thresh = float(gene_type.split('_ratio_')[1])
			detected_genes = max_min_ratio > max_min_thresh
		elif gene_type == 'cafeh':
			cafeh_pi_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_run) + '/simulation_new_' + str(simulation_run) + '_' + global_simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_sample_size) + '_cafeh_pi.npy'
			cafeh_pi = np.load(cafeh_pi_file)
			cafeh_p_active_file= simulated_learned_gene_models_dir + 'simulation_' + str(simulation_run) + '/simulation_new_' + str(simulation_run) + '_' + global_simulation_name_string + '_' + gene_name + '_eqtlss_' + str(eqtl_sample_size) + '_cafeh_p_active.npy'
			cafeh_p_active = np.load(cafeh_p_active_file)
			n_tiss = cafeh_p_active.shape[0]
			n_comp = cafeh_p_active.shape[1]
			detected_genes = np.asarray([False]*n_tiss)

			for kk in range(n_comp):
				max_min_ratio = np.max(cafeh_pi[kk,:])/np.min(cafeh_pi[kk,:])
				if max_min_ratio < 5.0:
					continue
				component_bools = cafeh_p_active[:,kk] > .2
				for ii, component_bool in enumerate(component_bools):
					if component_bool:
						detected_genes[ii] = True
		else:
			print('assumption eroror')
			pdb.set_trace()
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
	print(organized_number_detected_genes_output_file)
	t = open(organized_number_detected_genes_output_file,'w')
	t.write('simulation_number\teQTL_sample_size\tgene_model_type\tn_detected_heritable_genes\tn_heritable_genes\tn_detected_non_heritable_genes\tn_non_heritable_genes\n')
	for simulation_run in simulation_runs:
		print(simulation_run)
		# Within a simulation_run, the true causal heritable genes are the same across eqtl sample sizes
		# Extract dictionary list of heritable_genes and not heritable_genes
		simulation_gene_summary_file = simulated_gene_expression_dir + 'simulation_new_' + str(simulation_run) + '_' + global_simulation_name_string + '_causal_eqtl_effect_summary.txt'
		heritable_genes_dicti, not_heritable_genes_dicti, ordered_gene_names = extract_dictionary_lists_of_heritable_and_non_heritable_genes(simulation_gene_summary_file)

		# Loop through eqtl sample sizes
		for eqtl_sample_size in eqtl_sample_sizes:

			# At this sample size, compute the number of heritable genes with a model and the number of non-heritable genes with a model
			gene_type = 'max_min_ratio_2'
			n_heritable_genes_modeled, n_non_heritable_genes_modeled = compute_number_of_heritable_and_non_heritable_genes_with_a_model(ordered_gene_names, eqtl_sample_size, simulation_run, global_simulation_name_string, simulated_learned_gene_models_dir, heritable_genes_dicti, not_heritable_genes_dicti, gene_type)
			t.write(str(simulation_run) + '\t' + str(eqtl_sample_size) + '\t' + gene_type + '\t' + str(n_heritable_genes_modeled) + '\t' + str(len(heritable_genes_dicti)) + '\t' + str(n_non_heritable_genes_modeled) + '\t' + str(len(not_heritable_genes_dicti)) + '\n')
			gene_type = 'max_min_ratio_5'
			n_heritable_genes_modeled, n_non_heritable_genes_modeled = compute_number_of_heritable_and_non_heritable_genes_with_a_model(ordered_gene_names, eqtl_sample_size, simulation_run, global_simulation_name_string, simulated_learned_gene_models_dir, heritable_genes_dicti, not_heritable_genes_dicti, gene_type)
			t.write(str(simulation_run) + '\t' + str(eqtl_sample_size) + '\t' + gene_type + '\t' + str(n_heritable_genes_modeled) + '\t' + str(len(heritable_genes_dicti)) + '\t' + str(n_non_heritable_genes_modeled) + '\t' + str(len(not_heritable_genes_dicti)) + '\n')
			gene_type = 'max_min_ratio_10'
			n_heritable_genes_modeled, n_non_heritable_genes_modeled = compute_number_of_heritable_and_non_heritable_genes_with_a_model(ordered_gene_names, eqtl_sample_size, simulation_run, global_simulation_name_string, simulated_learned_gene_models_dir, heritable_genes_dicti, not_heritable_genes_dicti, gene_type)
			t.write(str(simulation_run) + '\t' + str(eqtl_sample_size) + '\t' + gene_type + '\t' + str(n_heritable_genes_modeled) + '\t' + str(len(heritable_genes_dicti)) + '\t' + str(n_non_heritable_genes_modeled) + '\t' + str(len(not_heritable_genes_dicti)) + '\n')
			gene_type = 'max_min_ratio_50'
			n_heritable_genes_modeled, n_non_heritable_genes_modeled = compute_number_of_heritable_and_non_heritable_genes_with_a_model(ordered_gene_names, eqtl_sample_size, simulation_run, global_simulation_name_string, simulated_learned_gene_models_dir, heritable_genes_dicti, not_heritable_genes_dicti, gene_type)
			t.write(str(simulation_run) + '\t' + str(eqtl_sample_size) + '\t' + gene_type + '\t' + str(n_heritable_genes_modeled) + '\t' + str(len(heritable_genes_dicti)) + '\t' + str(n_non_heritable_genes_modeled) + '\t' + str(len(not_heritable_genes_dicti)) + '\n')
			gene_type = 'max_min_ratio_100'
			n_heritable_genes_modeled, n_non_heritable_genes_modeled = compute_number_of_heritable_and_non_heritable_genes_with_a_model(ordered_gene_names, eqtl_sample_size, simulation_run, global_simulation_name_string, simulated_learned_gene_models_dir, heritable_genes_dicti, not_heritable_genes_dicti, gene_type)
			t.write(str(simulation_run) + '\t' + str(eqtl_sample_size) + '\t' + gene_type + '\t' + str(n_heritable_genes_modeled) + '\t' + str(len(heritable_genes_dicti)) + '\t' + str(n_non_heritable_genes_modeled) + '\t' + str(len(not_heritable_genes_dicti)) + '\n')
			gene_type = 'max_min_ratio_200'
			n_heritable_genes_modeled, n_non_heritable_genes_modeled = compute_number_of_heritable_and_non_heritable_genes_with_a_model(ordered_gene_names, eqtl_sample_size, simulation_run, global_simulation_name_string, simulated_learned_gene_models_dir, heritable_genes_dicti, not_heritable_genes_dicti, gene_type)
			t.write(str(simulation_run) + '\t' + str(eqtl_sample_size) + '\t' + gene_type + '\t' + str(n_heritable_genes_modeled) + '\t' + str(len(heritable_genes_dicti)) + '\t' + str(n_non_heritable_genes_modeled) + '\t' + str(len(not_heritable_genes_dicti)) + '\n')



			gene_type = 'component_gene'
			# At this sample size, compute the number of heritable genes with a model and the number of non-heritable genes with a model
			n_heritable_genes_modeled, n_non_heritable_genes_modeled = compute_number_of_heritable_and_non_heritable_genes_with_a_model(ordered_gene_names, eqtl_sample_size, simulation_run, global_simulation_name_string, simulated_learned_gene_models_dir, heritable_genes_dicti, not_heritable_genes_dicti, gene_type)
			t.write(str(simulation_run) + '\t' + str(eqtl_sample_size) + '\t' + gene_type + '\t' + str(n_heritable_genes_modeled) + '\t' + str(len(heritable_genes_dicti)) + '\t' + str(n_non_heritable_genes_modeled) + '\t' + str(len(not_heritable_genes_dicti)) + '\n')
			gene_type = 'all_non_zero_gene'
			# At this sample size, compute the number of heritable genes with a model and the number of non-heritable genes with a model
			n_heritable_genes_modeled, n_non_heritable_genes_modeled = compute_number_of_heritable_and_non_heritable_genes_with_a_model(ordered_gene_names, eqtl_sample_size, simulation_run, global_simulation_name_string, simulated_learned_gene_models_dir, heritable_genes_dicti, not_heritable_genes_dicti, gene_type)
			t.write(str(simulation_run) + '\t' + str(eqtl_sample_size) + '\t' + gene_type + '\t' + str(n_heritable_genes_modeled) + '\t' + str(len(heritable_genes_dicti)) + '\t' + str(n_non_heritable_genes_modeled) + '\t' + str(len(not_heritable_genes_dicti)) + '\n')
		t.flush()
	t.close()
	return

def create_file_containing_coloc_cs_calibration_per_high_pip_gene(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_coloc_results_dir, pip_threshold, coloc_cs_coverage_per_high_pip_snp_output_file):
	# Open output file handle
	t = open(coloc_cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tsimulation_number\twindow_name\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		for gene_id in [*causal_genes]:
			causal_genetic_elements[gene_id.split('_')[0]] = 1

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
				ensamble = gene_id.split('_')[0]
				if ensamble not in used_genes:
					used_genes[ensamble] = (pip, region)
				else:
					if pip > used_genes[ensamble][0]:
						used_genes[ensamble] = (pip, region)
			f.close()

			for gene_id in [*used_genes]:
				if used_genes[gene_id][0] < pip_threshold:
					continue
				gene_pip = used_genes[gene_id][0]
				region_name = used_genes[gene_id][1]
				if gene_id in causal_genetic_elements:
					causal_genetic_element_in_cs_boolean = 1
				else:
					causal_genetic_element_in_cs_boolean = 0
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_number) + '\t' + region_name + '\t' + gene_id + '\t' + 'gene' + '\t' + str(causal_genetic_element_in_cs_boolean) + '\t' + str(gene_pip) + '\n')
	t.close()
	return


def create_file_containing_two_step_coloc_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_two_step_tgfm_results_dir,simulated_coloc_results_dir, pip_threshold, coloc_cs_coverage_per_high_pip_snp_output_file):
	# Open output file handle
	t = open(coloc_cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tsimulation_number\twindow_name\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
			two_step_tissue_summary_file = simulated_two_step_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_pmces_uniform_two_step_tissues.txt'
			best_tissue_arr, sig_tissue_arr = extract_two_step_analyzed_tissues(two_step_tissue_summary_file)

			best_tissue_dicti = {}
			best_tissue_dicti[best_tissue_arr[0]] = 1


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
				cur_tissue_name = gene_id.split('_')[1]
				if cur_tissue_name not in best_tissue_dicti:
					continue
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
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_number) + '\t' + region_name + '\t' + gene_id + '\t' + 'gene' + '\t' + str(causal_genetic_element_in_cs_boolean) + '\t' + str(gene_pip) + '\n')
	t.close()
	return

def create_file_containing_coloc_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_coloc_results_dir, pip_threshold, coloc_cs_coverage_per_high_pip_snp_output_file):
	# Open output file handle
	t = open(coloc_cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tsimulation_number\twindow_name\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

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
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_number) + '\t' + region_name + '\t' + gene_id + '\t' + 'gene' + '\t' + str(causal_genetic_element_in_cs_boolean) + '\t' + str(gene_pip) + '\n')
	t.close()
	return

def create_file_containing_focus_tg_cs_calibration_per_high_pip_genes(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, pip_threshold, focus_cs_coverage_per_high_pip_snp_output_file, intercept=False):
	# Open output file handle
	t = open(focus_cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tsimulation_number\twindow_name\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		for gene_id in [*causal_genes]:
			causal_genetic_elements[gene_id.split('_')[0]] = 1

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
			focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_focus_res.focus.tsv'
			if intercept == True:
				focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_focus_w_intercept_res.focus.tsv'
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
				else:
					used_genes[gene_id] = (pip, region)
			f.close()

			# Add together pips across tissues
			used_genes2 = {}
			for gene_tissue_id in [*used_genes]:
				gene_id = gene_tissue_id.split('_')[0]
				if gene_id not in used_genes2:
					used_genes2[gene_id] = used_genes[gene_tissue_id]
				else:
					new_pip = used_genes2[gene_id][0] + used_genes[gene_tissue_id][0]
					if new_pip > 1.0:
						new_pip = 1.0
				
					used_genes2[gene_id] = (new_pip, used_genes[gene_tissue_id][1])


			for gene_id in [*used_genes2]:
				if used_genes2[gene_id][0] < pip_threshold:
					continue
				gene_pip = used_genes2[gene_id][0]
				region_name = used_genes2[gene_id][1]
				if gene_id in causal_genetic_elements:
					causal_genetic_element_in_cs_boolean = 1
				else:
					causal_genetic_element_in_cs_boolean = 0
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_number) + '\t' + region_name + '\t' + gene_id + '\t' + 'gene' + '\t' + str(causal_genetic_element_in_cs_boolean) + '\t' + str(gene_pip) + '\n')
	t.close()
	return




def create_file_containing_focus_tg_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, pip_threshold, focus_cs_coverage_per_high_pip_snp_output_file, intercept=False):
	# Open output file handle
	t = open(focus_cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tsimulation_number\twindow_name\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
			focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_focus_res.focus.tsv'
			if intercept:
				focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_focus_w_intercept_res.focus.tsv'
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
				else:
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
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_number) + '\t' + region_name + '\t' + gene_id + '\t' + 'gene' + '\t' + str(causal_genetic_element_in_cs_boolean) + '\t' + str(gene_pip) + '\n')
	t.close()
	return

def create_file_containing_focus_cs_calibration_per_high_pip_genes(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, pip_threshold, focus_cs_coverage_per_high_pip_snp_output_file, intercept=False):
	# Open output file handle
	t = open(focus_cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tsimulation_number\twindow_name\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		for gene_id in [*causal_genes]:
			ensamble = gene_id.split('_')[0]
			causal_genetic_elements[ensamble] = 1

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
			used_genes = {}
			for tissue_number in range(10):
				focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_number) + '_focus_res.focus.tsv'
				if intercept == True:
					focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_number) + '_focus_w_intercept_res.focus.tsv'
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
					ensamble = gene_id.split('_')[0]
					if ensamble in used_genes:
						if used_genes[ensamble][0] < pip:
							used_genes[ensamble] = (pip, region)
					else:
						used_genes[ensamble] = (pip, region)
				f.close()

			for gene_id in [*used_genes]:
				if used_genes[gene_id][0] < pip_threshold:
					continue
				gene_pip = used_genes[gene_id][0]
				region_name = used_genes[gene_id][1]
				if gene_id in causal_genetic_elements:
					causal_genetic_element_in_cs_boolean = 1
				else:
					causal_genetic_element_in_cs_boolean = 0
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_number) + '\t' + region_name + '\t' + gene_id + '\t' + 'gene' + '\t' + str(causal_genetic_element_in_cs_boolean) + '\t' + str(gene_pip) + '\n')
	t.close()
	return



def create_file_containing_two_step_focus_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_two_step_tgfm_results_dir, simulated_focus_results_dir, pip_threshold, focus_cs_coverage_per_high_pip_snp_output_file, intercept=False):
	# Open output file handle
	t = open(focus_cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tsimulation_number\twindow_name\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
			used_genes = {}

			two_step_tissue_summary_file = simulated_two_step_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_pmces_uniform_two_step_tissues.txt'
			best_tissue_arr, sig_tissue_arr = extract_two_step_analyzed_tissues(two_step_tissue_summary_file)

			for tissue_name in best_tissue_arr:
				tissue_number = int(tissue_name.split('ue')[1])

				focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_number) + '_focus_res.focus.tsv'
				if intercept == True:
					focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_number) + '_focus_w_intercept_res.focus.tsv'
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
					if gene_id in used_genes:
						if used_genes[gene_id][0] < pip:
							used_genes[gene_id] = (pip, region)
					else:
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
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_number) + '\t' + region_name + '\t' + gene_id + '\t' + 'gene' + '\t' + str(causal_genetic_element_in_cs_boolean) + '\t' + str(gene_pip) + '\n')
	t.close()
	return



def create_file_containing_focus_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, pip_threshold, focus_cs_coverage_per_high_pip_snp_output_file, intercept=False):
	# Open output file handle
	t = open(focus_cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tsimulation_number\twindow_name\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
			used_genes = {}
			for tissue_number in range(10):
				focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_number) + '_focus_res.focus.tsv'
				if intercept == True:
					focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_number) + '_focus_w_intercept_res.focus.tsv'
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
					if gene_id in used_genes:
						if used_genes[gene_id][0] < pip:
							used_genes[gene_id] = (pip, region)
					else:
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
				t.write(str(eqtl_sample_size) + '\t' + str(simulation_number) + '\t' + region_name + '\t' + gene_id + '\t' + 'gene' + '\t' + str(causal_genetic_element_in_cs_boolean) + '\t' + str(gene_pip) + '\n')
	t.close()
	return

def extract_gwas_eqtl_variants(gene_tissue_pairs, simulation_number, global_simulation_name_string, simulated_gene_expression_dir):
	dicti = {}
	for gene_tissue_pair in gene_tissue_pairs:
		gene_name = gene_tissue_pair.split('_')[0]
		tissue_num = int(gene_tissue_pair.split('issue')[1])
		causal_eqtl_effects_file = simulated_gene_expression_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_' + gene_name + '_' + 'causal_eqtl_effects.npy'
		causal_eqtl_effects = np.load(causal_eqtl_effects_file)
		window_snps_file = simulated_gene_expression_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_' + gene_name + '_' + 'cis_snpids.npy'
		window_snps = np.load(window_snps_file,allow_pickle=True)[:,1]
		causal_snp_indices = causal_eqtl_effects[:,tissue_num] != 0.0
		causal_snps = window_snps[causal_snp_indices]
		for causal_snp in causal_snps:
			dicti[causal_snp] =1
	return dicti

def extract_gwas_eqtl_variants_for_each_gene(gene_tissue_pairs, simulation_number, global_simulation_name_string, simulated_gene_expression_dir):
	dicti = {}
	for gene_tissue_pair in gene_tissue_pairs:
		gene_dicti = {}
		gene_name = gene_tissue_pair.split('_')[0]
		tissue_num = int(gene_tissue_pair.split('issue')[1])
		causal_eqtl_effects_file = simulated_gene_expression_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_' + gene_name + '_' + 'causal_eqtl_effects.npy'
		causal_eqtl_effects = np.load(causal_eqtl_effects_file)
		window_snps_file = simulated_gene_expression_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_' + gene_name + '_' + 'cis_snpids.npy'
		window_snps = np.load(window_snps_file,allow_pickle=True)[:,1]
		causal_snp_indices = causal_eqtl_effects[:,tissue_num] != 0.0
		causal_snps = window_snps[causal_snp_indices]
		for causal_snp in causal_snps:
			gene_dicti[causal_snp] =1
		dicti[gene_tissue_pair] = gene_dicti
	return dicti



def create_file_containing_tgfm_cs_calibration_per_high_pip_snp_non_mediated_variants_include_gene_variants_vary_ge_h2s(global_simulation_name_string, eqtl_sample_size, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_gene_expression_dir, ge_h2s):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('GWAS_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		for ge_h2 in ge_h2s:
			temp_simulation_name_string = global_simulation_name_string.split('ss_')[0] + 'ss_' + str(100000) + '_ge_h2_' + ge_h2
			# First extract dictionary list of causal genetic elements
			causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
			causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
			causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

			gene_based_variants = extract_gwas_eqtl_variants([*causal_genes], simulation_number, temp_simulation_name_string, simulated_gene_expression_dir)
			for gene_based_variant in [*gene_based_variants]:
				causal_genetic_elements[gene_based_variant] = 1

			for ln_pi_method in ln_pi_methods:
				for twas_method in twas_methods:
					if ln_pi_method.endswith('bootstrapped') and twas_method == 'susie_pmces':
						continue
					if ln_pi_method.startswith('tglr_bootstrapped') and twas_method == 'susie_pmces':
						continue
					if ln_pi_method == 'pmces_uniform_iterative_variant_gene_prior_pip_level_pmces' and twas_method == 'susie_pmces':
						continue					
					# Credible set file for this run
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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
						cs_genetic_element_pips = []
						for element_iter, cs_prob in enumerate(cs_probs):
							if cs_prob >= pip_threshold:
								cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
								cs_genetic_element_pips.append(cs_prob)
						cs_genetic_elements = np.asarray(cs_genetic_elements)
						if len(cs_genetic_elements) == 0:
							continue


						# Extract n_causal_components
						full_window_name = str(simulation_number) + '_' + component_window_name
						#n_causal_elements = full_window_name_to_n_causal_elements[full_window_name]

						for ii,genetic_element_name in enumerate(cs_genetic_elements):
							if genetic_element_name.startswith('ENSG'):
								class_name = 'gene'
							elif genetic_element_name.startswith('rs'):
								class_name = 'variant'
							else:
								print('assumptino eroror')
								pdb.set_trace()

							genetic_element_pip = cs_genetic_element_pips[ii]
							# Check if cs contains at least one causal genetic element
							causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)

							t.write(str(ge_h2) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')
					f.close()
	t.close()
	return



def create_file_containing_tgfm_cs_calibration_per_high_pip_snp_non_mediated_variants_include_gene_variants_vary_gwas_ss(global_simulation_name_string, eqtl_sample_size, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_gene_expression_dir, gwas_sample_sizes):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('GWAS_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		for gwas_sample_size in gwas_sample_sizes:
			temp_simulation_name_string = global_simulation_name_string.split('ss_')[0] + 'ss_' + str(gwas_sample_size) + '_ge_h2_075'
			# First extract dictionary list of causal genetic elements
			causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
			causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
			causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

			gene_based_variants = extract_gwas_eqtl_variants([*causal_genes], simulation_number, temp_simulation_name_string, simulated_gene_expression_dir)
			for gene_based_variant in [*gene_based_variants]:
				causal_genetic_elements[gene_based_variant] = 1

			for ln_pi_method in ln_pi_methods:
				for twas_method in twas_methods:
					if ln_pi_method.endswith('bootstrapped') and twas_method == 'susie_pmces':
						continue
					if ln_pi_method.startswith('tglr_bootstrapped') and twas_method == 'susie_pmces':
						continue
					if ln_pi_method == 'pmces_uniform_iterative_variant_gene_prior_pip_level_pmces' and twas_method == 'susie_pmces':
						continue					
					# Credible set file for this run
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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
						cs_genetic_element_pips = []
						for element_iter, cs_prob in enumerate(cs_probs):
							if cs_prob >= pip_threshold:
								cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
								cs_genetic_element_pips.append(cs_prob)
						cs_genetic_elements = np.asarray(cs_genetic_elements)
						if len(cs_genetic_elements) == 0:
							continue


						# Extract n_causal_components
						full_window_name = str(simulation_number) + '_' + component_window_name
						#n_causal_elements = full_window_name_to_n_causal_elements[full_window_name]

						for ii,genetic_element_name in enumerate(cs_genetic_elements):
							if genetic_element_name.startswith('ENSG'):
								class_name = 'gene'
							elif genetic_element_name.startswith('rs'):
								class_name = 'variant'
							else:
								print('assumptino eroror')
								pdb.set_trace()

							genetic_element_pip = cs_genetic_element_pips[ii]
							# Check if cs contains at least one causal genetic element
							causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)

							t.write(str(gwas_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')
					f.close()
	t.close()
	return



def create_file_containing_tgfm_cs_calibration_per_high_pip_snp_non_mediated_variants_include_gene_variants(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_gene_expression_dir, n_samp='100', gene_type='component_gene'):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		gene_based_variants = extract_gwas_eqtl_variants([*causal_genes], simulation_number, global_simulation_name_string, simulated_gene_expression_dir)
		for gene_based_variant in [*gene_based_variants]:
			causal_genetic_elements[gene_based_variant] = 1

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
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_nsamp_' + str(n_samp) + '_' + gene_type + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'

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
						cs_genetic_element_pips = []
						for element_iter, cs_prob in enumerate(cs_probs):
							if cs_prob >= pip_threshold:
								cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
								cs_genetic_element_pips.append(cs_prob)
						cs_genetic_elements = np.asarray(cs_genetic_elements)
						if len(cs_genetic_elements) == 0:
							continue


						# Extract n_causal_components
						full_window_name = str(simulation_number) + '_' + component_window_name
						#n_causal_elements = full_window_name_to_n_causal_elements[full_window_name]

						for ii,genetic_element_name in enumerate(cs_genetic_elements):
							if genetic_element_name.startswith('ENSG'):
								class_name = 'gene'
							elif genetic_element_name.startswith('rs'):
								class_name = 'variant'
							else:
								print('assumptino eroror')
								pdb.set_trace()

							genetic_element_pip = cs_genetic_element_pips[ii]
							# Check if cs contains at least one causal genetic element
							causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)

							t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')
					f.close()



	t.close()
	return


def create_file_containing_tgfm_cs_calibration_per_high_pip_snp_non_mediated_variants_include_gene_variants_mixed_ss(global_simulation_name_string, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_gene_expression_dir, simulated_learned_gene_models_dir, n_samp='100', gene_type='component_gene'):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	eqtl_sample_size = 'realistic'
	eqtl_ss_categories = ['low', 'high']
	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		gene_based_variants = extract_gwas_eqtl_variants([*causal_genes], simulation_number, global_simulation_name_string, simulated_gene_expression_dir)
		for gene_based_variant in [*gene_based_variants]:
			causal_genetic_elements[gene_based_variant] = 1

		# Extract causal eqtl sample size for this simulation
		causal_eqtl_ss_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_number) + '/simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_realistic_permuted_eqtl_ss.npy'
		cur_sim_sample_sizes = np.load(causal_eqtl_ss_file)


		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size_category in eqtl_ss_categories:
			for ln_pi_method in ln_pi_methods:
				for twas_method in twas_methods:
					if ln_pi_method.endswith('bootstrapped') and twas_method == 'susie_pmces':
						continue
					if ln_pi_method.startswith('tglr_bootstrapped') and twas_method == 'susie_pmces':
						continue
					if ln_pi_method == 'pmces_uniform_iterative_variant_gene_prior_pip_level_pmces' and twas_method == 'susie_pmces':
						continue					
					# Credible set file for this run
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_nsamp_' + str(n_samp) + '_' + gene_type + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'

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
						cs_genetic_element_pips = []
						for element_iter, cs_prob in enumerate(cs_probs):
							if cs_prob >= pip_threshold:
								cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
								cs_genetic_element_pips.append(cs_prob)
						cs_genetic_elements = np.asarray(cs_genetic_elements)
						if len(cs_genetic_elements) == 0:
							continue


						# Extract n_causal_components
						full_window_name = str(simulation_number) + '_' + component_window_name
						#n_causal_elements = full_window_name_to_n_causal_elements[full_window_name]

						for ii,genetic_element_name in enumerate(cs_genetic_elements):
							if genetic_element_name.startswith('ENSG'):
								class_name = 'gene'
							elif genetic_element_name.startswith('rs'):
								class_name = 'variant'
							else:
								print('assumptino eroror')
								pdb.set_trace()

							pass_booler = True
							if class_name == 'gene':
								tissue_num = int(genetic_element_name.split('ue')[1])
								if cur_sim_sample_sizes[tissue_num] < 200 and eqtl_sample_size_category == 'high':
									pass_booler = False
								if cur_sim_sample_sizes[tissue_num] > 200 and eqtl_sample_size_category == 'low':
									pass_booler = False
							if pass_booler == False:
								continue

							genetic_element_pip = cs_genetic_element_pips[ii]
							# Check if cs contains at least one causal genetic element
							causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)

							t.write(str(eqtl_sample_size_category) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')
					f.close()



	t.close()
	return



def compute_expected_pips_from_sampler_pis(gene_alphas):
	n_bs = gene_alphas[0].shape[0]
	n_genes = gene_alphas[0].shape[1]
	LL = len(gene_alphas)


	alpha_pips = np.ones((n_bs, n_genes))

	for component_iter in range(LL):
		alpha_pips = alpha_pips*(1.0 - gene_alphas[component_iter])

	alpha_pips = 1.0 - alpha_pips

	expected_alpha_pips = np.mean(alpha_pips,axis=0)
	return expected_alpha_pips


def extract_middle_gene_agg_pips_from_tgfm_obj( tgfm_res_file):

	g = open(tgfm_res_file, "rb")
	tgfm_res = pickle.load(g)
	g.close()


	if 'expected_alpha_pips' in tgfm_res:
		alphas = tgfm_res['alpha_phis']
		middle_genes_indices = tgfm_res['middle_gene_indices']
		middle_gene_tissues = tgfm_res['genes'][middle_genes_indices]

		gene_to_indices = {}
		for ii, middle_gene_tissue in enumerate(middle_gene_tissues):
			gene_name = middle_gene_tissue.split('_')[0]
			if gene_name not in gene_to_indices:
				gene_to_indices[gene_name] = []
			gene_to_indices[gene_name].append(ii)
		n_genes = len(gene_to_indices)
		ordered_genes = [*gene_to_indices]
		LL = len(alphas)
		gene_alphas = []
		for ll in range(LL):
			new_alpha = np.zeros((alphas[ll].shape[0], n_genes))
			for ii,gene_name in enumerate(ordered_genes):
				tmper = alphas[ll][:, middle_genes_indices]
				new_alpha[:, ii] = np.sum(tmper[:, gene_to_indices[gene_name]],axis=1)
			gene_alphas.append(new_alpha)
		expected_alpha_pips = compute_expected_pips_from_sampler_pis(gene_alphas)
	else:
		alphas = tgfm_res['alpha_pip']
		pdb.set_trace()


	return np.asarray(ordered_genes), np.asarray(expected_alpha_pips)

def create_file_containing_tgfm_cs_calibration_per_high_pip_genes_vary_ge_h2s(global_simulation_name_string, eqtl_sample_size, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_input_data_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, ge_h2s):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('GE_h2\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		for ge_h2 in ge_h2s:

			temp_simulation_name_string = global_simulation_name_string.split('ss_')[0] + 'ss_' + str(100000) + '_ge_h2_' + ge_h2

			# First extract dictionary list of causal genetic elements
			causal_variant_file = simulated_trait_dir + 'simulation_new_' + str(simulation_number) + '_' + temp_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
			causal_gene_file = simulated_trait_dir + 'simulation_new_' + str(simulation_number) + '_' + temp_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
			causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)
			for gene_tissue in [*causal_genes]:
				gene_name = gene_tissue.split('_')[0]
				causal_genetic_elements[gene_name] = 1


			# Now loop through eqtl sample sizes and ln_pi methods
			for ln_pi_method in ln_pi_methods:
				for twas_method in twas_methods:
					if ln_pi_method.endswith('bootstrapped') and twas_method == 'susie_pmces':
						continue
					if ln_pi_method.startswith('tglr_bootstrapped') and twas_method == 'susie_pmces':
						continue
					if ln_pi_method == 'pmces_uniform_iterative_variant_gene_prior_pip_level_pmces' and twas_method == 'susie_pmces':
						continue					
					# Credible set file for this run
					pdb.set_trace()
					cs_file = simulated_tgfm_results_dir + 'simulation_new_' + str(simulation_number) + '_' + temp_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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
						if len(data) == 3 and data[2] == 'NA':
							continue

						#tgfm_input_file = simulated_tgfm_input_data_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string +'_' + component_window_name + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_tgfm_input_data.pkl'
						tgfm_res_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_' + component_window_name + '_results.pkl'
						
						all_cs_genetic_elements, cs_probs = extract_middle_gene_agg_pips_from_tgfm_obj(tgfm_res_file)
						#all_cs_genetic_elements = data[1].split(';')
						#cs_probs = np.asarray(data[2].split(';')).astype(float)
						cs_genetic_elements = []
						cs_genetic_element_pips = []
						for element_iter, cs_prob in enumerate(cs_probs):
							if cs_prob >= pip_threshold:
								cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
								cs_genetic_element_pips.append(cs_prob)
						cs_genetic_elements = np.asarray(cs_genetic_elements)
						if len(cs_genetic_elements) == 0:
							continue


						# Extract n_causal_components
						full_window_name = str(simulation_number) + '_' + component_window_name
						#n_causal_elements = full_window_name_to_n_causal_elements[full_window_name]

						for ii,genetic_element_name in enumerate(cs_genetic_elements):
							if genetic_element_name.startswith('ENSG'):
								class_name = 'gene'
							elif genetic_element_name.startswith('rs'):
								class_name = 'variant'
							else:
								print('assumptino eroror')
								pdb.set_trace()

							genetic_element_pip = cs_genetic_element_pips[ii]
							# Check if cs contains at least one causal genetic element
							causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)

							t.write(str(ge_h2) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')
					f.close()



	t.close()
	return


def create_file_containing_tgfm_cs_calibration_per_high_pip_genes_vary_gwas_ss(global_simulation_name_string, eqtl_sample_size, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_input_data_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, gwas_sample_sizes):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('GWAS_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		for gwas_sample_size in gwas_sample_sizes:

			temp_simulation_name_string = global_simulation_name_string.split('ss_')[0] + 'ss_' + str(gwas_sample_size) + '_ge_h2_075'

			# First extract dictionary list of causal genetic elements
			causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
			causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
			causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)
			for gene_tissue in [*causal_genes]:
				gene_name = gene_tissue.split('_')[0]
				causal_genetic_elements[gene_name] = 1


			# Now loop through eqtl sample sizes and ln_pi methods
			for ln_pi_method in ln_pi_methods:
				for twas_method in twas_methods:
					if ln_pi_method.endswith('bootstrapped') and twas_method == 'susie_pmces':
						continue
					if ln_pi_method.startswith('tglr_bootstrapped') and twas_method == 'susie_pmces':
						continue
					if ln_pi_method == 'pmces_uniform_iterative_variant_gene_prior_pip_level_pmces' and twas_method == 'susie_pmces':
						continue					
					# Credible set file for this run
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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
						if len(data) == 3 and data[2] == 'NA':
							continue

						#tgfm_input_file = simulated_tgfm_input_data_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string +'_' + component_window_name + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_tgfm_input_data.pkl'
						tgfm_res_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_' + component_window_name + '_results.pkl'
						
						all_cs_genetic_elements, cs_probs = extract_middle_gene_agg_pips_from_tgfm_obj(tgfm_res_file)
						#all_cs_genetic_elements = data[1].split(';')
						#cs_probs = np.asarray(data[2].split(';')).astype(float)
						cs_genetic_elements = []
						cs_genetic_element_pips = []
						for element_iter, cs_prob in enumerate(cs_probs):
							if cs_prob >= pip_threshold:
								cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
								cs_genetic_element_pips.append(cs_prob)
						cs_genetic_elements = np.asarray(cs_genetic_elements)
						if len(cs_genetic_elements) == 0:
							continue


						# Extract n_causal_components
						full_window_name = str(simulation_number) + '_' + component_window_name
						#n_causal_elements = full_window_name_to_n_causal_elements[full_window_name]

						for ii,genetic_element_name in enumerate(cs_genetic_elements):
							if genetic_element_name.startswith('ENSG'):
								class_name = 'gene'
							elif genetic_element_name.startswith('rs'):
								class_name = 'variant'
							else:
								print('assumptino eroror')
								pdb.set_trace()

							genetic_element_pip = cs_genetic_element_pips[ii]
							# Check if cs contains at least one causal genetic element
							causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)

							t.write(str(gwas_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')
					f.close()



	t.close()
	return

def create_file_containing_tgfm_cs_calibration_per_high_pip_genes(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_input_data_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, gene_type='component_gene',n_samp='100'):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)
		for gene_tissue in [*causal_genes]:
			gene_name = gene_tissue.split('_')[0]
			causal_genetic_elements[gene_name] = 1


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
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_nsamp_' + str(n_samp) + '_' + gene_type + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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
						if len(data) == 3 and data[2] == 'NA':
							continue

						#tgfm_input_file = simulated_tgfm_input_data_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string +'_' + component_window_name + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_tgfm_input_data.pkl'
						tgfm_res_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string+'_nsamp_' + str(n_samp) + '_' + gene_type + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_' + component_window_name + '_results.pkl'
						
						all_cs_genetic_elements, cs_probs = extract_middle_gene_agg_pips_from_tgfm_obj(tgfm_res_file)
						#all_cs_genetic_elements = data[1].split(';')
						#cs_probs = np.asarray(data[2].split(';')).astype(float)
						cs_genetic_elements = []
						cs_genetic_element_pips = []
						for element_iter, cs_prob in enumerate(cs_probs):
							if cs_prob >= pip_threshold:
								cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
								cs_genetic_element_pips.append(cs_prob)
						cs_genetic_elements = np.asarray(cs_genetic_elements)
						if len(cs_genetic_elements) == 0:
							continue


						# Extract n_causal_components
						full_window_name = str(simulation_number) + '_' + component_window_name
						#n_causal_elements = full_window_name_to_n_causal_elements[full_window_name]

						for ii,genetic_element_name in enumerate(cs_genetic_elements):
							if genetic_element_name.startswith('ENSG'):
								class_name = 'gene'
							elif genetic_element_name.startswith('rs'):
								class_name = 'variant'
							else:
								print('assumptino eroror')
								pdb.set_trace()

							genetic_element_pip = cs_genetic_element_pips[ii]
							# Check if cs contains at least one causal genetic element
							causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)

							t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')
					f.close()



	t.close()
	return

def extract_dictionary_list_of_best_tagging_causal_gene_tissue_pairs(best_tagging_file, gene_level=False):
	f = open(best_tagging_file)
	dicti = {}
	dicti2 = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count+1
			continue
		if gene_level:
			dicti[data[1].split('_')[0]] = 1
		else:
			dicti[data[1]] = 1
		dicti2[data[0]] = data[1]
	f.close()
	return dicti, dicti2

def extract_dictionary_list_of_best_tagging_causal_gene_tissue_pairs2(best_tagging_file, gene_level=False):
	f = open(best_tagging_file)
	dicti = {}
	dicti2 = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count+1
			continue
		if gene_level:
			dicti[data[1].split('_')[0]] = 1
		else:
			dicti[data[1]] = 1
		dicti2[data[0].split('_')[0]] = (data[3], data[4])
	f.close()
	return dicti, dicti2



def create_file_containing_examples_of_tgfm_prioritized_gene_tissue_from_not_best_tagging_tissue(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, missingness_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_best_tagging_gt_dir, gene_level=False):
	twas_method = twas_methods[0]
	ln_pi_method = ln_pi_methods[0]
	missingness_methods = ['no_t0']
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tsimulation_number\twindow_name\tgenetic_element_name\tPIP\tcorrelation_with_causal_gene\tbest_correlation_with_causal_gene\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file, gene_level=gene_level)

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
			best_tagging_file = simulated_best_tagging_gt_dir+ 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_no_t0_eqtl_ss_' + str(eqtl_sample_size) + '_best_tagging_gt_pairs_from_sim_causal_effect.txt'
			tagging_causal_gt_pairs, dicti2 = extract_dictionary_list_of_best_tagging_causal_gene_tissue_pairs2(best_tagging_file, gene_level=gene_level)
			for missingness_method in missingness_methods:
				# Credible set file for this run
				cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_' + missingness_method + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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
					cs_genetic_element_pips = []
					for element_iter, cs_prob in enumerate(cs_probs):
						if cs_prob >= pip_threshold:
							cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
							cs_genetic_element_pips.append(cs_prob)
					cs_genetic_elements = np.asarray(cs_genetic_elements)
					if len(cs_genetic_elements) == 0:
						continue

					# Extract n_causal_components
					full_window_name = str(simulation_number) + '_' + component_window_name
					#n_causal_elements = full_window_name_to_n_causal_elements[full_window_name]

					for ii,genetic_element_name in enumerate(cs_genetic_elements):
						if genetic_element_name.startswith('ENSG'):
							class_name = 'gene'
							if gene_level:
								genetic_element_name = genetic_element_name.split('_')[0]
						elif genetic_element_name.startswith('rs'):
							class_name = 'variant'
						else:
							print('assumptino eroror')
							pdb.set_trace()

						if class_name != 'gene':
							continue

						genetic_element_pip = cs_genetic_element_pips[ii]

						# Limit to prioritized gene-tissue pairs where we got the gene right
						if genetic_element_name.split('_')[0] not in dicti2:
							continue

						gt_pairs = np.asarray(dicti2[genetic_element_name.split('_')[0]][0].split(';'))
						gt_abs_corrs = np.asarray(dicti2[genetic_element_name.split('_')[0]][1].split(';')).astype(float)
						index = np.where(gt_pairs == genetic_element_name)[0][0]

						corr1 = gt_abs_corrs[index]
						best_corr = np.max(gt_abs_corrs)

						t.write(str(eqtl_sample_size)  + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + genetic_element_name  + '\t' + str(genetic_element_pip) + '\t' + str(corr1) + '\t' + str(best_corr) + '\n')

				f.close()

	t.close()
	return



def create_file_containing_tgfm_cs_calibration_per_high_pip_snp_in_missing_causal_tissue_sim(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, missingness_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_best_tagging_gt_dir, gene_level=False):
	twas_method = twas_methods[0]
	ln_pi_method = ln_pi_methods[0]
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tmissingness\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file, gene_level=gene_level)

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
			best_tagging_file = simulated_best_tagging_gt_dir+ 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_no_t0_eqtl_ss_' + str(eqtl_sample_size) + '_best_tagging_gt_pairs_from_sim_causal_effect.txt'
			tagging_causal_gt_pairs, dicti2 = extract_dictionary_list_of_best_tagging_causal_gene_tissue_pairs(best_tagging_file, gene_level=gene_level)
			for missingness_method in missingness_methods:
				# Credible set file for this run
				cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_' + missingness_method + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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
					cs_genetic_element_pips = []
					for element_iter, cs_prob in enumerate(cs_probs):
						if cs_prob >= pip_threshold:
							cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
							cs_genetic_element_pips.append(cs_prob)
					cs_genetic_elements = np.asarray(cs_genetic_elements)
					if len(cs_genetic_elements) == 0:
						continue

					# Extract n_causal_components
					full_window_name = str(simulation_number) + '_' + component_window_name
					#n_causal_elements = full_window_name_to_n_causal_elements[full_window_name]

					for ii,genetic_element_name in enumerate(cs_genetic_elements):
						if genetic_element_name.startswith('ENSG'):
							class_name = 'gene'
							if gene_level:
								genetic_element_name = genetic_element_name.split('_')[0]
						elif genetic_element_name.startswith('rs'):
							class_name = 'variant'
						else:
							print('assumptino eroror')
							pdb.set_trace()

						genetic_element_pip = cs_genetic_element_pips[ii]
						# Check if cs contains at least one causal genetic element
						causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)


						if missingness_method == 'no_t0' and class_name == 'gene':
							if genetic_element_name in tagging_causal_gt_pairs:
								causal_genetic_element_in_cs_boolean = 1.0
								info = dicti2[genetic_element_name.split('_')[0]]

							else:
								if genetic_element_name.split('_')[0] in dicti2:
									info = dicti2[genetic_element_name.split('_')[0]]
									#print('#############')
									#print(eqtl_sample_size)
									#print(simulation_number)
									#print(full_window_name)
									#print(genetic_element_name)
									#print(info)

						t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + missingness_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')
				f.close()



	t.close()
	return



def create_file_containing_tgfm_cs_calibration_per_high_pip_snp_in_missing_causal_tissue_sim_stratefied_by_n_causal_genetic_elements_bins(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, missingness_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_best_tagging_gt_dir, bim_file, simulated_gene_position_file, gene_level=False):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	# create gene to position mapping
	ele_to_position_mapping = {}
	gene_to_position_mapping = {}
	for ii, gene in enumerate(all_genes):
		gene_to_position_mapping[gene] = all_genes_positions[ii]
		ele_to_position_mapping[gene] = all_genes_positions[ii]
	# Create variant to position mapping
	var_to_position_mapping = {}
	for ii, variant in enumerate(all_variants):
		var_to_position_mapping[variant] = all_variants_positions[ii]
		ele_to_position_mapping[variant] = all_variants_positions[ii]


	twas_method = twas_methods[0]
	ln_pi_method = ln_pi_methods[0]
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tmissingness\tn_causal_genetic_elements\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	aa = []
	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file, gene_level=gene_level)

		# Get positions of causal genetic elements
		causal_genetic_element_positions = []
		for ele in [*causal_genetic_elements]:
			if ele.startswith('ENSG'):
				causal_genetic_element_positions.append(ele_to_position_mapping[ele.split('_')[0]])
			else:
				causal_genetic_element_positions.append(ele_to_position_mapping[ele])
		causal_genetic_element_positions = np.asarray(causal_genetic_element_positions)

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
			best_tagging_file = simulated_best_tagging_gt_dir+ 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_no_t0_eqtl_ss_' + str(eqtl_sample_size) + '_best_tagging_gt_pairs_from_sim_causal_effect.txt'
			tagging_causal_gt_pairs, dicti2 = extract_dictionary_list_of_best_tagging_causal_gene_tissue_pairs(best_tagging_file, gene_level=gene_level)
			for missingness_method in missingness_methods:
				# Credible set file for this run
				cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_' + missingness_method + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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
					cs_genetic_element_pips = []
					for element_iter, cs_prob in enumerate(cs_probs):
						if cs_prob >= pip_threshold:
							cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
							cs_genetic_element_pips.append(cs_prob)
					cs_genetic_elements = np.asarray(cs_genetic_elements)
					if len(cs_genetic_elements) == 0:
						continue

					# Extract n_causal_components
					full_window_name = str(simulation_number) + '_' + component_window_name
					window_start_index = int(component_window_name.split('_')[1])
					window_end_index = int(component_window_name.split('_')[2])
					n_causal_elements = np.sum((causal_genetic_element_positions >= window_start_index) & (causal_genetic_element_positions < window_end_index))

					#n_causal_elements = full_window_name_to_n_causal_elements[full_window_name]

					for ii,genetic_element_name in enumerate(cs_genetic_elements):
						if genetic_element_name.startswith('ENSG'):
							class_name = 'gene'
							if gene_level:
								genetic_element_name = genetic_element_name.split('_')[0]
						elif genetic_element_name.startswith('rs'):
							class_name = 'variant'
						else:
							print('assumptino eroror')
							pdb.set_trace()

						genetic_element_pip = cs_genetic_element_pips[ii]
						# Check if cs contains at least one causal genetic element
						causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)


						if missingness_method == 'no_t0' and class_name == 'gene':
							if genetic_element_name in tagging_causal_gt_pairs:
								causal_genetic_element_in_cs_boolean = 1.0
								#info = dicti2[genetic_element_name.split('_')[0]]


						t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + missingness_method + '\t' + str(n_causal_elements) + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')
				f.close()


	t.close()
	return

def create_file_containing_tgfm_cs_calibration_only_correct_gene_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file,tiss_filter=False, gene_type='component_gene'):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_new_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_new_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		for gene in np.asarray([*causal_genes]):
			causal_genetic_elements[gene.split('_')[0]] =1


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
					cs_file = simulated_tgfm_results_dir + 'simulation_new_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
					if tiss_filter:
						cs_file = simulated_tgfm_results_dir + 'simulation_new_' + str(simulation_number) + '_' + global_simulation_name_string + '_all_t' +'_' + gene_type + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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
						cs_genetic_element_pips = []
						for element_iter, cs_prob in enumerate(cs_probs):
							if cs_prob >= pip_threshold:
								cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
								cs_genetic_element_pips.append(cs_prob)
						cs_genetic_elements = np.asarray(cs_genetic_elements)
						if len(cs_genetic_elements) == 0:
							continue


						# Extract n_causal_components
						full_window_name = str(simulation_number) + '_' + component_window_name
						#n_causal_elements = full_window_name_to_n_causal_elements[full_window_name]

						for ii,genetic_element_name in enumerate(cs_genetic_elements):
							if genetic_element_name.startswith('ENSG'):
								class_name = 'gene'
								tmp_genetic_element_name = genetic_element_name.split('_')[0]
							elif genetic_element_name.startswith('rs'):
								class_name = 'variant'
								tmp_genetic_element_name = genetic_element_name
							else:
								print('assumptino eroror')
								pdb.set_trace()

							genetic_element_pip = cs_genetic_element_pips[ii]
							# Check if cs contains at least one causal genetic element
							causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([tmp_genetic_element_name], causal_genetic_elements)

							t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')
					f.close()



	t.close()
	return


def create_file_containing_tgfm_cs_calibration_correct_gene_only_per_high_pip_snp_realistic_qtl_sample_size(global_simulation_name_string, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_learned_gene_models_dir, tiss_filter=False, gene_type='component_gene'):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')
	booler_vec = []
	twas_z_vec = []

	eqtl_sample_sizes = ['low_eqtl_ss', 'high_eqtl_ss']

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		tmp_causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		causal_genetic_elements = {}
		for gene_name in [*tmp_causal_genetic_elements]:
			if gene_name.startswith('ENSG'):
				causal_genetic_elements[gene_name.split('_')[0]] = 1
			else:
				causal_genetic_elements[gene_name] = 1

		# Extract causal eqtl sample size for this simulation
		causal_eqtl_ss_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_number) + '/simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_permuted_eqtl_ss.npy'
		causal_sample_size = np.load(causal_eqtl_ss_file)[0]

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
			for ln_pi_method in ln_pi_methods:
				for twas_method in twas_methods:
					if ln_pi_method.endswith('bootstrapped') and twas_method == 'susie_pmces':
						continue
					if ln_pi_method.endswith('uniform') and twas_method == 'susie_sampler':
						continue
					if ln_pi_method.startswith('tglr_bootstrapped') and twas_method == 'susie_pmces':
						continue
					if ln_pi_method == 'pmces_uniform_iterative_variant_gene_prior_pip_level_pmces' and twas_method == 'susie_pmces':
						continue
					# Credible set file for this run
					cs_file = simulated_tgfm_results_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string +  '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
					if tiss_filter:
						cs_file = simulated_tgfm_results_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_all_t' +'_' + gene_type + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
					if eqtl_sample_size == 'low_eqtl_ss' and causal_sample_size > 200:
						continue
					if eqtl_sample_size == 'high_eqtl_ss' and causal_sample_size < 200:
						continue
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
						cs_genetic_element_pips = []
						for element_iter, cs_prob in enumerate(cs_probs):
							if cs_prob >= pip_threshold:
								cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
								cs_genetic_element_pips.append(cs_prob)
						cs_genetic_elements = np.asarray(cs_genetic_elements)
						if len(cs_genetic_elements) == 0:
							continue


						# Extract n_causal_components
						full_window_name = str(simulation_number) + '_' + component_window_name
						#n_causal_elements = full_window_name_to_n_causal_elements[full_window_name]

						for ii,genetic_element_name in enumerate(cs_genetic_elements):
							if genetic_element_name.startswith('ENSG'):
								class_name = 'gene'
								tmp_genetic_element_name = genetic_element_name.split('_')[0]
							elif genetic_element_name.startswith('rs'):
								class_name = 'variant'
								tmp_genetic_element_name = genetic_element_name
							else:
								print('assumptino eroror')
								pdb.set_trace()

							genetic_element_pip = cs_genetic_element_pips[ii]
							# Check if cs contains at least one causal genetic element
							causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([tmp_genetic_element_name], causal_genetic_elements)

							t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')
					f.close()

	t.close()
	return



def create_file_containing_ctwas_tg_cs_calibration_per_high_pip_snp_realistic_qtl_sample_size(global_simulation_name_string, simulation_runs, simulated_trait_dir, simulated_ctwas_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_learned_gene_models_dir, gene_type='component_gene'):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')
	booler_vec = []
	twas_z_vec = []

	ln_pi_method='ctwas'
	twas_method='lasso'

	#eqtl_sample_sizes = ['low_eqtl_ss', 'high_eqtl_ss']
	# First loop through simulations
	for simulation_number in simulation_runs:

		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		# Extract causal eqtl sample size for this simulation
		causal_eqtl_ss_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_number) + '/simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_permuted_eqtl_ss.npy'
		causal_sample_size = np.load(causal_eqtl_ss_file)[0]

		for eqtl_sample_size in eqtl_sample_sizes:
			if eqtl_sample_size == 'low_eqtl_ss' and causal_sample_size > 200:
				continue
			if eqtl_sample_size == 'high_eqtl_ss' and causal_sample_size < 200:
				continue
			# Get list of used genetic elements
			used_ge = {}

			ctwas_result_file = simulated_ctwas_results_dir + 'simulation_new_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_realistic_qtl_ss/testing.susieIrss.txt'
			if os.path.isfile(ctwas_result_file) == False:
				print(simulation_number)
				continue
			# Stream ctwas results file
			head_count = 0
			f = open(ctwas_result_file)
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				genetic_element_name = data[1]
				genetic_element_pip = float(data[7])

				if genetic_element_name in used_ge:
					print('assumption eroror')
					pdb.set_trace()
				used_ge[genetic_element_name] = genetic_element_pip

				if genetic_element_pip < pip_threshold:
					continue

				if genetic_element_name.startswith('ENSG'):
					class_name = 'gene'
					# Quick fix from 'tissue_0' to 'tissue0' for exmample
					name_info = genetic_element_name.split('_')
					if len(name_info) != 3:
						print('assumtion erorro')
						pdb.set_trace()
					genetic_element_name = name_info[0] + '_' + name_info[1] + name_info[2]
				elif genetic_element_name.startswith('rs'):
					class_name = 'variant'
				else:
					print('assumptino eroror')
					pdb.set_trace()

				# Placeholder
				component_window_name = 'null'
				component_num = 'null'

				# Check if cs contains at least one causal genetic element
				causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)

				t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')

			f.close()

def create_file_containing_two_step_ctwas_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_two_step_tgfm_results_dir, simulated_ctwas_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_learned_gene_models_dir, gene_type='component_gene'):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')
	booler_vec = []
	twas_z_vec = []

	ln_pi_method='ctwas'
	twas_method='lasso'

	# First loop through simulations
	for simulation_number in simulation_runs:

		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		# Extract causal eqtl sample size for this simulation
		#causal_eqtl_ss_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_number) + '/simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_permuted_eqtl_ss.npy'
		#causal_sample_size = np.load(causal_eqtl_ss_file)[0]

		for eqtl_sample_size in eqtl_sample_sizes:
			# Get list of used genetic elements
			used_ge = {}

			two_step_tissue_summary_file = simulated_two_step_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_pmces_uniform_two_step_tissues.txt'

			best_tissue_arr, sig_tissue_arr = extract_two_step_analyzed_tissues(two_step_tissue_summary_file)

			for tissue_name in best_tissue_arr:
				tissue_number = tissue_name.split('ue')[1]
				ctwas_result_file = simulated_ctwas_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_number) + '/testing.susieIrss.txt'
				if os.path.isfile(ctwas_result_file) == False:
					continue
				# Stream ctwas results file
				head_count = 0
				f = open(ctwas_result_file)
				for line in f:
					line = line.rstrip()
					data = line.split('\t')
					if head_count == 0:
						head_count = head_count + 1
						continue
					genetic_element_name = data[1]
					genetic_element_pip = float(data[7])

					if genetic_element_name in used_ge and genetic_element_name.startswith('rs') == False:
						print('assumption eroror')
						pdb.set_trace()
					used_ge[genetic_element_name] = genetic_element_pip

					if genetic_element_pip < pip_threshold:
						continue

					if genetic_element_name.startswith('ENSG'):
						class_name = 'gene'
						# Quick fix from 'tissue_0' to 'tissue0' for exmample
						name_info = genetic_element_name.split('_')
						if len(name_info) != 3:
							print('assumtion erorro')
							pdb.set_trace()
						genetic_element_name = name_info[0] + '_' + name_info[1] + name_info[2]
					elif genetic_element_name.startswith('rs'):
						class_name = 'variant'
					else:
						print('assumptino eroror')
						pdb.set_trace()

					# Placeholder
					component_window_name = 'null'
					component_num = 'null'

					# Check if cs contains at least one causal genetic element
					causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)

					t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')

				f.close()



def create_file_containing_ctwas_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_ctwas_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_learned_gene_models_dir, gene_type='component_gene'):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')
	booler_vec = []
	twas_z_vec = []

	ln_pi_method='ctwas'
	twas_method='lasso'

	# First loop through simulations
	for simulation_number in simulation_runs:

		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		# Extract causal eqtl sample size for this simulation
		#causal_eqtl_ss_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_number) + '/simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_permuted_eqtl_ss.npy'
		#causal_sample_size = np.load(causal_eqtl_ss_file)[0]

		for eqtl_sample_size in eqtl_sample_sizes:
			# Get list of used genetic elements
			used_ge = {}

			for tissue_number in range(10):
				ctwas_result_file = simulated_ctwas_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_number) + '/testing.susieIrss.txt'
				if os.path.isfile(ctwas_result_file) == False:
					continue
				# Stream ctwas results file
				head_count = 0
				f = open(ctwas_result_file)
				for line in f:
					line = line.rstrip()
					data = line.split('\t')
					if head_count == 0:
						head_count = head_count + 1
						continue
					genetic_element_name = data[1]
					genetic_element_pip = float(data[7])

					if genetic_element_name in used_ge and genetic_element_name.startswith('rs') == False:
						print('assumption eroror')
						pdb.set_trace()
					used_ge[genetic_element_name] = genetic_element_pip

					if genetic_element_pip < pip_threshold:
						continue

					if genetic_element_name.startswith('ENSG'):
						class_name = 'gene'
						# Quick fix from 'tissue_0' to 'tissue0' for exmample
						name_info = genetic_element_name.split('_')
						if len(name_info) != 3:
							print('assumtion erorro')
							pdb.set_trace()
						genetic_element_name = name_info[0] + '_' + name_info[1] + name_info[2]
					elif genetic_element_name.startswith('rs'):
						class_name = 'variant'
					else:
						print('assumptino eroror')
						pdb.set_trace()

					# Placeholder
					component_window_name = 'null'
					component_num = 'null'

					# Check if cs contains at least one causal genetic element
					causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)

					t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')

				f.close()

def create_file_containing_ctwas_tg_cs_calibration_per_high_pip_gene(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_ctwas_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_learned_gene_models_dir, gene_type='component_gene'):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	booler_vec = []
	twas_z_vec = []

	ln_pi_method='ctwas'
	twas_method='lasso'

	# First loop through simulations
	for simulation_number in simulation_runs:

		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		for gene_id in [*causal_genes]:
			ensamble = gene_id.split('_')[0]
			causal_genetic_elements[ensamble] = 1

		# Extract causal eqtl sample size for this simulation
		#causal_eqtl_ss_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_number) + '/simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_permuted_eqtl_ss.npy'
		#causal_sample_size = np.load(causal_eqtl_ss_file)[0]

		for eqtl_sample_size in eqtl_sample_sizes:
			# Get list of used genetic elements
			used_ge = {}
			used_genes = {}

			ctwas_result_file = simulated_ctwas_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '/testing.susieIrss.txt'
			if os.path.isfile(ctwas_result_file) == False:
				continue
			# Stream ctwas results file
			head_count = 0
			f = open(ctwas_result_file)
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				genetic_element_name = data[1]
				pip = float(data[7])

				if genetic_element_name.startswith('ENSG'):
					if genetic_element_name in used_ge:
						print('assumption eororor')
					used_ge[genetic_element_name] = 1
					class_name = 'gene'
					# Quick fix from 'tissue_0' to 'tissue0' for exmample
					name_info = genetic_element_name.split('_')
					ensamble_id = name_info[0]
					if ensamble_id not in used_genes:
						used_genes[ensamble_id] = pip
					else:
						old_pip = used_genes[ensamble_id]
						new_pip = old_pip + pip
						if new_pip > 1.0:
							new_pip = 1.0
						used_genes[ensamble_id] = new_pip
				else:
					continue
			f.close()
			for gene_id in [*used_genes]:
				if used_genes[gene_id] < pip_threshold:
					continue
				gene_pip = used_genes[gene_id]
				if gene_id in causal_genetic_elements:
					causal_genetic_element_in_cs_boolean = 1
				else:
					causal_genetic_element_in_cs_boolean = 0
				region_name = 'null'
				component_number = 'null'
				t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + region_name  + '\t' + component_number + '\t' + gene_id + '\t' + 'gene' + '\t' + str(causal_genetic_element_in_cs_boolean) + '\t' + str(gene_pip) + '\n')
	t.close()
	return



def create_file_containing_ctwas_cs_calibration_per_high_pip_gene(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_ctwas_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_learned_gene_models_dir, gene_type='component_gene'):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')

	booler_vec = []
	twas_z_vec = []

	ln_pi_method='ctwas'
	twas_method='lasso'

	# First loop through simulations
	for simulation_number in simulation_runs:

		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		for gene_id in [*causal_genes]:
			ensamble = gene_id.split('_')[0]
			causal_genetic_elements[ensamble] = 1

		# Extract causal eqtl sample size for this simulation
		#causal_eqtl_ss_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_number) + '/simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_permuted_eqtl_ss.npy'
		#causal_sample_size = np.load(causal_eqtl_ss_file)[0]

		for eqtl_sample_size in eqtl_sample_sizes:
			# Get list of used genetic elements
			used_ge = {}
			used_genes = {}

			for tissue_number in range(10):
				ctwas_result_file = simulated_ctwas_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_number) + '/testing.susieIrss.txt'
				if os.path.isfile(ctwas_result_file) == False:
					continue
				# Stream ctwas results file
				head_count = 0
				f = open(ctwas_result_file)
				for line in f:
					line = line.rstrip()
					data = line.split('\t')
					if head_count == 0:
						head_count = head_count + 1
						continue
					genetic_element_name = data[1]
					pip = float(data[7])

					if genetic_element_name.startswith('ENSG'):
						class_name = 'gene'
						# Quick fix from 'tissue_0' to 'tissue0' for exmample
						name_info = genetic_element_name.split('_')
						ensamble_id = name_info[0]
						if ensamble_id not in used_genes:
							used_genes[ensamble_id] = pip
						else:
							old_pip = used_genes[ensamble_id]
							if pip > old_pip:
								used_genes[ensamble_id] = pip
					else:
						continue
				f.close()
			for gene_id in [*used_genes]:
				if used_genes[gene_id] < pip_threshold:
					continue
				gene_pip = used_genes[gene_id]
				if gene_id in causal_genetic_elements:
					causal_genetic_element_in_cs_boolean = 1
				else:
					causal_genetic_element_in_cs_boolean = 0
				region_name = 'null'
				component_number = 'null'
				t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + region_name  + '\t' + component_number + '\t' + gene_id + '\t' + 'gene' + '\t' + str(causal_genetic_element_in_cs_boolean) + '\t' + str(gene_pip) + '\n')
	t.close()
	return

def create_file_containing_ctwas_tg_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_ctwas_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_learned_gene_models_dir, gene_type='component_gene'):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')
	booler_vec = []
	twas_z_vec = []

	ln_pi_method='ctwas'
	twas_method='lasso'

	# First loop through simulations
	for simulation_number in simulation_runs:

		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		# Extract causal eqtl sample size for this simulation
		#causal_eqtl_ss_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_number) + '/simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_permuted_eqtl_ss.npy'
		#causal_sample_size = np.load(causal_eqtl_ss_file)[0]

		for eqtl_sample_size in eqtl_sample_sizes:
			# Get list of used genetic elements
			used_ge = {}

			ctwas_result_file = simulated_ctwas_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '/testing.susieIrss.txt'
			if os.path.isfile(ctwas_result_file) == False:
				continue
			# Stream ctwas results file
			head_count = 0
			f = open(ctwas_result_file)
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				genetic_element_name = data[1]
				genetic_element_pip = float(data[7])

				if genetic_element_name in used_ge:
					print('assumption eroror')
					pdb.set_trace()
				used_ge[genetic_element_name] = genetic_element_pip

				if genetic_element_pip < pip_threshold:
					continue

				if genetic_element_name.startswith('ENSG'):
					class_name = 'gene'
					# Quick fix from 'tissue_0' to 'tissue0' for exmample
					name_info = genetic_element_name.split('_')
					if len(name_info) != 3:
						print('assumtion erorro')
						pdb.set_trace()
					genetic_element_name = name_info[0] + '_' + name_info[1] + name_info[2]
				elif genetic_element_name.startswith('rs'):
					class_name = 'variant'
				else:
					print('assumptino eroror')
					pdb.set_trace()

				# Placeholder
				component_window_name = 'null'
				component_num = 'null'

				# Check if cs contains at least one causal genetic element
				causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)

				t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')

			f.close()

def create_file_containing_tgfm_cs_calibration_per_high_pip_snp_realistic_qtl_sample_size(global_simulation_name_string, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_learned_gene_models_dir, tiss_filter=False, gene_type='component_gene'):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')
	booler_vec = []
	twas_z_vec = []

	eqtl_sample_sizes = ['low_eqtl_ss', 'high_eqtl_ss']

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		# Extract causal eqtl sample size for this simulation
		causal_eqtl_ss_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_number) + '/simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_permuted_eqtl_ss.npy'
		causal_sample_size = np.load(causal_eqtl_ss_file)[0]

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
					cs_file = simulated_tgfm_results_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string +  '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
					if tiss_filter:
						cs_file = simulated_tgfm_results_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_all_t' +'_' + gene_type + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
					if eqtl_sample_size == 'low_eqtl_ss' and causal_sample_size > 200:
						continue
					if eqtl_sample_size == 'high_eqtl_ss' and causal_sample_size < 200:
						continue
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
						cs_genetic_element_pips = []
						for element_iter, cs_prob in enumerate(cs_probs):
							if cs_prob >= pip_threshold:
								cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
								cs_genetic_element_pips.append(cs_prob)
						cs_genetic_elements = np.asarray(cs_genetic_elements)
						if len(cs_genetic_elements) == 0:
							continue


						# Extract n_causal_components
						full_window_name = str(simulation_number) + '_' + component_window_name
						#n_causal_elements = full_window_name_to_n_causal_elements[full_window_name]

						for ii,genetic_element_name in enumerate(cs_genetic_elements):
							if genetic_element_name.startswith('ENSG'):
								class_name = 'gene'
							elif genetic_element_name.startswith('rs'):
								class_name = 'variant'
							else:
								print('assumptino eroror')
								pdb.set_trace()

							genetic_element_pip = cs_genetic_element_pips[ii]
							# Check if cs contains at least one causal genetic element
							causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)

							t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')
					f.close()

	t.close()
	return


def create_file_containing_tgfm_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, n_samp='100', gene_type='component_gene'):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')
	booler_vec = []
	twas_z_vec = []

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
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_nsamp_' + str(n_samp) + '_' + gene_type + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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
						cs_genetic_element_pips = []
						for element_iter, cs_prob in enumerate(cs_probs):
							if cs_prob >= pip_threshold:
								cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
								cs_genetic_element_pips.append(cs_prob)
						cs_genetic_elements = np.asarray(cs_genetic_elements)
						if len(cs_genetic_elements) == 0:
							continue


						# Extract n_causal_components
						full_window_name = str(simulation_number) + '_' + component_window_name
						#n_causal_elements = full_window_name_to_n_causal_elements[full_window_name]

						for ii,genetic_element_name in enumerate(cs_genetic_elements):
							if genetic_element_name.startswith('ENSG'):
								class_name = 'gene'
							elif genetic_element_name.startswith('rs'):
								class_name = 'variant'
							else:
								print('assumptino eroror')
								pdb.set_trace()

							genetic_element_pip = cs_genetic_element_pips[ii]
							# Check if cs contains at least one causal genetic element
							causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)

							t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')
					f.close()

	t.close()
	return

def extract_two_step_analyzed_tissues(two_step_tissue_summary_file,sig_thresh=0.05):
	f = open(two_step_tissue_summary_file)
	best_tissue_arr = []
	sig_tissue_arr = []
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		counter = counter + 1
		if counter == 1:
			best_tissue_arr.append(data[0])
		if float(data[1]) <= sig_thresh:
			sig_tissue_arr.append(data[0])
	f.close()
	best_tissue_arr = np.asarray(best_tissue_arr)
	sig_tissue_arr = np.asarray(sig_tissue_arr)

	return best_tissue_arr, sig_tissue_arr



def create_file_containing_two_step_tgfm_cs_calibration_per_high_pip_snp(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, n_samp='100', gene_type='component_gene'):
	# Open output file handle
	t = open(cs_coverage_per_high_pip_snp_output_file,'w')
	t.write('eQTL_sample_size\ttwo_step_tissue_method\tln_pi_method\ttwas_method\tsimulation_number\twindow_name\tcomponent_number\tgenetic_element_name\tgenetic_element_class\tcausal_genetic_element_in_cs\tPIP\n')
	booler_vec = []
	twas_z_vec = []

	ln_pi_method = 'susie_sampler'
	twas_method = 'iterative'
	two_step_tissue_methods = ['best_tissue', 'significant_tissues']

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
			two_step_tissue_summary_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_pmces_uniform_two_step_tissues.txt'

			best_tissue_arr, sig_tissue_arr = extract_two_step_analyzed_tissues(two_step_tissue_summary_file)
			
			for two_step_tissue_method in two_step_tissue_methods:
				if two_step_tissue_method == 'best_tissue':
					tissue_vec = np.copy(best_tissue_arr)
				elif two_step_tissue_method == 'significant_tissues':
					tissue_vec = np.copy(sig_tissue_arr)
				else:
					print('tissue method assumption eroror')
					pdb.set_trace()
				for tissue_name in tissue_vec:
					# Credible set file for this run
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_sampler_iterative_two_step_' + tissue_name + '_tgfm_pip_summary.txt'
					if os.path.isfile(cs_file) == False:
						print('skipped')
						continue
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
						cs_genetic_element_pips = []
						for element_iter, cs_prob in enumerate(cs_probs):
							if cs_prob >= pip_threshold:
								cs_genetic_elements.append(all_cs_genetic_elements[element_iter])
								cs_genetic_element_pips.append(cs_prob)
						cs_genetic_elements = np.asarray(cs_genetic_elements)
						if len(cs_genetic_elements) == 0:
							continue


						# Extract n_causal_components
						full_window_name = str(simulation_number) + '_' + component_window_name
						#n_causal_elements = full_window_name_to_n_causal_elements[full_window_name]

						for ii,genetic_element_name in enumerate(cs_genetic_elements):
							if genetic_element_name.startswith('ENSG'):
								class_name = 'gene'
							elif genetic_element_name.startswith('rs'):
								class_name = 'variant'
							else:
								print('assumptino eroror')
								pdb.set_trace()

							genetic_element_pip = cs_genetic_element_pips[ii]
							# Check if cs contains at least one causal genetic element
							causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs([genetic_element_name], causal_genetic_elements)

							t.write(str(eqtl_sample_size) + '\t' + two_step_tissue_method + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + genetic_element_name + '\t' + class_name  + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\t' + str(genetic_element_pip) + '\n')
					f.close()

	t.close()
	return



def create_file_containing_averaged_focus_high_pip_calibration(focus_cs_coverage_per_high_pip_snp_output_file, focus_cs_high_pip_coverage_output_file, eqtl_sample_sizes):
	per_component_coverage = np.loadtxt(focus_cs_coverage_per_high_pip_snp_output_file, dtype=str)[1:,:]
	t = open(focus_cs_high_pip_coverage_output_file,'w')
	t.write('eQTL_sample_size\tgenetic_element_class\tn_detected_elements\tcoverage\tcoverage_lb\tcoverage_ub\texpected_coverage\n')

	for eqtl_sample_size in eqtl_sample_sizes:

		subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size))
		components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
		prop = np.sum(components_covered)/len(components_covered)
		prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
		prop_lb = prop - (1.96*prop_se)
		prop_ub = prop + (1.96*prop_se)
		n_elements = np.sum(subset_indices)
		expected_coverage = np.mean((per_component_coverage[subset_indices,:][:,-1]).astype(float))
		t.write(str(eqtl_sample_size) + '\t' + 'gene' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')

	t.close()
	return

def create_file_containing_averaged_tgfm_high_pip_calibration_vary_ge_h2s(cs_coverage_per_component_output_file, cs_coverage_output_file, eqtl_sample_size, ge_h2s, ln_pi_methods, twas_methods):
	per_component_coverage = np.loadtxt(cs_coverage_per_component_output_file, dtype=str)[1:,:]
	t = open(cs_coverage_output_file,'w')
	t.write('eQTL_sample_size\tGE_h2\tln_pi_method\ttwas_method\tgenetic_element_class\tn_detected_elements\tcoverage\tcoverage_lb\tcoverage_ub\texpected_coverage\n')

	for ge_h2 in ge_h2s:
		for ln_pi_method in ln_pi_methods:
			for twas_method in twas_methods:
				subset_indices = (per_component_coverage[:,0] == str(ge_h2)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method)
				components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
				prop = np.sum(components_covered)/len(components_covered)
				prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
				prop_lb = prop - (1.96*prop_se)
				prop_ub = prop + (1.96*prop_se)
				n_elements = np.sum(subset_indices)
				expected_coverage = np.mean(per_component_coverage[subset_indices,9].astype(float))
				t.write(str(eqtl_sample_size) + '\t' + str(ge_h2) + '\t' + ln_pi_method +'\t' + twas_method  + '\t' + 'all' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')

				subset_indices = (per_component_coverage[:,0] == str(ge_h2)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method) & (per_component_coverage[:, 7] == 'variant')
				components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
				prop = np.sum(components_covered)/len(components_covered)
				prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
				prop_lb = prop - (1.96*prop_se)
				prop_ub = prop + (1.96*prop_se)
				n_elements = np.sum(subset_indices)
				expected_coverage = np.mean(per_component_coverage[subset_indices,9].astype(float))
				t.write(str(eqtl_sample_size)+ '\t' + str(ge_h2) + '\t' + ln_pi_method +'\t' + twas_method  + '\t' + 'variant' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')

				subset_indices = (per_component_coverage[:,0] == str(ge_h2)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method) & (per_component_coverage[:, 7] == 'gene')
				components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
				prop = np.sum(components_covered)/len(components_covered)
				prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
				prop_lb = prop - (1.96*prop_se)
				prop_ub = prop + (1.96*prop_se)
				n_elements = np.sum(subset_indices)	
				expected_coverage = np.mean(per_component_coverage[subset_indices,9].astype(float))		
				t.write(str(eqtl_sample_size)+ '\t' + str(ge_h2) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + 'gene' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')
	t.close()
	return


def create_file_containing_averaged_tgfm_high_pip_calibration_vary_gwas_ss(cs_coverage_per_component_output_file, cs_coverage_output_file, eqtl_sample_size, gwas_sample_sizes, ln_pi_methods, twas_methods):
	per_component_coverage = np.loadtxt(cs_coverage_per_component_output_file, dtype=str)[1:,:]
	t = open(cs_coverage_output_file,'w')
	t.write('eQTL_sample_size\tGWAS_sample_size\tln_pi_method\ttwas_method\tgenetic_element_class\tn_detected_elements\tcoverage\tcoverage_lb\tcoverage_ub\texpected_coverage\n')

	for gwas_sample in gwas_sample_sizes:
		for ln_pi_method in ln_pi_methods:
			for twas_method in twas_methods:
				subset_indices = (per_component_coverage[:,0] == str(gwas_sample)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method)
				components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
				prop = np.sum(components_covered)/len(components_covered)
				prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
				prop_lb = prop - (1.96*prop_se)
				prop_ub = prop + (1.96*prop_se)
				n_elements = np.sum(subset_indices)
				expected_coverage = np.mean(per_component_coverage[subset_indices,9].astype(float))
				t.write(str(eqtl_sample_size) + '\t' + str(gwas_sample) + '\t' + ln_pi_method +'\t' + twas_method  + '\t' + 'all' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')

				subset_indices = (per_component_coverage[:,0] == str(gwas_sample)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method) & (per_component_coverage[:, 7] == 'variant')
				components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
				prop = np.sum(components_covered)/len(components_covered)
				prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
				prop_lb = prop - (1.96*prop_se)
				prop_ub = prop + (1.96*prop_se)
				n_elements = np.sum(subset_indices)
				expected_coverage = np.mean(per_component_coverage[subset_indices,9].astype(float))
				t.write(str(eqtl_sample_size)+ '\t' + str(gwas_sample) + '\t' + ln_pi_method +'\t' + twas_method  + '\t' + 'variant' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')

				subset_indices = (per_component_coverage[:,0] == str(gwas_sample)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method) & (per_component_coverage[:, 7] == 'gene')
				components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
				prop = np.sum(components_covered)/len(components_covered)
				prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
				prop_lb = prop - (1.96*prop_se)
				prop_ub = prop + (1.96*prop_se)
				n_elements = np.sum(subset_indices)	
				expected_coverage = np.mean(per_component_coverage[subset_indices,9].astype(float))		
				t.write(str(eqtl_sample_size)+ '\t' + str(gwas_sample) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + 'gene' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')
	t.close()
	return

def create_file_containing_averaged_tgfm_high_pip_calibration_in_missing_causal_tissue_sim_stratefied_by_n_causal_genetic_elements_bins(cs_coverage_per_component_output_file, cs_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods, missingness_methods, n_causal_genetic_element_bins):
	per_component_coverage = np.loadtxt(cs_coverage_per_component_output_file, dtype=str)[1:,:]
	t = open(cs_coverage_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tmissingness_method\tn_causal_ele_bin\tgenetic_element_class\tn_detected_elements\tcoverage\tcoverage_lb\tcoverage_ub\texpected_coverage\n')

	for eqtl_sample_size in eqtl_sample_sizes:
		for ln_pi_method in ln_pi_methods:
			for twas_method in twas_methods:
				for missingness_method in missingness_methods:
					for n_causal_element_bin in n_causal_genetic_element_bins:
						n_causal_ele_lower = int(n_causal_element_bin.split('_')[0])
						n_causal_ele_upper = int(n_causal_element_bin.split('_')[1])

						subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method) & (per_component_coverage[:, 9] == 'variant')& (per_component_coverage[:,3] == missingness_method) & (per_component_coverage[:,4].astype(float) >= n_causal_ele_lower) & (per_component_coverage[:,4].astype(float) <= n_causal_ele_upper)
						components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
						prop = np.sum(components_covered)/len(components_covered)
						prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
						prop_lb = prop - (1.96*prop_se)
						prop_ub = prop + (1.96*prop_se)
						n_elements = np.sum(subset_indices)
						expected_coverage = np.mean(per_component_coverage[subset_indices,11].astype(float))
						t.write(str(eqtl_sample_size) + '\t' + ln_pi_method +'\t' + twas_method + '\t' + missingness_method + '\t' + n_causal_element_bin  + '\t' + 'variant' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')

						subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method) & (per_component_coverage[:, 9] == 'gene')& (per_component_coverage[:,3] == missingness_method) & (per_component_coverage[:,4].astype(float) >= n_causal_ele_lower) & (per_component_coverage[:,4].astype(float) <= n_causal_ele_upper)
						components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
						prop = np.sum(components_covered)/len(components_covered)
						prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
						prop_lb = prop - (1.96*prop_se)
						prop_ub = prop + (1.96*prop_se)
						n_elements = np.sum(subset_indices)	
						expected_coverage = np.mean(per_component_coverage[subset_indices,11].astype(float))		
						t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + missingness_method + '\t' + n_causal_element_bin + '\t' + 'gene' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')
	t.close()
	return



def create_file_containing_averaged_tgfm_high_pip_calibration_in_missing_causal_tissue_sim(cs_coverage_per_component_output_file, cs_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods, missingness_methods):
	per_component_coverage = np.loadtxt(cs_coverage_per_component_output_file, dtype=str)[1:,:]
	t = open(cs_coverage_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tmissingness_method\tgenetic_element_class\tn_detected_elements\tcoverage\tcoverage_lb\tcoverage_ub\texpected_coverage\n')

	for eqtl_sample_size in eqtl_sample_sizes:
		for ln_pi_method in ln_pi_methods:
			for twas_method in twas_methods:
				for missingness_method in missingness_methods:
					subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method) & (per_component_coverage[:,3] == missingness_method)
					components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
					prop = np.sum(components_covered)/len(components_covered)
					prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
					prop_lb = prop - (1.96*prop_se)
					prop_ub = prop + (1.96*prop_se)
					n_elements = np.sum(subset_indices)
					expected_coverage = np.mean(per_component_coverage[subset_indices,10].astype(float))
					#t.write(str(eqtl_sample_size) + '\t' + ln_pi_method +'\t' + twas_method  + '\t' + 'all' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')

					subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method) & (per_component_coverage[:, 8] == 'variant')& (per_component_coverage[:,3] == missingness_method)
					components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
					prop = np.sum(components_covered)/len(components_covered)
					prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
					prop_lb = prop - (1.96*prop_se)
					prop_ub = prop + (1.96*prop_se)
					n_elements = np.sum(subset_indices)
					expected_coverage = np.mean(per_component_coverage[subset_indices,10].astype(float))
					t.write(str(eqtl_sample_size) + '\t' + ln_pi_method +'\t' + twas_method + '\t' + missingness_method  + '\t' + 'variant' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')

					subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method) & (per_component_coverage[:, 8] == 'gene')& (per_component_coverage[:,3] == missingness_method)
					components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
					prop = np.sum(components_covered)/len(components_covered)
					prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
					prop_lb = prop - (1.96*prop_se)
					prop_ub = prop + (1.96*prop_se)
					n_elements = np.sum(subset_indices)	
					expected_coverage = np.mean(per_component_coverage[subset_indices,10].astype(float))		
					t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + missingness_method + '\t' + 'gene' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')
	t.close()
	return


def create_file_containing_averaged_tgfm_high_pip_calibration_two_step(cs_coverage_per_component_output_file, cs_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods, two_step_methods):
	per_component_coverage = np.loadtxt(cs_coverage_per_component_output_file, dtype=str)[1:,:]
	t = open(cs_coverage_output_file,'w')
	t.write('eQTL_sample_size\ttwo_step_tissue_method\tln_pi_method\ttwas_method\tgenetic_element_class\tn_detected_elements\tcoverage\tcoverage_lb\tcoverage_ub\texpected_coverage\n')

	for eqtl_sample_size in eqtl_sample_sizes:
		for ln_pi_method in ln_pi_methods:
			for twas_method in twas_methods:
				for two_step_method in two_step_methods:
					subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size)) & (per_component_coverage[:,1] == two_step_method) & (per_component_coverage[:,2] == ln_pi_method) & (per_component_coverage[:,3] == twas_method)
					components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
					prop = np.sum(components_covered)/len(components_covered)
					prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
					prop_lb = prop - (1.96*prop_se)
					prop_ub = prop + (1.96*prop_se)
					n_elements = np.sum(subset_indices)
					expected_coverage = np.mean(per_component_coverage[subset_indices,10].astype(float))
					t.write(str(eqtl_sample_size) + '\t' + two_step_method + '\t' + ln_pi_method +'\t' + twas_method  + '\t' + 'all' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')

					subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size))& (per_component_coverage[:,1] == two_step_method) & (per_component_coverage[:,2] == ln_pi_method) & (per_component_coverage[:,3] == twas_method) & (per_component_coverage[:, 8] == 'variant')
					components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
					prop = np.sum(components_covered)/len(components_covered)
					prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
					prop_lb = prop - (1.96*prop_se)
					prop_ub = prop + (1.96*prop_se)
					n_elements = np.sum(subset_indices)
					expected_coverage = np.mean(per_component_coverage[subset_indices,10].astype(float))
					t.write(str(eqtl_sample_size)+ '\t' + two_step_method  + '\t' + ln_pi_method +'\t' + twas_method  + '\t' + 'variant' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')

					subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size))& (per_component_coverage[:,1] == two_step_method) & (per_component_coverage[:,2] == ln_pi_method) & (per_component_coverage[:,3] == twas_method) & (per_component_coverage[:, 8] == 'gene')
					components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
					prop = np.sum(components_covered)/len(components_covered)
					prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
					prop_lb = prop - (1.96*prop_se)
					prop_ub = prop + (1.96*prop_se)
					n_elements = np.sum(subset_indices)	
					expected_coverage = np.mean(per_component_coverage[subset_indices,10].astype(float))		
					t.write(str(eqtl_sample_size) + '\t' + two_step_method + '\t' + ln_pi_method + '\t' + twas_method + '\t' + 'gene' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')	

	t.close()
	return




def create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_component_output_file, cs_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods, agg_eqtl_ss=False):
	per_component_coverage = np.loadtxt(cs_coverage_per_component_output_file, dtype=str)[1:,:]
	t = open(cs_coverage_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\ttwas_method\tgenetic_element_class\tn_detected_elements\tcoverage\tcoverage_lb\tcoverage_ub\texpected_coverage\n')

	for eqtl_sample_size in eqtl_sample_sizes:
		for ln_pi_method in ln_pi_methods:
			for twas_method in twas_methods:
				subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method)
				components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
				prop = np.sum(components_covered)/len(components_covered)
				prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
				prop_lb = prop - (1.96*prop_se)
				prop_ub = prop + (1.96*prop_se)
				n_elements = np.sum(subset_indices)
				expected_coverage = np.mean(per_component_coverage[subset_indices,9].astype(float))
				t.write(str(eqtl_sample_size) + '\t' + ln_pi_method +'\t' + twas_method  + '\t' + 'all' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')

				subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method) & (per_component_coverage[:, 7] == 'variant')
				components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
				prop = np.sum(components_covered)/len(components_covered)
				prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
				prop_lb = prop - (1.96*prop_se)
				prop_ub = prop + (1.96*prop_se)
				n_elements = np.sum(subset_indices)
				expected_coverage = np.mean(per_component_coverage[subset_indices,9].astype(float))
				t.write(str(eqtl_sample_size) + '\t' + ln_pi_method +'\t' + twas_method  + '\t' + 'variant' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')

				subset_indices = (per_component_coverage[:,0] == str(eqtl_sample_size)) & (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method) & (per_component_coverage[:, 7] == 'gene')
				components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
				prop = np.sum(components_covered)/len(components_covered)
				prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
				prop_lb = prop - (1.96*prop_se)
				prop_ub = prop + (1.96*prop_se)
				n_elements = np.sum(subset_indices)	
				expected_coverage = np.mean(per_component_coverage[subset_indices,9].astype(float))		
				t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + 'gene' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')
	if agg_eqtl_ss:
		for ln_pi_method in ln_pi_methods:
			for twas_method in twas_methods:
				subset_indices = (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method)
				components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
				prop = np.sum(components_covered)/len(components_covered)
				prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
				prop_lb = prop - (1.96*prop_se)
				prop_ub = prop + (1.96*prop_se)
				n_elements = np.sum(subset_indices)
				expected_coverage = np.mean(per_component_coverage[subset_indices,9].astype(float))
				t.write(str("Aggregate") + '\t' + ln_pi_method +'\t' + twas_method  + '\t' + 'all' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')

				subset_indices = (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method) & (per_component_coverage[:, 7] == 'variant')
				components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
				prop = np.sum(components_covered)/len(components_covered)
				prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
				prop_lb = prop - (1.96*prop_se)
				prop_ub = prop + (1.96*prop_se)
				n_elements = np.sum(subset_indices)
				expected_coverage = np.mean(per_component_coverage[subset_indices,9].astype(float))
				t.write(str("Aggregate") + '\t' + ln_pi_method +'\t' + twas_method  + '\t' + 'variant' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')

				subset_indices = (per_component_coverage[:,1] == ln_pi_method) & (per_component_coverage[:,2] == twas_method) & (per_component_coverage[:, 7] == 'gene')
				components_covered = (per_component_coverage[subset_indices,:][:,-2]).astype(float)
				prop = np.sum(components_covered)/len(components_covered)
				prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
				prop_lb = prop - (1.96*prop_se)
				prop_ub = prop + (1.96*prop_se)
				n_elements = np.sum(subset_indices)	
				expected_coverage = np.mean(per_component_coverage[subset_indices,9].astype(float))		
				t.write(str("Aggregate") + '\t' + ln_pi_method + '\t' + twas_method + '\t' + 'gene' + '\t' + str(n_elements) + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\t' + str(expected_coverage) + '\n')		

	t.close()
	return

def create_file_containing_focus_tg_high_pip_gene_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file, intercept=False):
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

		for gene_id in causal_gene:
			causal_genetic_elements[gene_id.split('_')[0]] = 1


		for eqtl_sample_size in eqtl_sample_sizes:
			# Create mapping from ln_pi_method to dictionary of elements in 95% cs for this simulation
			discovered_dicti = {}
			
			focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_focus_res.focus.tsv'
			if intercept == True:
				focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_focus_w_intercept_res.focus.tsv'
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
				if gene_id not in discovered_dicti:
					discovered_dicti[gene_id] = pip
				else:
					if pip > discovered_dicti[gene_id]:
						discovered_dicti[gene_id] = pip
			f.close()

			# get gene estimate by summing across gene-tissue pairs
			discovered_dicti2 = {}
			for gene_tissue in [*discovered_dicti]:
				gene_id = gene_tissue.split('_')[0]
				if gene_id not in discovered_dicti2:
					discovered_dicti2[gene_id] = discovered_dicti[gene_tissue]
				else:
					new_pip = discovered_dicti2[gene_id] + discovered_dicti[gene_tissue]
					if new_pip > 1:
						new_pip = 1.0
					discovered_dicti2[gene_id] = new_pip

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


				for gene_name in middle_genes:
					if gene_name not in causal_genetic_elements:
						continue
					# THis is a causal gene
					booler = 0.0
					if gene_name in discovered_dicti2:
						if discovered_dicti2[gene_name] >= pip_threshold:
							booler = 1.0
					t.write(str(eqtl_sample_size) + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gene_name + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return



def create_file_containing_focus_tg_high_pip_snp_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file, intercept=False):
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
			if intercept == True:
				focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_focus_w_intercept_res.focus.tsv'
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

def create_file_containing_focus_high_pip_gene_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file, intercept=False):
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

		for gene_id in [*causal_gene]:
			causal_genetic_elements[gene_id.split('_')[0]] = 1


		for eqtl_sample_size in eqtl_sample_sizes:
			# Create mapping from ln_pi_method to dictionary of elements in 95% cs for this simulation
			discovered_dicti = {}
			
			for tissue_iter in range(10):
				focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_iter) + '_focus_res.focus.tsv'
				if intercept == True:
					focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_iter) + '_focus_w_intercept_res.focus.tsv'
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
					ensamble = gene_id.split('_')[0]
					if pip >= pip_threshold:
						discovered_dicti[ensamble] = 1
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


				for gene_name in middle_genes:
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




def create_file_containing_focus_high_pip_snp_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file, intercept=False):
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
			
			for tissue_iter in range(10):
				focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_iter) + '_focus_res.focus.tsv'
				if intercept == True:
					focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_iter) + '_focus_w_intercept_res.focus.tsv'
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

def create_file_containing_two_step_focus_high_pip_snp_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_two_step_tgfm_results_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file, intercept=False):
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

			two_step_tissue_summary_file = simulated_two_step_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_pmces_uniform_two_step_tissues.txt'
			best_tissue_arr, sig_tissue_arr = extract_two_step_analyzed_tissues(two_step_tissue_summary_file)

			if len(best_tissue_arr) != 1:
				print('assumption eroror')
				pdb.set_trace()

			for tissue_name in best_tissue_arr:
				tissue_iter = int(tissue_name.split('ue')[1])

				focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_iter) + '_focus_res.focus.tsv'
				if intercept == True:
					focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_iter) + '_focus_w_intercept_res.focus.tsv'
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





def create_file_containing_coloc_high_pip_gene_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file):
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

		for causal_gene_name in [*causal_gene]:
			causal_genetic_elements[causal_gene_name.split('_')[0]] = 1


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
				gene_id = data[0].split('_')[0]
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


				for gene_name in middle_genes:
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

def create_file_containing_two_step_coloc_high_pip_snp_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_two_step_tgfm_results_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file):
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

			two_step_tissue_summary_file = simulated_two_step_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_pmces_uniform_two_step_tissues.txt'
			best_tissue_arr, sig_tissue_arr = extract_two_step_analyzed_tissues(two_step_tissue_summary_file)

			best_tissue_dicti = {}
			best_tissue_dicti[best_tissue_arr[0]] = 1

			
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
				cur_tissue_name = gene_id.split('_')[1]
				if cur_tissue_name not in best_tissue_dicti:
					continue
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

def create_file_containing_tgfm_high_pip_gene_power_per_component_vary_ge_h2s(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, ge_h2s, gene_type):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	# Open output file
	t = open(cs_power_per_component_output_file,'w')
	t.write('GE_h2\tln_pi_method\twas_method\tsimulation_number\twindow_name\tgenetic_element_class\tgenetic_element_name\tcausal_variant_in_cs\n')
	# First loop through simulations
	for simulation_number in simulation_runs:

		for ge_h2 in ge_h2s:

			temp_simulation_name_string = global_simulation_name_string.split('ss_')[0] + 'ss_' + str(100000) + '_ge_h2_' + ge_h2

			# Extract dictionary list of causal genetic elements
			causal_variant_file = simulated_trait_dir + 'simulation_new_' + str(simulation_number) + '_' + temp_simulation_name_string +'_all_t_' + gene_type + '_non_mediated_variant_causal_effect_sizes.txt'
			causal_gene_file = simulated_trait_dir + 'simulation_new_' + str(simulation_number) + '_' + temp_simulation_name_string +'_all_t_' + gene_type + '_expression_mediated_gene_causal_effect_sizes.txt'
			causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)
			for causal_gene_tissue in [*causal_gene]:
				causal_gene_name = causal_gene_tissue.split('_')[0]
				causal_genetic_elements[causal_gene_name] =1

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
					cs_file = simulated_tgfm_results_dir + 'simulation_new_' + str(simulation_number) + '_' + temp_simulation_name_string+'_all_t_' + gene_type  + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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

						if len(data) == 3 and data[2] == 'NA':
							continue
						tgfm_res_file = simulated_tgfm_results_dir + 'simulation_new_' + str(simulation_number) + '_' + temp_simulation_name_string+'_all_t_' + gene_type  + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_' + component_window_name + '_results.pkl'
						all_cs_genetic_elements, cs_probs = extract_middle_gene_agg_pips_from_tgfm_obj(tgfm_res_file)

						#all_cs_genetic_elements = data[1].split(';')
						#cs_probs = np.asarray(data[2].split(';')).astype(float)
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


				for gene_name in middle_genes:
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
							t.write(str(ge_h2) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gene_name + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return

def create_file_containing_tgfm_high_pip_gene_power_per_component_vary_gwas_ss(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, gwas_sample_sizes):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	# Open output file
	t = open(cs_power_per_component_output_file,'w')
	t.write('GWAS_sample_size\tln_pi_method\twas_method\tsimulation_number\twindow_name\tgenetic_element_class\tgenetic_element_name\tcausal_variant_in_cs\n')
	# First loop through simulations
	for simulation_number in simulation_runs:

		for gwas_sample_size in gwas_sample_sizes:


			temp_simulation_name_string = global_simulation_name_string.split('ss_')[0] + 'ss_' + str(gwas_sample_size) + '_ge_h2_075'

			# Extract dictionary list of causal genetic elements
			causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
			causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
			causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)
			for causal_gene_tissue in [*causal_gene]:
				causal_gene_name = causal_gene_tissue.split('_')[0]
				causal_genetic_elements[causal_gene_name] =1

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
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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

						if len(data) == 3 and data[2] == 'NA':
							continue
						tgfm_res_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_' + component_window_name + '_results.pkl'
						all_cs_genetic_elements, cs_probs = extract_middle_gene_agg_pips_from_tgfm_obj(tgfm_res_file)

						#all_cs_genetic_elements = data[1].split(';')
						#cs_probs = np.asarray(data[2].split(';')).astype(float)
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


				for gene_name in middle_genes:
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
							t.write(str(gwas_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gene_name + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return



def create_file_containing_tgfm_high_pip_gene_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, gene_type='component_gene', n_samp='100'):
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
		for causal_gene_tissue in [*causal_gene]:
			causal_gene_name = causal_gene_tissue.split('_')[0]
			causal_genetic_elements[causal_gene_name] =1

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
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_nsamp_' + str(n_samp) + '_' + gene_type+ '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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

						if len(data) == 3 and data[2] == 'NA':
							continue
						tgfm_res_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string+ '_nsamp_' + str(n_samp) + '_'+ gene_type + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_' + component_window_name + '_results.pkl'
						all_cs_genetic_elements, cs_probs = extract_middle_gene_agg_pips_from_tgfm_obj(tgfm_res_file)

						#all_cs_genetic_elements = data[1].split(';')
						#cs_probs = np.asarray(data[2].split(';')).astype(float)
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


				for gene_name in middle_genes:
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
			f.close()
	t.close()
	return


def create_file_containing_tgfm_high_pip_snp_power_per_component_vary_ge_h2s(global_simulation_name_string, eqtl_sample_size, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, ge_h2s):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	# Open output file
	t = open(cs_power_per_component_output_file,'w')
	t.write('GE_h2\tln_pi_method\twas_method\tsimulation_number\twindow_name\tgenetic_element_class\tgenetic_element_name\tcausal_variant_in_cs\n')
	# First loop through simulations
	for simulation_number in simulation_runs:
		for ge_h2 in ge_h2s:
			# Extract dictionary list of causal genetic elements
			temp_simulation_name_string = global_simulation_name_string.split('ss_')[0] + 'ss_' + str(100000) + '_ge_h2_' + ge_h2
			causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
			causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
			causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

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
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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
								t.write(str(ge_h2) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gene_name + '\t' + str(booler) + '\n')
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
							t.write(str(ge_h2) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'variant\t' + variant_name + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return




def create_file_containing_tgfm_high_pip_snp_power_per_component_vary_gwas_ss(global_simulation_name_string, eqtl_sample_size, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, gwas_sample_sizes):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	# Open output file
	t = open(cs_power_per_component_output_file,'w')
	t.write('GWAS_sample_size\tln_pi_method\twas_method\tsimulation_number\twindow_name\tgenetic_element_class\tgenetic_element_name\tcausal_variant_in_cs\n')
	# First loop through simulations
	for simulation_number in simulation_runs:
		for gwas_sample_size in gwas_sample_sizes:
			# Extract dictionary list of causal genetic elements
			temp_simulation_name_string = global_simulation_name_string.split('ss_')[0] + 'ss_' + str(gwas_sample_size) + '_ge_h2_075'
			causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
			causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
			causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

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
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + temp_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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
								t.write(str(gwas_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gene_name + '\t' + str(booler) + '\n')
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
							t.write(str(gwas_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'variant\t' + variant_name + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return

def create_file_containing_causal_gene_tissue_pair_stratification_in_missing_causal_tissue_sim(local_simulation_name_string, eqtl_sample_sizes, missingness_methods, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, simulated_best_tagging_gt_dir, simulated_gene_expression_dir):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	ln_pi_method = ln_pi_methods[0]
	twas_method = twas_methods[0]

	# Open output file
	t = open(cs_power_per_component_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\twas_method\tmissingness_method\tsimulation_number\twindow_name\tgenetic_element_class\tgenetic_element_name\tdiscovered_element_class\tdiscovered_element_boolean\n')
	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)


		gt_pair_to_causal_eqtl_var = extract_gwas_eqtl_variants_for_each_gene([*causal_gene], simulation_number, local_simulation_name_string, simulated_gene_expression_dir)

		for eqtl_sample_size in eqtl_sample_sizes:
			best_tagging_file = simulated_best_tagging_gt_dir+ 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_no_t0_eqtl_ss_' + str(eqtl_sample_size) + '_best_tagging_gt_pairs_from_sim_causal_effect.txt'
			tagging_causal_gt_pairs, causal_gt_pair_to_tagging_gt_pair = extract_dictionary_list_of_best_tagging_causal_gene_tissue_pairs(best_tagging_file)

			# Create mapping from ln_pi_method to dictionary of elements in 95% cs for this simulation
			discovered_dicti = {}
			for missingness_method in missingness_methods:
				discovered_pi_dicti = {}
				cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_' + missingness_method + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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
				discovered_dicti[missingness_method] = discovered_pi_dicti

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


				for missingness_method in missingness_methods:
					discovered_element_dicti = discovered_dicti[missingness_method]
					discovered_element_arr = np.asarray([*discovered_element_dicti])

					for gene_name_stem in middle_genes:
						for tissue_iter in range(10):
							gt_name = gene_name_stem + '_tissue' + str(tissue_iter)
							if gt_name not in causal_genetic_elements:
								continue
							# Step 1 causal gene-tissue pair category
							booler = 0.0
							booler2 = 0.0
							genetic_element_class = 'causal_gt'
							if gt_name in discovered_element_dicti:
								booler = 1.0
								booler2 = 1.0
							t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + missingness_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gt_name + '\t' + genetic_element_class + '\t' + str(booler) + '\n')
							# Step 2: best tagging gene-tissue pair
							best_tagging_gt_pair = causal_gt_pair_to_tagging_gt_pair[gt_name]
							booler = 0.0
							genetic_element_class = 'best_tagging_gt'
							if best_tagging_gt_pair in discovered_element_dicti:
								booler = 1.0
								booler2 = 1.0
							t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + missingness_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gt_name + '\t' + genetic_element_class + '\t' + str(booler) + '\n')
							# Step 3: other gene-tissue pair
							gene_name = gt_name.split('_')[0]
							genetic_element_class = 'other_gt'
							booler = 0.0
							if booler2 == 0.0:  # This means this causal gt was not discovered by the causal gt or the best tagging gt
								for ele in discovered_element_arr:
									if ele.startswith('ENSG') == False:
										continue
									if ele.split('_')[0] == gene_name:
										booler = 1.0
							t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + missingness_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gt_name + '\t' + genetic_element_class + '\t' + str(booler) + '\n')
							# Step 4: NM-variant
							gt_nm_var = gt_pair_to_causal_eqtl_var[gt_name]
							booler = 0.0
							for ele in discovered_element_arr:
								if ele in gt_nm_var:
									booler = 1.0
							genetic_element_class = 'nm_variant'
							t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + missingness_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gt_name + '\t' + genetic_element_class + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return




def create_file_containing_tgfm_high_pip_snp_power_per_component_in_missing_causal_tissue_sim(local_simulation_name_string, eqtl_sample_sizes, missingness_methods, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, simulated_best_tagging_gt_dir):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	ln_pi_method = ln_pi_methods[0]
	twas_method = twas_methods[0]

	# Open output file
	t = open(cs_power_per_component_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\twas_method\tmissingness_method\tsimulation_number\twindow_name\tgenetic_element_class\tgenetic_element_name\tcausal_variant_in_cs\n')
	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		for eqtl_sample_size in eqtl_sample_sizes:
			best_tagging_file = simulated_best_tagging_gt_dir+ 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_no_t0_eqtl_ss_' + str(eqtl_sample_size) + '_best_tagging_gt_pairs_from_sim_causal_effect.txt'
			tagging_causal_gt_pairs, dicti2 = extract_dictionary_list_of_best_tagging_causal_gene_tissue_pairs(best_tagging_file)

			# Create mapping from ln_pi_method to dictionary of elements in 95% cs for this simulation
			discovered_dicti = {}
			for missingness_method in missingness_methods:
				discovered_pi_dicti = {}
				cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + local_simulation_name_string + '_' + missingness_method + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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
				discovered_dicti[missingness_method] = discovered_pi_dicti

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
					for missingness_method in missingness_methods:
						for tissue_iter in range(10):
							gene_name = gene_name_stem + '_tissue' + str(tissue_iter)

							if missingness_method == 'all_t' and gene_name not in causal_genetic_elements:
								continue
							if missingness_method == 'no_t0' and gene_name not in tagging_causal_gt_pairs:
								continue
							# THis is a causal gene
							booler = 0.0
							if gene_name in discovered_dicti[missingness_method]:
								booler = 1.0
							t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + missingness_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gene_name + '\t' + str(booler) + '\n')
				for variant_name in middle_variants:
					if variant_name not in causal_genetic_elements:
						continue
					# This is a causal variant
					for missingness_method in missingness_methods:
						booler = 0.0
						if variant_name in discovered_dicti[missingness_method]:
							booler = 1.0
						t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + missingness_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'variant\t' + variant_name + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return

def create_file_containing_ctwas_tg_gene_tissue_fdr_power_curve_data_realistic_qtl_ss(global_simulation_name_string, eqtl_sample_size, simulation_runs, simulated_trait_dir, simulated_ctwas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file, fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)


	# Open raw data output file handle
	t = open(fdr_power_curve_raw_data_output_file,'w')
	t.write('simulation_number\twindow\tgene_name\tTGFM_PIP\tcausal_gene\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)


		# Extract causal eqtl sample size for this simulation
		causal_eqtl_ss_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_number) + '/simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_permuted_eqtl_ss.npy'
		causal_sample_size = np.load(causal_eqtl_ss_file)[0]

		if eqtl_sample_size == 'low_eqtl_ss' and causal_sample_size > 200:
			continue
		if eqtl_sample_size == 'high_eqtl_ss' and causal_sample_size < 200:
			continue			

		discovered_dicti = {}
		ctwas_result_file = simulated_ctwas_results_dir + 'simulation_new_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_realistic_qtl_ss/testing.susieIrss.txt'
		if os.path.isfile(ctwas_result_file) == False:
			continue
		f = open(ctwas_result_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			genetic_element_name = data[1]
			genetic_element_pip = float(data[7])

			if genetic_element_name.startswith('ENSG') == False:
				continue
			name_info = genetic_element_name.split('_')
			genetic_element_name = name_info[0] + '_' + name_info[1] + name_info[2]

			discovered_dicti[genetic_element_name] = genetic_element_pip
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
					causal_gener = 'null'
					if gene_name in causal_genetic_elements:
						causal_gener = 'causal'
					ctwas_pip = 0.0
					if gene_name in discovered_dicti:
						ctwas_pip = discovered_dicti[gene_name]
					t.write(str(simulation_number) + '\t' + window_name + '\t' + gene_name + '\t' + str(ctwas_pip) + '\t' + causal_gener + '\n')
		f.close()
	t.close()

	# Get vector of pips
	# Get vector of labels
	pips = []
	labels = []
	f = open(fdr_power_curve_raw_data_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pips.append(float(data[3]))
		labels.append(data[4])
	f.close()
	pips = np.asarray(pips)
	labels = np.asarray(labels)


	# Open output handle to processed fdr power curve data
	t = open(fdr_power_curve_data_output_file,'w')
	# Header
	t.write('pip_threshold\tfdr\tpower\n')
	prev_fdr = 0.0
	prev_power = 0.0
	for pip_threshold in np.arange(.05,1,.01):
		if np.sum(pips > pip_threshold) == 0:
			continue
		fdr = 1.0 - np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(pips > pip_threshold)
		power = np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(labels=='causal')
		if prev_fdr == fdr and prev_power == power:
			continue
		prev_fdr = fdr
		prev_power = power
		t.write(str(pip_threshold) + '\t' + str(fdr) + '\t' + str(power) + '\n')
	t.close()

	print(fdr_power_curve_data_output_file)
	return




def create_file_containing_ctwas_tg_gene_tissue_fdr_power_curve_data(global_simulation_name_string, eqtl_sample_size, simulation_runs, simulated_trait_dir, simulated_ctwas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file, fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)


	# Open raw data output file handle
	t = open(fdr_power_curve_raw_data_output_file,'w')
	t.write('simulation_number\twindow\tgene_name\tTGFM_PIP\tcausal_gene\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)



		discovered_dicti = {}
		ctwas_result_file = simulated_ctwas_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '/testing.susieIrss.txt'
		if os.path.isfile(ctwas_result_file) == False:
			continue
		f = open(ctwas_result_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			genetic_element_name = data[1]
			genetic_element_pip = float(data[7])

			if genetic_element_name.startswith('ENSG') == False:
				continue
			name_info = genetic_element_name.split('_')
			genetic_element_name = name_info[0] + '_' + name_info[1] + name_info[2]

			discovered_dicti[genetic_element_name] = genetic_element_pip
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
					causal_gener = 'null'
					if gene_name in causal_genetic_elements:
						causal_gener = 'causal'
					ctwas_pip = 0.0
					if gene_name in discovered_dicti:
						ctwas_pip = discovered_dicti[gene_name]
					t.write(str(simulation_number) + '\t' + window_name + '\t' + gene_name + '\t' + str(ctwas_pip) + '\t' + causal_gener + '\n')
		f.close()
	t.close()

	# Get vector of pips
	# Get vector of labels
	pips = []
	labels = []
	f = open(fdr_power_curve_raw_data_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pips.append(float(data[3]))
		labels.append(data[4])
	f.close()
	pips = np.asarray(pips)
	labels = np.asarray(labels)


	# Open output handle to processed fdr power curve data
	t = open(fdr_power_curve_data_output_file,'w')
	# Header
	t.write('pip_threshold\tfdr\tpower\n')
	prev_fdr = 0.0
	prev_power = 0.0
	for pip_threshold in np.arange(.05,1,.01):
		if np.sum(pips > pip_threshold) == 0:
			continue
		fdr = 1.0 - np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(pips > pip_threshold)
		power = np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(labels=='causal')
		if prev_fdr == fdr and prev_power == power:
			continue
		prev_fdr = fdr
		prev_power = power
		t.write(str(pip_threshold) + '\t' + str(fdr) + '\t' + str(power) + '\n')
	t.close()

	return


def create_file_containing_focus_tg_gene_tissue_fdr_power_curve_data(global_simulation_name_string, eqtl_sample_size, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file, fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir, intercept=False):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)


	# Open raw data output file handle
	t = open(fdr_power_curve_raw_data_output_file,'w')
	t.write('simulation_number\twindow\tgene_name\tTGFM_PIP\tcausal_gene\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)



		discovered_dicti = {}

		focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_focus_res.focus.tsv'
		if intercept == True:
			focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_focus_w_intercept_res.focus.tsv'
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
			discovered_dicti[gene_id] = pip
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
					causal_gener = 'null'
					if gene_name in causal_genetic_elements:
						causal_gener = 'causal'
					ctwas_pip = 0.0
					if gene_name in discovered_dicti:
						ctwas_pip = discovered_dicti[gene_name]
					t.write(str(simulation_number) + '\t' + window_name + '\t' + gene_name + '\t' + str(ctwas_pip) + '\t' + causal_gener + '\n')
		f.close()
	t.close()

	# Get vector of pips
	# Get vector of labels
	pips = []
	labels = []
	f = open(fdr_power_curve_raw_data_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pips.append(float(data[3]))
		labels.append(data[4])
	f.close()
	pips = np.asarray(pips)
	labels = np.asarray(labels)


	# Open output handle to processed fdr power curve data
	t = open(fdr_power_curve_data_output_file,'w')
	# Header
	t.write('pip_threshold\tfdr\tpower\n')
	prev_fdr = 0.0
	prev_power = 0.0
	for pip_threshold in np.arange(.05,1,.01):
		if np.sum(pips > pip_threshold) == 0:
			continue
		fdr = 1.0 - np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(pips > pip_threshold)
		power = np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(labels=='causal')
		if prev_fdr == fdr and prev_power == power:
			continue
		prev_fdr = fdr
		prev_power = power
		t.write(str(pip_threshold) + '\t' + str(fdr) + '\t' + str(power) + '\n')
	t.close()

	return




def create_file_containing_two_step_focus_gene_tissue_fdr_power_curve_data(global_simulation_name_string, eqtl_sample_size, simulation_runs, simulated_trait_dir,simulated_two_step_tgfm_results_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file, fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir, intercept=False):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)


	# Open raw data output file handle
	t = open(fdr_power_curve_raw_data_output_file,'w')
	t.write('simulation_number\twindow\tgene_name\tTGFM_PIP\tcausal_gene\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)



		discovered_dicti = {}


		two_step_tissue_summary_file = simulated_two_step_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_pmces_uniform_two_step_tissues.txt'
		best_tissue_arr, sig_tissue_arr = extract_two_step_analyzed_tissues(two_step_tissue_summary_file)

		if len(best_tissue_arr) != 1:
			print('assumption eroror')
			pdb.set_trace()


		for tissue_name in best_tissue_arr:
			tissue_iter = int(tissue_name.split('ue')[1])
			focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_iter) + '_focus_res.focus.tsv'
			if intercept == True:
				focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_iter) + '_focus_w_intercept_res.focus.tsv'
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
				discovered_dicti[gene_id] = pip
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
					causal_gener = 'null'
					if gene_name in causal_genetic_elements:
						causal_gener = 'causal'
					ctwas_pip = 0.0
					if gene_name in discovered_dicti:
						ctwas_pip = discovered_dicti[gene_name]
					t.write(str(simulation_number) + '\t' + window_name + '\t' + gene_name + '\t' + str(ctwas_pip) + '\t' + causal_gener + '\n')
		f.close()
	t.close()

	# Get vector of pips
	# Get vector of labels
	pips = []
	labels = []
	f = open(fdr_power_curve_raw_data_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pips.append(float(data[3]))
		labels.append(data[4])
	f.close()
	pips = np.asarray(pips)
	labels = np.asarray(labels)


	# Open output handle to processed fdr power curve data
	t = open(fdr_power_curve_data_output_file,'w')
	# Header
	t.write('pip_threshold\tfdr\tpower\n')
	prev_fdr = 0.0
	prev_power = 0.0
	for pip_threshold in np.arange(.05,1,.01):
		if np.sum(pips > pip_threshold) == 0:
			continue
		fdr = 1.0 - np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(pips > pip_threshold)
		power = np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(labels=='causal')
		if prev_fdr == fdr and prev_power == power:
			continue
		prev_fdr = fdr
		prev_power = power
		t.write(str(pip_threshold) + '\t' + str(fdr) + '\t' + str(power) + '\n')
	t.close()

	return



def create_file_containing_focus_gene_tissue_fdr_power_curve_data(global_simulation_name_string, eqtl_sample_size, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file, fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir, intercept=False):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)


	# Open raw data output file handle
	t = open(fdr_power_curve_raw_data_output_file,'w')
	t.write('simulation_number\twindow\tgene_name\tTGFM_PIP\tcausal_gene\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)



		discovered_dicti = {}
		for tissue_iter in range(10):
			focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_iter) + '_focus_res.focus.tsv'
			if intercept == True:
				focus_results_file = simulated_focus_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_iter) + '_focus_w_intercept_res.focus.tsv'
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
				discovered_dicti[gene_id] = pip
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
					causal_gener = 'null'
					if gene_name in causal_genetic_elements:
						causal_gener = 'causal'
					ctwas_pip = 0.0
					if gene_name in discovered_dicti:
						ctwas_pip = discovered_dicti[gene_name]
					t.write(str(simulation_number) + '\t' + window_name + '\t' + gene_name + '\t' + str(ctwas_pip) + '\t' + causal_gener + '\n')
		f.close()
	t.close()

	# Get vector of pips
	# Get vector of labels
	pips = []
	labels = []
	f = open(fdr_power_curve_raw_data_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pips.append(float(data[3]))
		labels.append(data[4])
	f.close()
	pips = np.asarray(pips)
	labels = np.asarray(labels)


	# Open output handle to processed fdr power curve data
	t = open(fdr_power_curve_data_output_file,'w')
	# Header
	t.write('pip_threshold\tfdr\tpower\n')
	prev_fdr = 0.0
	prev_power = 0.0
	for pip_threshold in np.arange(.05,1,.01):
		if np.sum(pips > pip_threshold) == 0:
			continue
		fdr = 1.0 - np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(pips > pip_threshold)
		power = np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(labels=='causal')
		if prev_fdr == fdr and prev_power == power:
			continue
		prev_fdr = fdr
		prev_power = power
		t.write(str(pip_threshold) + '\t' + str(fdr) + '\t' + str(power) + '\n')
	t.close()

	return

def create_file_containing_two_step_ctwas_gene_tissue_fdr_power_curve_data(global_simulation_name_string, eqtl_sample_size, simulation_runs, simulated_trait_dir, simulated_two_step_tgfm_results_dir, simulated_ctwas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file, fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)


	# Open raw data output file handle
	t = open(fdr_power_curve_raw_data_output_file,'w')
	t.write('simulation_number\twindow\tgene_name\tTGFM_PIP\tcausal_gene\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)



		discovered_dicti = {}
		missing=False


		two_step_tissue_summary_file = simulated_two_step_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_pmces_uniform_two_step_tissues.txt'

		best_tissue_arr, sig_tissue_arr = extract_two_step_analyzed_tissues(two_step_tissue_summary_file)

		for tissue_name in best_tissue_arr:
			tissue_number = tissue_name.split('ue')[1]
			ctwas_result_file = simulated_ctwas_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_number) + '/testing.susieIrss.txt'
			if os.path.isfile(ctwas_result_file) == False:
				missing=True
				continue
			f = open(ctwas_result_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				genetic_element_name = data[1]
				genetic_element_pip = float(data[7])

				if genetic_element_name.startswith('ENSG') == False:
					continue
				name_info = genetic_element_name.split('_')
				genetic_element_name = name_info[0] + '_' + name_info[1] + name_info[2]

				if genetic_element_name in discovered_dicti:
					print('assumption oeroror')
					pdb.set_trace()

				discovered_dicti[genetic_element_name] = genetic_element_pip
			f.close()
		if missing == True:
			continue


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
					causal_gener = 'null'
					if gene_name in causal_genetic_elements:
						causal_gener = 'causal'
					ctwas_pip = 0.0
					if gene_name in discovered_dicti:
						ctwas_pip = discovered_dicti[gene_name]
					t.write(str(simulation_number) + '\t' + window_name + '\t' + gene_name + '\t' + str(ctwas_pip) + '\t' + causal_gener + '\n')
		f.close()
	t.close()

	# Get vector of pips
	# Get vector of labels
	pips = []
	labels = []
	f = open(fdr_power_curve_raw_data_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pips.append(float(data[3]))
		labels.append(data[4])
	f.close()
	pips = np.asarray(pips)
	labels = np.asarray(labels)


	# Open output handle to processed fdr power curve data
	t = open(fdr_power_curve_data_output_file,'w')
	# Header
	t.write('pip_threshold\tfdr\tpower\n')
	prev_fdr = 0.0
	prev_power = 0.0
	for pip_threshold in np.arange(.05,1,.01):
		if np.sum(pips > pip_threshold) == 0:
			continue
		fdr = 1.0 - np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(pips > pip_threshold)
		power = np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(labels=='causal')
		if prev_fdr == fdr and prev_power == power:
			continue
		prev_fdr = fdr
		prev_power = power
		t.write(str(pip_threshold) + '\t' + str(fdr) + '\t' + str(power) + '\n')
	t.close()

	return


def create_file_containing_ctwas_gene_tissue_fdr_power_curve_data(global_simulation_name_string, eqtl_sample_size, simulation_runs, simulated_trait_dir, simulated_ctwas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file, fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)


	# Open raw data output file handle
	t = open(fdr_power_curve_raw_data_output_file,'w')
	t.write('simulation_number\twindow\tgene_name\tTGFM_PIP\tcausal_gene\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)



		discovered_dicti = {}
		missing=False
		for tissue_number in range(10):
			ctwas_result_file = simulated_ctwas_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_number)+ '/testing.susieIrss.txt'
			if os.path.isfile(ctwas_result_file) == False:
				missing=True
				continue
			f = open(ctwas_result_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				genetic_element_name = data[1]
				genetic_element_pip = float(data[7])

				if genetic_element_name.startswith('ENSG') == False:
					continue
				name_info = genetic_element_name.split('_')
				genetic_element_name = name_info[0] + '_' + name_info[1] + name_info[2]

				if genetic_element_name in discovered_dicti:
					print('assumption oeroror')
					pdb.set_trace()

				discovered_dicti[genetic_element_name] = genetic_element_pip
			f.close()
		if missing == True:
			continue


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
					causal_gener = 'null'
					if gene_name in causal_genetic_elements:
						causal_gener = 'causal'
					ctwas_pip = 0.0
					if gene_name in discovered_dicti:
						ctwas_pip = discovered_dicti[gene_name]
					t.write(str(simulation_number) + '\t' + window_name + '\t' + gene_name + '\t' + str(ctwas_pip) + '\t' + causal_gener + '\n')
		f.close()
	t.close()

	# Get vector of pips
	# Get vector of labels
	pips = []
	labels = []
	f = open(fdr_power_curve_raw_data_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pips.append(float(data[3]))
		labels.append(data[4])
	f.close()
	pips = np.asarray(pips)
	labels = np.asarray(labels)


	# Open output handle to processed fdr power curve data
	t = open(fdr_power_curve_data_output_file,'w')
	# Header
	t.write('pip_threshold\tfdr\tpower\n')
	prev_fdr = 0.0
	prev_power = 0.0
	for pip_threshold in np.arange(.05,1,.01):
		if np.sum(pips > pip_threshold) == 0:
			continue
		fdr = 1.0 - np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(pips > pip_threshold)
		power = np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(labels=='causal')
		if prev_fdr == fdr and prev_power == power:
			continue
		prev_fdr = fdr
		prev_power = power
		t.write(str(pip_threshold) + '\t' + str(fdr) + '\t' + str(power) + '\n')
	t.close()

	return


def create_file_containing_two_step_coloc_gene_tissue_fdr_power_curve_data(global_simulation_name_string, eqtl_sample_size, simulation_runs, simulated_trait_dir, simulated_two_step_tgfm_results_dir, simulated_coloc_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file, fdr_power_curve_data_output_file, global_window_file, simulated_learned_gene_models_dir):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)


	# Open raw data output file handle
	t = open(fdr_power_curve_raw_data_output_file,'w')
	t.write('simulation_number\twindow\tgene_name\tTGFM_PIP\tcausal_gene\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		two_step_tissue_summary_file = simulated_two_step_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_pmces_uniform_two_step_tissues.txt'
		best_tissue_arr, sig_tissue_arr = extract_two_step_analyzed_tissues(two_step_tissue_summary_file)

		best_tissue_dicti = {}
		best_tissue_dicti[best_tissue_arr[0]] = 1


		discovered_dicti = {}
		coloc_results_file = simulated_coloc_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_coloc_results.txt'
		f = open(coloc_results_file)
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
			cur_tissue_name = gene_id.split('_')[1]
			if cur_tissue_name not in best_tissue_dicti:
				continue

			if gene_id in discovered_dicti:
				print('assumption eororor')
				pdb.set_trace()
			discovered_dicti[gene_id] = pip
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
					causal_gener = 'null'
					if gene_name in causal_genetic_elements:
						causal_gener = 'causal'
					ctwas_pip = 0.0
					if gene_name in discovered_dicti:
						ctwas_pip = discovered_dicti[gene_name]
					t.write(str(simulation_number) + '\t' + window_name + '\t' + gene_name + '\t' + str(ctwas_pip) + '\t' + causal_gener + '\n')
		f.close()
	t.close()

	# Get vector of pips
	# Get vector of labels
	pips = []
	labels = []
	f = open(fdr_power_curve_raw_data_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pips.append(float(data[3]))
		labels.append(data[4])
	f.close()
	pips = np.asarray(pips)
	labels = np.asarray(labels)


	# Open output handle to processed fdr power curve data
	t = open(fdr_power_curve_data_output_file,'w')
	# Header
	t.write('pip_threshold\tfdr\tpower\n')
	prev_fdr = 0.0
	prev_power = 0.0
	for pip_threshold in np.arange(.05,1,.01):
		if np.sum(pips > pip_threshold) == 0:
			continue
		fdr = 1.0 - np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(pips > pip_threshold)
		power = np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(labels=='causal')
		if prev_fdr == fdr and prev_power == power:
			continue
		prev_fdr = fdr
		prev_power = power
		t.write(str(pip_threshold) + '\t' + str(fdr) + '\t' + str(power) + '\n')
	t.close()

	return




def create_file_containing_coloc_gene_tissue_fdr_power_curve_data(global_simulation_name_string, eqtl_sample_size, simulation_runs, simulated_trait_dir, simulated_coloc_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file, fdr_power_curve_data_output_file, global_window_file, simulated_learned_gene_models_dir):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)


	# Open raw data output file handle
	t = open(fdr_power_curve_raw_data_output_file,'w')
	t.write('simulation_number\twindow\tgene_name\tTGFM_PIP\tcausal_gene\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)



		discovered_dicti = {}
		coloc_results_file = simulated_coloc_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_coloc_results.txt'
		f = open(coloc_results_file)
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
			if gene_id in discovered_dicti:
				print('assumption eororor')
				pdb.set_trace()
			discovered_dicti[gene_id] = pip
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
					causal_gener = 'null'
					if gene_name in causal_genetic_elements:
						causal_gener = 'causal'
					ctwas_pip = 0.0
					if gene_name in discovered_dicti:
						ctwas_pip = discovered_dicti[gene_name]
					t.write(str(simulation_number) + '\t' + window_name + '\t' + gene_name + '\t' + str(ctwas_pip) + '\t' + causal_gener + '\n')
		f.close()
	t.close()

	# Get vector of pips
	# Get vector of labels
	pips = []
	labels = []
	f = open(fdr_power_curve_raw_data_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pips.append(float(data[3]))
		labels.append(data[4])
	f.close()
	pips = np.asarray(pips)
	labels = np.asarray(labels)


	# Open output handle to processed fdr power curve data
	t = open(fdr_power_curve_data_output_file,'w')
	# Header
	t.write('pip_threshold\tfdr\tpower\n')
	prev_fdr = 0.0
	prev_power = 0.0
	for pip_threshold in np.arange(.05,1,.01):
		if np.sum(pips > pip_threshold) == 0:
			continue
		fdr = 1.0 - np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(pips > pip_threshold)
		power = np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(labels=='causal')
		if prev_fdr == fdr and prev_power == power:
			continue
		prev_fdr = fdr
		prev_power = power
		t.write(str(pip_threshold) + '\t' + str(fdr) + '\t' + str(power) + '\n')
	t.close()

	return




def create_file_containing_tgfm_gene_tissue_fdr_power_curve_data_realistic_qtl_ss(global_simulation_name_string, eqtl_sample_size, simulation_runs, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file, fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir, gene_type='component_gene'):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	if gene_type.startswith('max_min_ratio_'):
		simulation_runs = simulation_runs + 100

	# Open raw data output file handle
	t = open(fdr_power_curve_raw_data_output_file,'w')
	t.write('simulation_number\twindow\tgene_name\tTGFM_PIP\tcausal_gene\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)


		# Extract causal eqtl sample size for this simulation
		causal_eqtl_ss_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_number) + '/simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_permuted_eqtl_ss.npy'
		causal_sample_size = np.load(causal_eqtl_ss_file)[0]


		if eqtl_sample_size == 'low_eqtl_ss' and causal_sample_size > 200:
			continue
		if eqtl_sample_size == 'high_eqtl_ss' and causal_sample_size < 200:
			continue			


		discovered_dicti = {}
		cs_file = simulated_tgfm_results_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_all_t_' + gene_type + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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
				ele_name = all_cs_genetic_elements[element_iter]
				if ele_name.startswith('ENSG') == False:
					continue
				discovered_dicti[ele_name] = cs_prob
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
					causal_gener = 'null'
					if gene_name in causal_genetic_elements:
						causal_gener = 'causal'
					tgfm_pip = 0.0
					if gene_name in discovered_dicti:
						tgfm_pip = discovered_dicti[gene_name]
					t.write(str(simulation_number) + '\t' + window_name + '\t' + gene_name + '\t' + str(tgfm_pip) + '\t' + causal_gener + '\n')
		f.close()
	t.close()

	# Get vector of pips
	# Get vector of labels
	pips = []
	labels = []
	f = open(fdr_power_curve_raw_data_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pips.append(float(data[3]))
		labels.append(data[4])
	f.close()
	pips = np.asarray(pips)
	labels = np.asarray(labels)


	# Open output handle to processed fdr power curve data
	t = open(fdr_power_curve_data_output_file,'w')
	# Header
	t.write('pip_threshold\tfdr\tpower\n')
	prev_fdr = 0.0
	prev_power = 0.0
	for pip_threshold in np.arange(.05,1,.01):
		if np.sum(pips > pip_threshold) == 0:
			continue
		fdr = 1.0 - np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(pips > pip_threshold)
		power = np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(labels=='causal')
		if prev_fdr == fdr and prev_power == power:
			continue
		prev_fdr = fdr
		prev_power = power
		t.write(str(pip_threshold) + '\t' + str(fdr) + '\t' + str(power) + '\n')
	t.close()

	print(fdr_power_curve_data_output_file)
	return


def create_file_containing_tgfm_gene_tissue_fdr_power_curve_data(global_simulation_name_string, simulation_runs, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file, fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file,eqtl_sample_size, gene_type='component_gene', n_samp='100'):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)


	# Open raw data output file handle
	t = open(fdr_power_curve_raw_data_output_file,'w')
	t.write('simulation_number\twindow\tgene_name\tTGFM_PIP\tcausal_gene\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)


		discovered_dicti = {}
		cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_nsamp_' + str(n_samp) + '_' + gene_type + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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
				ele_name = all_cs_genetic_elements[element_iter]
				if ele_name.startswith('ENSG') == False:
					continue
				discovered_dicti[ele_name] = cs_prob
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
					causal_gener = 'null'
					if gene_name in causal_genetic_elements:
						causal_gener = 'causal'
					tgfm_pip = 0.0
					if gene_name in discovered_dicti:
						tgfm_pip = discovered_dicti[gene_name]
					t.write(str(simulation_number) + '\t' + window_name + '\t' + gene_name + '\t' + str(tgfm_pip) + '\t' + causal_gener + '\n')
		f.close()
	t.close()

	# Get vector of pips
	# Get vector of labels
	pips = []
	labels = []
	f = open(fdr_power_curve_raw_data_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pips.append(float(data[3]))
		labels.append(data[4])
	f.close()
	pips = np.asarray(pips)
	labels = np.asarray(labels)


	# Open output handle to processed fdr power curve data
	t = open(fdr_power_curve_data_output_file,'w')
	# Header
	t.write('pip_threshold\tfdr\tpower\n')
	prev_fdr = 0.0
	prev_power = 0.0
	for pip_threshold in np.arange(.05,1,.01):
		if np.sum(pips > pip_threshold) == 0:
			continue
		fdr = 1.0 - np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(pips > pip_threshold)
		power = np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(labels=='causal')
		if prev_fdr == fdr and prev_power == power:
			continue
		prev_fdr = fdr
		prev_power = power
		t.write(str(pip_threshold) + '\t' + str(fdr) + '\t' + str(power) + '\n')
	t.close()
	return




def create_file_containing_two_step_tgfm_gene_tissue_fdr_power_curve_data(global_simulation_name_string, simulation_runs, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file, fdr_power_curve_data_output_file,two_step_method, global_window_file,eqtl_sample_size):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)


	# Open raw data output file handle
	t = open(fdr_power_curve_raw_data_output_file,'w')
	t.write('simulation_number\twindow\tgene_name\tTGFM_PIP\tcausal_gene\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)



		two_step_tissue_summary_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_pmces_uniform_two_step_tissues.txt'
		best_tissue_arr, sig_tissue_arr = extract_two_step_analyzed_tissues(two_step_tissue_summary_file)
			
		if two_step_method == 'best_tissue':
			tissue_vec = np.copy(best_tissue_arr)
		elif two_step_method == 'significant_tissues':
			tissue_vec = np.copy(sig_tissue_arr)
		else:
			print('tissue method assumption eroror')
			pdb.set_trace()
		discovered_dicti = {}
		skipped = False
		for tissue_name in tissue_vec:
			# Credible set file for this run
			cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_sampler_iterative_two_step_' + tissue_name + '_tgfm_pip_summary.txt'
			if os.path.isfile(cs_file) == False:
				skipped = True
				continue
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
					ele_name = all_cs_genetic_elements[element_iter]
					if ele_name.startswith('ENSG') == False:
						continue
					discovered_dicti[ele_name] = cs_prob
			f.close()

		if skipped:
			continue

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
					causal_gener = 'null'
					if gene_name in causal_genetic_elements:
						causal_gener = 'causal'
					tgfm_pip = 0.0
					if gene_name in discovered_dicti:
						tgfm_pip = discovered_dicti[gene_name]
					t.write(str(simulation_number) + '\t' + window_name + '\t' + gene_name + '\t' + str(tgfm_pip) + '\t' + causal_gener + '\n')
		f.close()
	t.close()

	# Get vector of pips
	# Get vector of labels
	pips = []
	labels = []
	f = open(fdr_power_curve_raw_data_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pips.append(float(data[3]))
		labels.append(data[4])
	f.close()
	pips = np.asarray(pips)
	labels = np.asarray(labels)


	# Open output handle to processed fdr power curve data
	t = open(fdr_power_curve_data_output_file,'w')
	# Header
	t.write('pip_threshold\tfdr\tpower\n')
	prev_fdr = 0.0
	prev_power = 0.0
	for pip_threshold in np.arange(.05,1,.01):
		if np.sum(pips > pip_threshold) == 0:
			continue
		fdr = 1.0 - np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(pips > pip_threshold)
		power = np.sum(labels[pips > pip_threshold] == 'causal')/np.sum(labels=='causal')
		if prev_fdr == fdr and prev_power == power:
			continue
		prev_fdr = fdr
		prev_power = power
		t.write(str(pip_threshold) + '\t' + str(fdr) + '\t' + str(power) + '\n')
	t.close()
	return




def create_file_containing_ctwas_tg_high_pip_snp_power_per_component_realistic_qtl_sample_size(global_simulation_name_string, simulation_runs, simulated_trait_dir, simulated_ctwas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, global_window_file, simulated_learned_gene_models_dir,tiss_filter=False, gene_type='component_gene'):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	eqtl_sample_sizes = ['low_eqtl_ss', 'high_eqtl_ss']


	ln_pi_method='ctwas'
	twas_method='lasso'



	# Open output file
	t = open(cs_power_per_component_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\twas_method\tsimulation_number\twindow_name\tgenetic_element_class\tgenetic_element_name\tcausal_variant_in_cs\n')
	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)


		# Extract causal eqtl sample size for this simulation
		causal_eqtl_ss_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_number) + '/simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_permuted_eqtl_ss.npy'
		causal_sample_size = np.load(causal_eqtl_ss_file)[0]

		discovered_dicti = {}
		for eqtl_sample_size in eqtl_sample_sizes:
			# Create mapping from ln_pi_method to dictionary of elements in 95% cs for this simulation
			if eqtl_sample_size == 'low_eqtl_ss' and causal_sample_size > 200:
				continue
			if eqtl_sample_size == 'high_eqtl_ss' and causal_sample_size < 200:
				continue			

			ctwas_result_file = simulated_ctwas_results_dir + 'simulation_new_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_realistic_qtl_ss/testing.susieIrss.txt'
			if os.path.isfile(ctwas_result_file) == False:
				continue
			discovered_pi_dicti = {}
			f = open(ctwas_result_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				component_window_name = data[4]
				component_num = 'NA'
				genetic_element_name = data[1]
				genetic_element_pip = float(data[7])
				if genetic_element_pip < pip_threshold:
					continue
				if genetic_element_name.startswith('ENSG'):
					name_info = genetic_element_name.split('_')
					if len(name_info) != 3:
						print('assumtion erorro')
						pdb.set_trace()
					genetic_element_name = name_info[0] + '_' + name_info[1] + name_info[2]
				discovered_pi_dicti[genetic_element_name] = 1
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
						booler = 0.0
						if gene_name in discovered_dicti[ln_pi_method + '_' + twas_method]:
							booler = 1.0
						t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gene_name + '\t' + str(booler) + '\n')
				for variant_name in middle_variants:
					if variant_name not in causal_genetic_elements:
						continue
					# This is a causal variant
					booler = 0.0
					if variant_name in discovered_dicti[ln_pi_method + '_' + twas_method]:
						booler = 1.0
					t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'variant\t' + variant_name + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return


def create_file_containing_ctwas_high_pip_snp_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_ctwas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, global_window_file, simulated_learned_gene_models_dir, gene_type='component_gene'):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	ln_pi_method='ctwas'
	twas_method='lasso'

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
			discovered_dicti = {}
			missing = False
			discovered_pi_dicti = {}

			for tissue_number in range(10):
				ctwas_result_file = simulated_ctwas_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_number) + '/testing.susieIrss.txt'
				if os.path.isfile(ctwas_result_file) == False:
					missing = True
					continue

				f = open(ctwas_result_file)
				head_count = 0
				for line in f:
					line = line.rstrip()
					data = line.split('\t')
					if head_count == 0:
						head_count = head_count + 1
						continue
					component_window_name = data[4]
					component_num = 'NA'
					genetic_element_name = data[1]
					genetic_element_pip = float(data[7])
					if genetic_element_pip < pip_threshold:
						continue
					if genetic_element_name.startswith('ENSG'):
						name_info = genetic_element_name.split('_')
						if len(name_info) != 3:
							print('assumtion erorro')
							pdb.set_trace()
						genetic_element_name = name_info[0] + '_' + name_info[1] + name_info[2]
					discovered_pi_dicti[genetic_element_name] = 1
				f.close()
				discovered_dicti[ln_pi_method + '_' + twas_method] = discovered_pi_dicti
			if missing == True:
				continue

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
						if gene_name in discovered_dicti[ln_pi_method + '_' + twas_method]:
							booler = 1.0
						t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gene_name + '\t' + str(booler) + '\n')
				for variant_name in middle_variants:
					if variant_name not in causal_genetic_elements:
						continue
					# This is a causal variant
					booler = 0.0
					if variant_name in discovered_dicti[ln_pi_method + '_' + twas_method]:
						booler = 1.0
					t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'variant\t' + variant_name + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return


def create_file_containing_two_step_ctwas_high_pip_snp_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_two_step_tgfm_results_dir, simulated_ctwas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, global_window_file, simulated_learned_gene_models_dir, gene_type='component_gene'):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	ln_pi_method='ctwas'
	twas_method='lasso'

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
			discovered_dicti = {}
			missing = False
			discovered_pi_dicti = {}


			two_step_tissue_summary_file = simulated_two_step_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_pmces_uniform_two_step_tissues.txt'
			best_tissue_arr, sig_tissue_arr = extract_two_step_analyzed_tissues(two_step_tissue_summary_file)

			for tissue_name in best_tissue_arr:
				tissue_number = tissue_name.split('ue')[1]
				ctwas_result_file = simulated_ctwas_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_number) + '/testing.susieIrss.txt'
				if os.path.isfile(ctwas_result_file) == False:
					missing = True
					continue

				f = open(ctwas_result_file)
				head_count = 0
				for line in f:
					line = line.rstrip()
					data = line.split('\t')
					if head_count == 0:
						head_count = head_count + 1
						continue
					component_window_name = data[4]
					component_num = 'NA'
					genetic_element_name = data[1]
					genetic_element_pip = float(data[7])
					if genetic_element_pip < pip_threshold:
						continue
					if genetic_element_name.startswith('ENSG'):
						name_info = genetic_element_name.split('_')
						if len(name_info) != 3:
							print('assumtion erorro')
							pdb.set_trace()
						genetic_element_name = name_info[0] + '_' + name_info[1] + name_info[2]
					discovered_pi_dicti[genetic_element_name] = 1
				f.close()
				discovered_dicti[ln_pi_method + '_' + twas_method] = discovered_pi_dicti
			if missing == True:
				continue

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
						if gene_name in discovered_dicti[ln_pi_method + '_' + twas_method]:
							booler = 1.0
						t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gene_name + '\t' + str(booler) + '\n')
				for variant_name in middle_variants:
					if variant_name not in causal_genetic_elements:
						continue
					# This is a causal variant
					booler = 0.0
					if variant_name in discovered_dicti[ln_pi_method + '_' + twas_method]:
						booler = 1.0
					t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'variant\t' + variant_name + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return


def create_file_containing_ctwas_tg_high_pip_gene_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_ctwas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, global_window_file, simulated_learned_gene_models_dir, gene_type='component_gene'):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	ln_pi_method='ctwas'
	twas_method='lasso'

	# Open output file
	t = open(cs_power_per_component_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\twas_method\tsimulation_number\twindow_name\tgenetic_element_class\tgenetic_element_name\tcausal_variant_in_cs\n')
	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		for gene_id in [*causal_genes]:
			ensamble = gene_id.split('_')[0]
			causal_genetic_elements[ensamble] = 1

		for eqtl_sample_size in eqtl_sample_sizes:
			discovered_dicti = {}
			missing = False
			discovered_pi_dicti = {}

			ctwas_result_file = simulated_ctwas_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '/testing.susieIrss.txt'
			if os.path.isfile(ctwas_result_file) == False:
				continue
			# Stream ctwas results file
			head_count = 0
			f = open(ctwas_result_file)

			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				component_window_name = data[4]
				component_num = 'NA'
				genetic_element_name = data[1]
				if genetic_element_name.startswith('ENSG') == False:
					continue
				genetic_element_pip = float(data[7])
				name_info = genetic_element_name.split('_')
				if len(name_info) != 3:
					print('assumtion erorro')
					pdb.set_trace()
				genetic_element_name = name_info[0] + '_' + name_info[1] + name_info[2]
				ensamble_id = name_info[0]
				if ensamble_id not in discovered_pi_dicti:
					discovered_pi_dicti[ensamble_id] = genetic_element_pip
				else:
					# Take sum
					old_pip = discovered_pi_dicti[ensamble_id]
					new_pip = old_pip + genetic_element_pip
					if new_pip > 1.0:
						new_pip = 1.0
					discovered_pi_dicti[ensamble_id] = new_pip
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


				for gene_name in middle_genes:
					if gene_name not in causal_genetic_elements:
						continue
					# THis is a causal gene
					booler = 0.0
					if gene_name in discovered_dicti[ln_pi_method + '_' + twas_method]:
						gene_pip = discovered_dicti[ln_pi_method + '_' + twas_method][gene_name]
						if gene_pip >= pip_threshold:
							booler = 1.0
					t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gene_name + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return







def create_file_containing_ctwas_high_pip_gene_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_ctwas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, global_window_file, simulated_learned_gene_models_dir, gene_type='component_gene'):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	ln_pi_method='ctwas'
	twas_method='lasso'

	# Open output file
	t = open(cs_power_per_component_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\twas_method\tsimulation_number\twindow_name\tgenetic_element_class\tgenetic_element_name\tcausal_variant_in_cs\n')
	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		for gene_id in [*causal_genes]:
			ensamble = gene_id.split('_')[0]
			causal_genetic_elements[ensamble] = 1

		for eqtl_sample_size in eqtl_sample_sizes:
			discovered_dicti = {}
			missing = False
			discovered_pi_dicti = {}

			for tissue_number in range(10):
				ctwas_result_file = simulated_ctwas_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_tissue' + str(tissue_number) + '/testing.susieIrss.txt'
				if os.path.isfile(ctwas_result_file) == False:
					missing = True
					continue

				f = open(ctwas_result_file)
				head_count = 0
				for line in f:
					line = line.rstrip()
					data = line.split('\t')
					if head_count == 0:
						head_count = head_count + 1
						continue
					component_window_name = data[4]
					component_num = 'NA'
					genetic_element_name = data[1]
					if genetic_element_name.startswith('ENSG') == False:
						continue
					genetic_element_pip = float(data[7])
					name_info = genetic_element_name.split('_')
					if len(name_info) != 3:
						print('assumtion erorro')
						pdb.set_trace()
					genetic_element_name = name_info[0] + '_' + name_info[1] + name_info[2]
					ensamble_id = name_info[0]
					if ensamble_id not in discovered_pi_dicti:
						discovered_pi_dicti[ensamble_id] = genetic_element_pip
					else:
						# Take max
						old_pip = discovered_pi_dicti[ensamble_id]
						if genetic_element_pip > old_pip:
							discovered_pi_dicti[ensamble_id] = genetic_element_pip
				f.close()
				discovered_dicti[ln_pi_method + '_' + twas_method] = discovered_pi_dicti
			if missing == True:
				continue

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


				for gene_name in middle_genes:
					if gene_name not in causal_genetic_elements:
						continue
					# THis is a causal gene
					booler = 0.0
					if gene_name in discovered_dicti[ln_pi_method + '_' + twas_method]:
						gene_pip = discovered_dicti[ln_pi_method + '_' + twas_method][gene_name]
						if gene_pip >= pip_threshold:
							booler = 1.0
					t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gene_name + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return




def create_file_containing_ctwas_tg_high_pip_snp_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_ctwas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, global_window_file, simulated_learned_gene_models_dir, gene_type='component_gene'):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	ln_pi_method='ctwas'
	twas_method='lasso'

	# Open output file
	t = open(cs_power_per_component_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\twas_method\tsimulation_number\twindow_name\tgenetic_element_class\tgenetic_element_name\tcausal_variant_in_cs\n')
	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		discovered_dicti = {}
		for eqtl_sample_size in eqtl_sample_sizes:
			ctwas_result_file = simulated_ctwas_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '/testing.susieIrss.txt'
			if os.path.isfile(ctwas_result_file) == False:
				print(simulation_number)
				continue

			discovered_pi_dicti = {}
			f = open(ctwas_result_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				component_window_name = data[4]
				component_num = 'NA'
				genetic_element_name = data[1]
				genetic_element_pip = float(data[7])
				if genetic_element_pip < pip_threshold:
					continue
				if genetic_element_name.startswith('ENSG'):
					name_info = genetic_element_name.split('_')
					if len(name_info) != 3:
						print('assumtion erorro')
						pdb.set_trace()
					genetic_element_name = name_info[0] + '_' + name_info[1] + name_info[2]
				discovered_pi_dicti[genetic_element_name] = 1
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
						booler = 0.0
						if gene_name in discovered_dicti[ln_pi_method + '_' + twas_method]:
							booler = 1.0
						t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'gene\t' + gene_name + '\t' + str(booler) + '\n')
				for variant_name in middle_variants:
					if variant_name not in causal_genetic_elements:
						continue
					# This is a causal variant
					booler = 0.0
					if variant_name in discovered_dicti[ln_pi_method + '_' + twas_method]:
						booler = 1.0
					t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + str(simulation_number) + '\t' + window_name + '\t' + 'variant\t' + variant_name + '\t' + str(booler) + '\n')
			f.close()
	t.close()
	return


def create_file_containing_tgfm_high_pip_snp_power_per_component_realistic_qtl_sample_size(global_simulation_name_string, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, simulated_learned_gene_models_dir,tiss_filter=False, gene_type='component_gene'):
	all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)
	all_genes, all_genes_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)

	eqtl_sample_sizes = ['low_eqtl_ss', 'high_eqtl_ss']


	# Open output file
	t = open(cs_power_per_component_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\twas_method\tsimulation_number\twindow_name\tgenetic_element_class\tgenetic_element_name\tcausal_variant_in_cs\n')
	# First loop through simulations
	for simulation_number in simulation_runs:
		# Extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_gene = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)


		# Extract causal eqtl sample size for this simulation
		causal_eqtl_ss_file = simulated_learned_gene_models_dir + 'simulation_' + str(simulation_number) + '/simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_permuted_eqtl_ss.npy'
		causal_sample_size = np.load(causal_eqtl_ss_file)[0]


		for eqtl_sample_size in eqtl_sample_sizes:
			# Create mapping from ln_pi_method to dictionary of elements in 95% cs for this simulation
			if eqtl_sample_size == 'low_eqtl_ss' and causal_sample_size > 200:
				continue
			if eqtl_sample_size == 'high_eqtl_ss' and causal_sample_size < 200:
				continue			

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
					cs_file = simulated_tgfm_results_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
					if tiss_filter:
						cs_file = simulated_tgfm_results_dir + 'simulation_realistic_qtl_ss_' + str(simulation_number) + '_' + global_simulation_name_string + '_all_t_' + gene_type  + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'
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


def create_file_containing_two_step_tgfm_high_pip_snp_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, n_samp='100', gene_type='component_gene'):
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
			two_step_tissue_summary_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_pmces_uniform_two_step_tissues.txt'
			best_tissue_arr, sig_tissue_arr = extract_two_step_analyzed_tissues(two_step_tissue_summary_file)
			tissue_vec = np.copy(best_tissue_arr)
			tissue_name = tissue_vec[0]

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
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_sampler_iterative_two_step_' + tissue_name + '_tgfm_pip_summary.txt'

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


def create_file_containing_tgfm_high_pip_snp_power_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, n_samp='100', gene_type='component_gene'):
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
					cs_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_nsamp_' + str(n_samp) + '_' + gene_type + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_tgfm_pip_summary.txt'

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
def create_file_containing_averaged_tgfm_cs_power_vary_ge_h2s(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_size, ln_pi_methods, twas_methods, ge_h2s):
	per_component_power = np.loadtxt(cs_power_per_component_output_file, dtype=str)[1:,:]
	t = open(cs_power_output_file,'w')
	t.write(('eQTL_sample_size\tGE_h2\tln_pi_method\ttwas_methods\tgenetic_element_class\tpower\tpower_lb\tpower_ub\n'))
	genetic_element_classes = np.asarray(['gene', 'variant'])
	for ge_h2 in ge_h2s:
		for ln_pi_method in ln_pi_methods:
			for twas_method in twas_methods:
				for element_class in genetic_element_classes:

					subset_indices = (per_component_power[:,0] == str(ge_h2)) & (per_component_power[:,1] == ln_pi_method) & (per_component_power[:,2] == twas_method) & ( per_component_power[:,5] == element_class)
					components_covered = (per_component_power[subset_indices,:][:,-1]).astype(float)
					prop = np.sum(components_covered)/len(components_covered)
					prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
					prop_lb = prop - (1.96*prop_se)
					prop_ub = prop + (1.96*prop_se)

					t.write(str(eqtl_sample_size) + '\t' + str(ge_h2) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + element_class + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\n')

	t.close()
	return

def create_file_containing_averaged_tgfm_cs_power_vary_gwas_ss(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_size, ln_pi_methods, twas_methods, gwas_sample_sizes):
	per_component_power = np.loadtxt(cs_power_per_component_output_file, dtype=str)[1:,:]
	t = open(cs_power_output_file,'w')
	t.write(('eQTL_sample_size\tGWAS_sample_size\tln_pi_method\ttwas_methods\tgenetic_element_class\tpower\tpower_lb\tpower_ub\n'))
	genetic_element_classes = np.asarray(['gene', 'variant'])
	for gwas_sample_size in gwas_sample_sizes:
		for ln_pi_method in ln_pi_methods:
			for twas_method in twas_methods:
				for element_class in genetic_element_classes:

					subset_indices = (per_component_power[:,0] == str(gwas_sample_size)) & (per_component_power[:,1] == ln_pi_method) & (per_component_power[:,2] == twas_method) & ( per_component_power[:,5] == element_class)
					components_covered = (per_component_power[subset_indices,:][:,-1]).astype(float)
					prop = np.sum(components_covered)/len(components_covered)
					prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
					prop_lb = prop - (1.96*prop_se)
					prop_ub = prop + (1.96*prop_se)

					t.write(str(eqtl_sample_size) + '\t' + str(gwas_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + element_class + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\n')

	t.close()
	return

def create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods, agg_eqtl_ss=False):
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

	if agg_eqtl_ss == True:
		for ln_pi_method in ln_pi_methods:
			for twas_method in twas_methods:
				for element_class in genetic_element_classes:
					subset_indices = (per_component_power[:,1] == ln_pi_method) & (per_component_power[:,2] == twas_method) & ( per_component_power[:,5] == element_class)
					components_covered = (per_component_power[subset_indices,:][:,-1]).astype(float)
					prop = np.sum(components_covered)/len(components_covered)
					prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
					prop_lb = prop - (1.96*prop_se)
					prop_ub = prop + (1.96*prop_se)

					t.write(str("Aggregate") + '\t' + ln_pi_method + '\t' + twas_method + '\t' + element_class + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\n')

	t.close()
	return

def create_file_containing_averaged_tgfm_cs_power_in_missing_causal_tissue_sim(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods, missingness_methods):
	per_component_power = np.loadtxt(cs_power_per_component_output_file, dtype=str)[1:,:]
	t = open(cs_power_output_file,'w')
	t.write(('eQTL_sample_size\tln_pi_method\ttwas_methods\tmissingness_method\tgenetic_element_class\tpower\tpower_lb\tpower_ub\n'))
	genetic_element_classes = np.asarray(['gene', 'variant'])
	for eqtl_sample_size in eqtl_sample_sizes:
		for ln_pi_method in ln_pi_methods:
			for twas_method in twas_methods:
				for missingness_method in missingness_methods:
					for element_class in genetic_element_classes:

						subset_indices = (per_component_power[:,0] == str(eqtl_sample_size)) & (per_component_power[:,1] == ln_pi_method) & (per_component_power[:,2] == twas_method)& (per_component_power[:,3] == missingness_method) & ( per_component_power[:,6] == element_class)
						components_covered = (per_component_power[subset_indices,:][:,-1]).astype(float)
						prop = np.sum(components_covered)/len(components_covered)
						prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
						prop_lb = prop - (1.96*prop_se)
						prop_ub = prop + (1.96*prop_se)

						t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + twas_method + '\t' + missingness_method + '\t' + element_class + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\n')

	t.close()
	return

def create_file_containing_averaged_causal_gene_tissue_pair_stratification_in_missing_causal_tissue_sim(causal_gene_tissue_pair_stratification_per_component_output_file, causal_gene_tissue_pair_stratification_output_file, eqtl_sample_sizes, missingness_methods):
	genetic_element_classes = ['causal_gt', 'best_tagging_gt', 'other_gt', 'nm_variant']
	per_component_power = np.loadtxt(causal_gene_tissue_pair_stratification_per_component_output_file, dtype=str)[1:,:]
	t = open(causal_gene_tissue_pair_stratification_output_file,'w')
	t.write(('eQTL_sample_size\tmissingness_method\tgenetic_element_class\tpower\tpower_lb\tpower_ub\n'))
	for eqtl_sample_size in eqtl_sample_sizes:
		for missingness_method in missingness_methods:
			for element_class in genetic_element_classes:
				subset_indices = (per_component_power[:,0] == str(eqtl_sample_size)) & (per_component_power[:,3] == missingness_method) & (per_component_power[:,8] == element_class)
				components_covered = (per_component_power[subset_indices,:][:,-1]).astype(float)
				prop = np.sum(components_covered)/len(components_covered)
				prop_se = np.sqrt((prop)*(1.0-prop))/np.sqrt(len(components_covered))
				prop_lb = prop - (1.96*prop_se)
				prop_ub = prop + (1.96*prop_se)

				t.write(str(eqtl_sample_size) + '\t' + missingness_method + '\t' + element_class + '\t' + str(prop) + '\t' + str(prop_lb) + '\t' + str(prop_ub) + '\n')

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
				tgfm_results_pkl = simulated_tgfm_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_nsamp_100_component_gene' + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method + '_' + window_name + '_results.pkl'

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
			tgfm_results_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_nsamp_100_component_gene' + '_eqtl_ss_' + str(eqtl_sample_size) + '_' + twas_method + '_' + ln_pi_method  + '_tgfm_pip_summary.txt'
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


def get_simulation_runs(simulated_tgfm_results_dir):
	valid_sim_runs = []

	for sim_run in range(1,101):
		eqtl_ss = '300'
		filer = simulated_tgfm_results_dir + 'simulation_' + str(sim_run) + '_chrom1_cis_window_100000_ss_100000_ge_h2_075_eqtl_ss_' + eqtl_ss + '_susie_pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped.txt'
		if os.path.isfile(filer) == False:
			continue
		eqtl_ss = '500'
		filer = simulated_tgfm_results_dir + 'simulation_' + str(sim_run) + '_chrom1_cis_window_100000_ss_100000_ge_h2_075_eqtl_ss_' + eqtl_ss + '_susie_pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped.txt'
		if os.path.isfile(filer) == False:
			continue
		eqtl_ss = '1000'
		filer = simulated_tgfm_results_dir + 'simulation_' + str(sim_run) + '_chrom1_cis_window_100000_ss_100000_ge_h2_075_eqtl_ss_' + eqtl_ss + '_susie_pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped.txt'
		if os.path.isfile(filer) == False:
			continue
		valid_sim_runs.append(sim_run)



	return np.asarray(valid_sim_runs)


def extract_window_variant_pips_for_comparison(susie_window_summary_file, tgfm_window_summary_file):
	# Load in TGFM results file
	g = open(tgfm_window_summary_file, "rb")
	tgfm_res = pickle.load(g)
	g.close()

	# Load in susie results file
	g = open(susie_window_summary_file, "rb")
	susie_res = pickle.load(g)
	g.close()


	tgfm_pips = tgfm_res['expected_beta_pips'][tgfm_res['middle_variant_indices']]
	susie_pips = susie_res['pips'][susie_res['middle_variant_indices']]
	# abs correlation thresh: .25, gene pip thresh .05
	valid_indices = susie_res['valid_indices_mat'][0,:][susie_res['middle_variant_indices']]


	return susie_pips, tgfm_pips, valid_indices


def create_file_summarizing_tgfm_variant_pip_and_susie_variant_pip_comparison(eqtl_sample_size, simulation_runs, simulated_tgfm_results_dir, global_simulation_name_string, variant_comparison_file):
	tgfm_pip_arr = []
	susie_pip_arr = []
	valid_index_arr = []

	for simulation_run in simulation_runs:
		susie_summary_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_variant_only_uniform_tgfm_pip_summary.txt'
		f = open(susie_summary_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			window_name = data[0]

			susie_window_summary_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_variant_only_uniform_' + window_name + '_results.pkl'
			tgfm_window_summary_file = simulated_tgfm_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_susie_sampler_pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped_' + window_name + '_results.pkl'

			if os.path.isfile(susie_window_summary_file) == False or os.path.isfile(tgfm_window_summary_file) == False:
				continue

			window_susie_pips, window_tgfm_pips, valid_snps = extract_window_variant_pips_for_comparison(susie_window_summary_file, tgfm_window_summary_file)

			tgfm_pip_arr.append(window_tgfm_pips)
			susie_pip_arr.append(window_susie_pips)
			valid_index_arr.append(valid_snps)

		f.close()

	tgfm_pip_arr = np.hstack(tgfm_pip_arr)
	susie_pip_arr = np.hstack(susie_pip_arr)
	valid_index_arr = np.hstack(valid_index_arr)
	
	print(np.corrcoef(tgfm_pip_arr, susie_pip_arr)[0,1])
	print(np.corrcoef(tgfm_pip_arr[valid_index_arr], susie_pip_arr[valid_index_arr])[0,1])
	print(np.corrcoef(tgfm_pip_arr[valid_index_arr==False], susie_pip_arr[valid_index_arr==False])[0,1])


	reg = LinearRegression().fit(tgfm_pip_arr.reshape((len(tgfm_pip_arr),1)), susie_pip_arr)
	print(reg.coef_)
	reg = LinearRegression().fit(tgfm_pip_arr[valid_index_arr].reshape((len(tgfm_pip_arr[valid_index_arr]),1)), susie_pip_arr[valid_index_arr])
	print(reg.coef_)
	reg = LinearRegression().fit(tgfm_pip_arr[valid_index_arr==False].reshape((len(tgfm_pip_arr[valid_index_arr==False]),1)), susie_pip_arr[valid_index_arr==False])
	print(reg.coef_)


	# Print to output file
	t = open(variant_comparison_file,'w')
	t.write('susie_pip\ttgfm_pip\tcorrelated_with_causal_gene\n')

	# Distinguish low and high pip snps
	low_pip_snps = (tgfm_pip_arr < .01) & (susie_pip_arr < .01)
	high_pip_snps = low_pip_snps == False

	# First print high pip snps
	tgfm_pip_arr_high = tgfm_pip_arr[high_pip_snps]
	susie_pip_arr_high = susie_pip_arr[high_pip_snps]
	valid_index_arr_high = valid_index_arr[high_pip_snps]
	for ii in range(len(tgfm_pip_arr_high)):
		if valid_index_arr_high[ii] == True:
			stringer = 'False'
		else:
			stringer = 'True'
		t.write(str(susie_pip_arr_high[ii]) + '\t' + str(tgfm_pip_arr_high[ii]) + '\t' + str(stringer) + '\n')

	# Now print subsampled low pip snps
	tgfm_pip_arr_low = tgfm_pip_arr[low_pip_snps]
	susie_pip_arr_low = susie_pip_arr[low_pip_snps]
	valid_index_arr_low = valid_index_arr[low_pip_snps]
	# Subsample
	n_low_pip_snps = len(tgfm_pip_arr_low)
	desired_n_low_pip_snps = 100000
	subsamp_indices = np.random.choice(np.arange(n_low_pip_snps), size=desired_n_low_pip_snps, replace=False)
	tgfm_pip_arr_low_sub = tgfm_pip_arr_low[subsamp_indices]
	susie_pip_arr_low_sub = susie_pip_arr_low[subsamp_indices]
	valid_index_arr_low_sub = valid_index_arr_low[subsamp_indices]

	for ii in range(len(tgfm_pip_arr_low_sub)):
		if valid_index_arr_low_sub[ii] == True:
			stringer = 'False'
		else:
			stringer = 'True'
		t.write(str(susie_pip_arr_low_sub[ii]) + '\t' + str(tgfm_pip_arr_low_sub[ii]) + '\t' + str(stringer) + '\n')

	t.close()
	return


def aggregate_average_prior_value_by_tissue(average_prior_value_by_tissue_output_file, agg_average_prior_value_by_tissue_output_file, eqtl_sample_sizes):
	aa = np.loadtxt(average_prior_value_by_tissue_output_file, dtype=str,delimiter='\t')
	aa = aa[1:,:]

	t = open(agg_average_prior_value_by_tissue_output_file,'w')
	t.write('eQTL_sample_size\ttissue_number\taverage_prior_value\tse_average_prior_value\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		for tissue_number in range(10):
			indices = (aa[:,0] == str(eqtl_sample_size)) & (aa[:,2] == str(tissue_number))
			priors =  aa[indices,3].astype(float)
			se_mean = np.std(priors)/np.sqrt(len(priors))
			t.write(str(eqtl_sample_size) + '\t' + str(tissue_number) + '\t' + str(np.mean(priors)) + '\t' + str(se_mean) + '\n')
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
simulated_organized_results_dir = sys.argv[7]
simulated_tgfm_results_dir = sys.argv[8]
simulated_trait_dir = sys.argv[9]
simulated_gene_expression_dir = sys.argv[10]
simulated_learned_gene_models_dir = sys.argv[11]
simulated_tgfm_input_data_dir = sys.argv[12]
simulated_gene_position_file = sys.argv[13]
processed_genotype_data_dir = sys.argv[14]
simulated_focus_results_dir = sys.argv[15]
simulated_coloc_results_dir = sys.argv[16]
simulated_causal_twas_results_dir = sys.argv[17]
simulated_two_step_tgfm_results_dir = sys.argv[18]

processed_genotype_data_dir = processed_genotype_data_dir + 'gwas_sample_size_' + str(n_gwas_individuals) + '/'



#############################################################
# Genome-wide analysis
#############################################################
# Simulation architecture
gene_trait_architecture = '2_caus_t'
eqtl_architecture = 'default'
local_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture

#Bim file
bim_file = processed_genotype_data_dir + 'simulated_gwas_data_' + chrom_num + '.bim'



# Used eQTL sample sizes
eqtl_sample_sizes = np.asarray(['realistic', '100','300', '500' ,'1000'])

# Simulation runs
# Currently hacky because had some failed simulations
# Simulation runs
# Currently hacky because had some failed simulations
tmp_simulation_runs = np.arange(1,101)
bads = {}
bads[55] = 1
bads[73] = 1
simulation_runs = []
for sim_run in tmp_simulation_runs:
	if sim_run not in bads:
		simulation_runs.append(sim_run)
simulation_runs = np.asarray(simulation_runs)
thresholds = [1e-1, 1e-2, 1e-3, 5e-4, 2.5e-4, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-14, 1e-16, 1e-18]

model_names = ['susie_pmces_uniform_iterative_variant_gene_prior_pip_level']

for model_name in model_names:
	##########################
	# Power and Type 1 Error
	##########################
	# Create file showing p-value in causal tissues and non-causal tissues
	mediated_iterative_sampler_pvalue_by_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + model_name + '_mediated_h2_pvalue_by_tissue_est_across_thresholds.txt'
	create_file_containing_iterative_bootstrapped_sampler_prior_mediated_h2_pvalue_in_causal_and_non_causal_tissues_across_thresholds(simulated_tgfm_results_dir, local_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_iterative_sampler_pvalue_by_tissue_output_file, thresholds, model_name)
	# Create file showing power to detect causal tissues
	power_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + model_name + '_mediated_h2_power_across_thresholds.txt'
	create_file_containing_mediated_h2_power_to_detect_causal_tissues_across_thresholds(mediated_iterative_sampler_pvalue_by_tissue_output_file, power_output_file, eqtl_sample_sizes, thresholds)
	# Create file showing type 1 error for null tissues
	type_1_error_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + model_name + '_h2_type_1_error_across_thresholds.txt'
	create_file_containing_mediated_h2_type_1_error_across_thresholds(mediated_iterative_sampler_pvalue_by_tissue_output_file, type_1_error_output_file, eqtl_sample_sizes, thresholds)

	##########################
	# Power and Type 1 Error based on gaussian approximation
	##########################
	# Create file showing p-value in causal tissues and non-causal tissues
	mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + model_name + '_mediated_h2_gaussian_approximation_pvalue_by_tissue_est.txt'
	create_file_containing_iterative_bootstrapped_sampler_prior_mediated_h2_gaussian_approximation_pvalue_in_causal_and_non_causal_tissues(simulated_tgfm_results_dir, local_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file, model_name)
	# Create file showing type 1 error for null tissues using gaussian approximation
	type_1_error_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + model_name + '_h2_type_1_error_gaussian_approximation.txt'
	create_file_containing_mediated_h2_type_1_error_gaussian_approximation(mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file, type_1_error_output_file, eqtl_sample_sizes)
	# Create file showing power for null tissues using gaussian approximation
	power_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + model_name + '_mediated_h2_power_gaussian_approximation.txt'
	create_file_containing_mediated_h2_power_to_detect_causal_tissues_gaussian_approximation(mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file, power_output_file, eqtl_sample_sizes)

	##########################
	# Bias in fraction of causal components going through gene expression
	##########################
	# Create file showing mediated h2 in causal tissues and non-causal tissues
	fraction_expr_med_disease_components_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + model_name + '_expected_fraction_expression_mediated_disease_components.txt'
	#create_file_containing_expected_fraction_of_expression_mediated_disease_components(simulated_tgfm_results_dir, local_simulation_name_string, eqtl_sample_sizes, simulation_runs, model_name, simulated_ld_scores_dir, fraction_expr_med_disease_components_output_file)
	# Average estimates across simulations
	organized_avg_fraction_expr_med_disease_components_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + model_name + '_avg_expected_fraction_expression_mediated_disease_components.txt'
	#create_file_containing_avg_fraction_h2_across_simulation_runs(fraction_expr_med_disease_components_output_file, organized_avg_fraction_expr_med_disease_components_output_file, eqtl_sample_sizes)
	#print(organized_avg_fraction_expr_med_disease_components_output_file)

	##########################
	# Bias in fraction of causal genes per tissue
	##########################
	# Create file showing mediated h2 in causal tissues and non-causal tissues
	fraction_causal_by_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + model_name + '_fraction_causal_by_tissue.txt'
	create_file_containing_fraction_causal_in_causal_and_non_causal_tissues(simulated_tgfm_results_dir, local_simulation_name_string, eqtl_sample_sizes, simulation_runs, model_name, fraction_causal_by_tissue_output_file)
	# Average estimates across simulations
	organized_avg_fraction_causal_by_tissue_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + model_name + '_avg_fraction_causal_by_tissue.txt'
	create_file_containing_avg_fraction_causal_by_tissue_across_simulation_runs(fraction_causal_by_tissue_output_file, organized_avg_fraction_causal_by_tissue_output_file, eqtl_sample_sizes)
	# Average estimates across simulations
	organized_avg_fraction_causal_by_tissue_causal_status_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + model_name + '_avg_fraction_causal_by_tissue_causal_status.txt'
	create_file_containing_avg_fraction_causal_by_tissue_causal_status_across_simulation_runs(fraction_causal_by_tissue_output_file, organized_avg_fraction_causal_by_tissue_causal_status_output_file, eqtl_sample_sizes)
	print(organized_avg_fraction_causal_by_tissue_causal_status_output_file)


#############################################################
# Fraction of TGFM elements mediated by gene expression
#############################################################
# Used eQTL sample sizes
eqtl_sample_sizes = np.asarray(['realistic', '100','300', '500' ,'1000'])


# Thresholds to consider
pip_thresholds=np.asarray([.1, .3, .5, .7, .9])

# Model specification
twas_method = 'susie_sampler'
ln_pi_method = 'pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'

# Windows
global_window_file = processed_genotype_data_dir + 'chromosome_' + str(chrom_num) + '_windows_3_mb.txt'

# Create file containing fraction of high pip elememnts mediated by gene expresssion
fraction_high_pip_elements_mediated_by_gene_expression_summary_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + twas_method + '_' + ln_pi_method + '_fraction_high_pip_elements_mediated_by_gene_expresssion.txt'
create_file_containing_fraction_of_high_pip_tgfm_elements_mediated_by_gene_expression(global_window_file, eqtl_sample_sizes, simulation_runs, pip_thresholds, twas_method, ln_pi_method, local_simulation_name_string, simulated_tgfm_results_dir, fraction_high_pip_elements_mediated_by_gene_expression_summary_file)

# Create file containing expected fraction of elememnts mediated by gene expresssion
expected_fraction_of_elements_mediated_by_gene_expression_summary_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + twas_method + '_' + ln_pi_method + '_expected_fraction_elements_mediated_by_gene_expresssion.txt'
create_file_containing_expected_fraction_of_tgfm_elements_mediated_by_gene_expression(global_window_file, eqtl_sample_sizes, simulation_runs, twas_method, ln_pi_method, local_simulation_name_string, simulated_tgfm_results_dir, expected_fraction_of_elements_mediated_by_gene_expression_summary_file)


#############################################################
# Default Fine-mapping simulation
#############################################################
gene_trait_architecture = '2_caus_t'
eqtl_architecture = 'default'
local_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture


# Used eQTL sample sizes
eqtl_sample_sizes = np.asarray(['realistic', '100','300', '500' ,'1000'])

# Simulation runs
# Currently hacky because had some failed simulations
tmp_simulation_runs = np.arange(1,101)
bads = {}
bads[55] = 1
bads[73] = 1
simulation_runs = []
for sim_run in tmp_simulation_runs:
	if sim_run not in bads:
		simulation_runs.append(sim_run)
simulation_runs = np.asarray(simulation_runs)


# ln_pi methods used
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
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods)
	print(cs_high_pip_coverage_output_file)

##################################
# Coverage/Calibration to detect snps with PIP > threshold
# Non-mediated variants also include expression-mediated variants
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_nm_variant_alt_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp_non_mediated_variants_include_gene_variants(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_gene_expression_dir)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_nm_variant_alt_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods)


##################################
# Coverage/Calibration to detect snps with PIP > threshold
# Non-mediated variants also include expression-mediated variants
# Limit to realistic simulation and stratify into low and high eqtl sample size
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_mixed_100_300_ss_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp_non_mediated_variants_include_gene_variants_mixed_ss(local_simulation_name_string, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_gene_expression_dir, simulated_learned_gene_models_dir)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_mixed_100_300_ss_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, ['low','high'], ln_pi_methods, twas_methods)
	print(cs_high_pip_coverage_output_file)

##################################
# Power to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]

for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_tgfm_high_pip_snp_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods,twas_methods)

##################################
# Coverage/Calibration to detect genes (not gene-tissue pairs)
##################################
# ln_pi methods used
ln_pi_methods = np.asarray(['pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])

# twas method
twas_methods = np.asarray(['susie_sampler'])
pip_thresholds = [.3, .5, .7, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_gene_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_genes(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_input_data_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_gene_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods)
	print(cs_high_pip_coverage_output_file)

##################################
# Power to detect genes (not gene-tissue pairs)
##################################
# ln_pi methods used
ln_pi_methods = np.asarray(['pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])

# twas method
twas_methods = np.asarray(['susie_sampler'])
pip_thresholds = [.3, .5, .7, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_gene_power_per_component.txt'
	create_file_containing_tgfm_high_pip_gene_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_gene_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods,twas_methods)


##################################
# TGFM FDR-Power Curve
##################################
twas_method='susie_sampler'
ln_pi_method='pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'
gene_type='component_gene'
for eqtl_ss in eqtl_sample_sizes:
	print(eqtl_ss)
	fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_tgfm_gene_tissue_fdr_power_curve_raw_data.txt'
	fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_tgfm_gene_tissue_fdr_power_curve_data.txt'
	create_file_containing_tgfm_gene_tissue_fdr_power_curve_data(local_simulation_name_string, simulation_runs, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, eqtl_ss, gene_type=gene_type)



tmp_simulation_runs = np.arange(1,101)
bads = {}
bads[55] = 1
bads[73] = 1
simulation_runs = []
for sim_run in tmp_simulation_runs:
	if sim_run not in bads:
		simulation_runs.append(sim_run)
simulation_runs = np.asarray(simulation_runs)



# ln_pi methods used
ln_pi_methods = np.asarray(['pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])

# twas method
twas_methods = np.asarray(['susie_sampler'])

##################################
# Two step TGFM Coverage/Calibration to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_two_step_tgfm_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	#create_file_containing_two_step_tgfm_cs_calibration_per_high_pip_snp(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_two_step_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_two_step_tgfm_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration_two_step(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ['susie_sampler'], ['iterative'], ['best_tissue', 'significant_tissues'])
	print(cs_high_pip_coverage_output_file)

##################################
# Two step TGFM power to detect gene-tissue pairs with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_two_step_tgfm_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_two_step_tgfm_high_pip_snp_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_two_step_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_two_step_tgfm_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods,twas_methods)
	print(cs_power_output_file)




##################################
# Two step TGFM FDR-power curve
##################################
two_step_tissue_methods = ['best_tissue', 'significant_tissues']
for eqtl_ss in eqtl_sample_sizes:
	for two_step_method in two_step_tissue_methods:
		fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + two_step_method + '_' + eqtl_ss + '_two_step_tgfm_gene_tissue_fdr_power_curve_raw_data.txt'
		fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + two_step_method + '_' + eqtl_ss + '_two_step_tgfm_gene_tissue_fdr_power_curve_data.txt'
		create_file_containing_two_step_tgfm_gene_tissue_fdr_power_curve_data(local_simulation_name_string, simulation_runs, simulated_trait_dir, simulated_two_step_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, two_step_method, global_window_file, eqtl_ss)
		print(fdr_power_curve_data_output_file)

print('cTWAS')
# Simulation runs
# Currently hacky because had some failed simulations
tmp_simulation_runs = np.arange(1,101)
bads = {}
bads[55] = 1
bads[73] = 1
simulation_runs = []
for sim_run in tmp_simulation_runs:
	if sim_run not in bads:
		simulation_runs.append(sim_run)
simulation_runs = np.asarray(simulation_runs)
##################################
# cTWAS Coverage/Calibration to detect snps with PIP > threshold
##################################
gene_type='component_gene'
pip_thresholds = [.3, .5, .7, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_ctwas_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	#create_file_containing_ctwas_cs_calibration_per_high_pip_snp(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_causal_twas_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_learned_gene_models_dir)
	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_ctwas_pip_' + str(pip_threshold) + '_calibration.txt'
	#create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ['ctwas'], ['lasso'], agg_eqtl_ss=False)

	# gene fine-mapping
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_ctwas_pip_' + str(pip_threshold) + '_gene_calibration_per_component.txt'
	#create_file_containing_ctwas_cs_calibration_per_high_pip_gene(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_causal_twas_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_learned_gene_models_dir)
	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_ctwas_pip_' + str(pip_threshold) + '_gene_calibration.txt'
	#create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ['ctwas'], ['lasso'], agg_eqtl_ss=False)

##################################
# cTWAS Power to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
gene_type='component_gene'
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_ctwas_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_ctwas_high_pip_snp_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_causal_twas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, global_window_file,simulated_learned_gene_models_dir)
	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_ctwas_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ['ctwas'],['lasso'], agg_eqtl_ss=False)

	# gene power
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_ctwas_pip_' + str(pip_threshold) + '_gene_power_per_component.txt'
	create_file_containing_ctwas_high_pip_gene_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_causal_twas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, global_window_file,simulated_learned_gene_models_dir)
	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_ctwas_pip_' + str(pip_threshold) + '_gene_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ['ctwas'],['lasso'], agg_eqtl_ss=False)
	print(cs_power_output_file)

##################################
# cTWAS FDR-power curve
##################################
# ln_pi methods used
ln_pi_method = 'ctwas'
# twas method
twas_method = 'lasso'
for eqtl_ss in eqtl_sample_sizes:
	fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'ctwas' + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_fdr_power_curve_raw_data.txt'
	fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'ctwas' + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_fdr_power_curve_data.txt'
	create_file_containing_ctwas_gene_tissue_fdr_power_curve_data(local_simulation_name_string, eqtl_ss, simulation_runs, simulated_trait_dir, simulated_causal_twas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir)
	print(fdr_power_curve_data_output_file)


##################################
# cTWAS-TG Coverage/Calibration to detect snps with PIP > threshold
##################################
gene_type='component_gene'
pip_thresholds = [.3, .5, .7, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_ctwas_tg_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	#create_file_containing_ctwas_tg_cs_calibration_per_high_pip_snp(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_causal_twas_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_learned_gene_models_dir)
	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_ctwas_tg_pip_' + str(pip_threshold) + '_calibration.txt'
	#create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ['ctwas'], ['lasso'], agg_eqtl_ss=False)

	# Gene Calibration
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_ctwas_tg_pip_' + str(pip_threshold) + '_gene_calibration_per_component.txt'
	create_file_containing_ctwas_tg_cs_calibration_per_high_pip_gene(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_causal_twas_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_learned_gene_models_dir)
	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_ctwas_tg_pip_' + str(pip_threshold) + '_gene_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ['ctwas'], ['lasso'], agg_eqtl_ss=False)
	print(cs_high_pip_coverage_output_file)

##################################
# cTWAS-TG Power to detect snps with PIP > threshold
##################################

pip_thresholds = [.3, .5, .7, .9, .95, .99]
gene_type='component_gene'
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_ctwas_tg_pip_' + str(pip_threshold) + '_power_per_component.txt'
	#create_file_containing_ctwas_tg_high_pip_snp_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_causal_twas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, global_window_file,simulated_learned_gene_models_dir)
	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_ctwas_tg_pip_' + str(pip_threshold) + '_power.txt'
	#create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ['ctwas'],['lasso'], agg_eqtl_ss=False)

	# cTWAS-TG gene power
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_ctwas_tg_pip_' + str(pip_threshold) + '_gene_power_per_component.txt'
	create_file_containing_ctwas_tg_high_pip_gene_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_causal_twas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, global_window_file,simulated_learned_gene_models_dir)
	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_ctwas_tg_pip_' + str(pip_threshold) + '_gene_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ['ctwas'],['lasso'], agg_eqtl_ss=False)

##################################
# cTWAS-TG FDR-power curve
##################################
# ln_pi methods used
ln_pi_method = 'ctwas'
# twas method
twas_method = 'lasso'
for eqtl_ss in eqtl_sample_sizes:
	fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'ctwas_tg' + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_fdr_power_curve_raw_data.txt'
	fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'ctwas_tg' + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_fdr_power_curve_data.txt'
	create_file_containing_ctwas_tg_gene_tissue_fdr_power_curve_data(local_simulation_name_string, eqtl_ss, simulation_runs, simulated_trait_dir, simulated_causal_twas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir)
	print(fdr_power_curve_data_output_file)

##################################
# Two-step cTWAS Coverage/Calibration to detect snps with PIP > threshold
##################################
tmp_simulation_runs = np.arange(1,101)
bads = {}
bads[55] = 1
bads[73] = 1
simulation_runs = []
for sim_run in tmp_simulation_runs:
	if sim_run not in bads:
		simulation_runs.append(sim_run)
simulation_runs = np.asarray(simulation_runs)
gene_type='component_gene'
pip_thresholds = [.3, .5, .7, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_two_step_ctwas_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_two_step_ctwas_cs_calibration_per_high_pip_snp(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_two_step_tgfm_results_dir, simulated_causal_twas_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_learned_gene_models_dir)
	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_two_step_ctwas_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ['ctwas'], ['lasso'], agg_eqtl_ss=False)

##################################
# Two-step cTWAS Power to detect gene-tissue pairs with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
gene_type='component_gene'
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_two_step_ctwas_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_two_step_ctwas_high_pip_snp_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_two_step_tgfm_results_dir, simulated_causal_twas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, global_window_file,simulated_learned_gene_models_dir)
	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_two_step_ctwas_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ['ctwas'],['lasso'], agg_eqtl_ss=False)
	print(cs_power_output_file)


##################################
# Two-step cTWAS FDR-power curve
##################################
# ln_pi methods used
ln_pi_method = 'ctwas'
# twas method
twas_method = 'lasso'
for eqtl_ss in eqtl_sample_sizes:
	fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'two_step_ctwas' + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_fdr_power_curve_raw_data.txt'
	fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'two_step_ctwas' + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_fdr_power_curve_data.txt'
	create_file_containing_two_step_ctwas_gene_tissue_fdr_power_curve_data(local_simulation_name_string, eqtl_ss, simulation_runs, simulated_trait_dir,simulated_two_step_tgfm_results_dir, simulated_causal_twas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir)
	print(fdr_power_curve_data_output_file)




print('FOCUS')

##################################
# FOCUS Coverage/Calibration to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Calibration to detect causal gene-tissue pairs (no genotype intercept)
	focus_cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_focus_cs_calibration_per_high_pip_snp(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, pip_threshold, focus_cs_coverage_per_high_pip_snp_output_file)
	focus_cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_focus_high_pip_calibration(focus_cs_coverage_per_high_pip_snp_output_file, focus_cs_high_pip_coverage_output_file, eqtl_sample_sizes)
	
	# Calibration to detect causal gene-tissue pairs (genotype intercept)
	focus_cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_intercept_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_focus_cs_calibration_per_high_pip_snp(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, pip_threshold, focus_cs_coverage_per_high_pip_snp_output_file, intercept=True)
	focus_cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_intercept_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_focus_high_pip_calibration(focus_cs_coverage_per_high_pip_snp_output_file, focus_cs_high_pip_coverage_output_file, eqtl_sample_sizes)

	# Calibration to detect causal genes (not gene tissue pairs) (no genotype intercept)
	focus_cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_pip_' + str(pip_threshold) + '_gene_calibration_per_component.txt'
	create_file_containing_focus_cs_calibration_per_high_pip_genes(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, pip_threshold, focus_cs_coverage_per_high_pip_snp_output_file)
	focus_cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_pip_' + str(pip_threshold) + '_gene_calibration.txt'
	create_file_containing_averaged_focus_high_pip_calibration(focus_cs_coverage_per_high_pip_snp_output_file, focus_cs_high_pip_coverage_output_file, eqtl_sample_sizes)

	# Calibration to detect causal genes (not gene tissue pairs) (genotype intercept)
	focus_cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_intercept_pip_' + str(pip_threshold) + '_gene_calibration_per_component.txt'
	create_file_containing_focus_cs_calibration_per_high_pip_genes(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, pip_threshold, focus_cs_coverage_per_high_pip_snp_output_file, intercept=True)
	focus_cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_intercept_pip_' + str(pip_threshold) + '_gene_calibration.txt'
	create_file_containing_averaged_focus_high_pip_calibration(focus_cs_coverage_per_high_pip_snp_output_file, focus_cs_high_pip_coverage_output_file, eqtl_sample_sizes)


##################################
# FOCUS Power to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Power to detect causal gene-tissue pairs (no genotype intercept)
	focus_cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_focus_high_pip_snp_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file)
	focus_cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_focus_cs_power(focus_cs_power_per_component_output_file, focus_cs_power_output_file, eqtl_sample_sizes)

	# Power to detect causal gene-tissue pairs (genotype intercept)
	focus_cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_intercept_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_focus_high_pip_snp_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file, intercept=True)
	focus_cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_intercept_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_focus_cs_power(focus_cs_power_per_component_output_file, focus_cs_power_output_file, eqtl_sample_sizes)

	# Power to detect causal genes (not gene tissue pairs) (no genotype intercept)
	focus_cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_pip_' + str(pip_threshold) + '_gene_power_per_component.txt'
	create_file_containing_focus_high_pip_gene_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file)
	focus_cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_pip_' + str(pip_threshold) + '_gene_power.txt'
	create_file_containing_averaged_focus_cs_power(focus_cs_power_per_component_output_file, focus_cs_power_output_file, eqtl_sample_sizes)

	# Power to detect causal genes (not gene tissue pairs) (genotype intercept)
	focus_cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_intercept_pip_' + str(pip_threshold) + '_gene_power_per_component.txt'
	create_file_containing_focus_high_pip_gene_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file, intercept=True)
	focus_cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_intercept_pip_' + str(pip_threshold) + '_gene_power.txt'
	create_file_containing_averaged_focus_cs_power(focus_cs_power_per_component_output_file, focus_cs_power_output_file, eqtl_sample_sizes)


##################################
# FOCUS FDR-power curve
##################################
for eqtl_ss in eqtl_sample_sizes:
	# W/O intercept
	fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'focus' + '_' + eqtl_ss + '_fdr_power_curve_raw_data.txt'
	fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'focus' + '_' + eqtl_ss + '_fdr_power_curve_data.txt'
	create_file_containing_focus_gene_tissue_fdr_power_curve_data(local_simulation_name_string, eqtl_ss, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir)
	print(fdr_power_curve_data_output_file)

	# W intercept
	fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'focus_intercept' + '_' + eqtl_ss + '_fdr_power_curve_raw_data.txt'
	fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'focus_intercept' + '_' + eqtl_ss + '_fdr_power_curve_data.txt'
	create_file_containing_focus_gene_tissue_fdr_power_curve_data(local_simulation_name_string, eqtl_ss, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir, intercept=True)
	print(fdr_power_curve_data_output_file)


##################################
# FOCUS-TG Coverage/Calibration to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]

# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Calibration to detect causal gene-tissue pairs (no genotype intercept)
	# Create file with one line per cs (columns: eQTL_sample_size, simulation_num, window_num, boolean_causal_variant_in_cs)
	focus_cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_tg_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_focus_tg_cs_calibration_per_high_pip_snp(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, pip_threshold, focus_cs_coverage_per_high_pip_snp_output_file)
	focus_cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_tg_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_focus_high_pip_calibration(focus_cs_coverage_per_high_pip_snp_output_file, focus_cs_high_pip_coverage_output_file, eqtl_sample_sizes)

	# Calibration to detect causal gene-tissue pairs (genotype intercept)
	# Create file with one line per cs (columns: eQTL_sample_size, simulation_num, window_num, boolean_causal_variant_in_cs)
	focus_cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_tg_intercept_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_focus_tg_cs_calibration_per_high_pip_snp(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, pip_threshold, focus_cs_coverage_per_high_pip_snp_output_file, intercept=True)
	focus_cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_tg_intercept_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_focus_high_pip_calibration(focus_cs_coverage_per_high_pip_snp_output_file, focus_cs_high_pip_coverage_output_file, eqtl_sample_sizes)

	# Calibration to detect causal genes (not gene tissue pairs; no genotype intercept)
	focus_cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_tg_pip_' + str(pip_threshold) + '_gene_calibration_per_component.txt'
	create_file_containing_focus_tg_cs_calibration_per_high_pip_genes(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, pip_threshold, focus_cs_coverage_per_high_pip_snp_output_file)
	focus_cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_tg_pip_' + str(pip_threshold) + '_gene_calibration.txt'
	create_file_containing_averaged_focus_high_pip_calibration(focus_cs_coverage_per_high_pip_snp_output_file, focus_cs_high_pip_coverage_output_file, eqtl_sample_sizes)

	# Calibration to detect causal genes (not gene tissue pairs; genotype intercept)
	focus_cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_tg_intercept_pip_' + str(pip_threshold) + '_gene_calibration_per_component.txt'
	create_file_containing_focus_tg_cs_calibration_per_high_pip_genes(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, pip_threshold, focus_cs_coverage_per_high_pip_snp_output_file, intercept=True)
	focus_cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_tg_intercept_pip_' + str(pip_threshold) + '_gene_calibration.txt'
	create_file_containing_averaged_focus_high_pip_calibration(focus_cs_coverage_per_high_pip_snp_output_file, focus_cs_high_pip_coverage_output_file, eqtl_sample_sizes)



##################################
# FOCUS-TG Power to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Power to detect causal gene-tissue pairs (no genotype intercept)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	focus_cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_tg_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_focus_tg_high_pip_snp_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file)
	focus_cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_tg_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_focus_cs_power(focus_cs_power_per_component_output_file, focus_cs_power_output_file, eqtl_sample_sizes)

	# Power to detect causal gene-tissue pairs (genotype intercept)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	focus_cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_tg_intercept_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_focus_tg_high_pip_snp_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file, intercept=True)
	focus_cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_tg_intercept_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_focus_cs_power(focus_cs_power_per_component_output_file, focus_cs_power_output_file, eqtl_sample_sizes)

	# Power to detect causal genes (not gene tissue pairs) (no genotype intercept)
	focus_cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_tg_pip_' + str(pip_threshold) + '_gene_power_per_component.txt'
	create_file_containing_focus_tg_high_pip_gene_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file)
	focus_cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_tg_pip_' + str(pip_threshold) + '_gene_power.txt'
	create_file_containing_averaged_focus_cs_power(focus_cs_power_per_component_output_file, focus_cs_power_output_file, eqtl_sample_sizes)

	# Power to detect causal genes (not gene tissue pairs) (genotype intercept)
	focus_cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_tg_intercept_pip_' + str(pip_threshold) + '_gene_power_per_component.txt'
	create_file_containing_focus_tg_high_pip_gene_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file, intercept=True)
	focus_cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_focus_tg_intercept_pip_' + str(pip_threshold) + '_gene_power.txt'
	create_file_containing_averaged_focus_cs_power(focus_cs_power_per_component_output_file, focus_cs_power_output_file, eqtl_sample_sizes)

##################################
# FOCUS-TG FDR-power curve
##################################
for eqtl_ss in eqtl_sample_sizes:
	# W/O intercept
	fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'focus_tg' + '_' + eqtl_ss + '_fdr_power_curve_raw_data.txt'
	fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'focus_tg' + '_' + eqtl_ss + '_fdr_power_curve_data.txt'
	create_file_containing_focus_tg_gene_tissue_fdr_power_curve_data(local_simulation_name_string, eqtl_ss, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir)

	# W intercept
	fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'focus_tg_intercept' + '_' + eqtl_ss + '_fdr_power_curve_raw_data.txt'
	fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'focus_tg_intercept' + '_' + eqtl_ss + '_fdr_power_curve_data.txt'
	create_file_containing_focus_tg_gene_tissue_fdr_power_curve_data(local_simulation_name_string, eqtl_ss, simulation_runs, simulated_trait_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir, intercept=True)

##################################
# Two-step FOCUS Coverage/Calibration to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Calibration to detect causal gene-tissue pairs (no genotype intercept)
	focus_cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_two_step_focus_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_two_step_focus_cs_calibration_per_high_pip_snp(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir,simulated_two_step_tgfm_results_dir, simulated_focus_results_dir, pip_threshold, focus_cs_coverage_per_high_pip_snp_output_file)
	focus_cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_two_step_focus_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_focus_high_pip_calibration(focus_cs_coverage_per_high_pip_snp_output_file, focus_cs_high_pip_coverage_output_file, eqtl_sample_sizes)
	print(focus_cs_high_pip_coverage_output_file)

##################################
# Two-step FOCUS Power to detect gene-tissue pairs with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Power to detect causal gene-tissue pairs (no genotype intercept)
	focus_cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_two_step_focus_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_two_step_focus_high_pip_snp_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir,simulated_two_step_tgfm_results_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, focus_cs_power_per_component_output_file, pip_threshold, global_window_file)
	focus_cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_two_step_focus_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_focus_cs_power(focus_cs_power_per_component_output_file, focus_cs_power_output_file, eqtl_sample_sizes)
	print(focus_cs_power_output_file)

##################################
# Two-step FOCUS FDR-power curve
##################################
for eqtl_ss in eqtl_sample_sizes:
	# W/O intercept
	fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'two_step_focus' + '_' + eqtl_ss + '_fdr_power_curve_raw_data.txt'
	fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'two_step_focus' + '_' + eqtl_ss + '_fdr_power_curve_data.txt'
	create_file_containing_two_step_focus_gene_tissue_fdr_power_curve_data(local_simulation_name_string, eqtl_ss, simulation_runs, simulated_trait_dir, simulated_two_step_tgfm_results_dir, simulated_focus_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir)




print('coloc')

##################################
# coloc Coverage/Calibration to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, simulation_num, window_num, boolean_causal_variant_in_cs)
	coloc_cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_coloc_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_coloc_cs_calibration_per_high_pip_snp(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_coloc_results_dir, pip_threshold, coloc_cs_coverage_per_high_pip_snp_output_file)

	coloc_cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_coloc_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_focus_high_pip_calibration(coloc_cs_coverage_per_high_pip_snp_output_file, coloc_cs_high_pip_coverage_output_file, eqtl_sample_sizes)

	# Create file with one line per cs (columns: eQTL_sample_size, simulation_num, window_num, boolean_causal_variant_in_cs)
	coloc_cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_coloc_pip_' + str(pip_threshold) + '_gene_calibration_per_component.txt'
	create_file_containing_coloc_cs_calibration_per_high_pip_gene(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_coloc_results_dir, pip_threshold, coloc_cs_coverage_per_high_pip_snp_output_file)

	coloc_cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_coloc_pip_' + str(pip_threshold) + '_gene_calibration.txt'
	create_file_containing_averaged_focus_high_pip_calibration(coloc_cs_coverage_per_high_pip_snp_output_file, coloc_cs_high_pip_coverage_output_file, eqtl_sample_sizes)


##################################
# coloc Power to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	coloc_cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_coloc_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_coloc_high_pip_snp_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_coloc_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, coloc_cs_power_per_component_output_file, pip_threshold, global_window_file)
	coloc_cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_coloc_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_focus_cs_power(coloc_cs_power_per_component_output_file, coloc_cs_power_output_file, eqtl_sample_sizes)

	# Power to detect causal genes (not gene tissue pairs)
	coloc_cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_coloc_pip_' + str(pip_threshold) + '_gene_power_per_component.txt'
	create_file_containing_coloc_high_pip_gene_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_coloc_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, coloc_cs_power_per_component_output_file, pip_threshold, global_window_file)
	coloc_cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_coloc_pip_' + str(pip_threshold) + '_gene_power.txt'
	create_file_containing_averaged_focus_cs_power(coloc_cs_power_per_component_output_file, coloc_cs_power_output_file, eqtl_sample_sizes)

##################################
# coloc FDR-power curve
##################################
for eqtl_ss in eqtl_sample_sizes:
	fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'coloc' + '_' + eqtl_ss + '_fdr_power_curve_raw_data.txt'
	fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'coloc' + '_' + eqtl_ss + '_fdr_power_curve_data.txt'
	create_file_containing_coloc_gene_tissue_fdr_power_curve_data(local_simulation_name_string, eqtl_ss, simulation_runs, simulated_trait_dir, simulated_coloc_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, global_window_file, simulated_learned_gene_models_dir)
	print(fdr_power_curve_data_output_file)



##################################
# Two-step coloc Coverage/Calibration to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, simulation_num, window_num, boolean_causal_variant_in_cs)
	coloc_cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_two_step_coloc_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_two_step_coloc_cs_calibration_per_high_pip_snp(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_two_step_tgfm_results_dir,simulated_coloc_results_dir, pip_threshold, coloc_cs_coverage_per_high_pip_snp_output_file)
	coloc_cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_two_step_coloc_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_focus_high_pip_calibration(coloc_cs_coverage_per_high_pip_snp_output_file, coloc_cs_high_pip_coverage_output_file, eqtl_sample_sizes)
	print(coloc_cs_high_pip_coverage_output_file)

##################################
# Two-step coloc Power to detect gene-tissue pairs with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	coloc_cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_two_step_coloc_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_two_step_coloc_high_pip_snp_power_per_component(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_two_step_tgfm_results_dir, simulated_coloc_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, coloc_cs_power_per_component_output_file, pip_threshold, global_window_file)
	coloc_cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_two_step_coloc_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_focus_cs_power(coloc_cs_power_per_component_output_file, coloc_cs_power_output_file, eqtl_sample_sizes)
	print(coloc_cs_power_output_file)

##################################
# Two-step coloc FDR-power curve
##################################
for eqtl_ss in eqtl_sample_sizes:
	fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'two_step_coloc' + '_' + eqtl_ss + '_fdr_power_curve_raw_data.txt'
	fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_' + 'two_step_coloc' + '_' + eqtl_ss + '_fdr_power_curve_data.txt'
	create_file_containing_two_step_coloc_gene_tissue_fdr_power_curve_data(local_simulation_name_string, eqtl_ss, simulation_runs, simulated_trait_dir, simulated_two_step_tgfm_results_dir, simulated_coloc_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, global_window_file, simulated_learned_gene_models_dir)
	print(fdr_power_curve_data_output_file)




'''
#############################################################
# Compare TGFM (variant) PIP with SuSiE variant PIPs
#############################################################
eqtl_sample_size = 500
simulation_runs = np.arange(1,21)
variant_comparison_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_variant_pip_susie_variant_pip_comparison.txt'
create_file_summarizing_tgfm_variant_pip_and_susie_variant_pip_comparison(eqtl_sample_size, simulation_runs, simulated_tgfm_results_dir, global_simulation_name_string, variant_comparison_file)
'''











































































'''
#############################################################
# Fine-mapping evaluation metrics varying gwas sample sizes
#############################################################
# Used eQTL sample sizes
eqtl_sample_sizes = np.asarray([500])
eqtl_sample_size = 500


# Simulation runs
# Currently hacky because had some failed simulations
simulation_runs = np.arange(1,101)
simulation_runs = np.arange(1,51)
#simulation_runs = np.delete(simulation_runs,11)
simulation_runs = np.delete(simulation_runs,[11])

gwas_sample_sizes = [50000, 100000, 200000]
ge_h2s = ['05', '075', '1']



# ln_pi methods used
ln_pi_methods = np.asarray(['uniform', 'pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])
#ln_pi_methods = np.asarray(['uniform'])


# twas method
twas_methods = np.asarray(['susie_pmces', 'susie_sampler'])


#Bim file
bim_file = processed_genotype_data_dir + 'simulated_gwas_data_' + chrom_num + '.bim'

# Windows
global_window_file = processed_genotype_data_dir + 'chromosome_' + str(chrom_num) + '_windows_3_mb.txt'


##################################
# Coverage/Calibration to detect snps with PIP > threshold WHILE VARYING GWAS SAMPLE SIZES
# Non-mediated variants also include expression-mediated variants
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
pip_thresholds = [.5, .9]

# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_vary_gwas_ss_nm_variant_alt_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp_non_mediated_variants_include_gene_variants_vary_gwas_ss(global_simulation_name_string, eqtl_sample_size, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_gene_expression_dir, gwas_sample_sizes)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_vary_gwas_ss_nm_variant_alt_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration_vary_gwas_ss(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_size, gwas_sample_sizes, ln_pi_methods, twas_methods)




##################################
# Coverage/Calibration to detect snps with PIP > threshold WHILE VARYING GE H2s
# Non-mediated variants also include expression-mediated variants
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
pip_thresholds = [.5, .9]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_vary_ge_h2s_nm_variant_alt_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp_non_mediated_variants_include_gene_variants_vary_ge_h2s(global_simulation_name_string, eqtl_sample_size, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_gene_expression_dir, ge_h2s)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_vary_ge_h2s_nm_variant_alt_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration_vary_ge_h2s(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_size, ge_h2s, ln_pi_methods, twas_methods)


##################################
# Power to detect snps with PIP > threshold WHILE VARYING GWAS SAMPLE SIZES
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
pip_thresholds = [.5, .9]

for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_vary_gwas_ss_power_per_component.txt'
	create_file_containing_tgfm_high_pip_snp_power_per_component_vary_gwas_ss(global_simulation_name_string, eqtl_sample_size, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, gwas_sample_sizes)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_vary_gwas_ss_power.txt'
	create_file_containing_averaged_tgfm_cs_power_vary_gwas_ss(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_size, ln_pi_methods,twas_methods, gwas_sample_sizes)

##################################
# Power to detect snps with PIP > threshold WHILE VARYING GE H2
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
pip_thresholds = [.5, .9]

for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_vary_ge_h2s_power_per_component.txt'
	create_file_containing_tgfm_high_pip_snp_power_per_component_vary_ge_h2s(global_simulation_name_string, eqtl_sample_size, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, ge_h2s)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_vary_ge_h2s_power.txt'
	create_file_containing_averaged_tgfm_cs_power_vary_ge_h2s(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_size, ln_pi_methods,twas_methods, ge_h2s)


##################################
# Coverage/Calibration to detect genes (not gene-tissue pairs) WHILE VARYING GWAS SS
##################################
# ln_pi methods used
ln_pi_methods = np.asarray(['pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])

# twas method
twas_methods = np.asarray(['susie_sampler'])
pip_thresholds = [.3, .5, .7, .9, .95, .99]
pip_thresholds = [.5, .9]

# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_vary_gwas_ss_gene_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_genes_vary_gwas_ss(global_simulation_name_string, eqtl_sample_size, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_input_data_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, gwas_sample_sizes)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_vary_gwas_ss_gene_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration_vary_gwas_ss(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_size, gwas_sample_sizes, ln_pi_methods, twas_methods)

##################################
# Coverage/Calibration to detect genes (not gene-tissue pairs) WHILE VARYING GE H2s
##################################
# ln_pi methods used
ln_pi_methods = np.asarray(['pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])

# twas method
twas_methods = np.asarray(['susie_sampler'])
pip_thresholds = [.3, .5, .7, .9, .95, .99]
pip_thresholds = [.5, .9]

# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_vary_ge_h2s_gene_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_genes_vary_ge_h2s(global_simulation_name_string, eqtl_sample_size, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_input_data_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, ge_h2s)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_vary_ge_h2s_gene_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration_vary_ge_h2s(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_size, ge_h2s, ln_pi_methods, twas_methods)

##################################
# Power to detect genes (not gene-tissue pairs) WHILE VARYING GWAS SS
##################################
# ln_pi methods used
ln_pi_methods = np.asarray(['pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])

# twas method
twas_methods = np.asarray(['susie_sampler'])
pip_thresholds = [.3, .5, .7, .9, .95, .99]
pip_thresholds = [.5, .9]

# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_vary_gwas_ss_gene_power_per_component.txt'
	create_file_containing_tgfm_high_pip_gene_power_per_component_vary_gwas_ss(global_simulation_name_string, eqtl_sample_size, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, gwas_sample_sizes)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_vary_gwas_ss_gene_power.txt'
	create_file_containing_averaged_tgfm_cs_power_vary_gwas_ss(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_size, ln_pi_methods,twas_methods, gwas_sample_sizes)


##################################
# Power to detect genes (not gene-tissue pairs) WHILE VARYING GE H2s
##################################
# ln_pi methods used
ln_pi_methods = np.asarray(['pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])

# twas method
twas_methods = np.asarray(['susie_sampler'])
pip_thresholds = [.3, .5, .7, .9, .95, .99]
pip_thresholds = [.5, .9]

# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_vary_ge_h2s_gene_power_per_component.txt'
	create_file_containing_tgfm_high_pip_gene_power_per_component_vary_ge_h2s(global_simulation_name_string, eqtl_sample_size, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, ge_h2s)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_vary_ge_h2s_gene_power.txt'
	create_file_containing_averaged_tgfm_cs_power_vary_ge_h2s(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_size, ln_pi_methods,twas_methods, ge_h2s)

'''


'''
#############################################################
# Fine-mapping evaluation metrics (for realistic eQTL sample size simulation settings)
#############################################################
# Simulation runs
# Currently hacky because had some failed simulations
tmp_simulation_runs = np.arange(101,200)
bads = {}
bads[134] = 1
simulation_runs = []
for sim_run in tmp_simulation_runs:
	if sim_run not in bads:
		simulation_runs.append(sim_run)
simulation_runs = np.asarray(simulation_runs)


# ln_pi methods used
ln_pi_methods = np.asarray(['uniform','pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])
# twas method
twas_methods = np.asarray(['susie_pmces', 'susie_sampler'])

#Bim file
bim_file = processed_genotype_data_dir + 'simulated_gwas_data_' + chrom_num + '.bim'
# Windows
global_window_file = processed_genotype_data_dir + 'chromosome_' + str(chrom_num) + '_windows_3_mb.txt'

eqtl_architecture = 'default'
gene_trait_architecture = '1_caus_t'
temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture
'''


##################################
# Coverage/Calibration to detect snps with PIP > threshold
##################################
'''
gene_type='component_gene'
pip_thresholds = [.5, .9]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_realistic_eqtl_ss_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp_realistic_qtl_sample_size(temp_global_simulation_name_string, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_learned_gene_models_dir,tiss_filter=True,gene_type=gene_type)
	print(cs_coverage_per_high_pip_snp_output_file)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_realistic_eqtl_ss_simulation_' + temp_global_simulation_name_string  + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, ['low_eqtl_ss', 'high_eqtl_ss'], ln_pi_methods, twas_methods, agg_eqtl_ss=True)
	print(cs_high_pip_coverage_output_file)
'''


##################################
# Power to detect snps with PIP > threshold
##################################
'''
pip_thresholds = [.5,.9]
gene_type='component_gene'
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_realistic_eqtl_ss_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_tgfm_high_pip_snp_power_per_component_realistic_qtl_sample_size(temp_global_simulation_name_string, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file,simulated_learned_gene_models_dir,tiss_filter=True, gene_type=gene_type)

	cs_power_output_file = simulated_organized_results_dir + 'organized_realistic_eqtl_ss_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, ['low_eqtl_ss', 'high_eqtl_ss'], ln_pi_methods,twas_methods, agg_eqtl_ss=True)
	print(cs_power_output_file)
'''
'''
##################################
# Make FDR-power curve
##################################
gene_types = ['component_gene']
eqtl_sss = ['low_eqtl_ss', 'high_eqtl_ss', 'aggregate']
# ln_pi methods used
ln_pi_method = 'pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'
# twas method
twas_method = 'susie_sampler'

for gene_type in gene_types:
	for eqtl_ss in eqtl_sss:
		fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_realistic_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_fdr_power_curve_raw_data.txt'
		fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_realistic_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_fdr_power_curve_data.txt'
		create_file_containing_tgfm_gene_tissue_fdr_power_curve_data_realistic_qtl_ss(temp_global_simulation_name_string, eqtl_ss, simulation_runs, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir, gene_type=gene_type)
# ln_pi
ln_pi_method = 'uniform'
# twas method
twas_method ='susie_pmces'
for gene_type in gene_types:
	for eqtl_ss in eqtl_sss:
		fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_realistic_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + 'fdr_power_curve_raw_data.txt'
		fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_realistic_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_fdr_power_curve_data.txt'
		create_file_containing_tgfm_gene_tissue_fdr_power_curve_data_realistic_qtl_ss(temp_global_simulation_name_string, eqtl_ss, simulation_runs, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir, gene_type=gene_type)

# ln_pi
ln_pi_method = 'uniform'
# twas method
twas_method ='susie_sampler'
for gene_type in gene_types:
	for eqtl_ss in eqtl_sss:
		fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_realistic_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + 'fdr_power_curve_raw_data.txt'
		fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_realistic_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_fdr_power_curve_data.txt'
		create_file_containing_tgfm_gene_tissue_fdr_power_curve_data_realistic_qtl_ss(temp_global_simulation_name_string, eqtl_ss, simulation_runs, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir, gene_type=gene_type)
'''




tmp_simulation_runs = np.arange(101,200)
bads = {}
bads[134] = 1
simulation_runs = []
for sim_run in tmp_simulation_runs:
	if sim_run not in bads:
		simulation_runs.append(sim_run)
simulation_runs = np.asarray(simulation_runs)

##################################
# Coverage/Calibration to detect snps with PIP > threshold using causal twas
##################################
'''
gene_type='component_gene'
pip_thresholds = [.5, .9]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_realistic_eqtl_ss_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_ctwas_tg_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_ctwas_tg_cs_calibration_per_high_pip_snp_realistic_qtl_sample_size(temp_global_simulation_name_string, simulation_runs, simulated_trait_dir, simulated_causal_twas_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_learned_gene_models_dir,tiss_filter=True,gene_type=gene_type)
	print(cs_coverage_per_high_pip_snp_output_file)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_realistic_eqtl_ss_simulation_' + temp_global_simulation_name_string  + '_' + gene_type + '_ctwas_tg_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, ['low_eqtl_ss', 'high_eqtl_ss'], ['ctwas'], ['lasso'], agg_eqtl_ss=True)
	print(cs_high_pip_coverage_output_file)
'''

'''
##################################
# Power to detect snps with PIP > threshold
##################################
pip_thresholds = [.5,.9]
gene_type='component_gene'
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_realistic_eqtl_ss_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_ctwas_tg_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_ctwas_tg_high_pip_snp_power_per_component_realistic_qtl_sample_size(temp_global_simulation_name_string, simulation_runs, simulated_trait_dir, simulated_causal_twas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, global_window_file,simulated_learned_gene_models_dir,tiss_filter=True, gene_type=gene_type)
	print(cs_power_per_component_output_file)

	cs_power_output_file = simulated_organized_results_dir + 'organized_realistic_eqtl_ss_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_ctwas_tg_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, ['low_eqtl_ss', 'high_eqtl_ss'], ['ctwas'],['lasso'], agg_eqtl_ss=True)
	print(cs_power_output_file)
'''
##################################
# Make FDR-power curve
##################################
'''
eqtl_sss = ['low_eqtl_ss', 'high_eqtl_ss', 'aggregate']
# ln_pi methods used
ln_pi_method = 'ctwas'
# twas method
twas_method = 'lasso'

for eqtl_ss in eqtl_sss:
	fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_realistic_simulation_' + temp_global_simulation_name_string + '_' + 'ctwas_tg' + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_fdr_power_curve_raw_data.txt'
	fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_realistic_simulation_' + temp_global_simulation_name_string + '_' + 'ctwas_tg' + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_fdr_power_curve_data.txt'
	create_file_containing_ctwas_tg_gene_tissue_fdr_power_curve_data_realistic_qtl_ss(temp_global_simulation_name_string, eqtl_ss, simulation_runs, simulated_trait_dir, simulated_causal_twas_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, simulated_learned_gene_models_dir)
'''



#############################################################
# Fine-mapping evaluation metrics (for alternative simulation settings)
#############################################################
# Used eQTL sample sizes
eqtl_sample_sizes = np.asarray([100,1000])
eqtl_sample_sizes = np.asarray([100])


# Simulation runs
# Currently hacky because had some failed simulations
tmp_simulation_runs = np.arange(1,101)
bads = {}
bads[60] = 1
bads[81] = 1
simulation_runs = []
for sim_run in tmp_simulation_runs:
	if sim_run not in bads:
		simulation_runs.append(sim_run)
simulation_runs = np.asarray(simulation_runs)




# ln_pi methods used

ln_pi_methods = np.asarray(['pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])

# twas method
twas_methods = np.asarray(['susie_sampler'])


#Bim file
bim_file = processed_genotype_data_dir + 'simulated_gwas_data_' + chrom_num + '.bim'

# Windows
global_window_file = processed_genotype_data_dir + 'chromosome_' + str(chrom_num) + '_windows_3_mb.txt'



eqtl_architecture = 'default'

gene_trait_architecture = '2_caus_t'

temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture

##################################
# Coverage/Calibration to detect snps with PIP > threshold
##################################
gene_type='component_gene'
'''
pip_thresholds = [.5, .9]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file,tiss_filter=True,gene_type=gene_type)
	print(cs_coverage_per_high_pip_snp_output_file)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string  + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods)
	print(cs_high_pip_coverage_output_file)
'''
##################################
# Coverage/Calibration to detect snps with PIP > threshold
##################################
gene_type='all_non_zero_gene'

pip_thresholds = [.5, .9]
#simulation_runs = np.arange(1,100)

# Vary ln_pi_method
'''
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file,tiss_filter=True,gene_type=gene_type)
	print(cs_coverage_per_high_pip_snp_output_file)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string  + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods)
	print(cs_high_pip_coverage_output_file)
'''



##################################
# Coverage/Calibration to detect snps with PIP > threshold
##################################
gene_type='max_min_ratio_5'
eqtl_sample_sizes = np.asarray([100])

simulation_runs = np.arange(101,200)

pip_thresholds = [.5, .9]
'''
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file,tiss_filter=True,gene_type=gene_type)
	print(cs_coverage_per_high_pip_snp_output_file)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string  + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods)
	print(cs_high_pip_coverage_output_file)
'''
gene_type='max_min_ratio_50'
eqtl_sample_sizes = np.asarray([100])

simulation_runs = np.arange(101,200)
'''
pip_thresholds = [.5, .9]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file,tiss_filter=True,gene_type=gene_type)
	print(cs_coverage_per_high_pip_snp_output_file)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string  + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods)
	print(cs_high_pip_coverage_output_file)
'''


##################################
# Coverage/Calibration to detect snps with PIP > threshold
##################################
'''
gene_type='component_gene'
pip_thresholds = [.3, .5, .7, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_calibration_only_correct_gene_per_component.txt'
	create_file_containing_tgfm_cs_calibration_only_correct_gene_per_high_pip_snp(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file,tiss_filter=True,gene_type=gene_type)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string  + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_calibration_only_correct_gene.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods)
	print(cs_high_pip_coverage_output_file)

gene_type='all_non_zero_gene'
pip_thresholds = [.3, .5, .7, .9, .95, .99]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_calibration_only_correct_gene_per_component.txt'
	create_file_containing_tgfm_cs_calibration_only_correct_gene_per_high_pip_snp(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file,tiss_filter=True,gene_type=gene_type)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string  + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_calibration_only_correct_gene.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods)
	print(cs_high_pip_coverage_output_file)
'''

##################################
# Power to detect snps with PIP > threshold
##################################

# Used eQTL sample sizes
eqtl_sample_sizes = np.asarray([100,1000])
eqtl_sample_sizes = np.asarray([100])


# Simulation runs
# Currently hacky because had some failed simulations
tmp_simulation_runs = np.arange(1,101)
bads = {}
bads[60] = 1
bads[81] = 1
simulation_runs = []
for sim_run in tmp_simulation_runs:
	if sim_run not in bads:
		simulation_runs.append(sim_run)
simulation_runs = np.asarray(simulation_runs)

'''
pip_thresholds = [.5,.9]
eqtl_architecture = 'default'
gene_trait_architecture = '2_caus_t'
temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture
gene_type='component_gene'
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_tgfm_high_pip_snp_power_per_component(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file,tiss_filter=True, gene_type=gene_type)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods,twas_methods)
	print(cs_power_output_file)

gene_type='all_non_zero_gene'
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_tgfm_high_pip_snp_power_per_component(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file,tiss_filter=True, gene_type=gene_type)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods,twas_methods)
	print(cs_power_output_file)
'''

'''
gene_type='max_min_ratio_5'
eqtl_sample_sizes = np.asarray([100])
simulation_runs = np.arange(101,200)
pip_thresholds = [.5, .9]
eqtl_architecture = 'default'
gene_trait_architecture = '2_caus_t'
temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_tgfm_high_pip_snp_power_per_component(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file,tiss_filter=True, gene_type=gene_type)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods,twas_methods)
	print(cs_power_output_file)



gene_type='max_min_ratio_50'
eqtl_sample_sizes = np.asarray([100])
simulation_runs = np.arange(101,200)
pip_thresholds = [.5, .9]
eqtl_architecture = 'default'
gene_trait_architecture = '2_caus_t'
temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_tgfm_high_pip_snp_power_per_component(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file,tiss_filter=True, gene_type=gene_type)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods,twas_methods)
	print(cs_power_output_file)
'''

'''
##################################
# Coverage/Calibration to detect genes (not gene-tissue pairs)
##################################
# ln_pi methods used
ln_pi_methods = np.asarray(['pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])
# twas method
twas_methods = np.asarray(['susie_sampler'])


eqtl_architecture = 'default'

gene_trait_architecture = '2_caus_t'

temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture
'''
'''
pip_thresholds = [.5,.9]
gene_type='component_gene'
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_gene_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_genes(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_input_data_dir, simulated_tgfm_results_dir, pip_threshold,gene_type, cs_coverage_per_high_pip_snp_output_file)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_gene_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods)
	print(cs_high_pip_coverage_output_file)

pip_thresholds = [.5,.9]
gene_type='all_non_zero_gene'
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_gene_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_genes(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_input_data_dir, simulated_tgfm_results_dir, pip_threshold,gene_type, cs_coverage_per_high_pip_snp_output_file)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_gene_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods)
	print(cs_high_pip_coverage_output_file)
'''


'''
##################################
# Power to detect genes (not gene-tissue pairs)
##################################
# ln_pi methods used
ln_pi_methods = np.asarray(['pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])

# twas method
twas_methods = np.asarray(['susie_sampler'])
pip_thresholds = [.5, .9]
gene_type='component_gene'
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string +'_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_gene_power_per_component.txt'
	#create_file_containing_tgfm_high_pip_gene_power_per_component(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, gene_type)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_gene_power.txt'
	#create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods,twas_methods)
	print(cs_power_output_file)

gene_type='all_non_zero_gene'
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string +'_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_gene_power_per_component.txt'
	create_file_containing_tgfm_high_pip_gene_power_per_component(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, gene_type)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_tgfm_pip_' + str(pip_threshold) + '_gene_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods,twas_methods)
	print(cs_power_output_file)
'''




##################################
# Make FDR-Power curve data
##################################
gene_trait_architecture = '2_caus_t'
temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture
# ln_pi
ln_pi_method = 'pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'
# twas method
twas_method ='susie_sampler'

tmp_simulation_runs = np.arange(1,101)
bads = {}
bads[60] = 1
bads[81] = 1
simulation_runs = []
for sim_run in tmp_simulation_runs:
	if sim_run not in bads:
		simulation_runs.append(sim_run)
simulation_runs = np.asarray(simulation_runs)
'''
gene_types = ['component_gene', 'all_non_zero_gene', 'max_min_ratio_5', 'max_min_ratio_50']
eqtl_sss = ['100']
for gene_type in gene_types:
	for eqtl_ss in eqtl_sss:
		fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_tgfm_gene_tissue_fdr_power_curve_raw_data.txt'
		fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_tgfm_gene_tissue_fdr_power_curve_data.txt'
		create_file_containing_tgfm_gene_tissue_fdr_power_curve_data(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, eqtl_ss, gene_type=gene_type)
'''
'''
# ln_pi
ln_pi_method = 'uniform'
# twas method
twas_method ='susie_pmces'
for gene_type in gene_types:
	for eqtl_ss in eqtl_sss:
		fdr_power_curve_raw_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_tgfm_gene_tissue_fdr_power_curve_raw_data.txt'
		fdr_power_curve_data_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + gene_type + '_' + twas_method + '_' + ln_pi_method + '_' + eqtl_ss + '_tgfm_gene_tissue_fdr_power_curve_data.txt'
		create_file_containing_tgfm_gene_tissue_fdr_power_curve_data(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, fdr_power_curve_raw_data_output_file,fdr_power_curve_data_output_file, twas_method, ln_pi_method, global_window_file, eqtl_ss, gene_type=gene_type)
'''

#############################################################
# Genome-wide analysis
#############################################################
#Bim file
bim_file = processed_genotype_data_dir + 'simulated_gwas_data_' + chrom_num + '.bim'

# Used eQTL sample sizes
eqtl_sample_sizes = np.asarray([100])

simulation_runs = np.arange(1,100)


'''
model_name = 'susie_pmces_uniform_iterative_variant_gene_prior_w_prior_pip_level'
gene_types = ['component_gene', 'all_non_zero_gene']
gene_types = [ 'max_min_ratio_5']

temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture
for gene_type in gene_types:
	##########################
	# Power and Type 1 Error based on gaussian approximation
	##########################
	# Create file showing p-value in causal tissues and non-causal tissues
	mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + model_name + '_' + gene_type + '_mediated_h2_gaussian_approximation_pvalue_by_tissue_est.txt'
	create_file_containing_iterative_bootstrapped_sampler_prior_mediated_h2_gaussian_approximation_pvalue_in_causal_and_non_causal_tissues(simulated_tgfm_results_dir, temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file, model_name, gene_type)
	# Create file showing type 1 error for null tissues using gaussian approximation
	type_1_error_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + model_name+ '_' + gene_type + '_h2_type_1_error_gaussian_approximation.txt'
	create_file_containing_mediated_h2_type_1_error_gaussian_approximation(mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file, type_1_error_output_file, eqtl_sample_sizes)
	# Create file showing power for null tissues using gaussian approximation
	power_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + model_name + '_' + gene_type+ '_mediated_h2_power_gaussian_approximation.txt'
	create_file_containing_mediated_h2_power_to_detect_causal_tissues_gaussian_approximation(mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file, power_output_file, eqtl_sample_sizes)



model_name = 'susie_pmces_uniform_iterative_variant_gene_prior_pip_level'
gene_types = ['component_gene', 'all_non_zero_gene']
gene_types = ['all_non_zero_gene','component_gene', 'max_min_ratio_5', 'max_min_ratio_50']

temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture
for gene_type in gene_types:
	##########################
	# Power and Type 1 Error based on gaussian approximation
	##########################
	# Create file showing p-value in causal tissues and non-causal tissues
	mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + model_name + '_' + gene_type + '_mediated_h2_gaussian_approximation_pvalue_by_tissue_est.txt'
	create_file_containing_iterative_bootstrapped_sampler_prior_mediated_h2_gaussian_approximation_pvalue_in_causal_and_non_causal_tissues(simulated_tgfm_results_dir, temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file, model_name, gene_type)
	# Create file showing type 1 error for null tissues using gaussian approximation
	type_1_error_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + model_name+ '_' + gene_type + '_h2_type_1_error_gaussian_approximation.txt'
	create_file_containing_mediated_h2_type_1_error_gaussian_approximation(mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file, type_1_error_output_file, eqtl_sample_sizes)
	# Create file showing power for null tissues using gaussian approximation
	power_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + model_name + '_' + gene_type+ '_mediated_h2_power_gaussian_approximation.txt'
	create_file_containing_mediated_h2_power_to_detect_causal_tissues_gaussian_approximation(mediated_iterative_sampler_gaussian_approximation_pvalue_by_tissue_output_file, power_output_file, eqtl_sample_sizes)
'''

'''
# Create file containing expected fraction of elememnts mediated by gene expresssion
# Model specification
twas_method = 'susie_sampler'
ln_pi_method = 'pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'
for gene_type in gene_types:
	expected_fraction_of_elements_mediated_by_gene_expression_summary_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + twas_method + '_' + ln_pi_method + '_' + gene_type + '_expected_fraction.txt'
	create_file_containing_expected_fraction_of_tgfm_elements_mediated_by_gene_expression(global_window_file, eqtl_sample_sizes, simulation_runs, twas_method, ln_pi_method, temp_global_simulation_name_string, simulated_tgfm_results_dir, gene_type, expected_fraction_of_elements_mediated_by_gene_expression_summary_file)
'''

'''
model_name = 'susie_pmces_uniform_iterative_variant_gene_prior_pip_level'

gene_types = ['all_non_zero_gene','component_gene', 'max_min_ratio_5', 'max_min_ratio_50']
for gene_type in gene_types:
	##########################
	# Average prior value
	##########################
	# Create file showing p-value in causal tissues and non-causal tissues
	average_prior_value_by_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + model_name + '_' + gene_type + '_average_prior_by_tissue.txt'
	create_file_containing_average_prior_value_in_causal_and_non_causal_tissues(simulated_tgfm_results_dir, temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, average_prior_value_by_tissue_output_file, model_name, gene_type)

	agg_average_prior_value_by_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + model_name + '_' + gene_type + '_agg_average_prior_by_tissue.txt'
	aggregate_average_prior_value_by_tissue(average_prior_value_by_tissue_output_file, agg_average_prior_value_by_tissue_output_file, eqtl_sample_sizes)
'''


'''
model_name = 'susie_pmces_uniform_iterative_variant_gene_prior_w_prior_pip_level'

gene_types = ['max_min_ratio_5']
for gene_type in gene_types:
	##########################
	# Average prior value
	##########################
	# Create file showing p-value in causal tissues and non-causal tissues
	average_prior_value_by_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + model_name + '_' + gene_type + '_average_prior_by_tissue.txt'
	create_file_containing_average_prior_value_in_causal_and_non_causal_tissues(simulated_tgfm_results_dir, temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, average_prior_value_by_tissue_output_file, model_name, gene_type)

	agg_average_prior_value_by_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_' + model_name + '_' + gene_type + '_agg_average_prior_by_tissue.txt'
	aggregate_average_prior_value_by_tissue(average_prior_value_by_tissue_output_file, agg_average_prior_value_by_tissue_output_file, eqtl_sample_sizes)
'''


# Used eQTL sample sizes
eqtl_sample_sizes = np.asarray([100])
model_name = 'susie_pmces_uniform_iterative_variant_gene_prior_pip_level'
temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture
organized_number_detected_genes_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_number_of_genes_with_gene_model.txt'
#create_file_containing_number_of_detected_genes(temp_global_simulation_name_string, np.arange(101,108),eqtl_sample_sizes, simulated_gene_expression_dir, simulated_learned_gene_models_dir, organized_number_detected_genes_output_file)




















'''
#############################################################
# Fine-mapping evaluation metrics (for alternative simulation settings)
#############################################################
# Used eQTL sample sizes
eqtl_sample_sizes = np.asarray([500])


# Simulation runs
# Currently hacky because had some failed simulations
simulation_runs = np.arange(1,51)



# ln_pi methods used

ln_pi_methods = np.asarray(['pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])

# twas method
twas_methods = np.asarray(['susie_sampler'])


#Bim file
bim_file = processed_genotype_data_dir + 'simulated_gwas_data_' + chrom_num + '.bim'

# Windows
global_window_file = processed_genotype_data_dir + 'chromosome_' + str(chrom_num) + '_windows_3_mb.txt'



eqtl_architecture = 'selection_1'

gene_trait_architecture = '2_caus_t'

temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture



##################################
# Power to detect snps with PIP > threshold
##################################
pip_thresholds = [.5,.9]

eqtl_architecture = 'selection_1'

gene_trait_architecture = '2_caus_t'

temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture

for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_all_genes_tgfm_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_tgfm_high_pip_snp_power_per_component(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file,tiss_filter=True, all_genes_bool=True)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_all_genes_tgfm_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods,twas_methods)

##################################
# Coverage/Calibration to detect snps with PIP > threshold
##################################
pip_thresholds = [.5,.7,.9, .95]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_all_genes' + '_tgfm_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file,tiss_filter=True,all_genes=True)
	print(cs_coverage_per_high_pip_snp_output_file)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string  + '_all_genes'+ '_tgfm_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods)
	print(cs_high_pip_coverage_output_file)


eqtl_sample_sizes = np.asarray([300,500])
temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture

##################################
# Coverage/Calibration to detect snps with PIP > threshold
##################################
pip_thresholds = [.5,.7,.9, .95]

# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file,tiss_filter=True)
	print(cs_coverage_per_high_pip_snp_output_file)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods)
	print(cs_high_pip_coverage_output_file)

'''
'''
eqtl_architecture = 'selection_125'
gene_trait_architecture = '2_caus_t'

temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture

##################################
# Coverage/Calibration to detect snps with PIP > threshold
##################################
pip_thresholds = [.5,.9]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file,tiss_filter=True)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods)
	print(cs_high_pip_coverage_output_file)


eqtl_architecture = 'random_Neqtl'
gene_trait_architecture = '2_caus_t'

temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture

##################################
# Coverage/Calibration to detect snps with PIP > threshold
##################################
pip_thresholds = [.5,.9]
# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file,tiss_filter=True)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods)
	print(cs_high_pip_coverage_output_file)
'''



##################################
# Power to detect snps with PIP > threshold
##################################
'''
pip_thresholds = [.5,.9]

eqtl_architecture = 'selection_1'

gene_trait_architecture = '2_caus_t'

temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture

for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_tgfm_high_pip_snp_power_per_component(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file,tiss_filter=True)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods,twas_methods)

eqtl_architecture = 'selection_125'

gene_trait_architecture = '2_caus_t'

temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture

for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_tgfm_high_pip_snp_power_per_component(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file,tiss_filter=True)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods,twas_methods)


eqtl_architecture = 'random_Neqtl'

gene_trait_architecture = '2_caus_t'

temp_global_simulation_name_string = global_simulation_name_string + '_gt_arch_' + gene_trait_architecture + '_qtl_arch_' + eqtl_architecture

for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_power_per_component.txt'
	create_file_containing_tgfm_high_pip_snp_power_per_component(temp_global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file,tiss_filter=True)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + temp_global_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_power.txt'
	create_file_containing_averaged_tgfm_cs_power(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods,twas_methods)






'''



























































'''
#############################################################
# Fine-mapping evaluation metrics for missing causal tissue
#############################################################
# Used eQTL sample sizes
eqtl_sample_sizes = np.asarray([300,500,1000])


# Simulation runs
simulation_runs = np.arange(1,31)



# ln_pi methods used

ln_pi_methods = np.asarray(['pmces_uniform_iterative_variant_gene_prior_pip_level_bootstrapped'])


# twas method
twas_methods = np.asarray(['susie_sampler'])

# Missingness version
missingness_methods = np.asarray(['all_t', 'no_t0'])


#Bim file
bim_file = processed_genotype_data_dir + 'simulated_gwas_data_' + chrom_num + '.bim'

# Windows
global_window_file = processed_genotype_data_dir + 'chromosome_' + str(chrom_num) + '_windows_3_mb.txt'

local_simulation_name_string = global_simulation_name_string + '_gt_arch_' + '1_caus_t'

##################################
# Coverage/Calibration to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]

# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_missing_causal_tissue_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp_in_missing_causal_tissue_sim(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, missingness_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_best_tagging_gt_dir)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_missing_causal_tissue_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration_in_missing_causal_tissue_sim(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods, missingness_methods)

##################################
# Investigate cases where TGFM doesn't get best tagging gene-tissue pair
##################################
pip_thresholds = [.5]


for pip_threshold in pip_thresholds:
	not_best_tagging_gene_tissue_pairs_examples_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_not_best_tagging_gene_tissue_pairs_when_missing_causal_tissue.txt'
	create_file_containing_examples_of_tgfm_prioritized_gene_tissue_from_not_best_tagging_tissue(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, missingness_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, not_best_tagging_gene_tissue_pairs_examples_output_file, simulated_best_tagging_gt_dir)
	print(not_best_tagging_gene_tissue_pairs_examples_output_file)


##################################
# Coverage/Calibration to detect snps with PIP > threshold varying by number causal genetic elements in the region
##################################
pip_thresholds = [.5, .9]
n_causal_genetic_element_bins = ['0_5', '6_10', '11_20']

# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_missing_causal_tissue_n_causal_genetic_element_stratefied_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp_in_missing_causal_tissue_sim_stratefied_by_n_causal_genetic_elements_bins(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, missingness_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_best_tagging_gt_dir, bim_file, simulated_gene_position_file,)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_missing_causal_tissue_n_causal_genetic_element_stratefied_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration_in_missing_causal_tissue_sim_stratefied_by_n_causal_genetic_elements_bins(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods, missingness_methods,n_causal_genetic_element_bins)


##################################
# Power to detect snps with PIP > threshold
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]
pip_thresholds = [.5, .9]

for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per causal element (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, gene_or_variant, causal_element_name, boolean_causal_element_in_cs)
	cs_power_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_missing_causal_tissue_power_per_component.txt'
	create_file_containing_tgfm_high_pip_snp_power_per_component_in_missing_causal_tissue_sim(local_simulation_name_string, eqtl_sample_sizes, missingness_methods, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, cs_power_per_component_output_file, pip_threshold, twas_methods, global_window_file, simulated_best_tagging_gt_dir)

	cs_power_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_missing_causal_tissue_power.txt'
	create_file_containing_averaged_tgfm_cs_power_in_missing_causal_tissue_sim(cs_power_per_component_output_file, cs_power_output_file, eqtl_sample_sizes, ln_pi_methods,twas_methods, missingness_methods)
	print(cs_power_output_file)

##################################
# Coverage/Calibration to detect snps with PIP > threshold at gene level
##################################
pip_thresholds = [.3, .5, .7, .9, .95, .99]

# Vary ln_pi_method
for pip_threshold in pip_thresholds:
	print(pip_threshold)
	# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
	cs_coverage_per_high_pip_snp_output_file =  simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_missing_causal_tissue_gene_calibration_per_component.txt'
	create_file_containing_tgfm_cs_calibration_per_high_pip_snp_in_missing_causal_tissue_sim(local_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, twas_methods, missingness_methods, simulated_trait_dir, simulated_tgfm_results_dir, pip_threshold, cs_coverage_per_high_pip_snp_output_file, simulated_best_tagging_gt_dir, gene_level=True)

	cs_high_pip_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_missing_causal_tissue_gene_calibration.txt'
	create_file_containing_averaged_tgfm_high_pip_calibration_in_missing_causal_tissue_sim(cs_coverage_per_high_pip_snp_output_file, cs_high_pip_coverage_output_file, eqtl_sample_sizes, ln_pi_methods, twas_methods, missingness_methods)


##################################
# Fraction of causal gene-tissue pairs discovered by each class of genetic element (causal gene-tissue, best tagging gene-tissue, other gene-tissue, nm variant)
##################################
pip_thresholds = [.5, .9]
for pip_threshold in pip_thresholds:
	causal_gene_tissue_pair_stratification_per_component_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_missing_causal_tissue_causal_gene_tissue_pair_stratification_per_component.txt'
	create_file_containing_causal_gene_tissue_pair_stratification_in_missing_causal_tissue_sim(local_simulation_name_string, eqtl_sample_sizes, missingness_methods, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, simulated_tgfm_input_data_dir, bim_file, simulated_gene_position_file, causal_gene_tissue_pair_stratification_per_component_output_file, pip_threshold, twas_methods, global_window_file, simulated_best_tagging_gt_dir, simulated_gene_expression_dir)

	causal_gene_tissue_pair_stratification_output_file = simulated_organized_results_dir + 'organized_simulation_' + local_simulation_name_string + '_tgfm_pip_' + str(pip_threshold) + '_missing_causal_tissue_causal_gene_tissue_pair_stratification.txt'
	create_file_containing_averaged_causal_gene_tissue_pair_stratification_in_missing_causal_tissue_sim(causal_gene_tissue_pair_stratification_per_component_output_file, causal_gene_tissue_pair_stratification_output_file, eqtl_sample_sizes, missingness_methods)
	print(causal_gene_tissue_pair_stratification_output_file)
'''



























































