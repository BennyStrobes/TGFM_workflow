import numpy as np 
import os
import sys
import pdb
import scipy.stats

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
			h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_sldsc_results_organized_mediated_h2_5_50.txt'
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
			frac_h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_sldsc_results_h2_5_50_med.txt'
			frac_h2_raw = np.loadtxt(frac_h2_file,dtype=str)
			frac_h2 = frac_h2_raw[1,0]

			arr.append(float(frac_h2))

			# print to output file
			t.write(str(eqtl_sample_size) + '\t' + str(simulation_run) + '\t' + str(frac_h2) + '\t' + str(fraction_expression_mediated_heritability) + '\n')
		print(np.mean(arr))
	# Close file handle
	t.close()
	return

def create_file_containing_mediated_h2_pvalue_in_causal_and_non_causal_tissues(simulated_sldsc_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_pvalue_by_tissue_output_file):
	# Open output file and print header
	t = open(mediated_pvalue_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\ttissue_number\tz_score\tpvalue\tcausal_status\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_sldsc_results_organized_res.txt'
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

def create_file_containing_mediated_h2_in_causal_and_non_causal_tissues_sparse_estimate(simulated_sldsc_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_h2_by_tissue_output_file):
	# Open output file and print header
	t = open(mediated_h2_by_tissue_output_file,'w')
	t.write('eqtl_sample_size\tsimulation_number\ttissue_number\testimated_med_h2\tsimulated_mediated_h2\n')
	for eqtl_sample_size in eqtl_sample_sizes:
		arr = []
		for simulation_run in simulation_runs:
			# Extract per annotation mediated heritability file for this run
			h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_sldsc_results_organized_0.5_sparse_ard_eqtl_coefficients_mv_update_res.txt'
			# Load in per anno med h2
			h2_data_raw = np.loadtxt(h2_file, dtype=str, delimiter='\t')
			med_tau = h2_data_raw[-10:,1].astype(float)

			gene_m_file = '/n/scratch3/users/b/bes710/causal_eqtl_gwas/simulation/simulated_ld_scores/' + 'simulation_9_chrom21_cis_window_100000_gene_weighted_ld_scores_eqtlss_' + str(eqtl_sample_size) + '_M.txt'
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
			h2_file = simulated_sldsc_results_dir + 'simulation_' + str(simulation_run) + '_' + global_simulation_name_string + '_eqtl_ss_' + str(eqtl_sample_size) + '_sldsc_results_organized_mediated_h2_5_50.txt'
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
		valid_sim = False 
		for file_name in os.listdir(simulated_tgfm_results_dir):

			if file_name.startswith('simulation_' + str(sim_iter) + '_' + global_simulation_name_string + '_eqtl_ss_500_ln_pi_distribution_estimate_1e-10'):
				if file_name.endswith('_results.pkl'):
					valid_sim = True
		if valid_sim:
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

# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
def create_file_containing_tgfm_cs_calibration_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, cs_coverage_per_component_output_file):
	# Open output file handle
	t = open(cs_coverage_per_component_output_file,'w')
	t.write('eQTL_sample_size\tln_pi_method\tsimulation_number\twindow_name\tcomponent_number\tcausal_genetic_element_in_cs\n')

	# First loop through simulations
	for simulation_number in simulation_runs:
		# First extract dictionary list of causal genetic elements
		causal_variant_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
		causal_gene_file = simulated_trait_dir + 'simulation_' + str(simulation_number) + '_' + global_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
		causal_genetic_elements, causal_variants, causal_genes = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)

		# Now loop through eqtl sample sizes and ln_pi methods
		for eqtl_sample_size in eqtl_sample_sizes:
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
					cs_genetic_elements = data[3].split(';')
					cs_probs = np.asarray(data[4].split(';')).astype(float)
					# Quick error check
					if np.sum(cs_probs) < .95:
						print('assumption error')
						pdb.set_trace()
					# Check if cs contains at least one causal genetic element
					causal_genetic_element_in_cs_boolean = check_if_at_least_one_causal_genetic_element_is_in_cs(cs_genetic_elements, causal_genetic_elements)
					t.write(str(eqtl_sample_size) + '\t' + ln_pi_method + '\t' + str(simulation_number) + '\t' + component_window_name + '\t' + component_num + '\t' + str(int(causal_genetic_element_in_cs_boolean)) + '\n')
				f.close()
	t.close()
	return

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


# Used eQTL sample sizes
eqtl_sample_sizes = np.asarray([100,200,300,500,1000])

# ln_pi methods used
ln_pi_methods = np.asarray(['point_estimate_1e-10', 'sparse_estimate_1e-10', 'distribution_estimate_1e-10', 'point_estimate_1e-30', 'sparse_estimate_1e-30', 'distribution_estimate_1e-30'])

# Simulation runs
# Currently hacky because had some failed simulations
simulation_runs = extract_successfully_completed_simulation_runs(simulated_tgfm_results_dir, global_simulation_name_string)


'''
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
# Power and type 1 error in TGFM-SLDSC
##############################
# Create file showing p-value in causal tissues and non-causal tissues
mediated_pvalue_by_tissue_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_mediated_h2_pvalue_by_tissue_est.txt'
create_file_containing_mediated_h2_pvalue_in_causal_and_non_causal_tissues(simulated_sldsc_results_dir, global_simulation_name_string, eqtl_sample_sizes, simulation_runs, mediated_pvalue_by_tissue_output_file)

# Create file showing power to detect causal tissues
power_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_mediated_h2_power.txt'
create_file_containing_mediated_h2_power_to_detect_causal_tissues(mediated_pvalue_by_tissue_output_file, power_output_file, eqtl_sample_sizes)

# Create file showing type 1 error for null tissues
type_1_error_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_mediated_h2_type_1_error.txt'
create_file_containing_mediated_h2_type_1_error(mediated_pvalue_by_tissue_output_file, type_1_error_output_file, eqtl_sample_sizes)
'''


##################################
# Fine-mapping evaluation metrics
##################################
# (i) Coverage/Calibration: the proportion of credible sets that include at least one true causal genetic element across simulation replicates;
# (ii) Power: the number of true causal variants identified (i.e., covered by a credible set)
# (iii) Resolution: the size of credible sets and the number of fine-mapped variants with high confidence (e.g., PIP >95%);
# These should be done for both genes and variants


##################################
# Coverage/Calibration
##################################
# Create file with one line per cs (columns: eQTL_sample_size, ln_pvalue_method, simulation_num, window_num, component_num, boolean_causal_variant_in_cs)
cs_coverage_per_component_output_file =  simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_cs_calibration_per_component.txt'
create_file_containing_tgfm_cs_calibration_per_component(global_simulation_name_string, eqtl_sample_sizes, simulation_runs, ln_pi_methods, simulated_trait_dir, simulated_tgfm_results_dir, cs_coverage_per_component_output_file)

cs_coverage_output_file = simulated_organized_results_dir + 'organized_simulation_' + global_simulation_name_string + '_tgfm_cs_calibration.txt'
create_file_containing_averaged_tgfm_cs_calibration(cs_coverage_per_component_output_file, cs_coverage_output_file, eqtl_sample_sizes, ln_pi_methods)







