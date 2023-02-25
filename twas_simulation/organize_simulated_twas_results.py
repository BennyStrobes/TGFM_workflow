import numpy as np 
import os
import sys
import pdb

def extract_valid_simulation_numbers(simulated_twas_dir, simulation_name_string):
	valid_sims = []
	for sim_iter in range(200):
		sim_file = simulated_twas_dir + 'simulation_' + str(sim_iter) + '_' + simulation_name_string + '_simualated_twas_results.txt'
		if os.path.isfile(sim_file):
			valid_sims.append(sim_iter)


	return np.asarray(valid_sims)


def create_file_summarizing_effect_of_modeling_gene_distribution(eqtl_sample_size, model_type, sim_numbers, simulated_twas_dir, output_file):
	pmces_dicti = {}
	distr_dicti = {}
	ratio_dicti = {}
	used = []
	for sim_number in sim_numbers:
		sim_file = simulated_twas_dir + 'simulation_' + str(sim_number) + '_' + simulation_name_string + '_simualated_twas_results.txt'
		f = open(sim_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			if data[1] != eqtl_sample_size:
				continue
			gene_name = data[0]
			twas_model = data[2]
			if data[5] == 'nan':
				continue
			if twas_model == model_type + '_pmces':
				if str(sim_number) + '_' + gene_name in pmces_dicti:
					print('assumption eroror')
					pdb.set_trace()
				pmces_dicti[str(sim_number) + '_' + gene_name] = data[7]
				used.append(str(sim_number) + '_' + gene_name)
			elif twas_model == model_type + '_distr':
				if str(sim_number) + '_' + gene_name in distr_dicti:
					print('assumption eroror')
					pdb.set_trace()
				distr_dicti[str(sim_number) + '_' + gene_name] = data[7]
				ratio_dicti[str(sim_number) + '_' + gene_name] = data[8]
				used.append(str(sim_number) + '_' + gene_name)

		f.close()
	used = np.unique(used)

	t = open(output_file,'w')
	t.write('pmces_p_value\tdistribution_p_value\tvariance_ratio\n')
	for ele in used:
		t.write(pmces_dicti[ele] + '\t' + distr_dicti[ele] + '\t' + ratio_dicti[ele] + '\n')
	t.close()
	return


def run_fdr_analysis(p_value_threshold, parameters, sim_numbers, simulated_twas_dir, simulation_name_string, output_file):
	t = open(output_file,'w')
	t.write('twas_model\teqtl_sample_size\tfdr\tfdr_lb\tfdr_ub\n')

	for parameter_tuple in parameters:
		twas_model = parameter_tuple[0]
		eqtl_ss = parameter_tuple[1]
		fdr_numerator = 0.0
		fdr_denominator = 0.0
		for sim_number in sim_numbers:
			sim_file = simulated_twas_dir + 'simulation_' + str(sim_number) + '_' + simulation_name_string + '_simualated_twas_results.txt'
			f = open(sim_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				if data[1] != eqtl_ss:
					continue
				if data[2] != twas_model:
					continue
				pvalue = float(data[-2])
				if np.isnan(pvalue):
					continue
				if pvalue > p_value_threshold:
					continue
				fdr_denominator = fdr_denominator + 1
				if data[3] == 'False':
					fdr_numerator = fdr_numerator + 1
			f.close()
		fdr = fdr_numerator/fdr_denominator
		se = np.sqrt((fdr)*(1.0-fdr))/np.sqrt(fdr_denominator)
		t.write(twas_model + '\t' + str(eqtl_ss) + '\t' + str(fdr) + '\t' + str(fdr-(se*1.96)) + '\t' + str(fdr+(se*1.96)) + '\n')	
	t.close()
	return



def run_power_analysis(p_value_threshold, parameters, sim_numbers, simulated_twas_dir, simulation_name_string, output_file):
	t = open(output_file,'w')
	t.write('twas_model\teqtl_sample_size\tpower\tpower_lb\tpower_ub\n')

	for parameter_tuple in parameters:
		twas_model = parameter_tuple[0]
		eqtl_ss = parameter_tuple[1]
		power_numerator = 0.0
		power_denominator = 0.0
		for sim_number in sim_numbers:
			sim_file = simulated_twas_dir + 'simulation_' + str(sim_number) + '_' + simulation_name_string + '_simualated_twas_results.txt'
			f = open(sim_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				if data[1] != eqtl_ss:
					continue
				if data[2] != twas_model:
					continue
				if data[3] != 'True':
					continue
				pvalue = float(data[-2])
				power_denominator = power_denominator + 1
				if np.isnan(pvalue) == False and pvalue < p_value_threshold:
					power_numerator = power_numerator + 1
			f.close()
		power = power_numerator/power_denominator
		se = np.sqrt((power)*(1.0-power))/np.sqrt(power_denominator)
		t.write(twas_model + '\t' + str(eqtl_ss) + '\t' + str(power) + '\t' + str(power-(se*1.96)) + '\t' + str(power+(se*1.96)) + '\n')	
	t.close()
	return

def make_power_false_discovery_curve_input_data(parameters, sim_numbers, simulated_twas_dir, simulation_name_string, output_file):
	t = open(output_file,'w')
	t.write('twas_model\teqtl_sample_size\tp_value_thresh\tpower\tfdr\n')

	for parameter_tuple in parameters:
		twas_model = parameter_tuple[0]
		eqtl_ss = parameter_tuple[1]
		pvalues = []
		labels = []
		for sim_number in sim_numbers:
			sim_file = simulated_twas_dir + 'simulation_' + str(sim_number) + '_' + simulation_name_string + '_simualated_twas_results.txt'
			f = open(sim_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				if data[1] != eqtl_ss:
					continue
				if data[2] != twas_model:
					continue
				pvalue = float(data[-2])
				if np.isnan(pvalue):
					pvalue = 1.0
				label = data[3]
				pvalues.append(pvalue)
				labels.append(label)
			f.close()
		pvalues = np.asarray(pvalues)
		labels = np.asarray(labels)
		positive_labels = np.where(labels=='True')[0]

		thresholds = [.99, .98, .97, .96, .95, .94, .93, .92, .91, .9, .95, .8, .7, .6, .5, .4, .3, .2, .1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-20, 1e-30, 1e-40, 1e-50]
		for threshold in thresholds:
			# Compute power
			positive_pvalues = pvalues[positive_labels]
			power = np.sum(positive_pvalues <= threshold)/len(positive_pvalues)
			# Compute fdr
			identified_positives = np.where(pvalues <= threshold)[0]
			identified_positive_labels = labels[identified_positives]
			fdr = np.sum(identified_positive_labels == 'False')/len(identified_positive_labels)
			# Print
			t.write(twas_model + '\t' + eqtl_ss + '\t' + str(threshold) + '\t' + str(power) + '\t' + str(fdr) + '\n')
	t.close()
	return



simulated_twas_dir = sys.argv[1]
simulated_organized_results_dir = sys.argv[2]
simulation_name_string = sys.argv[3]



sim_numbers = extract_valid_simulation_numbers(simulated_twas_dir, simulation_name_string)

eqtl_sample_sizes = ['100', '200', '300', '500', '1000']
model_types = ['susie_pmces', 'susie_distr', 'marginal_pmces', 'marginal_distr', 'fusion_lasso_pmces']
parameters = []
for eqtl_sample_size in eqtl_sample_sizes:
	for model_type in model_types:
		parameters.append((model_type, eqtl_sample_size))
parameters.append(('true_causal_effects', 'inf'))


#############################################
# Create file comparing pmces with distribution
#############################################
for eqtl_sample_size in eqtl_sample_sizes:
	for model_type in model_types:
		output_file = simulated_organized_results_dir + simulation_name_string + '_organized_bf_comparison_' + model_type + '_' + eqtl_sample_size + '.txt'
		create_file_summarizing_effect_of_modeling_gene_distribution(eqtl_sample_size, model_type, sim_numbers, simulated_twas_dir, output_file)
#############################################
# Power analysis
#############################################
p_value_thresholds = [.05, .01]
for p_value_threshold in p_value_thresholds:
	output_file = simulated_organized_results_dir + simulation_name_string + '_statistical_power_p_' + str(p_value_threshold) + '.txt'
	run_power_analysis(p_value_threshold, parameters, sim_numbers, simulated_twas_dir, simulation_name_string, output_file)
#############################################
# FDR analysis
#############################################
p_value_thresholds = [.05, .01]
for p_value_threshold in p_value_thresholds:
	output_file = simulated_organized_results_dir + simulation_name_string + '_statistical_fdr_p_' + str(p_value_threshold) + '.txt'
	run_fdr_analysis(p_value_threshold, parameters, sim_numbers, simulated_twas_dir, simulation_name_string, output_file)

#############################################
# Power-False discovery curve
#############################################
output_file = simulated_organized_results_dir + simulation_name_string + '_power_false_discovery_curve_input.txt'
make_power_false_discovery_curve_input_data(parameters, sim_numbers, simulated_twas_dir, simulation_name_string, output_file)


