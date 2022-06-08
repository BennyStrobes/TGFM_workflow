import numpy as np 
import os
import sys
import pdb





def extract_study_names_from_window_file(window_file):
	data = np.loadtxt(window_file,dtype=str,delimiter='\t')
	arbitrary_study_file = data[1,-3]  # Could be any row
	study_names = np.loadtxt(arbitrary_study_file,dtype=str,delimiter='\t')
	return study_names

def extract_window_info(window_id):
	window_info = window_id.split(':')
	chrom_num = window_info[0]
	start = int(window_info[1])
	end = int(window_info[2])
	return chrom_num, start, end

def windows_overlap(window_id_1, window_id_2):
	chrom_num1, start1, end1 = extract_window_info(window_id_1)
	chrom_num2, start2, end2 = extract_window_info(window_id_2)
	if window_id_1 == window_id_2:
		return 'identical'
	if chrom_num1 != chrom_num2:
		return 'no_overlap'
	# Windows are not identical, but are on the same chromosome
	if start2 >= start1 and start2 <= end1:
		return 'overlap'
	elif end2 >= start1 and end2 <= end1:
		return 'overlap'
	else:
		return 'no_overlap'

def get_ordered_intersection_variants(variants, pmcefs, intersection_variants_dicti, num_variants):
	new_arr = np.zeros(num_variants) - 50000.0
	for var_index, var_name in enumerate(variants):
		pmcef = pmcefs[var_index]
		if var_name in intersection_variants_dicti:
			new_arr[intersection_variants_dicti[var_name]] = pmcef

	if sum(new_arr == -50000.0)  > 0:
		print('assumption erroro')
		pdb.set_trace()
	return new_arr


def get_component_genotype(window_variant_ids, window_genotype, component_variant_id_dicti):
	indices = []
	for index, val in enumerate(window_variant_ids):
		if val in component_variant_id_dicti:
			indices.append(index)
	indices = np.asarray(indices)
	if len(indices) != len(component_variant_id_dicti):
		print('assumptoinoeieroerieo)')
		pdb.set_trace()

	return window_genotype[:, indices]

def get_mean_abs_corr_between_two_genotype_matrices(geno1, geno2):
	col1 = geno1.shape[1]
	col2 = geno2.shape[1]

	corry = []
	for ii in range(col1):
		for jj in range(col2):
			corrz = np.corrcoef(geno1[:,ii], geno2[:,jj])[0,1]
			corry.append(corrz)
	return np.mean(np.abs(corry))


def filter_component_file(study_name, study_component_file, window_to_genotype_file, study_filtered_component_file):
	f = open(study_component_file)
	t = open(study_filtered_component_file,'w')

	prior_windows = {}
	head_count= 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		counter = counter + 1
		# Extract relevent fields
		lead_variant = data[2]
		window_id = data[4]
		component_id = lead_variant + '_' + window_id

		window_genotype_file = window_to_genotype_file[window_id][0]
		window_variant_id_file = window_to_genotype_file[window_id][1]
		window_variant_ids = np.loadtxt(window_variant_id_file,dtype=str,delimiter='\t')
		window_genotype = np.loadtxt(window_genotype_file, delimiter='\t')

		component_variant_ids = data[5].split(',')
		component_variant_id_dicti = {}
		for component_variant_id in component_variant_ids:
			component_variant_id_dicti[component_variant_id] = 1

		component_genotype = get_component_genotype(window_variant_ids, window_genotype, component_variant_id_dicti)

		# Loop through old components
		prev_window_ids = [*prior_windows]

		pass_filter = True

		for prev_window_id in prev_window_ids:
			# Check if current window id overlaps prev window id
			if windows_overlap(window_id, prev_window_id) == 'overlap':
				#print('windows overlap')
				prev_window_component_arr = prior_windows[prev_window_id]  # This is a vector of length number of components in the window. each element is a genotype matrix corresponding to that componetn
				for prev_window_component_genotype in prev_window_component_arr:
					mean_abs_corr = get_mean_abs_corr_between_two_genotype_matrices(prev_window_component_genotype, component_genotype)
					if mean_abs_corr > .5:
						pass_filter = False
			elif windows_overlap(window_id, prev_window_id) == 'no_overlap':
				#print('windows dont overlap')
				prior_windows.pop(prev_window_id, None)

		if window_id not in prior_windows:
			prior_windows[window_id] = []
		prior_windows[window_id].append(component_genotype)

		if pass_filter == True:
			t.write(line + '\n')
		else:
			print('throwing out ' + component_id)

	f.close()
	t.close()


def extract_study_component_correlations(study_name, study_component_file, study_component_correlation_file, ukbb_genome_wide_susie_organized_results_dir):
	f = open(study_component_file)
	prior_windows = {}
	head_count= 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent fields
		lead_variant = data[2]
		window_id = data[4]
		component_id = lead_variant + '_' + window_id

		pmcef_file = ukbb_genome_wide_susie_organized_results_dir + study_name + '_' + component_id + '_component_posterior_mean_causal_effect_sizes.txt'
		pmcef_data = np.loadtxt(pmcef_file,dtype=str,delimiter='\t')
		component_variant_ids = pmcef_data[1:,0]
		component_pmcef = pmcef_data[1:,1].astype(float)

		# Loop through old components
		for prev_window_id, prev_window_arr in prior_windows.items():
			# Check if current window id overlaps prev window id
			if windows_overlap(window_id, prev_window_id) == 'overlap':
				print('windows overlap')
				prev_window_component_arr = prior_windows[prev_window_id]
				for prev_window_component in prev_window_component_arr:
					prev_window_variants = prev_window_component[0]
					prev_window_pmcef = prev_window_component[1]
					intersection_variants = list(set(prev_window_variants).intersection(component_variant_ids))
					if len(intersection_variants) < 50:
						continue
					intersection_variants_dicti = {}
					for i, val in enumerate(intersection_variants):
						intersection_variants_dicti[val] = i 

					prior_intersection_variants_pmcef = get_ordered_intersection_variants(prev_window_variants, prev_window_pmcef, intersection_variants_dicti, len(intersection_variants))
					curr_intersection_variants_pmcef = get_ordered_intersection_variants(component_variant_ids, component_pmcef, intersection_variants_dicti, len(intersection_variants))
					pdb.set_trace()

			elif windows_overlap(window_id, prev_window_id) == 'no_overlap':
				print('windows dont overlap')
				pdb.set_trace()
			elif windows_overlap(window_id, prev_window_id) == 'identical':
				print('identical windows')
				pdb.set_trace()



		if window_id not in prior_windows:
			prior_windows[window_id] = []
			prior_windows[window_id].append((component_variant_ids, component_pmcef))

	f.close()


def create_mapping_from_window_to_genotype_file(window_file):
	f = open(window_file)
	head_count = 0
	mapping = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		genotype_file = data[-2]
		variant_id_file = data[-4]
		window_id = ':'.join(data[:3])
		if window_id in mapping:
			print('assumption oerororo')
			pdb.set_trace()
		mapping[window_id] = (genotype_file, variant_id_file)
	f.close()

	return mapping



window_file = sys.argv[1]  # only use this file to get study names
ukbb_genome_wide_susie_organized_results_dir = sys.argv[2]



study_names = extract_study_names_from_window_file(window_file)


window_to_genotype_file = create_mapping_from_window_to_genotype_file(window_file)



global_output_file = ukbb_genome_wide_susie_organized_results_dir + 'organized_susie_components_files.txt'
#t = open(global_output_file,'w')
#t.write('trait\ttrait_component_file\n')


# Perform analysis for each study seperately
for study_name in study_names:
	# File containing all components for this study
	study_component_file = ukbb_genome_wide_susie_organized_results_dir + study_name + '_organized_susie_components.txt'  # Input file
	study_filtered_component_file = ukbb_genome_wide_susie_organized_results_dir + study_name + '_organized_susie_filtered_components.txt'  # Output file
	filter_component_file(study_name, study_component_file, window_to_genotype_file, study_filtered_component_file)
	#t.write(study_name + '\t' + study_component_file + '\n')
	#extract_study_component_max_correlation_among_ci

#t.close()



