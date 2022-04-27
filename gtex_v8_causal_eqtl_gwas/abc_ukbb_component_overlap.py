import numpy as np 
import os
import sys
import pdb
import gzip


def get_epimap_cell_types_and_file_names(epimap_input_dir):
	cell_types = []
	for file_name in os.listdir(epimap_input_dir):
		if file_name.endswith('_abc_links_by_single_cell_type_hg38_liftover.tsv') == False:
			continue
		cell_type = file_name.split('_abc_links_by_single_cell_type_hg38_liftover')[0]
		cell_types.append(cell_type)
	cell_types = np.sort(cell_types)
	file_names = []
	for cell_type in cell_types:
		file_name = epimap_input_dir + cell_type + '_abc_links_by_single_cell_type_hg38_liftover.tsv'
		file_names.append(file_name)
	return cell_types, np.asarray(file_names)


def generate_epimap_chromosome(epimap_cell_types, epimap_file_names, correlation_threshold, chrom_string):
	# Initialize chromosome
	chromosome = ['NULL']*300000000

	for itera, epimap_cell_type in enumerate(epimap_cell_types):
		print(epimap_cell_type)
		epimap_file_name = epimap_file_names[itera]
		f = open(epimap_file_name)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			if len(data) != 8:
				print('assumptino eorroor')
				pdb.set_trace()
			line_chrom_string = data[4]
			if line_chrom_string != chrom_string:
				continue
			correlation = float(data[2])
			if correlation < 0.0:
				print('assumption eroror')
				pdb.set_trace()
			if correlation < correlation_threshold:
				continue
			start_pos = int(data[5])
			end_pos = int(data[6])
			if start_pos > end_pos:
				print('assumption eororr')
				pdb.set_trace()

			for genome_pos in range(start_pos, end_pos):
				if chromosome[genome_pos] == 'NULL':
					chromosome[genome_pos] = epimap_cell_type
				else:
					chromosome[genome_pos] = chromosome[genome_pos] + ',' + epimap_cell_type
		f.close()
	return chromosome

def extract_ukbb_studies_and_component_files(global_ukbb_component_file):
	f = open(global_ukbb_component_file)
	studies = []
	files = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		studies.append(data[0])
		files.append(data[1])
	f.close()
	return np.asarray(studies), np.asarray(files)

def variant_id_arr_to_position_arr(variant_ids_arr):
	pos_arr = []
	for variant_id in variant_ids_arr:
		variant_info = variant_id.split('_')
		if len(variant_info) != 4:
			print('assuumptoin oroeoror')
			pdb.set_trace()
		variant_pos = int(variant_info[1])
		pos_arr.append(variant_pos)

	return np.asarray(pos_arr)

def raw_overlap_ukbb_components_with_epimap_in_a_chromosome(ukbb_study, ukbb_study_component_file, chrom_string, epimap_chromosome, epimap_cell_types, t2_dicti):
	# first create mapping from epimap cell type to position
	epimap_cell_type_to_position = {}
	for itera, epimap_cell_type in enumerate(epimap_cell_types):
		epimap_cell_type_to_position[epimap_cell_type] = itera

	# Print header if chrom_1
	if chrom_string == 'chr1':
		t2_dicti[ukbb_study].write('chrom_num\tlead_variant_position\tlead_variant_name\tlead_variant_pi\twindow_id\tcs_variants\tcs_pi')
		for ct in epimap_cell_types:
			t2_dicti[ukbb_study].write('\t' + ct)
		t2_dicti[ukbb_study].write('\n')


	# Stream component file
	f = open(ukbb_study_component_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Each line is a component
		# Extract relevent fields
		line_chrom_string = 'chr' + data[0]
		# ignore components not on current chromosome
		if line_chrom_string != chrom_string:
			continue
		# extract relevent fields
		lead_variant_id = data[2]
		window_id = data[4]
		variant_ids_cs = data[5].split(',')
		variant_positions = variant_id_arr_to_position_arr(variant_ids_cs)
		pi_cs = np.asarray(data[6].split(',')).astype(float)

		# Now extract assignment of this component to each cell type
		un_normalized_assignments = np.zeros(len(epimap_cell_types))

		for variant_iter, variant_position in enumerate(variant_positions):
			variant_pi = pi_cs[variant_iter]
			variant_cell_types_str = epimap_chromosome[variant_position]
			if variant_cell_types_str == 'NULL':
				continue
			variant_cell_types = np.unique(variant_cell_types_str.split(','))
			for variant_cell_type in variant_cell_types:
				un_normalized_assignments[epimap_cell_type_to_position[variant_cell_type]] = un_normalized_assignments[epimap_cell_type_to_position[variant_cell_type]] + 1
		t2_dicti[ukbb_study].write(line + '\t' + '\t'.join(un_normalized_assignments.astype(str)) + '\n')
	f.close()
	return t2_dicti


def overlap_ukbb_components_with_epimap_in_a_chromosome(ukbb_study, ukbb_study_component_file, chrom_string, epimap_chromosome, epimap_cell_types, t_dicti):
	# first create mapping from epimap cell type to position
	epimap_cell_type_to_position = {}
	for itera, epimap_cell_type in enumerate(epimap_cell_types):
		epimap_cell_type_to_position[epimap_cell_type] = itera

	# Print header if chrom_1
	if chrom_string == 'chr1':
		t_dicti[ukbb_study].write('chrom_num\tlead_variant_position\tlead_variant_name\tlead_variant_pi\twindow_id\tcs_variants\tcs_pi\tsum_un_normalized_cell_type_weight')
		for ct in epimap_cell_types:
			t_dicti[ukbb_study].write('\t' + ct)
		t_dicti[ukbb_study].write('\n')


	# Stream component file
	f = open(ukbb_study_component_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Each line is a component
		# Extract relevent fields
		line_chrom_string = 'chr' + data[0]
		# ignore components not on current chromosome
		if line_chrom_string != chrom_string:
			continue
		# extract relevent fields
		lead_variant_id = data[2]
		window_id = data[4]
		variant_ids_cs = data[5].split(',')
		variant_positions = variant_id_arr_to_position_arr(variant_ids_cs)
		pi_cs = np.asarray(data[6].split(',')).astype(float)

		# Now extract assignment of this component to each cell type
		un_normalized_assignments = np.zeros(len(epimap_cell_types))

		for variant_iter, variant_position in enumerate(variant_positions):
			variant_pi = pi_cs[variant_iter]
			variant_cell_types_str = epimap_chromosome[variant_position]
			if variant_cell_types_str == 'NULL':
				continue
			variant_cell_types = np.unique(variant_cell_types_str.split(','))
			for variant_cell_type in variant_cell_types:
				un_normalized_assignments[epimap_cell_type_to_position[variant_cell_type]] = un_normalized_assignments[epimap_cell_type_to_position[variant_cell_type]] + variant_pi
		if np.sum(un_normalized_assignments) > 0:
			normalized_assignments = un_normalized_assignments/np.sum(pi_cs)
			x_cell_type_normalized_assignments = normalized_assignments/np.sum(normalized_assignments)
			max_cell_type_weight = np.sum(normalized_assignments)
		else:
			x_cell_type_normalized_assignments = un_normalized_assignments + (1.0/len(epimap_cell_types))
			max_cell_type_weight = 0.0
		t_dicti[ukbb_study].write(line + '\t' + str(max_cell_type_weight) + '\t' + '\t'.join(x_cell_type_normalized_assignments.astype(str)) + '\n')
	f.close()
	return t_dicti


global_ukbb_component_file = sys.argv[1]
epimap_ukbb_overlap_dir = sys.argv[2]
correlation_threshold = float(sys.argv[3])

ukbb_studies, ukbb_study_component_files = extract_ukbb_studies_and_component_files(global_ukbb_component_file)

t_dicti = {}
t2_dicti = {}
for ukbb_study in ukbb_studies:
	study_output_file = epimap_ukbb_overlap_dir + ukbb_study + '_component_overlap_with_abc_' + str(correlation_threshold) + '.txt'
	t_dicti[ukbb_study] = open(study_output_file,'w')
	study_output2_file = epimap_ukbb_overlap_dir + ukbb_study + '_raw_component_overlap_with_abc_' + str(correlation_threshold) + '.txt'
	t2_dicti[ukbb_study] = open(study_output2_file,'w')


epimap_cell_types, epimap_file_names = get_epimap_cell_types_and_file_names(epimap_ukbb_overlap_dir)

for chrom_num in range(1,23):
	print(chrom_num)
	# extract data structure where element is a position and values are epimap annotations at that position 
	epimap_chromosome = generate_epimap_chromosome(epimap_cell_types, epimap_file_names, correlation_threshold, 'chr' + str(chrom_num))

	for ukbb_study_index, ukbb_study in enumerate(ukbb_studies):
		ukbb_study_component_file = ukbb_study_component_files[ukbb_study_index]
		t_dicti = overlap_ukbb_components_with_epimap_in_a_chromosome(ukbb_study, ukbb_study_component_file, 'chr' + str(chrom_num), epimap_chromosome, epimap_cell_types, t_dicti)
		t2_dicti = raw_overlap_ukbb_components_with_epimap_in_a_chromosome(ukbb_study, ukbb_study_component_file, 'chr' + str(chrom_num), epimap_chromosome, epimap_cell_types, t2_dicti)

for ukbb_study in ukbb_studies:
	t_dicti[ukbb_study].close()
	t2_dicti[ukbb_study].close()
