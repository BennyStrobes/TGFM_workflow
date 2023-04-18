import numpy as np 
import os
import sys
import pdb




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
			causal_variants[rsid] = causal_effect_size
			causal_genetic_elements[rsid] = causal_effect_size
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

				causal_genes[gene_tissue_name] = causal_effect_size
				causal_genetic_elements[gene_tissue_name] = causal_effect_size
	f.close()

	return causal_genetic_elements, causal_variants, causal_genes


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

########################
# Command line args
########################
global_window_file = sys.argv[1]
simulation_window_list_file = sys.argv[2]
simulation_number = int(sys.argv[3])
n_windows_per_sim = int(sys.argv[4])
n_causal_genetic_elements_str = sys.argv[5]
simulated_trait_dir = sys.argv[6]
simulation_name_string = sys.argv[7]
simulated_gene_position_file = sys.argv[8]
processed_genotype_data_dir = sys.argv[9]
chrom_num = sys.argv[10]


# Set set
np.random.seed(simulation_number)


#Bim file
bim_file = processed_genotype_data_dir + 'simulated_gwas_data_' + chrom_num + '.bim'


# Extract all genes and their positions
all_genes, all_gene_positions = extract_all_genes_and_their_positions(simulated_gene_position_file)
# Do the same for variants
all_variants, all_variants_positions = extract_all_variants_and_their_positions(bim_file)


# Extract causal genetic elements
causal_variant_file = simulated_trait_dir + simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
causal_gene_file = simulated_trait_dir + simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
causal_genetic_elements_gt, causal_variants, causal_gene_tissue_pairs = extract_dictionary_list_of_causal_genetic_elements(causal_gene_file, causal_variant_file)




# Load in all windows
all_windows_data_raw = np.loadtxt(global_window_file, dtype=str, delimiter='\t')
header = all_windows_data_raw[0,:]
all_windows_data = all_windows_data_raw[1:,:]

# Filter windows to those with significant number of middle var
num_middle_var = all_windows_data[:,-1].astype(float)
subset_windows_data = all_windows_data[num_middle_var > 100,:]

# Number of window subsets
total_windows = subset_windows_data.shape[0]


# Extract number of causal genetic elements per window
n_causal_genetic_elements_per_window = []
valid_windows = []
for window_iter in range(total_windows):
	window_start = int(subset_windows_data[window_iter,1])
	window_end = int(subset_windows_data[window_iter,4])

	n_causal_elements = 0
	# Count genes
	for gene_iter, gene_name in enumerate(all_genes):
		gene_position = all_gene_positions[gene_iter]
		for tissue_iter in range(10):
			full_gene_name = gene_name + '_tissue' + str(tissue_iter)
			if full_gene_name in causal_genetic_elements_gt and gene_position >= window_start and gene_position < window_end:
				n_causal_elements = n_causal_elements + 1
	# Count variants
	for variant_iter, variant_name in enumerate(all_variants):
		variant_position = all_variants_positions[variant_iter]
		if variant_name in causal_genetic_elements_gt and variant_position >= window_start and variant_position < window_end:
			n_causal_elements = n_causal_elements + 1
	n_causal_genetic_elements_per_window.append(n_causal_elements)

	if n_causal_elements >= int(n_causal_genetic_elements_str.split('_')[0]) and n_causal_elements <= int(n_causal_genetic_elements_str.split('_')[1]):
		valid_windows.append(window_iter)

n_causal_genetic_elements_per_window = np.asarray(n_causal_genetic_elements_per_window)
valid_windows = np.asarray(valid_windows)


if len(valid_windows) > n_windows_per_sim:
	randomly_selected_windows = np.random.choice(valid_windows, size=n_windows_per_sim, replace=False)
else:
	print("assumption eroror: not of windows in this bin")
	randomly_selected_windows = np.copy(valid_windows)




# Print to output file
t = open(simulation_window_list_file,'w')
t.write('\t'.join(header) + '\n')
for randomly_selected_window in np.sort(randomly_selected_windows):
	t.write('\t'.join(subset_windows_data[randomly_selected_window,:]) + '\n')
t.close()

