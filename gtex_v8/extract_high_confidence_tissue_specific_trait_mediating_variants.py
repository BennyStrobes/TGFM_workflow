import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pyreadr
import numpy as np 
import os
import sys
import pdb
import pickle
import time
import gzip



def extract_tissue_names(gtex_pseudotissue_file, ignore_testis=True):
	f = open(gtex_pseudotissue_file)
	head_count = 0
	tissue_names = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0] == 'Testis' and ignore_testis:
			continue
		tissue_names.append(data[0])
	f.close()
	return np.asarray(tissue_names)


def extract_trait_names(trait_list_file):
	f = open(trait_list_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0])
	f.close()
	return np.asarray(arr)

def extract_high_tgfm_pip_gene_tissue_pairs_from_specified_tissue(tgfm_results_summary_file, tissue_name, gene_trait_pip_thresh):
	f = open(tgfm_results_summary_file)
	head_count = 0
	gene_tissue_pairs = []
	gene_tissue_pair_pips = []

	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		line_tissue_name = data[2]
		if line_tissue_name != tissue_name:
			continue
		line_pip = float(data[-1])
		if line_pip < gene_trait_pip_thresh:
			continue
		gene_tissue_name = data[0]
		gene_tissue_pairs.append(gene_tissue_name)
		gene_tissue_pair_pips.append(line_pip)
	f.close()
	return np.asarray(gene_tissue_pairs), np.asarray(gene_tissue_pair_pips)


def compute_pips(alpha_mat):
	LL = alpha_mat.shape[0]
	n_elements = alpha_mat.shape[1]
	anti_pips = np.ones(n_elements)

	for component_iter in range(LL):
		anti_pips = anti_pips*(1.0 - alpha_mat[component_iter,:])
	pips = 1.0 - anti_pips
	return pips

def extract_tissue_specific_variants(tissue_name, tgfm_results_dir, gtex_susie_gene_models_dir, trait_names, tissue_variants_output_file, gene_trait_pip_thresh=.5, variant_gene_pip_thresh=.5):
	# Open output file handle
	t = open(tissue_variants_output_file,'w')
	# Print header
	#t.write('variant_name\tchrom_num\tvariant_position\ttrait_name\tvariant_gene_pip\tgene_trait_pip\n')


	# Loop through traits
	for trait_name in trait_names:
		# TGFM results summary file
		tgfm_results_summary_file = tgfm_results_dir + 'tgfm_results_' + trait_name + '_component_gene_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_tgfm_per_gene_tissue_pip_summary.txt'

		# Extract high pip gene-tissue pairs from specified tissue
		gene_tissue_pairs, gene_tissue_pair_pips = extract_high_tgfm_pip_gene_tissue_pairs_from_specified_tissue(tgfm_results_summary_file, tissue_name, gene_trait_pip_thresh)
		

		# Loop through gene-tissue pairs
		for ii, gene_tissue_pair in enumerate(gene_tissue_pairs):
			# Extract pip for this gene-tissue pair
			gene_tissue_pair_pip = gene_tissue_pair_pips[ii]

			# Get name of gene corresponding to this gene-tissue pair
			gene_name = gene_tissue_pair.split('_')[0]

			# Load in gene-model for this gene tissue pair
			gene_model_file = gtex_susie_gene_models_dir + tissue_name + '/' + tissue_name + '_' + gene_name + '_1KG_only_fusion_output.wgt.RDat'
			# Load gene-tissue model from gene-tissue weight file
			gene_tissue_model = pyreadr.read_r(gene_model_file)

			# Gene model variant names
			gene_model_variant_names = np.asarray(gene_tissue_model['variant_names'])[:,0]

			# Extract pips for variants on this gene-tissue pair
			alpha_mat = np.asarray(gene_tissue_model['susie_alpha'])	
			variant_gene_pips = compute_pips(alpha_mat)

			# Loop through variants
			for var_iter, variant_name in enumerate(gene_model_variant_names):
				# Only consider high pip variants for gene
				variant_pip = variant_gene_pips[var_iter]
				if variant_pip >= variant_gene_pip_thresh:
					# If pass pip threshold, print to output
					var_info = variant_name.split('_')
					var_chrom = var_info[0]
					var_pos = var_info[1]
					t.write(var_chrom + '\t' + var_pos + '\t' + str(int(var_pos) + 1) + '\t' + variant_name + '\t' + trait_name + '\t' + str(variant_pip) + '\t' + str(gene_tissue_pair_pip) + '\n')
	t.close()

	return

def extract_high_tgfm_pip_non_mediated_variants(tgfm_results_summary_file, variant_pip_thresh):
	f = open(tgfm_results_summary_file)
	variants = []
	variant_pips = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 3:
			continue
		ele_names = np.asarray(data[1].split(';'))
		ele_pips = np.asarray(data[2].split(';')).astype(float)
		for ii, ele_name in enumerate(ele_names):
			ele_pip = ele_pips[ii]
			if ele_name.startswith('ENSG') == False and ele_pip >= variant_pip_thresh:
				variants.append(ele_name)
				variant_pips.append(ele_pip)
	f.close()
	return np.asarray(variants), np.asarray(variant_pips)


def extract_non_mediated_variants(tgfm_results_dir, trait_names, non_mediated_variants_output_file, variant_pip_thresh=.5):
	# Open output file handle
	t = open(non_mediated_variants_output_file,'w')
	# Print header
	#t.write('variant_name\tchrom_num\tvariant_position\ttrait_name\tvariant_trait_pip\n')


	# Loop through traits
	for trait_name in trait_names:
		# TGFM results summary file
		tgfm_results_summary_file = tgfm_results_dir + 'tgfm_results_' + trait_name + '_component_gene_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_tgfm_pip_summary.txt'

		# Extract high pip gene-tissue pairs from specified tissue
		nm_variants, nm_variant_pips = extract_high_tgfm_pip_non_mediated_variants(tgfm_results_summary_file, variant_pip_thresh)

		for var_iter, variant_name in enumerate(nm_variants):
			variant_pip = nm_variant_pips[var_iter]
			# If pass pip threshold, print to output
			var_info = variant_name.split('_')
			var_chrom = var_info[0]
			var_pos = var_info[1]
			t.write(var_chrom + '\t' + var_pos + '\t' + str(int(var_pos) + 1) + '\t' + variant_name + '\t' + trait_name + '\t' + str(variant_pip) + '\n')
	t.close()
	return

def extract_null_variants_nearby_genes(gtex_susie_gene_models_dir, tissue_names, null_variants_nearby_genes_output_file, n_var_per_tissue=100):
	# Open output file handle
	t = open(null_variants_nearby_genes_output_file,'w')
	# Print header
	#t.write('variant_name\tchrom_num\tvariant_position\ttrait_name\tvariant_trait_pip\n')

	# Loop through tissues
	for tissue_name in tissue_names:
		print(tissue_name)
		tissue_gene_file = gtex_susie_gene_models_dir + tissue_name + '/' + tissue_name + '_component_gene_pos_file.txt'
		tissue_gene_data = np.loadtxt(tissue_gene_file,dtype=str,delimiter='\t')[1:,:]
		n_genes = tissue_gene_data.shape[0]
		randomly_selected_genes = np.random.choice(np.arange(n_genes), size=n_var_per_tissue, replace=False, p=None)
		subset_tissue_gene_data = tissue_gene_data[randomly_selected_genes, :]
		for row_iter in range(subset_tissue_gene_data.shape[0]):
			gene_model_file = subset_tissue_gene_data[row_iter,0]
			
			# Load gene-tissue model from gene-tissue weight file
			gene_tissue_model = pyreadr.read_r(gene_model_file)

			# Gene model variant names
			gene_model_variant_names = np.asarray(gene_tissue_model['variant_names'])[:,0]

			variant_name = np.random.choice(gene_model_variant_names)
			var_info = variant_name.split('_')
			var_chrom = var_info[0]
			var_pos = var_info[1]
			t.write(var_chrom + '\t' + var_pos + '\t' + str(int(var_pos) + 1) + '\t' + variant_name + '\t' + tissue_name + '\t' + 'NA' + '\n')
	t.close()

#####################
# Command line args
#####################
tgfm_results_dir = sys.argv[1]
gtex_susie_gene_models_dir = sys.argv[2]
gtex_pseudotissue_file = sys.argv[3]
trait_list_file = sys.argv[4]
output_dir = sys.argv[5]



#########################
# Extract tissue names
tissue_names = extract_tissue_names(gtex_pseudotissue_file)

#########################
# Extract trait names
trait_names = extract_trait_names(trait_list_file)


# Extract null variants nearby genes

# Output file
null_variants_nearby_genes_output_file = output_dir + 'tgfm_null_variants_nearby_genes.txt'
# Extract non-mediating variants across traits
extract_null_variants_nearby_genes(gtex_susie_gene_models_dir, tissue_names, null_variants_nearby_genes_output_file)
'''
#########################
# Set filters on what to call a trait mediating variant from a given tissue
gene_trait_pip_thresh = .5
variant_gene_pip_thresh = .5

# Loop through tissues
for tissue_name in tissue_names:
	print(tissue_name)
	################################################################
	# Extract trait-mediating variants corresponding to this tissue
	# Output file
	tissue_variants_output_file = output_dir + 'tgfm_' + tissue_name + '_mediating_variants_' + str(variant_gene_pip_thresh) + '_' + str(gene_trait_pip_thresh) + '.txt'
	# Extract variants across traits
	extract_tissue_specific_variants(tissue_name, tgfm_results_dir, gtex_susie_gene_models_dir, trait_names, tissue_variants_output_file, gene_trait_pip_thresh=gene_trait_pip_thresh, variant_gene_pip_thresh=variant_gene_pip_thresh)


#########################
# Set filters on what to call a trait mediating variant from a given tissue
gene_trait_pip_thresh = .25
variant_gene_pip_thresh = .25

# Loop through tissues
for tissue_name in tissue_names:
	print(tissue_name)
	################################################################
	# Extract trait-mediating variants corresponding to this tissue
	# Output file
	tissue_variants_output_file = output_dir + 'tgfm_' + tissue_name + '_mediating_variants_' + str(variant_gene_pip_thresh) + '_' + str(gene_trait_pip_thresh) + '.txt'
	# Extract variants across traits
	extract_tissue_specific_variants(tissue_name, tgfm_results_dir, gtex_susie_gene_models_dir, trait_names, tissue_variants_output_file, gene_trait_pip_thresh=gene_trait_pip_thresh, variant_gene_pip_thresh=variant_gene_pip_thresh)


#########################
# Set filters on what to call a trait mediating variant from a given tissue
gene_trait_pip_thresh = .5
variant_gene_pip_thresh = .5

# Output file
non_mediated_variants_output_file = output_dir + 'tgfm_non_mediating_variants_' + str(variant_gene_pip_thresh) + '_' + str(gene_trait_pip_thresh) + '.txt'

# Extract non-mediating variants across traits
extract_non_mediated_variants(tgfm_results_dir, trait_names, non_mediated_variants_output_file, variant_pip_thresh=variant_gene_pip_thresh)



#########################
# Set filters on what to call a trait mediating variant from a given tissue
gene_trait_pip_thresh = .25
variant_gene_pip_thresh = .25

# Output file
non_mediated_variants_output_file = output_dir + 'tgfm_non_mediating_variants_' + str(variant_gene_pip_thresh) + '_' + str(gene_trait_pip_thresh) + '.txt'
# Extract non-mediating variants across traits
extract_non_mediated_variants(tgfm_results_dir, trait_names, non_mediated_variants_output_file, variant_pip_thresh=variant_gene_pip_thresh)
'''







