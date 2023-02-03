import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os
import pdb
from pandas_plink import read_plink1_bin



def load_in_sldsc_tau_from_sldsc_results_file(real_sldsc_results_file):
	aa = np.loadtxt(real_sldsc_results_file,dtype=str)
	tau = aa[1:,-3].astype(float)
	return tau

def extract_expected_per_snp_heritabilities(sldsc_tau, variant_annotation_file):
	snp_h2 = []
	snp_names = []
	f = open(variant_annotation_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rs_id = data[2]
		snp_names.append(rs_id)
		snp_anno_vec = np.asarray(data[4:]).astype(float)
		snp_h2.append(np.dot(sldsc_tau, snp_anno_vec))
	f.close()
	return np.asarray(snp_h2), np.asarray(snp_names)


# Simulate non-mediated variant causal effect sizes
def simulate_non_mediated_variant_causal_effect_sizes(variant_annotation_file, real_sldsc_results_file, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, non_mediated_causal_effects_output):
	# Extract underlying heritability model
	sldsc_tau = load_in_sldsc_tau_from_sldsc_results_file(real_sldsc_results_file)

	# Extract expected per snp heritability of snps
	snp_expected_heritability, snp_names = extract_expected_per_snp_heritabilities(sldsc_tau, variant_annotation_file)
	# Set negative per snp heritabilities to zero
	snp_expected_heritability[snp_expected_heritability < 0.0] = 0.0

	# Convert whole thing to probability
	snp_probability = snp_expected_heritability/np.sum(snp_expected_heritability)

	# Calculate number of causal snps
	total_non_mediated_h2 = total_heritability - (fraction_expression_mediated_heritability*total_heritability)
	n_causal_snps = int(np.round(total_non_mediated_h2/per_element_heritability))

	# Total number of snps
	total_num_snps = len(snp_probability)

	# Extract indices of non-mediated causal snps
	non_med_causal_snp_indices = np.random.choice(np.arange(total_num_snps), size=n_causal_snps, replace=False, p=snp_probability)

	# Extract causal effect sizes
	non_mediated_causal_effect_sizes = np.zeros(total_num_snps)  # Initialization
	non_mediated_causal_effect_sizes[non_med_causal_snp_indices] = np.random.normal(loc=0.0, scale=np.sqrt(per_element_heritability),size=n_causal_snps)

	# Print to output file
	t = open(non_mediated_causal_effects_output,'w')
	for var_iter in range(total_num_snps):
		t.write(snp_names[var_iter] + '\t' + str(non_mediated_causal_effect_sizes[var_iter]) + '\n')
	t.close()

def extract_ordered_gene_names_from_gene_summary_file(simulated_expression_summary_file):
	f = open(simulated_expression_summary_file)
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

# Simulate mediated gene-expression causal effect sizes
def simulate_expression_mediated_gene_causal_effect_sizes(simulated_expression_summary_file, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, expression_mediated_causal_effects_output):
	# Extract list of ordered gene names
	ordered_gene_names = extract_ordered_gene_names_from_gene_summary_file(simulated_expression_summary_file)
	total_genes = len(ordered_gene_names)

	# Calculate number of causal genes
	total_mediated_h2 = (fraction_expression_mediated_heritability*total_heritability)
	total_n_causal_genes = int(np.round(total_mediated_h2/per_element_heritability))

	# We are going to assume 50% of gene-mediated h2 is coming from tissue0 and 50% of gene-mediated h2 is coming from tissue3
	# And as a reminder there are 10 tissues
	n_causal_genes_t0 = int(np.round(.5*total_mediated_h2/per_element_heritability))
	n_causal_genes_t3 = int(np.round(.5*total_mediated_h2/per_element_heritability))

	# Extract indices of non-mediated causal genes
	t0_med_causal_gene_indices = np.random.choice(np.arange(total_genes), size=n_causal_genes_t0, replace=False)
	t3_med_causal_gene_indices = np.random.choice(np.arange(total_genes), size=n_causal_genes_t3, replace=False)

	# Extract gene causal effect sizes
	gene_causal_effect_sizes = np.zeros((total_genes, 10))
	gene_causal_effect_sizes[t0_med_causal_gene_indices, 0] = np.random.normal(loc=0.0, scale=np.sqrt(per_element_heritability),size=n_causal_genes_t0)
	gene_causal_effect_sizes[t3_med_causal_gene_indices, 3] = np.random.normal(loc=0.0, scale=np.sqrt(per_element_heritability),size=n_causal_genes_t3)
	
	# Print to output file
	t = open(expression_mediated_causal_effects_output,'w')
	for gene_iter in range(total_genes):
		t.write(ordered_gene_names[gene_iter] + '\t' + '\t'.join(gene_causal_effect_sizes[gene_iter,:].astype(str)) + '\n')
	t.close()

def mean_impute_and_standardize_genotype(G_obj_geno):
	G_obj_geno_stand = np.copy(G_obj_geno)
	ncol = G_obj_geno_stand.shape[1]
	#n_missing = []
	for col_iter in range(ncol):
		nan_indices = np.isnan(G_obj_geno[:,col_iter])
		non_nan_mean = np.mean(G_obj_geno[nan_indices==False, col_iter])
		G_obj_geno_stand[nan_indices, col_iter] = non_nan_mean
		G_obj_geno_stand[:, col_iter] = (G_obj_geno_stand[:, col_iter] - np.mean(G_obj_geno_stand[:, col_iter]))/np.std(G_obj_geno_stand[:, col_iter])
		#n_missing.append(np.sum(nan_indices))
	#n_missing = np.asarray(n_missing)

	#G_obj_geno_stand = (G_obj_geno_stand -np.mean(G_obj_geno_stand,axis=0))/np.std(G_obj_geno_stand,axis=0)

	return G_obj_geno_stand

def compute_expression_mediated_trait_values_for_single_gene(genotype_obj, gene_trait_effect_size_vec, causal_eqtl_effect_sizes, eqtl_snp_indices, computation_version='v1'):
	# Total number of snps
	global_n_snps = len(eqtl_snp_indices)

	# Extract eqtl variant names
	eqtl_index_positions = np.arange(global_n_snps)[eqtl_snp_indices].astype(int).astype(str)
	eqtl_index_names = []
	for eqtl_index_position in eqtl_index_positions:
		eqtl_index_names.append('variant' + eqtl_index_position)
	eqtl_index_names = np.asarray(eqtl_index_names)

	# Extract genotype matrix for these snps
	eqtl_genotype = np.asarray(genotype_obj.sel(variant=eqtl_index_names))
	stand_eqtl_genotype = mean_impute_and_standardize_genotype(eqtl_genotype)

	# NOTE: THE FOLLOWING TWO VERSIONS GIVE EQUIVALENT RESULTS
	# IT WAS MORE PLACED HERE FOR  A TEACHING EXCERCISE
	if computation_version == 'v1':
		# Get genetically predicted gene expression (in each tissue)
		genetically_predicted_gene_expression = np.dot(stand_eqtl_genotype, causal_eqtl_effect_sizes)  # Matrix of number of samplesXnumber of tissues

		# Standardize genetically predicted gene expression
		standardized_genetically_predicted_gene_expression = genetically_predicted_gene_expression/np.std(genetically_predicted_gene_expression,axis=0)

		tissue_corr = np.corrcoef(np.transpose(standardized_genetically_predicted_gene_expression))

		# Get expression mediated trait values
		expr_med_trait_for_current_gene = np.dot(standardized_genetically_predicted_gene_expression, gene_trait_effect_size_vec)
	elif computation_version == 'v2':
		LD_mat = np.corrcoef(np.transpose(stand_eqtl_genotype))
		genetically_predicted_gene_expression_variance = np.diag(np.dot(np.dot(np.transpose(causal_eqtl_effect_sizes), LD_mat), causal_eqtl_effect_sizes))
		variant_trait_effect_sizes_per_tissue = (causal_eqtl_effect_sizes/np.sqrt(genetically_predicted_gene_expression_variance))*gene_trait_effect_size_vec # n_varXn_tiss
		cross_tissue_variant_trait_effect_sizes = np.sum(variant_trait_effect_sizes_per_tissue,axis=1)
		expr_med_trait_for_current_gene = np.dot(stand_eqtl_genotype,cross_tissue_variant_trait_effect_sizes)

	return expr_med_trait_for_current_gene




def compute_non_mediated_variant_mediated_trait_values(non_mediated_variant_causal_effects_file, gwas_plink_stem, variant_mediated_trait_values_output, n_gwas_individuals):
	# Initialize variant mediated trait values
	variant_med_trait = np.zeros(n_gwas_individuals)

	# Load in genotype object
	genotype_obj = read_plink1_bin(gwas_plink_stem + '.bed', gwas_plink_stem + '.bim', gwas_plink_stem + '.fam', verbose=False)

	# Extract non-mediated variant-trait effect sizes
	var_trait_effect_sizes_raw = np.loadtxt(non_mediated_variant_causal_effects_file, dtype=str, delimiter='\t')
	ordered_rsids = var_trait_effect_sizes_raw[:,0]
	var_trait_effect_sizes = var_trait_effect_sizes_raw[:,1].astype(float)


	# Loop through variants
	total_nvar = len(ordered_rsids)
	for var_iter in range(total_nvar):
		print(var_iter)
		# Only consider variants with non-zero variant-trait effect
		var_trait_effect_size = var_trait_effect_sizes[var_iter]
		if var_trait_effect_size == 0.0:
			continue

		# Ok, we are at a variant that has non-zero effect on the trait

		# Extract genotype for this variant
		variant_genotype = np.asarray(genotype_obj.sel(variant='variant' + str(var_iter)))
		std_variant_genotype = np.copy(variant_genotype)
		# Mean impute genotype
		nan_indices = np.isnan(variant_genotype)
		non_nan_mean = np.mean(variant_genotype[nan_indices==False])
		std_variant_genotype[nan_indices] = non_nan_mean
		# Standardize genotype
		std_variant_genotype = (std_variant_genotype - np.mean(std_variant_genotype))/np.std(std_variant_genotype)

		# Get predicted trait effects from this one variant
		trait_effects_from_single_variant = std_variant_genotype*var_trait_effect_size

		# Update global effects
		variant_med_trait = variant_med_trait + trait_effects_from_single_variant

	# Save to output
	np.savetxt(variant_mediated_trait_values_output, variant_med_trait, fmt="%s", delimiter='\n')

	return




def compute_expression_mediated_trait_values(simulated_expression_summary_file, expression_mediated_causal_effects_file, gwas_plink_stem, expression_mediated_trait_values_output, n_gwas_individuals):
	# First extract gene-trait effect sizes
	gene_trait_effect_sizes_raw = np.loadtxt(expression_mediated_causal_effects_file, dtype=str, delimiter='\t')
	ordered_gene_names = gene_trait_effect_sizes_raw[:,0]
	gene_trait_effect_sizes = gene_trait_effect_sizes_raw[:,1:].astype(float)

	# Initialize expression mediated trait values
	expr_med_trait = np.zeros(n_gwas_individuals)

	# Load in genotype object
	genotype_obj = read_plink1_bin(gwas_plink_stem + '.bed', gwas_plink_stem + '.bim', gwas_plink_stem + '.fam', verbose=False)

	#global_tissue_corr = np.zeros((10,10))
	#tissue_corr_counter = 0
	# Loop through genes in simulated_expression_summary_file (should be of same length as ordered_gene_names) and only stop at genes with non-zero gene-trait-effect sizes
	f = open(simulated_expression_summary_file)
	head_count = 0
	gene_counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		print(gene_counter)
		# This corresponds to a single gene
		# Extract relevent fields for this gene
		gene_name = data[0]
		causal_eqtl_effect_file = data[3]
		eqtl_snp_id_file = data[4]
		eqtl_snp_indices_file = data[5]

		# Quick error check
		if gene_name != ordered_gene_names[gene_counter]:
			print('assumption error')
			pdb.set_trace()

		# We can skip genes that have no gene to trait effect
		gene_trait_effect_size_vec = gene_trait_effect_sizes[gene_counter, :]
		if np.array_equal(gene_trait_effect_size_vec, np.zeros(len(gene_trait_effect_size_vec))):
			gene_counter = gene_counter + 1
			continue

		# We are at a gene that has a non-zero effect on the trait
		
		# Load in causal eqtl effect sizes for this gene
		causal_eqtl_effect_sizes = np.load(causal_eqtl_effect_file)
		# Load in "global" indices of snps defining causal eqtl effect sizes
		eqtl_snp_indices = np.load(eqtl_snp_indices_file)

		# Compute expression-mediated trait value coming from this single gene (across tissues)
		expr_med_trait_for_current_gene = compute_expression_mediated_trait_values_for_single_gene(genotype_obj, gene_trait_effect_size_vec, causal_eqtl_effect_sizes, eqtl_snp_indices, computation_version='v1')
		#expr_med_trait_for_current_gene_v2 = compute_expression_mediated_trait_values_for_single_gene(genotype_obj, gene_trait_effect_size_vec, causal_eqtl_effect_sizes, eqtl_snp_indices, computation_version='v2')
		#global_tissue_corr = global_tissue_corr + tissue_corr
		#tissue_corr_counter = tissue_corr_counter + 1

		# Now add current gene to total
		expr_med_trait = expr_med_trait + expr_med_trait_for_current_gene

		# Update gene counter
		gene_counter = gene_counter + 1
	f.close()

	# Save to output
	np.savetxt(expression_mediated_trait_values_output, expr_med_trait, fmt="%s", delimiter='\n')

	return


def simulate_trait_values(expression_mediated_trait_values_file, variant_mediated_trait_values_file, trait_values_output):
	# Load in genetically predicted trait values
	expr_med_trait = np.loadtxt(expression_mediated_trait_values_file)
	non_med_trait = np.loadtxt(variant_mediated_trait_values_file)

	# Calculate total genetic variance
	expr_med_var = np.var(expr_med_trait)
	non_med_var = np.var(non_med_trait)
	total_genetic_var = np.var(expr_med_trait + non_med_trait)

	# Residual variance
	residual_var = 1.0 - total_genetic_var

	# Draw trait values
	trait_values = np.random.normal(loc=(expr_med_trait + non_med_trait), scale=np.sqrt(residual_var))

	# Standardize trait values
	standardized_trait_values = (trait_values - np.mean(trait_values))/np.std(trait_values)

	# Save to output
	np.savetxt(trait_values_output, standardized_trait_values, fmt="%s", delimiter='\n')

	return


# Command line args
simulation_number = int(sys.argv[1])
chrom_num = sys.argv[2]
cis_window = int(sys.argv[3])
simulated_gene_expression_dir = sys.argv[4]
simulation_name_string = sys.argv[5]
processed_genotype_data_dir = sys.argv[6]
ldsc_real_data_results_dir = sys.argv[7]
per_element_heritability = float(sys.argv[8])
total_heritability = float(sys.argv[9])
fraction_expression_mediated_heritability = float(sys.argv[10])
simulated_trait_dir = sys.argv[11]  # Output dir
n_gwas_individuals = int(sys.argv[12])


# Set seed
np.random.seed(simulation_number)


####################################################
# Simulate mediated gene-expression causal effect sizes
####################################################
simulated_expression_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'  # This file contains a line for each gene, and we will use it to select which genes in which tissue are used
expression_mediated_causal_effects_output = simulated_trait_dir + simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
simulate_expression_mediated_gene_causal_effect_sizes(simulated_expression_summary_file, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, expression_mediated_causal_effects_output)


####################################################
# Simulate non-mediated variant causal effect sizes
####################################################
variant_annotation_file = processed_genotype_data_dir + 'baseline.' + chrom_num + '.annot'  # Matrix of annotations describing genetic variants
real_sldsc_results_file = ldsc_real_data_results_dir + 'blood_WHITE_COUNT_sldsc_baseline.results'  # Output of sldsc run on real trait (contains heritability model)
non_mediated_causal_effects_output = simulated_trait_dir + simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
simulate_non_mediated_variant_causal_effect_sizes(variant_annotation_file, real_sldsc_results_file, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, non_mediated_causal_effects_output)


####################################################
# Extract component of trait values coming from genetically predicted gene expression
####################################################
simulated_expression_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'  # This file contains a line for each gene, and also contains causal eqtl effect variant-to-gene effect sizes
expression_mediated_causal_effects_file = simulated_trait_dir + simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'  # File containing gene-to-trait effect sizes
gwas_plink_stem = processed_genotype_data_dir + 'simulated_gwas_data_' + str(chrom_num)  # Genotype directory
expression_mediated_trait_values_output = simulated_trait_dir + simulation_name_string + '_expression_mediated_trait_values.txt'  # Output file
compute_expression_mediated_trait_values(simulated_expression_summary_file, expression_mediated_causal_effects_file, gwas_plink_stem, expression_mediated_trait_values_output, n_gwas_individuals)

####################################################
# Extract component of trait values coming from non-mediated variant effects
####################################################
non_mediated_variant_causal_effects_file = simulated_trait_dir + simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
variant_mediated_trait_values_output = simulated_trait_dir + simulation_name_string + '_non_mediated_variant_mediated_trait_values.txt'  # Output file
compute_non_mediated_variant_mediated_trait_values(non_mediated_variant_causal_effects_file, gwas_plink_stem, variant_mediated_trait_values_output, n_gwas_individuals)




####################################################
# Simulate trait values
####################################################
expression_mediated_trait_values_file = simulated_trait_dir + simulation_name_string + '_expression_mediated_trait_values.txt'  # Now an input file
variant_mediated_trait_values_file = simulated_trait_dir + simulation_name_string + '_non_mediated_variant_mediated_trait_values.txt'  # Now an input file
trait_values_output = simulated_trait_dir + simulation_name_string + '_trait_values.txt'  # Output file
simulate_trait_values(expression_mediated_trait_values_file, variant_mediated_trait_values_file, trait_values_output)


