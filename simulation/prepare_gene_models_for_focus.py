import numpy as np 
import os
import sys
import pdb


def genotype_alignment_error(snp_ids, cis_genotype_data):
	n_snps = snp_ids.shape[0]
	error_bool = False
	if n_snps != cis_genotype_data.shape[0]:
		error_bool = True

	for snp_index in range(n_snps):
		if cis_genotype_data[snp_index,1] != snp_ids[snp_index,1]:
			error_bool = True

		genotype_data_snp_id = 'chr' + cis_genotype_data[snp_index,0] + '_' + cis_genotype_data[snp_index,3] + '_' + cis_genotype_data[snp_index,4] + '_' + cis_genotype_data[snp_index,5]
		
		if snp_ids[snp_index,0] != genotype_data_snp_id:
			error_bool = True
	return error_bool




simulation_number = sys.argv[1]
chrom_num = sys.argv[2]
simulation_name_string = sys.argv[3]
simulated_gene_expression_dir = sys.argv[4]
simulated_learned_gene_models_dir = sys.argv[5]
eqtl_sample_size = sys.argv[6]
processed_genotype_data_dir = sys.argv[7]
focus_gene_models_dir = sys.argv[8]
gene_type = sys.argv[9]


# Extract gene-tissue pairs for this run
gene_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'

# Extract genotype data
genotype_data_file = processed_genotype_data_dir + 'simulated_gwas_data_' + chrom_num + '.bim'
genotype_data = np.loadtxt(genotype_data_file, dtype=str, delimiter='\t')

# Output file
output_summary_file = focus_gene_models_dir + simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_gene_models_python_summary.txt'
t = open(output_summary_file,'w')
t.write('gene_tissue\tgene_id\ttissue_id\tchr\ttss\twgt.matrix_file\tsnps_file\n')


# Stream gene summary file
f = open(gene_summary_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Extract relevent fields
	ensamble_id = data[0]
	line_chrom_num = data[1]
	tss = data[2]
	causal_eqtl_effect_file = data[3]
	cis_snp_id_file = data[4]
	cis_snp_index_file = data[5]
	total_n_genome_snps = int(data[6])


	# Quick error-checking
	if line_chrom_num != chrom_num:
		print('assumption eroorrr')
		pdb.set_trace()

	# Load in gene models for gene across tissues
	fitted_gene_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_gene_model_pmces.npy'
	gene_model_mat = np.load(fitted_gene_file)

	# Load in snp indices
	snp_indices_raw = np.load(cis_snp_index_file)
	snp_indices = np.asarray([False]*total_n_genome_snps)
	snp_indices[snp_indices_raw] = True


	snp_ids = np.load(cis_snp_id_file, allow_pickle=True)
	n_cis_snps = snp_ids.shape[0]

	# Get cis genotype data
	cis_genotype_data = genotype_data[snp_indices, :]

	# Quick error check
	if genotype_alignment_error(snp_ids, cis_genotype_data):
		print('assumption eroror')
		pdb.set_trace()
	if gene_model_mat.shape[1] != snp_ids.shape[0]:
		print('assumptione rororo')
		pdb.set_trace()

	# Save snps file
	# This can be shared across all tissues
	snps_file = focus_gene_models_dir + simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_' + ensamble_id + '_snps_file.txt'
	np.savetxt(snps_file, cis_genotype_data, fmt="%s", delimiter='\t')

	# Total number of tissues for this gene
	n_tiss = gene_model_mat.shape[0]

	# Extract vector of length number of tissues corresponding to whether tissue has a valid component
	valid_components_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + ensamble_id + '_eqtlss_' + str(eqtl_sample_size) + '_gene_valid_susie_comp.npy'
	valid_components = np.load(valid_components_file)

	# Looop through tissues
	for tissue_iter in range(n_tiss):
		# Skip gene, tissue pairs with no gene models
		if np.array_equal(gene_model_mat[tissue_iter,:], np.zeros(n_cis_snps)):
			continue
		if gene_type == 'component_gene' and valid_components[tissue_iter] == 'False':
			continue
		# Save weights
		weights_file = focus_gene_models_dir + simulation_name_string + '_eqtlss_' + str(eqtl_sample_size) + '_' + ensamble_id + '_' + str(tissue_iter) + '_susie_weights.txt'
		np.savetxt(weights_file, gene_model_mat[tissue_iter,:], fmt="%s", delimiter='\n')
		t.write(ensamble_id + '_tissue' + str(tissue_iter) + '\t' + ensamble_id + '\t' + str(tissue_iter) + '\t' + line_chrom_num + '\t' + tss + '\t' + weights_file + '\t' + snps_file + '\n')

f.close()
t.close()