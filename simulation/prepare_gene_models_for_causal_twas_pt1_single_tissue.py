import numpy as np
import os
import sys
import pdb





input_causal_eqtl_effects_file = sys.argv[1]
simulated_learned_gene_models_stem = sys.argv[2]
genotype_bim_file = sys.argv[3]
tmp_pos_file = sys.argv[4]
gene_model_output_root = sys.argv[5]
eqtl_sample_size = sys.argv[6]
tissue_number = sys.argv[7]

tiss_iter = int(tissue_number)

# First load in BIM file as one big matrix
bim_file = np.loadtxt(genotype_bim_file, dtype=str, delimiter='\t')

t = open(tmp_pos_file,'w')
t.write('PANEL\tWGT\tID\tCHR\tP0\tP1\tN\n')


# Now loop through genes in gene file
f = open(input_causal_eqtl_effects_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Extract relevent fields
	ensamble_id = data[0]
	chrom_num = data[1]
	tss = data[2]
	causal_eqtl_effects_file = data[3]
	cis_snp_ids_file = data[4]
	cis_snp_indices_file = data[5]

	# Load in (ordered) cis snp ids for this geene
	cis_indices = np.load(cis_snp_indices_file)
	cis_snp_ids = np.load(cis_snp_ids_file, allow_pickle=True)

	# Quick error checking
	if np.array_equal(bim_file[cis_indices,:][:,1], cis_snp_ids[:,1]) == False:
		print('assumption erroror')
		pdb.set_trace()

	# Get gene bim file
	gene_bim = bim_file[cis_indices,:]

	# Load in Gene models PMCES
	gene_pmces_file = simulated_learned_gene_models_stem + '_' + ensamble_id + '_eqtlss_' + eqtl_sample_size + '_lasso_gene_model_pmces.npy'
	gene_pmces = np.load(gene_pmces_file)

	# Quick error checking
	if gene_pmces.shape[1] != gene_bim.shape[0]:
		print('assumption eroror')
		pdb.set_trace()

	# Now loop through tissues
	gt_vec = gene_pmces[tiss_iter,:]
	if np.var(gt_vec) == 0.0:
		continue
	# Gene has gene model
	# Extract snp indices that have non-zero effect on gene
	non_zero_indices = np.where(gt_vec != 0.0)[0]
	# Filter gene_bim and gt_vec to these indices
	gt_vec_sub = gt_vec[non_zero_indices]
	gene_bim_sub = gene_bim[non_zero_indices,:]
	# Get int vector of snp positions for gene
	snp_pos = gene_bim_sub[:,3].astype(int)
	start_pos = np.min(snp_pos)
	end_pos = np.max(snp_pos)

	# Save gt_vec_sub
	gt_vec_file = gene_model_output_root + '.' + ensamble_id + '_tissue_' + str(tiss_iter) + '_snp_gene_effects.txt'
	np.savetxt(gt_vec_file, gt_vec_sub, delimiter='\t', fmt="%s")
	# Save gene_bim_sub
	gt_bim_file = gene_model_output_root + '.' + ensamble_id + '_tissue_' + str(tiss_iter) + '_snp_bim.txt'
	np.savetxt(gt_bim_file, gene_bim_sub, delimiter='\t', fmt="%s")
	
	t.write('multitissue\tmultitissue/multitissue.' + ensamble_id + '_tissue_' + str(tiss_iter) + '\t' + ensamble_id + '_tissue_' + str(tiss_iter) + '\t' + chrom_num + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t100\n')


f.close()
t.close()