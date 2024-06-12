import numpy as np 
import os
import sys
import pdb
import scipy.stats as stats
import time


def create_mapping_from_pos_to_gwas_t_statistic(gwas_sumstat_file):
	f = open(gwas_sumstat_file)
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid= data[1]
		pos = data[2]
		a1 = data[3]
		a2 = data[4]
		t_stat = float(data[6])/np.sqrt(float(data[7]))

		if pos in mapping:
			print('assumption eororor')
			pdb.set_trace()
		mapping[pos] = (t_stat, a1, a2, rsid)
	f.close()
	return mapping

def make_temporary_jlim_gwas_ss_file(tmp_jlim_gwas_summary_stat_file, snp_ids, pos_to_gwas_t_statistic, gwas_ss):
	t = open(tmp_jlim_gwas_summary_stat_file,'w')
	t.write('CHR\tBP\tSNP\tP\n')
	pvalues = []

	for ii, snp_id in enumerate(snp_ids):
		pos = snp_id.split('_')[1]
		tmp_a1 = snp_id.split('_')[2]
		tmp_a2 = snp_id.split('_')[3]

		gwas_t_stat = pos_to_gwas_t_statistic[pos][0]
		a1 = pos_to_gwas_t_statistic[pos][1]
		a2 = pos_to_gwas_t_statistic[pos][2]

		if tmp_a1 != a1:
			print('assumption eroror')
			pdb.set_trace()
		if tmp_a2 != a2:
			print('assumption eororro')
			pdb.set_trace()

		# One-tailed p-value
		p_value_one_tailed = stats.t.sf(np.abs(gwas_t_stat), gwas_ss-2)
		# Two-tailed p-value
		p_value_two_tailed = p_value_one_tailed * 2

		if p_value_two_tailed > 1:
			print('assumption erroro')
			pdb.set_trace()

		t.write('1\t' + pos + '\t' + snp_id + '\t' + str(p_value_two_tailed) + '\n')
		pvalues.append(p_value_two_tailed)
	t.close()
	return np.asarray(pvalues)

def make_temporary_jlim_eqtl_ss_file(tmp_jlim_eqtl_summary_stat_file, snp_ids, t_stat, tissue_ss, gene_name):
	# Quick error check
	if len(snp_ids) != len(t_stat):
		print('assumption erroror')
		pdb.set_trace()
	# open output file handles
	t = open(tmp_jlim_eqtl_summary_stat_file,'w')
	t.write('CHR\tBP\tSNP\tP\tGene\n')

	pvalues = []
	for ii, snp_id in enumerate(snp_ids):
		pos = snp_id.split('_')[1]

		# One-tailed p-value
		p_value_one_tailed = stats.t.sf(np.abs(t_stat[ii]), tissue_ss-2)
		# Two-tailed p-value
		p_value_two_tailed = p_value_one_tailed * 2


		t.write('1\t' + pos + '\t' + snp_id + '\t' + str(p_value_two_tailed) + '\t' + gene_name + '\n')
		pvalues.append(p_value_two_tailed)

	t.close()
	return np.asarray(pvalues)


def make_temporary_jlim_eqtl_ss_file_cross_tissues(tmp_jlim_eqtl_summary_stat_file, snp_ids, t_stat, eqtl_sss, gene_id):
	# open output file handles
	t = open(tmp_jlim_eqtl_summary_stat_file,'w')
	t.write('CHR\tBP\tSNP\tP\tGene\n')

	for tissue_iter, tissue_ss in enumerate(eqtl_sss):
		for ii, snp_id in enumerate(snp_ids):
			pos = snp_id.split('_')[1]

			# One-tailed p-value
			p_value_one_tailed = stats.t.sf(np.abs(t_stat[tissue_iter,ii]), tissue_ss-2)
			# Two-tailed p-value
			p_value_two_tailed = p_value_one_tailed * 2

			if p_value_two_tailed > 1:
				print('asssumption eroror')
				pdb.set_trace()

			t.write('1\t' + pos + '\t' + snp_id + '\t' + str(p_value_two_tailed) + '\t' + gene_id + '_tissue' + str(tissue_iter) + '\n')

	t.close()
	return

def make_temporary_jlim_eqtl_ss_file_cross_tissues_with_eqtl_sample_size_equal_to_320(tmp_jlim_eqtl_summary_stat_file, snp_ids, t_stat, eqtl_sss, gene_id):
	# open output file handles
	t = open(tmp_jlim_eqtl_summary_stat_file,'w')
	t.write('CHR\tBP\tSNP\tP\tGene\n')

	for tissue_iter, tissue_ss in enumerate(eqtl_sss):
		if tissue_ss != 320:
			continue
		for ii, snp_id in enumerate(snp_ids):
			pos = snp_id.split('_')[1]

			# One-tailed p-value
			p_value_one_tailed = stats.t.sf(np.abs(t_stat[tissue_iter,ii]), tissue_ss-2)
			# Two-tailed p-value
			p_value_two_tailed = p_value_one_tailed * 2

			if p_value_two_tailed > 1:
				print('assumption eroror')
				pdb.set_trace()

			t.write('1\t' + pos + '\t' + snp_id + '\t' + str(p_value_two_tailed) + '\t' + gene_id + '_tissue' + str(tissue_iter) + '\n')

	t.close()
	return


def run_jlim_with_ref_genotype_data(gwas_summary_stat_file, eqtl_summary_stat_file, jlim_output_file, tissue_ss, gene_id, tss, jlim_ref_geno_dir, jlim_source_code_dir, index_snp_pos):
	# Get CMD string
	jlim_cmd_str = 'sh ' + jlim_source_code_dir + 'run_jlim.sh'
	jlim_cmd_str = jlim_cmd_str + ' --maintr-file ' + gwas_summary_stat_file + ' --sectr-file ' + eqtl_summary_stat_file
	jlim_cmd_str = jlim_cmd_str + ' --ref-ld ' + jlim_ref_geno_dir
	jlim_cmd_str = jlim_cmd_str + ' --index-snp 1:' + str(index_snp_pos)
	jlim_cmd_str = jlim_cmd_str + ' --output-file ' + jlim_output_file
	jlim_cmd_str = jlim_cmd_str + ' --sectr-sample-size ' + str(tissue_ss)
	#jlim_cmd_str = jlim_cmd_str + ' --window-size 500000'
	jlim_cmd_str = jlim_cmd_str + ' --manual-window-boundary ' + str(int(tss)-100000) + '-' + str(int(tss)+100000)

	# Run JLIM
	os.system(jlim_cmd_str)

	return


def print_genes_results_to_global_output_file(t_global, tmp_jlim_output_file):
	f = open(tmp_jlim_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_name = data[8]
		t_stat = data[2]
		pvalue = data[3]
		t_global.write(gene_name + '\t' + pvalue + '\t' + t_stat + '\n')
	f.close()
	return t_global


def create_gene_to_eqtl_info_mapping(eqtl_summary_stat_file):
	mapping = {}
	used_genes = []
	f = open(eqtl_summary_stat_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent fields
		gene_id = data[0]
		rs_id = data[1]
		snp_id = data[2]
		t_stats = np.asarray(data[3].split(';')).astype(float)

		# Gene_id not seen before
		if gene_id not in mapping:
			mapping[gene_id] = {}
			mapping[gene_id]['t'] = []
			mapping[gene_id]['snp_names'] = []
			used_genes.append(gene_id)
		mapping[gene_id]['t'].append(t_stats)
		mapping[gene_id]['snp_names'].append(snp_id)
	f.close()

	for gene_id in used_genes:
		mapping[gene_id]['t'] = np.transpose(np.asarray(mapping[gene_id]['t']))
		mapping[gene_id]['snp_names'] = np.asarray(mapping[gene_id]['snp_names'])

	return mapping



########################
# Command line args
########################
simulation_number = sys.argv[1]
chrom_num = sys.argv[2]
eqtl_sample_size = sys.argv[3]
simulation_name_string = sys.argv[4]
processed_genotype_data_dir = sys.argv[5]
processed_jlim_genotype_data_dir = sys.argv[6]
jlim_ref_geno_dir = sys.argv[7]
jlim_source_code_dir = sys.argv[8]
simulated_learned_gene_models_dir = sys.argv[9]
simulated_jlim_data_dir = sys.argv[10] # Results dir1
simulated_jlim_results_dir = sys.argv[11]  # Results dir2
simulated_gwas_dir = sys.argv[12]


# GWAS summary stat file
gwas_sumstat_file = simulated_gwas_dir + simulation_name_string + '_merged_gwas_summary_stats.txt'
pos_to_gwas_t_statistic = create_mapping_from_pos_to_gwas_t_statistic(gwas_sumstat_file)

# Gene summary file
gene_summary_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_eqtl_results_summary.txt'

# eQTL summary stat file
eqtl_summary_stat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_eqtl_summary_stats.txt'

# Create gene to eQTL info mapping
gene_to_eqtl_info = create_gene_to_eqtl_info_mapping(eqtl_summary_stat_file)

# Open up global results results
global_output_file = simulated_jlim_results_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_jlim_results.txt'
t_global = open(global_output_file,'w')
# Print header
t_global.write('gene_name\tjlim_p\tjlim_t\n')

# Loop through genes
head_count = 0
f = open(gene_summary_file)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue

	# Extract relevent fields
	gene_id = data[0]

	tss = data[1]
	eqtl_sss = np.asarray(data[2].split(';')).astype(int)

	# Extract eqtl info for this gnee
	t_mat = gene_to_eqtl_info[gene_id]['t']
	snp_ids = gene_to_eqtl_info[gene_id]['snp_names']

	# 1. Create temporary jlim gwas summary stat file
	tmp_jlim_gwas_summary_stat_file = simulated_jlim_data_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_tmp_gwas_summary_stat.txt'
	gwas_pvalues = make_temporary_jlim_gwas_ss_file(tmp_jlim_gwas_summary_stat_file, snp_ids, pos_to_gwas_t_statistic, 100000)


	# Simulation settings where sample size is constant across tissues
	if eqtl_sample_size != 'realistic':
		# 2. Create temporary jlim eqtl summary stat file for eqtls in this gene (across tissues)
		tmp_jlim_eqtl_summary_stat_file = simulated_jlim_data_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_tmp_eqtl_summary_stat.txt'
		make_temporary_jlim_eqtl_ss_file_cross_tissues(tmp_jlim_eqtl_summary_stat_file, snp_ids, t_mat, eqtl_sss, gene_id)

		# Get index snp
		index_snp_pos = snp_ids[np.argmax(np.abs(t_mat[0,:]))].split('_')[1]

		# 3. Run JLIM using reference LD
		tmp_jlim_output_file = simulated_jlim_results_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_tmp_jlim_res.txt'
		run_jlim_with_ref_genotype_data(tmp_jlim_gwas_summary_stat_file, tmp_jlim_eqtl_summary_stat_file, tmp_jlim_output_file, int(eqtl_sample_size), gene_id, tss, jlim_ref_geno_dir, jlim_source_code_dir, index_snp_pos)

		# Print gene's results to global output file
		t_global = print_genes_results_to_global_output_file(t_global, tmp_jlim_output_file)
		t_global.flush()

	elif eqtl_sample_size == 'realistic':
		# Get index snp
		index_snp_pos = snp_ids[np.argmax(np.abs(t_mat[0,:]))].split('_')[1]

		###########################################
		# Part A
		###########################################
		# For computational efficiency
		# First run together all tissues (5 of them) with eqtl sample size == 320
		# 2. Create temporary jlim eqtl summary stat file for eqtls in this gene (across tissues)
		tmp_jlim_eqtl_summary_stat_file = simulated_jlim_data_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_tmp_eqtl_summary_stat.txt'
		make_temporary_jlim_eqtl_ss_file_cross_tissues_with_eqtl_sample_size_equal_to_320(tmp_jlim_eqtl_summary_stat_file, snp_ids, t_mat, eqtl_sss, gene_id)

		# 3. Run JLIM using reference LD
		tmp_jlim_output_file = simulated_jlim_results_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_tmp_jlim_res_a.txt'
		run_jlim_with_ref_genotype_data(tmp_jlim_gwas_summary_stat_file, tmp_jlim_eqtl_summary_stat_file, tmp_jlim_output_file, 320, gene_id, tss, jlim_ref_geno_dir, jlim_source_code_dir, index_snp_pos)

		# Print gene's results to global output file
		t_global = print_genes_results_to_global_output_file(t_global, tmp_jlim_output_file)
		t_global.flush()

		###########################################
		# Part B
		###########################################
		# Now run seperately for each tissue with eqtl sample size != 320
		for tissue_iter, tissue_ss in enumerate(eqtl_sss):
			if tissue_ss == 320:
				continue
			# 2. Create temporary jlim eqtl summary stat file for eqtls in this tissue
			tmp_jlim_eqtl_summary_stat_file = simulated_jlim_data_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_tmp_eqtl_summary_stat.txt'
			t_stat = t_mat[tissue_iter,:]
			eqtl_pvalues = make_temporary_jlim_eqtl_ss_file(tmp_jlim_eqtl_summary_stat_file, snp_ids, t_stat, tissue_ss, gene_id + '_tissue' + str(tissue_iter))


			# 3. Run JLIM using reference LD
			tmp_jlim_output_file = simulated_jlim_results_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_tmp_jlim_res_b.txt'
			run_jlim_with_ref_genotype_data(tmp_jlim_gwas_summary_stat_file, tmp_jlim_eqtl_summary_stat_file, tmp_jlim_output_file, tissue_ss, gene_id, tss, jlim_ref_geno_dir, jlim_source_code_dir, index_snp_pos)

			# Print gene's results to global output file
			t_global = print_genes_results_to_global_output_file(t_global, tmp_jlim_output_file)
		t_global.flush()

	else:
		print('assumption erroror')
		pdb.set_trace()
	

f.close()
t_global.close()

