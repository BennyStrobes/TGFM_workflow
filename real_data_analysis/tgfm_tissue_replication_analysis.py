import numpy as np
import os
import sys
import pdb



def get_trait_names(trait_names_file):
	trait_names = []
	f = open(trait_names_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		trait_names.append(data[0])
	f.close()
	return np.asarray(trait_names)


def get_tissue_names(tissue_names_file):
	tissue_names = []
	f = open(tissue_names_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0] == 'Testis':
			continue
		tissue_names.append(data[0])
	f.close()
	return np.asarray(tissue_names)

def get_trait_gene_pairs_for_replication(trait_names, orig_tgfm_results_dir, removed_tissue_name, pip_threshold, suffix1):
	dicti = {}
	# Loop through traits
	for trait_name in trait_names:
		# get file containing list of all gene-tissue pairs for this trait
		trait_gt_pip_file = orig_tgfm_results_dir + '_' + trait_name + '_' + suffix1 + '_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_tgfm_per_gene_tissue_full_pip_summary.txt'
		f = open(trait_gt_pip_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			tissue_name = data[2]
			if tissue_name != removed_tissue_name:
				continue
			gene_name = data[1].split('.')[0]
			gt_pip = float(data[5])
			if gt_pip < pip_threshold:
				continue
			# Quick error check
			if trait_name + ':' + gene_name in dicti:
				print('assumption eororo')
				pdb.set_trace()
			# Add to list
			dicti[trait_name + ':' + gene_name] = gt_pip
		f.close()
	return dicti


def run_replication_analysis(trait_names, original_tissue_names, new_tissue_names, removed_tissue_name, orig_tgfm_results_dir, new_tgfm_results_dir, pip_threshold, suffix1, suffix2, replication_root):
	# First get all trait-gene pairs from removed tissue with PIP > pip threshold in original analysis
	trait_gene_pairs_for_replication = get_trait_gene_pairs_for_replication(trait_names, orig_tgfm_results_dir, removed_tissue_name, pip_threshold, suffix1)

	# Open raw replication output file handle
	raw_replication_output_file = replication_root + 'raw_replication_results.txt'
	t = open(raw_replication_output_file,'w')
	t.write('trait_name\tgene_name\tremoved_tissue\tremoved_tissue_pip\tnew_tissue\tnew_tissue_pip\n')

	# Loop through replication tissues
	for tmp_tissue_name in new_tissue_names:
		avg = []
		for trait_name in trait_names:
			trait_gt_pip_file = new_tgfm_results_dir + '_' + trait_name + '_' + suffix2 + '_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_tgfm_per_gene_tissue_full_pip_summary.txt'
			f = open(trait_gt_pip_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				tissue_name = data[2]
				if tissue_name != tmp_tissue_name:
					continue
				gene_name = data[1].split('.')[0]
				gt_pip = float(data[5])
				if trait_name + ':' + gene_name not in trait_gene_pairs_for_replication:
					continue
				orig_pip = trait_gene_pairs_for_replication[trait_name + ':' + gene_name]
				t.write(trait_name + '\t' + gene_name + '\t' + removed_tissue_name + '\t' + str(orig_pip) + '\t' + tissue_name + '\t' + str(gt_pip) + '\n')
				avg.append(gt_pip)
			f.close()

		print('#########')
		print(tmp_tissue_name)
		print(np.mean(avg))
		print(len(avg))
	t.close()
	return



#####################
# Command line args
#####################
trait_names_file = sys.argv[1]
orig_tgfm_results_dir = sys.argv[2]
new_tgfm_results_dir = sys.argv[3]
original_tissue_names_file = sys.argv[4]
new_tissue_names_file = sys.argv[5]
removed_tissue_name = sys.argv[6]
replication_output_root = sys.argv[7]
suffix1 = sys.argv[8]
suffix2 = sys.argv[9]



# Extract relevent info
trait_names = get_trait_names(trait_names_file)
original_tissue_names = get_tissue_names(original_tissue_names_file)
new_tissue_names = get_tissue_names(new_tissue_names_file)


pip_threshold = 0.5
# run_replication_analysis(trait_names, original_tissue_names, new_tissue_names, removed_tissue_name, orig_tgfm_results_dir, new_tgfm_results_dir, pip_threshold, suffix1, suffix2, replication_output_root + '_pip_' + str(pip_threshold) + '_')





