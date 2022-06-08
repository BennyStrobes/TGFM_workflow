import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pyreadr
import numpy as np 
import os
import sys
import pdb
import pandas as pd

def meta_analysis(effects, weights):
	# From Omer Weissbrod
	#assert method in ['fixed', 'random']
	d = effects
	variances = 1.0/weights

	wt = weights

	#compute summtest
	summ = wt.dot(d) / wt.sum()
	varsum = np.sum(wt*wt*variances) / (np.sum(wt)**2)
	###summtest = summ / np.sqrt(varsum)

	summary=summ
	se_summary=np.sqrt(varsum)

	return summary, se_summary


def meta_analyze_weights(weight_arr, samp_size_arr, snp_arr):
	num_studies = len(samp_size_arr)
	float_samp_size_arr = np.asarray(samp_size_arr)
	if num_studies == 1:
		return weight_arr[0], samp_size_arr
	else:
		for study_num in range(1, num_studies):
			# Error checking
			if np.array_equal(snp_arr[0], snp_arr[study_num]) == False:
				print('assumption errorr')
				pdb.set_trace()
		num_variants = len(weight_arr[0])
		meta_analyzed_weights = []
		for variant_num in range(num_variants):
			variant_weights = []
			for study_num in range(num_studies):
				variant_weights.append(weight_arr[study_num][variant_num])
			variant_weights = np.asarray(variant_weights)
			variant_me_weight, variant_me_weight_se = meta_analysis(variant_weights,float_samp_size_arr)
			meta_analyzed_weights.append(variant_me_weight)
		meta_analyzed_weights = np.asarray(meta_analyzed_weights)
		return meta_analyzed_weights, samp_size_arr



pseudotissue_name = sys.argv[1]
composit_tissue_string = sys.argv[2]
gtex_fusion_weights_data_dir = sys.argv[3]
gtex_fusion_weights_dir = sys.argv[4]
pseudotissue_gtex_fusion_weights_dir = sys.argv[5]


print(pseudotissue_name)


# Extract composit tissues from composit tissue string
composit_tissues = composit_tissue_string.split(',')

genes = {}
# Loop through composit tissues
for composit_tissue in composit_tissues:
	# Get gene summary file for composit tissue
	gene_summary_file = gtex_fusion_weights_data_dir + composit_tissue + '/' + composit_tissue + '_gene_summary.txt'
	head_count = 0
	counter = 0
	f = open(gene_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		chrom_num = data[1]
		tss = data[2]
		gene_weight_file = gtex_fusion_weights_dir + composit_tissue + '/' + composit_tissue + '_' + gene_id + '_1KG_only_fusion_output.wgt.RDat'
		if os.path.exists(gene_weight_file) == False:
			continue
		#result = pyreadr.read_r(gene_weight_file)
		#hsq = np.asarray(result['hsq'])[0][0]
		#hsq_p = np.asarray(result['hsq.pv'])[0,0]
		#if hsq < 0.0 or hsq_p > .01:
			#print('error')
			#pdb.set_trace()
			#continue
		output_line = gene_weight_file + '\t' + composit_tissue + '\t' + gene_id + '\t' + chrom_num + '\t' + tss  + '\n'
		if gene_id not in genes:
			genes[gene_id] = []
		genes[gene_id].append(output_line)
	f.close()


pos_file = pseudotissue_gtex_fusion_weights_dir + pseudotissue_name + '_cis_heritable_genes.pos'
t = open(pos_file,'w')
t.write('WGT\tSNP\tID\tCHR\tP0\tP1\thsq\thsq_p\n')

gene_names = np.asarray([*genes])


for gene_name in gene_names:
	gene_arr = genes[gene_name]
	# Quick error checking
	if len(gene_arr) > len(composit_tissues):
		print('assumption eroror')
		pdb.set_trace()
	weight_arr = []
	samp_size_arr = []
	snp_arr = []
	hsq_arr = []
	hsq_p_arr = []
	cv_r_squared_arr = []
	for gene_ele_string in gene_arr:
		gene_ele_info = gene_ele_string.rstrip().split('\t')
		gene_ele_tissue = gene_ele_info[1]
		# Quick error checking
		if gene_ele_info[2] != gene_name:
			print('assumption erorro')
			pdb.set_trace()
		chrom_num = gene_ele_info[3]
		tss = gene_ele_info[4]
		gene_weight_file = gene_ele_info[0]
		result = pyreadr.read_r(gene_weight_file)

		best_model_index = np.argmax(np.asarray(result['cv.performance'])[0,:])
		cv_r_squared_arr.append(np.asarray(result['cv.performance'])[0, best_model_index])

		snp_arr.append(np.asarray(result['snps']['V2']))
		samp_size_arr.append(np.asarray(result['N.tot']['N.tot'])[0])
		hsq_arr.append(np.asarray(result['hsq'])[1,0])
		hsq_p_arr.append(np.asarray(result['hsq.pv'])[0,0])

		best_weight_vector = np.asarray(result['wgt.matrix'])[:,best_model_index]

		if result['wgt.matrix'].keys()[best_model_index] == 'top1':
			#best_weight_vector[np.abs(best_weight_vector) != np.max(np.abs(best_weight_vector))]=0
			top_index = np.argmax(np.abs(best_weight_vector))
			top_value = best_weight_vector[top_index]
			best_weight_vector = np.zeros(len(best_weight_vector))
			best_weight_vector[top_index] = top_value

		weight_arr.append(best_weight_vector)
	meta_analyzed_weights, meta_analyzed_samp_size_arr = meta_analyze_weights(weight_arr, samp_size_arr, snp_arr)

	new_weight_matrix = pd.DataFrame(data=meta_analyzed_weights, columns=['meta_weights'], index=snp_arr[0])
	new_hsq_mat = pd.DataFrame(data=hsq_arr, columns=['hsq'])
	new_hsq_pv_mat = pd.DataFrame(data=hsq_p_arr, columns=['hsq.pv'])
	new_n_tot_mat = pd.DataFrame(data=meta_analyzed_samp_size_arr, columns=['N.tot'])
	new_cv_mat = pd.DataFrame(data=cv_r_squared_arr, columns=['rsq'])

	result['wgt.matrix'] = new_weight_matrix
	result['hsq'] = new_hsq_mat
	result['N.tot'] = new_n_tot_mat
	result['hsq.pv'] = new_hsq_pv_mat
	result['cv.performance'] = new_cv_mat

	# Save new result file to rds
	meta_analyzed_weight_file = pseudotissue_gtex_fusion_weights_dir + pseudotissue_name + '_' + gene_name + '_1KG_only_fusion_output_meta_analyzed.wgt.RDat'
	pyreadr.write_rdata(meta_analyzed_weight_file, result['wgt.matrix'])
	
	meta_analyzed_snp_file = pseudotissue_gtex_fusion_weights_dir + pseudotissue_name + '_' + gene_name + '_1KG_only_fusion_output_meta_analyzed_snps.RDat'
	pyreadr.write_rdata(meta_analyzed_snp_file, result['snps'])

	t.write(meta_analyzed_weight_file + '\t' + meta_analyzed_snp_file + '\t' + gene_name + '\t' + chrom_num + '\t' + tss + '\t' + tss + '\t' + str(np.mean(hsq_arr)) + '\t' + str(np.mean(hsq_p_arr)) + '\n')
t.close()

