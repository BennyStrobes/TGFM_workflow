import numpy as np 
import os
import sys
import pdb




def create_ukbb_trait_names_supp_table(output_file, trait_names_file):
	t = open(output_file,'w')
	f = open(trait_names_file)
	head_count = 0
	sample_sizes = []

	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		t.write(data[0] + '\t' + data[4] + '\t' + data[2] + '\t' + data[3] + '\n')
		if head_count > 0:
			sample_sizes.append(float(data[2]))
		head_count = head_count + 1
	print(np.mean(sample_sizes))
	f.close()
	t.close()
	return

def get_composit_tissue_ss(tissue_cov_file):
	aa = np.loadtxt(tissue_cov_file,dtype=str,delimiter='\t')
	samples = aa[0,1:]
	return len(samples)


def create_gtex_tissue_info_supp_table(output_file, gtex_pseudotissue_file, gtex_covariate_dir):
	t = open(output_file,'w')
	f = open(gtex_pseudotissue_file)
	sample_size_arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write('pseudotissue_name\tpseudotissue_sample_size\tcomposit_tissues\tcomposit_tissues_sample_size\n')
			continue
		pseudotissue_name = data[0]
		pseudotissue_sample_size = data[1]
		composit_tissue_string = data[3]
		if pseudotissue_name == 'Testis':
			continue
		composit_tissue_arr = composit_tissue_string.split(',')
		composit_tissue_ss_arr = []
		for composit_tissue in composit_tissue_arr:
			tissue_cov_file = gtex_covariate_dir + composit_tissue + '_covariates.txt'
			composit_tissue_ss = get_composit_tissue_ss(tissue_cov_file)
			composit_tissue_ss_arr.append(composit_tissue_ss)
		composit_tissue_ss_arr = np.asarray(composit_tissue_ss_arr)
		# Quick error check
		if float(pseudotissue_sample_size) != sum(composit_tissue_ss_arr):
			if pseudotissue_name != 'Artery_Coronary':
				print('assumption eroror')
				pdb.set_trace()

		sample_size_arr.append(sum(composit_tissue_ss_arr))
		
		# print to output
		t.write(pseudotissue_name + '\t' + str(int(np.sum(composit_tissue_ss_arr))) + '\t' + composit_tissue_string + '\t' + ','.join(composit_tissue_ss_arr.astype(str)) + '\n')
	f.close()
	t.close()
	print(output_file)
	print(np.mean(sample_size_arr))
	return


def create_supp_table_figure3_numerical(output_file, tgfm_organized_results_dir):
	independent_traits = np.asarray(["body_HEIGHTz", "blood_MEAN_PLATELET_VOL", "bmd_HEEL_TSCOREz", "blood_MEAN_CORPUSCULAR_HEMOGLOBIN", "blood_MONOCYTE_COUNT", "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT", "pigment_HAIR", "lung_FEV1FVCzSMOKE", "body_BALDING1", "biochemistry_Cholesterol", "bp_DIASTOLICadjMEDz", "lung_FVCzSMOKE", "repro_MENARCHE_AGE", "disease_ALLERGY_ECZEMA_DIAGNOSED", "other_MORNINGPERSON", "repro_NumberChildrenEverBorn_Pooled"])
	t = open(output_file,'w')
	t.write('trait_name\tgenetic_element_class\tTGFM_PIP\n')

	for trait_name in independent_traits:
		element_class = 'variant'
		variant_input_file = tgfm_organized_results_dir + 'tgfm_results_' + trait_name + '_component_gene_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_tgfm_n_causal_' + element_class + 's_cross_pip_threshold_sqrt_plot_input.txt'
		f = open(variant_input_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			if float(data[1]) < .2:
				continue
			t.write(trait_name + '\t' + element_class + '\t' + data[1] + '\n')
		f.close()

		element_class = 'gene_tissue'
		variant_input_file = tgfm_organized_results_dir + 'tgfm_results_' + trait_name + '_component_gene_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_tgfm_n_causal_' + element_class + '_pairs_cross_pip_threshold_sqrt_plot_input.txt'
		f = open(variant_input_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			if float(data[1]) < .2:
				continue
			t.write(trait_name + '\t' + element_class + '\t' + data[1] + '\n')
		f.close()

		element_class = 'gene'
		variant_input_file = tgfm_organized_results_dir + 'tgfm_results_' + trait_name + '_component_gene_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_tgfm_n_causal_' + element_class + 's_cross_pip_threshold_sqrt_plot_input.txt'
		f = open(variant_input_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			if float(data[1]) < .2:
				continue
			t.write(trait_name + '\t' + element_class + '\t' + data[1] + '\n')
		f.close()
	t.close()
	f = open(output_file)
	head_count = 0
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		namer = data[0] + ':' + data[1]
		if namer not in dicti:
			dicti[namer] =1 
		else:
			dicti[namer] = dicti[namer] + 1
	f.close()
	print(dicti)
	return




supp_data_dir = sys.argv[1]
trait_names_file = sys.argv[2]
gtex_pseudotissue_file = sys.argv[3]
gtex_covariate_dir = sys.argv[4]
tgfm_organized_results_dir = sys.argv[5]


#########################################
# Create suppTable_ukbb_trait_names.txt
#########################################
output_file = supp_data_dir + 'suppTable_ukbb_trait_names.txt'
create_ukbb_trait_names_supp_table(output_file, trait_names_file)


#########################################
# Create suppTable_gtex_tissue_info.txt
#########################################
output_file = supp_data_dir + 'suppTable_gtex_tissue_info.txt'
#create_gtex_tissue_info_supp_table(output_file, gtex_pseudotissue_file, gtex_covariate_dir)


#########################################
# Create suppTable_figure3_numerical.txt
#########################################
output_file = supp_data_dir + 'suppTable_figure3_numerical.txt'
#create_supp_table_figure3_numerical(output_file, tgfm_organized_results_dir)
