import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pyreadr
import numpy as np 
import os
import sys
import pdb

def fail_multi_tissue_heritability_filter(hsq, hsq_p, pvalue_thresh):
	ntiss = hsq.shape[0]

	booler = True

	for tiss_num in range(ntiss):
		if hsq[tiss_num, 0] > 0 and hsq_p[tiss_num, 0] < pvalue_thresh:
			booler = False
	return booler 

def generate_all_gene_pos_file(gene_summary_file, cis_heritable_pos_file, pseudotissue_gene_model_dir, pseudotissue_name):
	# Open output file handle
	t = open(cis_heritable_pos_file,'w')
	t.write('WGT\tID\tCHR\tP0\tP1\thsq\thsq_p\tcomponent_boolean\tcomponent_indices\n')

	# Loop through genes in gene summary file
	head_count = 0
	counter = 0
	f = open(gene_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Extract relevent fields corresponding to this gene
		gene_id = data[0]
		chrom_num = data[1]
		tss = data[2]
		counter = counter + 1

		
		# Get name of file corresponding to gene model
		gene_weight_file = pseudotissue_gene_model_dir + pseudotissue_name + '_' + gene_id + '_1KG_only_fusion_output.wgt.RDat'
		if os.path.exists(gene_weight_file) == False:
			continue
			
		result = pyreadr.read_r(gene_weight_file)

		'''
		if np.var(np.asarray(result['susie_mu'])) <= 0.0:
			continue
		if np.sum(np.sum(np.asarray(result['susie_V']))) <= 0.0:
			continue
		if np.asarray(result['susie_converged']['susie_converged'])[0] == False:
			continue
		'''
		hsq = np.asarray(result['combined_heritabilities'])
		hsq_p = np.asarray(result['combined_heritabilities_p'])

		#if fail_multi_tissue_heritability_filter(hsq, hsq_p, .01):
			#continue

		if 'component_bool' not in result:
			print('skipped')
			continue

		component_boolean = np.asarray(result['component_bool']['component_bool'])[0]
		if component_boolean == True:
			cs_indices = np.asarray(result['susie_cs'])[:,0]
			cs_index_string = ','.join(cs_indices.astype(str))
			component_boolean_string = 'True'
		else:
			cs_index_string = 'NULL'
			component_boolean_string = 'False'

		output_line = gene_weight_file + '\t' + gene_id + '\t' + chrom_num + '\t' + tss + '\t' + tss + '\t' + str(np.mean(hsq)) + '\t' + str(np.mean(hsq_p)) + '\t' + component_boolean_string + '\t' + cs_index_string + '\n'
		t.write(output_line)
	f.close()
	t.close()
 
def generate_all_nonzero_gene_pos_file(gene_summary_file, cis_heritable_pos_file, pseudotissue_gene_model_dir, pseudotissue_name):
	# Open output file handle
	t = open(cis_heritable_pos_file,'w')
	t.write('WGT\tID\tCHR\tP0\tP1\thsq\thsq_p\tcomponent_boolean\tcomponent_indices\n')

	# Loop through genes in gene summary file
	head_count = 0
	counter = 0
	f = open(gene_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Extract relevent fields corresponding to this gene
		gene_id = data[0]
		chrom_num = data[1]
		tss = data[2]
		counter = counter + 1

		
		# Get name of file corresponding to gene model
		gene_weight_file = pseudotissue_gene_model_dir + pseudotissue_name + '_' + gene_id + '_1KG_only_fusion_output.wgt.RDat'
		if os.path.exists(gene_weight_file) == False:
			continue
			
		result = pyreadr.read_r(gene_weight_file)

		if np.var(np.asarray(result['susie_mu'])) <= 0.0:
			continue
		if np.sum(np.sum(np.asarray(result['susie_V']))) <= 0.0:
			continue
		if np.asarray(result['susie_converged']['susie_converged'])[0] == False:
			continue
		hsq = np.asarray(result['combined_heritabilities'])
		hsq_p = np.asarray(result['combined_heritabilities_p'])

		#if fail_multi_tissue_heritability_filter(hsq, hsq_p, .01):
			#continue

		if 'component_bool' not in result:
			print('skipped')
			continue

		component_boolean = np.asarray(result['component_bool']['component_bool'])[0]
		if component_boolean == True:
			cs_indices = np.asarray(result['susie_cs'])[:,0]
			cs_index_string = ','.join(cs_indices.astype(str))
			component_boolean_string = 'True'
		else:
			cs_index_string = 'NULL'
			component_boolean_string = 'False'

		output_line = gene_weight_file + '\t' + gene_id + '\t' + chrom_num + '\t' + tss + '\t' + tss + '\t' + str(np.mean(hsq)) + '\t' + str(np.mean(hsq_p)) + '\t' + component_boolean_string + '\t' + cs_index_string + '\n'
		t.write(output_line)
	f.close()
	t.close()

def generate_cis_heritable_gene_pos_file(gene_summary_file, cis_heritable_pos_file, pseudotissue_gene_model_dir, pseudotissue_name):
	# Open output file handle
	t = open(cis_heritable_pos_file,'w')
	t.write('WGT\tID\tCHR\tP0\tP1\thsq\thsq_p\tcomponent_boolean\tcomponent_indices\n')

	# Loop through genes in gene summary file
	head_count = 0
	counter = 0
	f = open(gene_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Extract relevent fields corresponding to this gene
		gene_id = data[0]
		chrom_num = data[1]
		tss = data[2]
		counter = counter + 1

		
		# Get name of file corresponding to gene model
		gene_weight_file = pseudotissue_gene_model_dir + pseudotissue_name + '_' + gene_id + '_1KG_only_fusion_output.wgt.RDat'
		if os.path.exists(gene_weight_file) == False:
			continue
			
		result = pyreadr.read_r(gene_weight_file)
		if np.var(np.asarray(result['susie_mu'])) <= 0.0:
			continue
		if np.sum(np.sum(np.asarray(result['susie_V']))) <= 0.0:
			continue
		if np.asarray(result['susie_converged']['susie_converged'])[0] == False:
			continue
		hsq = np.asarray(result['combined_heritabilities'])
		hsq_p = np.asarray(result['combined_heritabilities_p'])

		if fail_multi_tissue_heritability_filter(hsq, hsq_p, .01):
			continue

		if 'component_bool' not in result:
			print('skipped')
			continue

		component_boolean = np.asarray(result['component_bool']['component_bool'])[0]
		if component_boolean == True:
			cs_indices = np.asarray(result['susie_cs'])[:,0]
			cs_index_string = ','.join(cs_indices.astype(str))
			component_boolean_string = 'True'
		else:
			cs_index_string = 'NULL'
			component_boolean_string = 'False'

		output_line = gene_weight_file + '\t' + gene_id + '\t' + chrom_num + '\t' + tss + '\t' + tss + '\t' + str(np.mean(hsq)) + '\t' + str(np.mean(hsq_p)) + '\t' + component_boolean_string + '\t' + cs_index_string + '\n'
		t.write(output_line)
	f.close()
	t.close()

def generate_component_gene_pos_file(gene_summary_file, cis_heritable_pos_file, pseudotissue_gene_model_dir, pseudotissue_name):
	# Open output file handle
	t = open(cis_heritable_pos_file,'w')
	t.write('WGT\tID\tCHR\tP0\tP1\thsq\thsq_p\tcomponent_boolean\tcomponent_indices\n')

	# Loop through genes in gene summary file
	head_count = 0
	counter = 0
	f = open(gene_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Extract relevent fields corresponding to this gene
		gene_id = data[0]
		chrom_num = data[1]
		tss = data[2]
		counter = counter + 1

		
		# Get name of file corresponding to gene model
		gene_weight_file = pseudotissue_gene_model_dir + pseudotissue_name + '_' + gene_id + '_1KG_only_fusion_output.wgt.RDat'
		if os.path.exists(gene_weight_file) == False:
			continue
			
		result = pyreadr.read_r(gene_weight_file)
		if np.var(np.asarray(result['susie_mu'])) <= 0.0:
			continue
		if np.sum(np.sum(np.asarray(result['susie_V']))) <= 0.0:
			continue
		if np.asarray(result['susie_converged']['susie_converged'])[0] == False:
			continue
		hsq = np.asarray(result['combined_heritabilities'])
		hsq_p = np.asarray(result['combined_heritabilities_p'])

		if 'component_bool' not in result:
			print('skipped')
			continue
		component_boolean = np.asarray(result['component_bool']['component_bool'])[0]
		if component_boolean == True:
			cs_indices = np.asarray(result['susie_cs'])[:,0]
			cs_index_string = ','.join(cs_indices.astype(str))
			component_boolean_string = 'True'
		else:
			cs_index_string = 'NULL'
			component_boolean_string = 'False'

		if component_boolean == True:
			output_line = gene_weight_file + '\t' + gene_id + '\t' + chrom_num + '\t' + tss + '\t' + tss + '\t' + str(np.mean(hsq)) + '\t' + str(np.mean(hsq_p)) + '\t' + component_boolean_string + '\t' + cs_index_string + '\n'
			t.write(output_line)
	f.close()
	t.close()





pseudotissue_name = sys.argv[1]
gtex_pseudotissue_gene_model_input_dir = sys.argv[2]
gtex_susie_gene_models_dir = sys.argv[3]

# Directory containing pseudotissue gene models
pseudotissue_gene_model_dir = gtex_susie_gene_models_dir + pseudotissue_name + '/'


# gene_summary_file is a file containing a list of genes that we attempted to make gene models for in this pseudotissue
gene_summary_file = gtex_pseudotissue_gene_model_input_dir + pseudotissue_name + '_gene_summary.txt'


# Generate pos file with only component genes
gene_component_pos_file = pseudotissue_gene_model_dir + pseudotissue_name + '_component_gene_pos_file.txt'
generate_component_gene_pos_file(gene_summary_file, gene_component_pos_file, pseudotissue_gene_model_dir, pseudotissue_name)

# Generate pos file with only cis-heritable genes
cis_heritable_pos_file = pseudotissue_gene_model_dir + pseudotissue_name + '_cis_heritable_gene_pos_file.txt'
generate_cis_heritable_gene_pos_file(gene_summary_file, cis_heritable_pos_file, pseudotissue_gene_model_dir, pseudotissue_name)


# Generate pos file with only all genes
all_nonzero_pos_file = pseudotissue_gene_model_dir + pseudotissue_name + '_all_non_zero_gene_pos_file.txt'
#generate_all_nonzero_gene_pos_file(gene_summary_file, all_nonzero_pos_file, pseudotissue_gene_model_dir, pseudotissue_name)

all_pos_file = pseudotissue_gene_model_dir + pseudotissue_name + '_all_gene_pos_file.txt'
#generate_all_gene_pos_file(gene_summary_file, all_pos_file, pseudotissue_gene_model_dir, pseudotissue_name)



