import numpy as np
import os
import sys
import pdb




def get_list_of_gene_models_tested_by_tgfm(gtex_susie_gene_models_dir, gene_type):
	used_genes = {}
	for tissue_name in os.listdir(gtex_susie_gene_models_dir):
		if tissue_name.endswith('.py'):
			continue
		if tissue_name.endswith('.sh'):
			continue
		if tissue_name.startswith('temp'):
			continue
		if tissue_name == 'Testis':
			continue
		tissue_gene_list = gtex_susie_gene_models_dir + tissue_name + '/' + tissue_name + '_' + gene_type + '_pos_file.txt'
		head_count = 0
		f = open(tissue_gene_list)
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count = head_count + 1
				continue
			ensamble_id = data[1].split('.')[0]
			used_genes[ensamble_id] = 1
		f.close()
	return used_genes


def extract_silver_standard_and_background_genes(ldl_silver_standard_gene_set_file, gene_symbol_to_ensamble_id_mapping):
	silver_standard_set = {}
	bgrd_set= {}

	f = open(ldl_silver_standard_gene_set_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split(',')
		if head_count == 0:
			head_count = head_count + 1
			continue
		ens_id = data[1].split('.')[0]
		gene_id = data[0]

		mapped_ens_id_vec = gene_symbol_to_ensamble_id_mapping[gene_id]
		if len(mapped_ens_id_vec) != 1:
			print('assumption eroeoro')
			pdb.set_trace()
		mapped_ens_id = mapped_ens_id_vec[0]


		if ens_id == 'NA':
			ens_id = mapped_ens_id
		else:
			if mapped_ens_id != ens_id:
				pritn('assumption oeroror')
				pdb.set_trace()

		ctwas_pip = data[6]
		if data[-1] == 'known':
			silver_standard_set[ens_id] = ctwas_pip
		else:
			bgrd_set[ens_id] = ctwas_pip
	f.close()
	return silver_standard_set, bgrd_set

def create_mapping_from_gene_name_to_tgfm_gene_pip(tgfm_gene_pip_summary_file):
	mapping = {}
	f = open(tgfm_gene_pip_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_name = data[0]
		pip = float(data[2])
		if gene_name in mapping:
			print('assumptione rororoeoer')
			pdb.set_trace()
		mapping[gene_name] = pip
	f.close()


	return mapping

def get_ensamble_id_and_gene_name(info_str):
	info = info_str.split('"')
	passed1 = False
	passed2 = False
	for ii, ele in enumerate(info):
		if ele == 'gene_id ':
			ensamble_id = info[(ii+1)]
			passed1 = True
		if ele == '; gene_name ':
			gene_name = info[(ii+1)]
			passed2 = True
	if passed1 == False or passed2 == False:
		print('assumption oerororr')
		pdb.set_trace()
	return ensamble_id, gene_name

def create_mapping_from_gene_symbol_id_to_ensamble_id(gene_annotation_file, ):
	f = open(gene_annotation_file)
	gene_name_to_ensg = {}
	for line in f:
		if line.startswith('##'):
			continue
		line = line.rstrip()
		data = line.split('\t')
		if data[2] != 'gene':
			continue
		info_str = data[8]
		ensamble_id, gene_name = get_ensamble_id_and_gene_name(info_str)
		ensamble_id = ensamble_id.split('.')[0]
		if gene_name not in gene_name_to_ensg:
			gene_name_to_ensg[gene_name] = []
			gene_name_to_ensg[gene_name].append(ensamble_id)
		else:
			gene_name_to_ensg[gene_name].append(ensamble_id)
	f.close()
	return gene_name_to_ensg


######################
# Command line args
######################
tgfm_organized_results_dir = sys.argv[1]
gtex_susie_gene_models_dir = sys.argv[2]
gene_annotation_file = sys.argv[3]
ldl_silver_standard_gene_set_file = sys.argv[4]
output_dir = sys.argv[5]


# Trait name (as run in tgfm)
tgfm_trait_name = "biochemistry_LDLdirect"
# Type of gene model used by TGFM
gene_type = "component_gene"

# Get list of genes tested by TGFM
tgfm_gene_model_genes = get_list_of_gene_models_tested_by_tgfm(gtex_susie_gene_models_dir,gene_type)

# Create mapping from gene name to ensamble id
gene_symbol_to_ensamble_id_mapping = create_mapping_from_gene_symbol_id_to_ensamble_id(gene_annotation_file)

# Extract dictionary list of silver standard genes and background genes
silver_standard_genes, background_genes = extract_silver_standard_and_background_genes(ldl_silver_standard_gene_set_file,gene_symbol_to_ensamble_id_mapping)

# Create mapping from gene name to TGFM gene pip
tgfm_gene_pip_summary_file = tgfm_organized_results_dir + 'tgfm_results_' + tgfm_trait_name + '_component_gene_susie_sampler_uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler_tgfm_per_gene_pip_summary.txt'
gene_name_to_tgfm_gene_pip = create_mapping_from_gene_name_to_tgfm_gene_pip(tgfm_gene_pip_summary_file)


# Open outut summary fiel
output_file = output_dir + 'ldl_silver_standard_enrichment_summary.txt'
t = open(output_file,'w')
t.write('gene_name\tsilver_standard\tgene_pip\n')
silver_gene_pips = []
for silver_standard_gene in [*silver_standard_genes]:
	gene_pip = np.nan
	if silver_standard_gene in tgfm_gene_model_genes:
		gene_pip = 0.0
		if silver_standard_gene in gene_name_to_tgfm_gene_pip:
			gene_pip = gene_name_to_tgfm_gene_pip[silver_standard_gene]
	silver_gene_pips.append(gene_pip)
	if np.isnan(gene_pip) == False:
		t.write(silver_standard_gene + '\t' + 'silver standard\t' + str(gene_pip) + '\n')
silver_gene_pips = np.asarray(silver_gene_pips)
background_gene_pips = []
for bgrd_gene in [*background_genes]:
	gene_pip = np.nan
	if bgrd_gene in tgfm_gene_model_genes:
		gene_pip = 0.0
		if bgrd_gene in gene_name_to_tgfm_gene_pip:
			gene_pip = gene_name_to_tgfm_gene_pip[bgrd_gene]
	background_gene_pips.append(gene_pip)
	if np.isnan(gene_pip) == False:
		t.write(silver_standard_gene + '\t' + 'background\t' + str(gene_pip) + '\n')
background_gene_pips = np.asarray(background_gene_pips)
t.close()
print(output_file)


# Open outut summary file for ctwas
output_file = output_dir + 'ldl_silver_standard_ctwas_enrichment_summary.txt'
t = open(output_file,'w')
t.write('gene_name\tsilver_standard\tgene_pip\n')
silver_gene_pips = []
for silver_standard_gene in [*silver_standard_genes]:
	gene_pip = silver_standard_genes[silver_standard_gene]
	if gene_pip == 'NA':
		continue
	gene_pip = float(gene_pip)
	silver_gene_pips.append(gene_pip)
	t.write(silver_standard_gene + '\t' + 'silver standard\t' + str(gene_pip) + '\n')
silver_gene_pips = np.asarray(silver_gene_pips)
background_gene_pips = []
for bgrd_gene in [*background_genes]:
	gene_pip = background_genes[bgrd_gene]
	if gene_pip == 'NA':
		continue
	gene_pip = float(gene_pip)
	background_gene_pips.append(gene_pip)
	t.write(silver_standard_gene + '\t' + 'background\t' + str(gene_pip) + '\n')
background_gene_pips = np.asarray(background_gene_pips)
t.close()



