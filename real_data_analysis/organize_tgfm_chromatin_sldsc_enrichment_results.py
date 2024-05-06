import numpy as np
import os 
import sys
import pdb
import scipy.stats




def get_trait_names(trait_file):
	f = open(trait_file)
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

def get_cell_type_z_score_from_sldsc_file(sldsc_res_file):
	aa = np.loadtxt(sldsc_res_file,dtype=str,delimiter='\t')
	return float(aa[1,-1])


def get_mapping_from_trait_cell_type_pair_to_sldsc_z_score(trait_names, chromatin_cell_type_group_ldsc_dir):
	mapping = {}
	for trait_name in trait_names:
		for cell_type_num in range(1,11):
			sldsc_res_file = chromatin_cell_type_group_ldsc_dir + trait_name + '_sldsc_cell_type_group_' + str(cell_type_num) + '_baselineLD.results'
			z_score = get_cell_type_z_score_from_sldsc_file(sldsc_res_file)
			mapping[trait_name + ':' + str(cell_type_num)] = z_score
	return mapping

def create_tgfm_sldsc_chromatin_overlap_summary_file(tgfm_trait_tissue_significance_file, tissue_name_to_cell_type_group, cell_type_group_to_num, tgfm_fdr,output_file, independent_traits, p_thresh=.05):
	# Open output file handle
	t = open(output_file,'w')
	# Print header to output file
	t.write('trait_name\ttrait_type\ttissue_name\tcell_type_group_name\tcell_type_group_num\ttgfm_sig\ttgfm_p\tsldsc_z\tsldsc_p\n')

	aa = 0
	bb = 0
	# Stream tgfm trait tissue significance file
	f = open(tgfm_trait_tissue_significance_file)
	used_trait_tissues = {}
	valid_traits = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if tgfm_fdr == 0.05 and data[-1] != '**':
			continue
		if tgfm_fdr == 0.2:
			if data[-1] != '**' and data[-1] != '*':
				continue
		tissue_name = data[0]
		trait_name = data[1]
		tgfm_sig = float(data[2])
		cell_type_group = tissue_name_to_cell_type_group[tissue_name]
		cell_type_group_num = cell_type_group_to_num[cell_type_group]
		if cell_type_group == 'Other':
			continue
		trait_type = 'non_independent'
		if trait_name in independent_traits:
			trait_type = 'independent'
		valid_traits[trait_name] = 1
		sldsc_z = trait_cell_type_to_sldsc_z[trait_name + ':' + str(cell_type_group_num)]
		sldsc_p = scipy.stats.norm.sf(sldsc_z) #onesided
		t.write(trait_name + '\t' + trait_type + '\t' + tissue_name + '\t' + cell_type_group + '\t' + str(cell_type_group_num) + '\t' + 'True' + '\t' + str(tgfm_sig) + '\t' + str(sldsc_z) + '\t' + str(sldsc_p) + '\n')
		used_trait_tissues[trait_name + ':' + tissue_name] =1
		if trait_type == 'non_independent':
			continue
		if sldsc_p < p_thresh and sldsc_z > 0:
			aa = aa + 1
		else:
			bb = bb + 1
	f.close()

	cc = 0
	dd = 0
	# Stream tgfm trait tissue significance file
	f = open(tgfm_trait_tissue_significance_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissue_name = data[0]
		trait_name = data[1]
		tgfm_sig = float(data[2])
		cell_type_group = tissue_name_to_cell_type_group[tissue_name]
		cell_type_group_num = cell_type_group_to_num[cell_type_group]
		if trait_name not in valid_traits:
			continue
		if cell_type_group == 'Other':
			continue
		trait_type = 'non_independent'
		if trait_name in independent_traits:
			trait_type = 'independent'
		valid_traits[trait_name] = 1
		sldsc_z = trait_cell_type_to_sldsc_z[trait_name + ':' + str(cell_type_group_num)]
		sldsc_p = scipy.stats.norm.sf(abs(sldsc_z))*2 #twosided
		if trait_name + ':' + tissue_name in used_trait_tissues:
			continue
		t.write(trait_name + '\t' + trait_type + '\t' + tissue_name + '\t' + cell_type_group + '\t' + str(cell_type_group_num) + '\t' + 'False' + '\t' + str(tgfm_sig) + '\t' + str(sldsc_z) + '\t' + str(sldsc_p) + '\n')

		if trait_type == 'non_independent':
			continue

		if sldsc_p < p_thresh and sldsc_z > 0:
			cc = cc + 1
		else:
			dd = dd + 1
	f.close()

	print(aa)
	print(bb)
	print(cc)
	print(dd)

	table = np.asarray([[aa, bb], [cc,dd]])
	print(scipy.stats.fisher_exact(table, alternative='two-sided'))

	t.close()

	print(output_file)

	return


#####################
# Command line args
#####################
tgfm_trait_tissue_significance_file = sys.argv[1]
trait_file = sys.argv[2]
chromatin_cell_type_group_ldsc_dir = sys.argv[3]
gtex_pseudotissue_and_cell_type_group_file = sys.argv[4]



# First get trait names
trait_names = get_trait_names(trait_file)

# Get mapping from cell type group to cell type group number
cell_type_group_to_num = {}
cell_type_group_to_num['Adrenal_Pancreas'] = 1
cell_type_group_to_num['Cardiovascular'] = 2
cell_type_group_to_num['Brain'] = 3
cell_type_group_to_num['Connective_Bone'] = 4
cell_type_group_to_num['Digestive'] = 5
cell_type_group_to_num['Blood_Immune'] = 6
cell_type_group_to_num['Kidney'] = 7
cell_type_group_to_num['Liver'] = 8
cell_type_group_to_num['Other'] = 9
cell_type_group_to_num['Muscle'] = 10

# Get list of independent traits
independent_traits = {}
independent_traits['body_HEIGHTz'] = 1
independent_traits['blood_MONOCYTE_COUNT'] = 1
independent_traits['lung_FEV1FVCzSMOKE'] = 1
independent_traits['lung_FVCzSMOKE'] = 1
independent_traits['blood_MEAN_PLATELET_VOL'] = 1
independent_traits['bp_DIASTOLICadjMEDz'] = 1
independent_traits['blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT'] = 1
independent_traits['biochemistry_Cholesterol'] = 1
independent_traits['blood_MEAN_CORPUSCULAR_HEMOGLOBIN'] = 1
independent_traits['bmd_HEEL_TSCOREz'] = 1
independent_traits['pigment_HAIR'] = 1
independent_traits['disease_ALLERGY_ECZEMA_DIAGNOSED'] = 1
independent_traits['repro_MENARCHE_AGE'] = 1
independent_traits['body_BALDING1'] = 1
independent_traits['repro_NumberChildrenEverBorn_Pooled'] = 1
independent_traits['other_MORNINGPERSON'] = 1


# Get mapping gtex tissue name to cell type group
tissue_name_to_cell_type_group = {}
f = open(gtex_pseudotissue_and_cell_type_group_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	tissue_name = data[0]
	cell_type_group = data[4]
	if cell_type_group not in cell_type_group_to_num:
		print('assumptione roror')
		pdb.set_trace()
	tissue_name_to_cell_type_group[tissue_name] = cell_type_group
f.close()

# Get mapping from trait_cell-type pair to sldsc z score
trait_cell_type_to_sldsc_z =  get_mapping_from_trait_cell_type_pair_to_sldsc_z_score(trait_names, chromatin_cell_type_group_ldsc_dir)


# Extract TGFM trait-tissue pairs with corresponding cell type p-value

print(tgfm_trait_tissue_significance_file)
tgfm_fdr  = .05
output_file = chromatin_cell_type_group_ldsc_dir + 'tgfm_sldsc_chromatin_overlap_summary_' + str(tgfm_fdr) + '.txt'
create_tgfm_sldsc_chromatin_overlap_summary_file(tgfm_trait_tissue_significance_file, tissue_name_to_cell_type_group, cell_type_group_to_num, tgfm_fdr,output_file, independent_traits)
tgfm_fdr  = .2
output_file = chromatin_cell_type_group_ldsc_dir + 'tgfm_sldsc_chromatin_overlap_summary_' + str(tgfm_fdr) + '.txt'
create_tgfm_sldsc_chromatin_overlap_summary_file(tgfm_trait_tissue_significance_file, tissue_name_to_cell_type_group, cell_type_group_to_num, tgfm_fdr,output_file, independent_traits)
tgfm_fdr  = .5
output_file = chromatin_cell_type_group_ldsc_dir + 'tgfm_sldsc_chromatin_overlap_summary_' + str(tgfm_fdr) + '.txt'
create_tgfm_sldsc_chromatin_overlap_summary_file(tgfm_trait_tissue_significance_file, tissue_name_to_cell_type_group, cell_type_group_to_num, tgfm_fdr,output_file, independent_traits)
tgfm_fdr  = .75
output_file = chromatin_cell_type_group_ldsc_dir + 'tgfm_sldsc_chromatin_overlap_summary_' + str(tgfm_fdr) + '.txt'
create_tgfm_sldsc_chromatin_overlap_summary_file(tgfm_trait_tissue_significance_file, tissue_name_to_cell_type_group, cell_type_group_to_num, tgfm_fdr,output_file, independent_traits)
tgfm_fdr  = .9
output_file = chromatin_cell_type_group_ldsc_dir + 'tgfm_sldsc_chromatin_overlap_summary_' + str(tgfm_fdr) + '.txt'
create_tgfm_sldsc_chromatin_overlap_summary_file(tgfm_trait_tissue_significance_file, tissue_name_to_cell_type_group, cell_type_group_to_num, tgfm_fdr,output_file, independent_traits)
tgfm_fdr  = .95
output_file = chromatin_cell_type_group_ldsc_dir + 'tgfm_sldsc_chromatin_overlap_summary_' + str(tgfm_fdr) + '.txt'
create_tgfm_sldsc_chromatin_overlap_summary_file(tgfm_trait_tissue_significance_file, tissue_name_to_cell_type_group, cell_type_group_to_num, tgfm_fdr,output_file, independent_traits)


