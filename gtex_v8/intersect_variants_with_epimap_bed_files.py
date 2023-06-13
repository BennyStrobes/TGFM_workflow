import numpy as np 
import os
import sys
import pdb
import scipy.stats




def extract_tissue_names(gtex_pseudotissue_file, ignore_testis=True):
	f = open(gtex_pseudotissue_file)
	head_count = 0
	tissue_names = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0] == 'Testis' and ignore_testis:
			continue
		tissue_names.append(data[0])
	f.close()
	return np.asarray(tissue_names)



def extract_track_names(epimap_track_file):
	f = open(epimap_track_file)
	head_count = 0
	track_names = []
	track_cts = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		track_names.append(data[0])
		track_cts.append(data[3])
	f.close()

	return np.asarray(track_names), np.asarray(track_cts)

def get_number_of_overlapping_variants(overlap_file):
	f = open(overlap_file)
	counter = 0
	for line in f:
		counter = counter + 1
	f.close()
	return counter


def intersect_nm_snps_with_epimap_tracks(track_names, epimap_data_dir, epimap_enrichment_dir, pip_thresh):
	nm_variant_file = epimap_enrichment_dir + 'tgfm_non_mediating_variants_' + str(pip_thresh) + '_' + str(pip_thresh) + '.txt'
	n_nm_var = np.loadtxt(nm_variant_file, dtype=str).shape[0]

	mapping = {}
	# Loop through tracks
	for track_name in track_names:
		track_file = epimap_data_dir + 'Epimap.' + track_name + '.hg38.bed'
		overlap_file = epimap_data_dir + 'nm_variant_' + str(pip_thresh) + '_' + str(pip_thresh) + '_' + track_name + '_overlap.txt'

		# Bedtools command
		bedtools_cmd = 'bedtools intersect -a ' + nm_variant_file + ' -b ' + track_file + ' > ' + overlap_file
		os.system(bedtools_cmd)

		# Extract number of overlapping variants
		n_overlapping_nm_variants = get_number_of_overlapping_variants(overlap_file)

		mapping[track_name] = (n_overlapping_nm_variants, n_nm_var)

	return mapping


def intersect_tissue_snps_with_epimap_tracks(tissue_name, track_names, track_cts, epimap_data_dir, epimap_enrichment_dir, pip_thresh, tissue_track_enrichment_summary_file, epi_map_track_name_to_nm_overlap_numbers):
	t = open(tissue_track_enrichment_summary_file,'w')
	t.write('track_name\ttrack_cell_type\taa\tbb\tcc\tdd\torat\torat_p\n')

	tissue_variant_file = epimap_enrichment_dir + 'tgfm_' + tissue_name + '_mediating_variants_' + str(pip_thresh) + '_' + str(pip_thresh) + '.txt'
	n_tissue_var = get_number_of_overlapping_variants(tissue_variant_file)
	
	# Loop through tracks
	for ii, track_name in enumerate(track_names):
		track_ct = track_cts[ii]
		track_file = epimap_data_dir + 'Epimap.' + track_name + '.hg38.bed'
		overlap_file = epimap_data_dir + tissue_name + '_variant_' + str(pip_thresh) + '_' + str(pip_thresh) + '_' + track_name + '_overlap.txt'

		# Bedtools command
		bedtools_cmd = 'bedtools intersect -a ' + tissue_variant_file + ' -b ' + track_file + ' > ' + overlap_file
		os.system(bedtools_cmd)

		# Extract number of overlapping variants
		n_overlapping_nm_variants = get_number_of_overlapping_variants(overlap_file)

		aa = n_overlapping_nm_variants
		bb = n_tissue_var
		cc, dd = epi_map_track_name_to_nm_overlap_numbers[track_name]

		orat, orat_p = scipy.stats.fisher_exact(np.asarray([[aa,bb],[cc,dd]]))
		t.write(track_name + '\t' + track_ct + '\t' + str(aa) + '\t' + str(bb) + '\t' + str(cc) + '\t' + str(dd) + '\t' + str(orat) + '\t' + str(orat_p) + '\n')
	t.close()
	return

def order_summary_file_by_odds_ratios(tissue_track_enrichment_summary_file, tissue_track_enrichment_summary_or_ordered_file):
	aa = np.loadtxt(tissue_track_enrichment_summary_file, dtype=str, delimiter='\t')
	header = aa[0,:]
	raw_data = aa[1:,:]
	orat = raw_data[:,-2].astype(float)
	raw_data_no_nan = raw_data[np.isnan(orat) == False,:]
	orat_no_nan = raw_data_no_nan[:,-2].astype(float)
	ordering = np.argsort(-orat_no_nan)

	t = open(tissue_track_enrichment_summary_or_ordered_file,'w')
	t.write('\t'.join(header) + '\n')
	for indexer in ordering:
		t.write('\t'.join(raw_data_no_nan[indexer, :]) + '\n')
	t.close()
	return

gtex_pseudotissue_file = sys.argv[1]
epimap_track_file = sys.argv[2]
epimap_data_dir = sys.argv[3]
epimap_enrichment_dir = sys.argv[4]


pip_thresh = .25

#########################
# Extract tissue names
tissue_names = extract_tissue_names(gtex_pseudotissue_file)

#########################
# Extract track names
track_names, track_cts = extract_track_names(epimap_track_file)
track_names = track_names

#########################
# intersect nm snps with epimap tracks
#epi_map_track_name_to_nm_overlap_numbers = intersect_nm_snps_with_epimap_tracks(track_names, epimap_data_dir, epimap_enrichment_dir, pip_thresh)

#########################
# intersect tissue snps with epimap tracks
for tissue_name in tissue_names:
	print(tissue_name)
	tissue_track_enrichment_summary_file = epimap_enrichment_dir + tissue_name + '_' + str(pip_thresh) + '_epimap_enrichment_summary.txt'
	#intersect_tissue_snps_with_epimap_tracks(tissue_name, track_names, track_cts, epimap_data_dir, epimap_enrichment_dir, pip_thresh, tissue_track_enrichment_summary_file, epi_map_track_name_to_nm_overlap_numbers)


	# Order output file by odds ratios
	tissue_track_enrichment_summary_or_ordered_file = epimap_enrichment_dir + tissue_name + '_' + str(pip_thresh) + '_epimap_enrichment_summary_or_ordered.txt'
	order_summary_file_by_odds_ratios(tissue_track_enrichment_summary_file, tissue_track_enrichment_summary_or_ordered_file)






