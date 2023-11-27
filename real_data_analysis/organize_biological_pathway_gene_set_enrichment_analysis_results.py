import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import sys
import pdb
import statsmodels.stats.multitest
from scipy.stats import rankdata
import scipy.stats


def compute_fdr(p_vals):
	ranked_p_values = rankdata(p_vals)
	fdr = p_vals * len(p_vals) / ranked_p_values
	fdr[fdr > 1] = 1
	return fdr

def load_in_gsea_results(trait_gsea_results_file):
	f = open(trait_gsea_results_file)
	res = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		aa = float(data[2])
		bb = float(data[3])
		cc = float(data[4])
		dd = float(data[5])
		fe_orat, fe_pvalue = scipy.stats.fisher_exact(np.asarray([[aa,bb],[cc,dd]]), alternative='greater')
		if aa + bb < 20:
			continue
		if data[0].startswith('MP') == False and data[0].startswith('GO') == False:
			continue
		orat = float(data[6])
		pvalue = fe_pvalue
		data[7] = str(fe_pvalue)
		line = '\t'.join(data)
		res.append((pvalue,orat,line))
	f.close()
	return res


def get_fdr_significant_results(gsea_results, fdr_thresh):
	pvalz = []
	for gsea_result in gsea_results:
		pvalz.append(gsea_result[0])
	pvalz = np.asarray(pvalz)

	best_val = 0
	n_tests = len(pvalz)
	for ii, pvaler in enumerate(pvalz):
		counter = ii + 1
		#if pvaler <= .005:
		#if pvaler <= (counter/n_tests)*fdr_thresh:
		if pvaler <= .005:
			best_val = counter
	return gsea_results[:best_val]







######################
# Command line args
######################
biological_pathway_enrichment_dir = sys.argv[1]



independent_traits= np.asarray(["body_HEIGHTz", "blood_MEAN_PLATELET_VOL", "bmd_HEEL_TSCOREz", "blood_MEAN_CORPUSCULAR_HEMOGLOBIN", "blood_MONOCYTE_COUNT", "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT", "pigment_HAIR", "lung_FEV1FVCzSMOKE", "body_BALDING1", "biochemistry_Cholesterol", "bp_DIASTOLICadjMEDz", "lung_FVCzSMOKE", "repro_MENARCHE_AGE", "disease_ALLERGY_ECZEMA_DIAGNOSED", "other_MORNINGPERSON", "repro_NumberChildrenEverBorn_Pooled"])

pip_thresh = .5

for independent_trait in independent_traits:
	# trait gene set enrichment results file
	trait_gsea_results_file = biological_pathway_enrichment_dir + 'biologial_pathway_enrichment_analysis_2_' + independent_trait + '_' + str(pip_thresh) + '.txt'
	gsea_results = load_in_gsea_results(trait_gsea_results_file)
	# Sort gsea results
	gsea_results.sort(key=lambda a: a[0])
	#print(len(gsea_results))

	fdrs = [.1]
	for fdr_thresh in fdrs:
		fdr_sig_gsea_results = get_fdr_significant_results(gsea_results, fdr_thresh)
		if len(fdr_sig_gsea_results) > 0:
			fdr_sig_gsea_results.sort(key=lambda a: a[1],reverse=True)
			print(independent_trait)
			print(fdr_thresh)
			print(fdr_sig_gsea_results)


