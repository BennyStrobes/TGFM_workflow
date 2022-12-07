import numpy as np 
import os
import sys
import pdb
import gzip




def filter_anno_file(input_anno_file, output_anno_file):
	f = gzip.open(input_anno_file)
	t = open(output_anno_file,'w')
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		t.write('\t'.join(data[:5]) + '\t' + data[76] + '\n')
	f.close()
	t.close()

def filter_anno_file2(input_anno_file, output_anno_file):
	f = gzip.open(input_anno_file)
	t = open(output_anno_file,'w')
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		t.write('\t'.join(data[:77]) + '\t' + '\t'.join(data[80:]) + '\n')
	f.close()
	t.close()


def filter_ld_score_file(input_ld_score_file, output_ld_score_file):
	f = gzip.open(input_ld_score_file)
	t = open(output_ld_score_file,'w')
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		t.write('\t'.join(data[:4]) + '\t' + data[75] + '\n')
	f.close()
	t.close()

def filter_ld_score_file2(input_ld_score_file, output_ld_score_file):
	f = gzip.open(input_ld_score_file)
	t = open(output_ld_score_file,'w')
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		t.write('\t'.join(data[:76]) + '\t' + '\t'.join(data[79:]) + '\n')
	f.close()
	t.close()


#72
def filter_m_file(input_m_file, output_m_file):
	aa = np.loadtxt(input_m_file)
	#out = aa[:1]
	out = np.hstack((aa[0],aa[72]))
	#np.savetxt(output_m_file,out.astype(str), fmt="%s", delimiter='\t')
	t = open(output_m_file,'w')
	t.write(str(out[0]) + '\t' + str(out[1]) + '\n')
	t.close()

def filter_m_file2(input_m_file, output_m_file):
	aa = np.loadtxt(input_m_file)
	#out = aa[:1]
	out = np.hstack((aa[:73], aa[76:]))
	t = open(output_m_file,'w')
	t.write('\t'.join(out.astype(str)) + '\n')
	t.close()


baseline_ld_dir = sys.argv[1]
output_root = sys.argv[2]



for chrom_num in range(1,23):
	print(chrom_num)

	# Annotation file
	input_anno_file = baseline_ld_dir + 'baselineLD.' + str(chrom_num) + '.annot.gz'
	output_anno_file = output_root + 'baselineld_non_eqtl_cpp_removed.' + str(chrom_num) + '.annot'
	filter_anno_file2(input_anno_file, output_anno_file)

	# LD Score file
	input_ld_score_file = baseline_ld_dir + 'baselineLD.' + str(chrom_num) + '.l2.ldscore.gz'
	output_ld_score_file = output_root + 'baselineld_non_eqtl_cpp_removed.' + str(chrom_num) + '.l2.ldscore'
	filter_ld_score_file2(input_ld_score_file, output_ld_score_file)

	# M file
	input_m_file = baseline_ld_dir + 'baselineLD.' + str(chrom_num) + '.l2.M'
	output_m_file = output_root + 'baselineld_non_eqtl_cpp_removed.' + str(chrom_num) + '.l2.M'
	filter_m_file2(input_m_file, output_m_file)

	input_m_5_50_file = baseline_ld_dir + 'baselineLD.' + str(chrom_num) + '.l2.M_5_50'
	output_m_5_50_file = output_root + 'baselineld_non_eqtl_cpp_removed.' + str(chrom_num) + '.l2.M_5_50'
	filter_m_file2(input_m_5_50_file, output_m_5_50_file)

	'''
	# Annotation file
	input_anno_file = baseline_ld_dir + 'baselineLD.' + str(chrom_num) + '.annot.gz'
	output_anno_file = output_root + 'baseline_intercept_only.' + str(chrom_num) + '.annot'
	filter_anno_file(input_anno_file, output_anno_file)

	# LD Score file
	input_ld_score_file = baseline_ld_dir + 'baselineLD.' + str(chrom_num) + '.l2.ldscore.gz'
	output_ld_score_file = output_root + 'baseline_intercept_only.' + str(chrom_num) + '.l2.ldscore'
	filter_ld_score_file(input_ld_score_file, output_ld_score_file)

	# M file
	input_m_file = baseline_ld_dir + 'baselineLD.' + str(chrom_num) + '.l2.M'
	output_m_file = output_root + 'baseline_intercept_only.' + str(chrom_num) + '.l2.M'
	filter_m_file(input_m_file, output_m_file)

	input_m_5_50_file = baseline_ld_dir + 'baselineLD.' + str(chrom_num) + '.l2.M_5_50'
	output_m_5_50_file = output_root + 'baseline_intercept_only.' + str(chrom_num) + '.l2.M_5_50'
	filter_m_file(input_m_5_50_file, output_m_5_50_file)
	'''
