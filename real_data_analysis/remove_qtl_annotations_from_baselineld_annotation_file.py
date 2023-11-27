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
		t.write('\t'.join(data[:76]) + '\t' + '\t'.join(data[80:]) + '\n')
	f.close()
	t.close()


def filter_m_file(input_m_file, output_m_file):
	aa = np.loadtxt(input_m_file)
	#out = aa[:1]
	out = np.hstack((aa[:72], aa[76:]))
	t = open(output_m_file,'w')
	t.write('\t'.join(out.astype(str)) + '\n')
	t.close()



input_stem = sys.argv[1]
output_stem = sys.argv[2]


# Annotation file
input_anno_file = input_stem + '.annot.gz'
output_anno_file = output_stem + '.annot'
filter_anno_file(input_anno_file, output_anno_file)

