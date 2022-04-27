import numpy as np 
import os
import sys
import pdb
import gzip
import scipy.stats




gtex_pseudotissue_file = sys.argv[1]
output_dir = sys.argv[2]
tiffany_dir = sys.argv[3]



tissue_names_df = np.loadtxt(gtex_pseudotissue_file, dtype=str, delimiter='\t')
tissues = tissue_names_df[1:,0]

output_file = output_dir + 'UKB_460K.blood_WHITE_COUNT_marginal_twas_pvalues.txt'
t = open(output_file,'w')
t.write('tissue\tnominal_pvalue\n')

for tissue in tissues:
	print(tissue)
	tiffany_file = tiffany_dir + 'Marginal_alphas_UKB_460K.blood_WHITE_COUNT_' + tissue + '.txt.gz'
	if os.path.exists(tiffany_file) == False:
		pdb.set_trace()
	f = gzip.open(tiffany_file)
	for line in f:
		line = line.rstrip().decode('utf-8')
		if line == 'NA':
			continue
		data = line.split('\t')
		zscore = float(line)
		pvalue = scipy.stats.norm.sf(abs(zscore))
		t.write(tissue + '\t' + str(pvalue) + '\n')
	f.close()
t.close()