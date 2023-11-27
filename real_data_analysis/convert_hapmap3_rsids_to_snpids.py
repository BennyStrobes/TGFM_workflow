import numpy as np 
import os
import sys
import pdb




def get_rs_ids(rsid_file):
	arr = []
	f = open(rsid_file)
	for line in f:
		line = line.rstrip()
		arr.append(line)
	f.close()
	return np.asarray(arr)

def create_rsid_to_snp_id_mapping(ukbb_hg38_genotype_dir):
	mapping = {}
	for chrom_num in range(1,23):
		chrom_bim_file = ukbb_hg38_genotype_dir + '1000G.EUR.hg38.' + str(chrom_num) + '.bim'
		f = open(chrom_bim_file)
		for line in f:
			line = line.rstrip()
			data = line.split('\t')

			rsid = data[1]
			snp_id = 'chr' + data[0] + '_' + data[3] + '_' + data[4] + '_' + data[5]

			if rsid in mapping:
				if snp_id != mapping[rsid]:
					print('assumption eroror')
					pdb.set_trace()
			mapping[rsid] = snp_id


		f.close()

	return mapping





rsid_file = sys.argv[1]
ukbb_hg38_genotype_dir = sys.argv[2]
snpid_file = sys.argv[3]

rs_ids = get_rs_ids(rsid_file)

rsid_to_snpid_mapping = create_rsid_to_snp_id_mapping(ukbb_hg38_genotype_dir)

t = open(snpid_file,'w')

for rs_id in rs_ids:
	if rs_id not in rsid_to_snpid_mapping:
		print('missing')
		continue

	snp_id = rsid_to_snpid_mapping[rs_id]
	snp_info = snp_id.split('_')


	snp_id1 = snp_info[0] + '_' + snp_info[1] + '_' + snp_info[2] + '_' + snp_info[3]
	snp_id2 = snp_info[0] + '_' + snp_info[1] + '_' + snp_info[3] + '_' + snp_info[2]

	t.write(snp_id1 + '\n')
	t.write(snp_id2 + '\n')



t.close()