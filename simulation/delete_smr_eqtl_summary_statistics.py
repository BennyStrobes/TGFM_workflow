import numpy as np
import os
import sys
import pdb








eqtl_flist_file = sys.argv[1]
besd_file_stem = sys.argv[2]


# First, delete each esd file in flist
f = open(eqtl_flist_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	# skip header
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Error checking
	if len(data) != 7:
		print('assumption eroror')
		pdb.set_trace()

	# Delete esd file corresponding to this line
	esd_file = data[6]
	os.system('rm ' + esd_file)
f.close()
# Remove flist file
os.system('rm ' + eqtl_flist_file)


# Remove besd files
os.system('rm ' + besd_file_stem + '.besd')
os.system('rm ' + besd_file_stem + '.epi')
os.system('rm ' + besd_file_stem + '.esi')
os.system('rm ' + besd_file_stem + '.summary')
