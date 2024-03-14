import numpy as np 
import os
import sys
import pdb






filer1 = sys.argv[1]
filer2 = sys.argv[2]



aa = np.loadtxt(filer1, dtype=str)
bb = np.loadtxt(filer2, dtype=str)

if np.array_equal(aa[:,1], bb[:,1]) == False:
	print('assumption erororor!')
	pdb.set_trace()
