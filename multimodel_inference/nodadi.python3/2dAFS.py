#!/usr/bin/env python

import matplotlib
matplotlib.use('PDF')
import moments
import pylab
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
import sys
infile=sys.argv[1]
pop_ids=[sys.argv[2],sys.argv[3]]
projections=[int(sys.argv[4]),int(sys.argv[5])]

fs = moments.Spectrum.from_file(infile)
data=fs.project(projections)
ns=data.sample_sizes
np.set_printoptions(precision=3)     

moments.Plotting.plot_single_2d_sfs(data, vmin=0.1)
plt.savefig('2dAFS_'+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]+"_"+sys.argv[4]+"_"+sys.argv[5]+'.pdf')
nalleles=data.S()
print( "N alleles: ",nalleles)
