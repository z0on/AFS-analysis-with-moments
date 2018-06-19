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
#pop_ids=[sys.argv[2],sys.argv[3]]
projections=[int(sys.argv[2]),int(sys.argv[3])]

#dd = Misc.make_data_dict(infile)
#data = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=True)

fs = moments.Spectrum.from_file(infile)
fs=fs.project(projections)
#fs = data.fold()
nalleles=fs.S()
print "N alleles: ",nalleles
#ns=fs.sample_sizes
#np.set_printoptions(precision=3)     
#pylab.figure()
moments.Plotting.plot_single_2d_sfs(fs, vmin=1)
#plt.show()
plt.savefig('2dAFS'+"_"+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]+'.pdf')
