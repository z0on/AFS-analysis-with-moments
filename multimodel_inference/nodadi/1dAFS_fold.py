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
pop_ids=[sys.argv[2]]
projections=[int(sys.argv[3])]

fs = moments.Spectrum.from_file(infile)
fs=fs.fold()
data=fs.project(projections)
ns=data.sample_sizes
np.set_printoptions(precision=3)     

moments.Plotting.plot_1d_fs(data)
plt.savefig('1dAFSf_'+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]+'.pdf')
nalleles=data.S()
print "N alleles: ",nalleles
np.savetxt('1dsfss', data[1:-1], fmt='%.2f',newline=" ")

