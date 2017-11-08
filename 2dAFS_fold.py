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

dd = Misc.make_data_dict(infile)
data = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=False)
ns=data.sample_sizes

'''
# masking singletons 
#data.mask[:,0]=True
#data.mask[0,:]=True
data.mask[:,1]=True
data.mask[1,:]=True
#data.mask[:,projections[1]]=True
#data.mask[projections[0],:]=True
data.mask[:,projections[1]-1]=True
data.mask[projections[0]-1,:]=True
'''
np.set_printoptions(precision=3)     

moments.Plotting.plot_single_2d_sfs(data, vmin=1)
plt.savefig('2dAFSf_'+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]+"_"+sys.argv[4]+"_"+sys.argv[5]+'.pdf')
