#!/usr/bin/env python
import matplotlib
matplotlib.use('PDF')
import moments
import random
import pylab
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
import sys
infile=sys.argv[1]
params=array([  float(sys.argv[2]),   float(sys.argv[3]),  float(sys.argv[4]),   float(sys.argv[5])])

import os
data = moments.Spectrum.from_file(infile)
nalleles=data.S()
print "N alleles: ",nalleles
ns=data.sample_sizes
np.set_printoptions(precision=3)     

#-------------------
# one growth period with jumps in between

def gj1(params , ns):
    nu01,nu1,T1, p_misid = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0])
    fs = moments.Spectrum(sts)
    nu_func1 = lambda t: [nu01 * np.exp(np.log(nu1/nu01) * t / T1)]
    fs.integrate(nu_func1, T1)
    return (1-p_misid)*fs + p_misid*moments.Numerics.reverse_array(fs)
#    return fs
 
func=gj1
upper_bound = [100,100,10,0.3]
lower_bound = [1e-5,1e-5, 1e-5,1e-5]
params = moments.Misc.perturb_params(params, fold=1, upper_bound=upper_bound,
						  lower_bound=lower_bound)

poptg = moments.Inference.optimize_log(params, data, func,
							   lower_bound=lower_bound,
							   upper_bound=upper_bound,
							   verbose=len(params), maxiter=20)
model = func(poptg, ns)
ll_model = moments.Inference.ll_multinom(model, data)
theta = moments.Inference.optimal_sfs_scaling(model, data)

# index for this replicate    
ind=str(random.randint(0,99999))

print "gj1",ind,sys.argv[1],' ll: ', ll_model,' p: ', poptg, " t: ",theta

moments.Plotting.plot_1d_comp_multinom(model, data)

plt.savefig("gj1_"+ind+"_"+sys.argv[1]+'.pdf')
