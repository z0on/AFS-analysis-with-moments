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
fs = moments.Spectrum.from_file(infile)
data=fs.fold()
fns=data.sample_sizes
np.set_printoptions(precision=3)     

#-------------------
# two growth periods without jumps

def g2(params , ns):
    nu1,nu2,T1, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0])
    fs = moments.Spectrum(sts)
    nu_func1 = lambda t: [np.exp(np.log(nu1) * t / T1)]
    fs.integrate(nu_func1, T1)
    nu_func2 = lambda t: [nu1 * np.exp(np.log(nu2/nu1) * t / T2)]
    fs.integrate(nu_func2, T2)
    return (1-p_misid)*fs + p_misid*moments.Numerics.reverse_array(fs)
 #   return fs
 
func=g2
upper_bound = [100,100,10,10]
lower_bound = [1e-3,1e-3, 1e-3,1e-3]
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
ind=str(random.randint(0,999999))

print "g2f",ind, sys.argv[1],' ll: ', ll_model,' p: ', poptg, " t: ",theta

moments.Plotting.plot_1d_comp_multinom(model, data)

plt.savefig("g2f_"+ind+"_"+sys.argv[1]+'.pdf')
