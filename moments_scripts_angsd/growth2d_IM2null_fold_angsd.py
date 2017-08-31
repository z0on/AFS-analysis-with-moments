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
#pop_ids=[sys.argv[2],sys.argv[3]]
#projections=[int(sys.argv[4]),int(sys.argv[5])]
params=[float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4])]

import os
fs = moments.Spectrum.from_file(infile)
data = fs.fold()
nalleles=data.S()
print "N alleles: ",nalleles
#dd = Misc.make_data_dict(infile)
#data = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=False)
ns=data.sample_sizes
np.set_printoptions(precision=3)     

#-------------------
# 2d growth model: null for IM2

def growth2d(params, ns):
    """
    size change then exponential growth, then sample twice

    params = (nu0,nu,T,p_misid)

    nu0: initial size after change
    nu: Ratio of contemporary to initial population size
    T: Time in the past at which growth began (in units of 2*Na 
       generations) 
    p_misid: proportion of misidentified ancestral states
    """
    nu0, nu, T = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    nu_func = lambda t: [nu0*np.exp(np.log(nu/nu0) * t / T)]
    fs.integrate(nu_func, T, 0.01)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    return fs
 
func=growth2d
upper_bound = [100,100,10]
lower_bound = [1e-5,1e-5,1e-5]
params = moments.Misc.perturb_params(params, fold=1, upper_bound=upper_bound,
						  lower_bound=lower_bound)

poptg = moments.Inference.optimize_log(params, data, func,
							   lower_bound=lower_bound,
							   upper_bound=upper_bound,
							   verbose=False,maxiter=10)
# extracting model predictions, likelihood and theta
model = func(poptg, ns)
ll_model = moments.Inference.ll_multinom(model, data)
theta = moments.Inference.optimal_sfs_scaling(model, data)

# random index for this replicate
ind=str(random.randint(0,999999))

# printing parameters 
print "gr2dResF",ind,sys.argv[1],' ll: ', ll_model,' p: ', poptg, " t: ",theta, 
moments.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3)
                                                            
#plotting quad-panel figure wit AFS, model, residuals:
plt.savefig("gr2df"+ind+"_"+sys.argv[1]+'.pdf')
