#!/usr/bin/env python

import matplotlib
matplotlib.use('PDF')
import moments
import pylab
import random
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
import sys
infile=sys.argv[1]
#pop_ids=[sys.argv[2],sys.argv[3]]
#projections=[int(sys.argv[4]),int(sys.argv[5])]
params=[float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]),float(sys.argv[6])]

# mutation rate per sequenced portion of genome per generation
mu=0.018
# generation time, in thousand years
gtime=0.005 

fs = moments.Spectrum.from_file(infile)
data = fs.fold()
nalleles=data.S()
print "N alleles: ",nalleles
ns=data.sample_sizes
np.set_printoptions(precision=3)     

#-------------------
# split into unequal pop sizes with asymmetrical migration

def s2m(params , ns):
#    p_misid: proportion of misidentified ancestral states
    nu1, nu2, T, m12, m21 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T, m = np.array([[0, m12], [m21, 0]]))
    return fs
 
func=s2m
upper_bound = [100, 100, 100, 200,200]
lower_bound = [1e-3,1e-3, 1e-3,0.1,0.1]
params = moments.Misc.perturb_params(params, fold=2, upper_bound=upper_bound,
                              lower_bound=lower_bound)

poptg = moments.Inference.optimize_log(params, data, func,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(params), maxiter=30)

# extracting model predictions, likelihood and theta
model = func(poptg, ns)
ll_model = moments.Inference.ll_multinom(model, data)
theta = moments.Inference.optimal_sfs_scaling(model, data)

# random index for this replicate
ind=str(random.randint(0,999999))

# plotting demographic model
plot_mod = moments.ModelPlot.generate_model(func, poptg, ns)
moments.ModelPlot.plot_model(plot_mod, save_file="s2m_"+ind+".png", draw_scale=False, nref=theta/(4*mu), gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

# bootstrapping for SDs of params and theta
#all_boot=moments.Misc.bootstrap(dd,pop_ids,projections)
#uncert=moments.Godambe.GIM_uncert(s2m,all_boot,poptg,data)

# printing parameters and their SDs
print "s2mfRes",ind,sys.argv[1],' ll: ', ll_model,' p: ', poptg, " t: ",theta
moments.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3)
                                    
# plotting quad-panel figure wit AFS, model, residuals:
plt.savefig("s2mf_"+ind+"_"+sys.argv[1]+'.pdf')

