#!/usr/bin/env python

# growth before split, two epochs in each pop, symmetric migration in middle epoch
# genomic islands
# n(para): 12


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
pop_ids=[sys.argv[2],sys.argv[3]]
projections=[int(sys.argv[4]),int(sys.argv[5])]
if len(sys.argv)==9:
    params = np.loadtxt(sys.argv[8], delimiter=" ", unpack=False)
else:
    params=[1,1,1,1,1,1,1,1,1,1,0.5]

# mutation rate per sequenced portion of genome per generation: for A.millepora, 0.02
mu=float(sys.argv[6])
# generation time, in thousand years: 0.005  (5 years)
gtime=float(sys.argv[7]) 

# set Polarized=False below for folded AFS analysis
fs = moments.Spectrum.from_file(infile)
fs=fs.fold()
data=fs.project(projections)
ns=data.sample_sizes
np.set_printoptions(precision=3)     

#-------------------
# split into unequal pop sizes with asymmetrical migration

def sc12imsm(params , ns):
#    p_misid: proportion of misidentified ancestral states
    nu0, nu1_2,nu2_2,nu1_3,nu2_3,T1, T2, T3,m,mi, P = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu0], T1)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1_2, nu2_2], T2, m = np.array([[0, m], [m, 0]]))
    fs.integrate([nu1_3, nu2_3], T3, m = np.array([[0, 0], [0, 0]]))

    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi.integrate([nu0], T1)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate([nu1_2, nu2_2], T2, m = np.array([[0, mi], [mi, 0]]))
    fsi.integrate([nu1_3, nu2_3], T3, m = np.array([[0, 0], [0, 0]]))

    fs2=P*fsi+(1-P)*fs
    return fs2
 
func=sc12imsm
upper_bound = [100, 100,100,100, 100, 100, 100,100, 200,200,0.999]
lower_bound = [1e-3, 1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,0.001]
params = moments.Misc.perturb_params(params, fold=2, upper_bound=upper_bound,
                              lower_bound=lower_bound)

poptg = moments.Inference.optimize_log(params, data, func,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=False, maxiter=30)

# extracting model predictions, likelihood and theta
model = func(poptg, ns)
ll_model = moments.Inference.ll_multinom(model, data)
theta = moments.Inference.optimal_sfs_scaling(model, data)

# random index for this replicate
ind=str(random.randint(0,999999))

# plotting demographic model
plot_mod = moments.ModelPlot.generate_model(func, poptg, ns)
moments.ModelPlot.plot_model(plot_mod, save_file="sc12imsm_"+ind+".png", pop_labels=pop_ids, nref=theta/(4*mu), draw_scale=False, gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

# bootstrapping for SDs of params and theta

# printing parameters and their SDs
print( "RESULT","sc12imsm",ind,len(params),ll_model,sys.argv[1],sys.argv[2],sys.argv[3],poptg,theta)
                                    
# plotting quad-panel figure witt AFS, model, residuals:
moments.Plotting.plot_2d_comp_multinom(model, data, vmin=0.1, resid_range=3,
                                    pop_ids =pop_ids)
plt.savefig("sc12imsm_"+ind+"_"+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]+"_"+sys.argv[4]+"_"+sys.argv[5]+'.pdf')

