#!/usr/bin/env python

# split into two diff sizes, then secondary contact with instantaneous size change followed by growth, with symmetric migration
# n(para): 11

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
pop_ids=[sys.argv[2],sys.argv[3]]
projections=[int(sys.argv[4]),int(sys.argv[5])]
if len(sys.argv)==9:
  params = np.loadtxt(sys.argv[8], delimiter=" ", unpack=False)
else:
params=array([1,1,1,1,1,1,1,1,1,1,0.01])

# mutation rate per sequenced portion of genome per generation: for A.millepora, 0.02
mu=float(sys.argv[6])
# generation time, in thousand years: 0.005  (5 years)
gtime=float(sys.argv[7]) 

dd = Misc.make_data_dict(infile)
data = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=True)
ns=data.sample_sizes
np.set_printoptions(precision=3)     


#-------------------
# split with growth and asymmetrical migration; with genomic islands
def IM2iSC(params, ns):
    """
    Isolation-with-migration model with split into two arbtrary sizes
    p_misid: proportion of misidentified ancestral states
    
    """
    nua,nub,nu1_0,nu2_0,nu1,nu2,T,T0,m12,m21,p_misid = params
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/T)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T)
    nu_func = lambda t: [nu1_func(t), nu2_func(t)]

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nua, nub], T0, m = np.array([[0, 0], [0, 0]]))
    fs.integrate(nu_func, T, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))

    return (1-p_misid)*fs + p_misid*moments.Numerics.reverse_array(fs)

func=IM2iSC
upper_bound = [100,100,100,100,100, 100, 10,10, 200,200,0.25]
lower_bound = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3, 1e-3,1e-3,1e-5,1e-5,1e-5]
params = moments.Misc.perturb_params(params, fold=2, upper_bound=upper_bound,
                              lower_bound=lower_bound)

# fitting (poptg = optimal parameters):
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
moments.ModelPlot.plot_model(plot_mod, save_file="IMsc2_"+ind+"_"+sys.argv[1]+".png",pop_labels=pop_ids, nref=theta/(4*mu), draw_scale=False, gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

# bootstrapping for SDs of params and theta
all_boot=moments.Misc.bootstrap(dd,pop_ids,projections)
uncert=moments.Godambe.GIM_uncert(func,all_boot,poptg,data)

# printing parameters and their SDs
print "RESULT","IMsc2",ind,len(params),ll_model,sys.argv[1],sys.argv[2],sys.argv[3],poptg,theta,uncert
                                    
# plotting quad-panel figure witt AFS, model, residuals:
moments.Plotting.plot_2d_comp_multinom(model, data, vmin=0.1, resid_range=3,
                                    pop_ids =pop_ids)
plt.savefig("IMsc2_"+ind+"_"+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]+"_"+sys.argv[4]+"_"+sys.argv[5]+'.pdf')
