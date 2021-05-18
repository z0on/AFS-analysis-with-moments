#!/usr/bin/env python3

# split into two different sizes with growth, aymmetric migration scales with pop size
# genomic islands (lower migration)
# may take a LONG TIME to fit (3.5 hours)


import matplotlib
matplotlib.use('PDF')
import moments
import random
import pylab
import gadma
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
    params=array([ 1,1,1,1,1,1,1,0.01])

# mutation rate per sequenced portion of genome per generation: for A.millepora, 0.02
mu=float(sys.argv[6])
# generation time, in thousand years: 0.005  (5 years)
gtime=float(sys.argv[7]) 

#infile="5kA3_dadi.data"
#pop_ids=["W","K"]
#projections=[32,38]

fs = moments.Spectrum.from_file(infile)
fs=fs.fold()
data=fs.project(projections)
ns=data.sample_sizes
np.set_printoptions(precision=3)     


#-------------------
# split with growth and asymmetrical migration; with genomic islands
def IMi(params, ns):
    """
    Isolation-with-migration model with split into two arbtrary sizes
    Migration is asymmetric and scales with the size of the source population
    p_misid: proportion of misidentified ancestral states
    
    """
    nu1_0,nu2_0,nu1,nu2,T,m12,m21,Fi,P = params
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/T)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T)
    nu_func = lambda t: [nu1_func(t), nu2_func(t)]

    m21_func = lambda t: m21 * nu2_func(t)
    m12_func = lambda t: m12 * nu1_func(t)
    migs = lambda t: np.array([[0, m12_func(t)], [m21_func(t), 0]])
    migs.i = lambda t: np.array([[0, Fi*m12_func(t)], [Fi*m21_func(t), 0]])
 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, T, dt_fac=0.01, m=migs)
 
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate(nu_func, T, dt_fac=0.01, m=migs.i)

    fs2=P*fsi+(1-P)*fs
    return fs2

func=IMi
upper_bound = [100,100,100, 100, 10, 100,100,0.999,0.999]
lower_bound = [1e-3,1e-3,1e-3,1e-3, 1e-3,1e-5,1e-5,1e-5,1e-5]

# if starting parameterss are supplied, don't run GA; if not, run 150 generations of GA
if len(sys.argv)==9:
     params = np.loadtxt(sys.argv[8], delimiter=" ", unpack=False)
#     params = moments.Misc.perturb_params(params, fold=1.5, upper_bound=upper_bound, lower_bound=lower_bound)
     Xinit=[params]
     nGA=1
else:
     Xinit=None
     nGA=150

par_labels = ('nu1_0','nu2_0','nu1','nu2','T','m12','m21','F_isl','F_gen')

import timeit
# allowed fold-excess in evaluation time
X = 20
num_init=50
variables = gadma.cli.get_variables(par_labels, lower_bound, upper_bound)
X_init = [[var.resample() for var in variables] for _ in range(num_init)]
Y_init = []
def f():
    for x in X_init:
        Y_init.append(func(x, ns))
total_time = timeit.timeit(f, number=1)
mean_time = total_time / num_init

result = gadma.Inference.optimize_ga(data=data,
                                     model_func=func,
                                     verbose=0,
                                     X_init=Xinit,
                                     engine='moments',
                                     args=(),
                                     p_ids = par_labels,
                                     maxtime_per_eval = mean_time * X,
                                     lower_bound=lower_bound,
                                     upper_bound=upper_bound,
                                     local_optimizer='BFGS_log',
                                     ga_maxiter=nGA,
                                     ls_maxiter=1)
poptg=result.x  

# extracting model predictions, likelihood and theta
model = func(poptg, ns)
ll_model = moments.Inference.ll_multinom(model, data)
theta = moments.Inference.optimal_sfs_scaling(model, data)

# random index for this replicate
ind=str(random.randint(0,999999))

# plotting demographic model
plot_mod = moments.ModelPlot.generate_model(func, poptg, ns)
moments.ModelPlot.plot_model(plot_mod, save_file="IMimne_"+ind+".png", pop_labels=pop_ids, nref=theta/(4*mu), draw_scale=False, gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

# bootstrapping for SDs of params and theta

# printing parameters and their SDs
print( "RESULT","IMimne",ind,len(par_labels),ll_model,sys.argv[1],sys.argv[2],sys.argv[3],poptg,theta)
                                    
# plotting quad-panel figure witt AFS, model, residuals:
moments.Plotting.plot_2d_comp_multinom(model, data, vmin=0.1, resid_range=3,
                                    pop_ids =pop_ids)
plt.savefig("IMimne_"+ind+"_"+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]+"_"+sys.argv[4]+"_"+sys.argv[5]+'.pdf')

