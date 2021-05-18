#!/usr/bin/env python3

# split, two epochs in each pop, asymmetric migration only in the second epoch ("secondary contact").
# genomic islands (lower migration)


# does NOT use GADMA

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
#import gadma
infile=sys.argv[1]
pop_ids=[sys.argv[2],sys.argv[3]]
projections=[int(sys.argv[4]),int(sys.argv[5])]

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

def sc3ei(params , ns):
#    p_misid: proportion of misidentified ancestral states
# P: proportion of sites with lower migration
# Fi: factor of migration reduction (1e-5 - 0.999)
# migration is asymmetric, scales with pop size
    nu1_1,nu2_1,nu1_2,nu2_2,T1,T2,m12,m21,Fi,P = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1_1, nu2_1], T1, m = np.array([[0, 0], [0, 0]]))
    fs.integrate([nu1_2, nu2_2], T2, m = np.array([[0, m12], [m21, 0]]))

    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate([nu1_1, nu2_1], T1, m = np.array([[0, 0], [0, 0]]))
    fsi.integrate([nu1_2, nu2_2], T2, m = np.array([[0, m12*Fi], [m21*Fi, 0]]))

    fs2=P*fsi+(1-P)*fs
    return fs2
 
func=sc3ei

upper_bound = [100, 100, 100,100,100,100,200,200,0.999,0.999]
lower_bound = [1e-5,1e-5, 1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5]

'''
if len(sys.argv)==9:
     params = np.loadtxt(sys.argv[8], delimiter=" ", unpack=False)
#     params = moments.Misc.perturb_params(params, fold=1.5, upper_bound=upper_bound, lower_bound=lower_bound)
     Xinit=[params]
     nGA=1
else:
     Xinit=None
     nGA=150

par_labels = ('nu1_1','nu2_1','nu1_2','nu2_2','T1','T2','m12','m21','F_isl','F_gen')

# calculating time limit for GADMA evaluations (the generation will re-spawn if stuck for longer than that)

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
'''
if len(sys.argv)==9:
    params = np.loadtxt(sys.argv[8], delimiter=" ", unpack=False)
else:
    params=[[1]*(len(upper_bound)-1),0.01]

params = moments.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

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
moments.ModelPlot.plot_model(plot_mod, save_file="sci_"+ind+".png", pop_labels=pop_ids, nref=theta/(4*mu), draw_scale=False, gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

# bootstrapping for SDs of params and theta

# printing parameters and their SDs
print( "RESULT","sci",ind,len(par_labels),ll_model,sys.argv[1],sys.argv[2],sys.argv[3],poptg,theta)
                                    
# plotting quad-panel figure witt AFS, model, residuals:
moments.Plotting.plot_2d_comp_multinom(model, data, vmin=0.1, resid_range=3,
                                    pop_ids =pop_ids)
plt.savefig("sci_"+ind+"_"+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]+"_"+sys.argv[4]+"_"+sys.argv[5]+'.pdf')

