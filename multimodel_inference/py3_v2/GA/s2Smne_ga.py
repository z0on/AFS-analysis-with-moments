#!/usr/bin/env python3

# split, two epochs in each pop, asymmetric migration.
# lower Ne in a fractiion of genome (background selection)
# migration scales with size of source population


# uses genetic algorithm from GADMA for optimization

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
import gadma
infile=sys.argv[1]
pop_ids=[sys.argv[2],sys.argv[3]]
projections=[int(sys.argv[4]),int(sys.argv[5])]

# mutation rate per sequenced portion of genome per generation: for A.millepora, 0.02
mu=float(sys.argv[6])
# generation time, in thousand years: 0.005  (5 years)
gtime=float(sys.argv[7]) 

# set Polarized=False below for folded AFS analysis
fs = moments.Spectrum.from_file(infile)
data=fs.project(projections)
ns=data.sample_sizes
np.set_printoptions(precision=3)     

#-------------------
# split into unequal pop sizes with asymmetrical migration

def sc3ei(params , ns):
#    p_misid: proportion of misidentified ancestral states
# Ps: proportion of sites with lower Ne
# Fs: factor of Ne reduction (0.1 - 0.999)
    nu1_1,nu2_1,nu1_2,nu2_2,T1,T2,m12,m21,Fs,Ps,p_misid = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1_1, nu2_1], T1, m = np.array([[0, nu2_1*m12], [nu1_1*m21, 0]]))
    fs.integrate([nu1_2, nu2_2], T2, m = np.array([[0, nu2_2*m12], [nu1_2*m21, 0]]))

    stss = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fss = moments.Spectrum(stss)
    fss = moments.Manips.split_1D_to_2D(fss, ns[0], ns[1])
    fss.integrate([nu1_1*Fs, nu2_1*Fs], T1, m = np.array([[0, nu2_1*m12*Fs], [nu1_1*m21*Fs, 0]]))
    fss.integrate([nu1_2*Fs, nu2_2*Fs], T2, m = np.array([[0, nu2_2*m12*Fs], [nu1_2*m21*Fs, 0]]))

    fs2=P*fss+(1-P)*fs
    return (1-p_misid)*fs2 + p_misid*moments.Numerics.reverse_array(fs2)
 
func=sc3ei

upper_bound = [100, 100, 100,100,100,100,100,100,0.999,0.999,0.25]
lower_bound = [1e-5,1e-5, 1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-1,1e-3,1e-5]
if len(sys.argv)==9:
     params = np.loadtxt(sys.argv[8], delimiter=" ", unpack=False)
#     params = moments.Misc.perturb_params(params, fold=1.5, upper_bound=upper_bound, lower_bound=lower_bound)
     Xinit=[params]
     nGA=1
else:
     Xinit=None
     nGA=150

par_labels = ('nu1_1','nu2_1','nu1_2','nu2_2','T1','T2','m12','m21','Fs','F_gs','f_misid')

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
                                     
# extracting model predictions, likelihood and theta
model = func(poptg, ns)
ll_model = moments.Inference.ll_multinom(model, data)
theta = moments.Inference.optimal_sfs_scaling(model, data)

# random index for this replicate
ind=str(random.randint(0,999999))

# plotting demographic model
plot_mod = moments.ModelPlot.generate_model(func, poptg, ns)
moments.ModelPlot.plot_model(plot_mod, save_file="s2Simne_"+ind+".png", pop_labels=pop_ids, nref=theta/(4*mu), draw_scale=False, gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

# bootstrapping for SDs of params and theta

# printing parameters and their SDs
print( "RESULT","s2Simne",ind,len(par_labels),ll_model,sys.argv[1],sys.argv[2],sys.argv[3],poptg,theta)
                                    
# plotting quad-panel figure witt AFS, model, residuals:
moments.Plotting.plot_2d_comp_multinom(model, data, vmin=0.1, resid_range=3,
                                    pop_ids =pop_ids)
plt.savefig("s2Simne_"+ind+"_"+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]+"_"+sys.argv[4]+"_"+sys.argv[5]+'.pdf')

