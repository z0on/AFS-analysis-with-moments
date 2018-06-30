Moments is AFS analysis method superficially similar to dadi but operating
on different math (ordinary differential equations rather than diffusion approximation)
It is considerably faster than dadi, handles up to five populations simultaneously,
has easy functions for bootstrapping and estimating parameter uncertainties,
and plots cartoons of inferred demographic scenarios.

This collection of command-line scrips is to run basic pairwise analysis of 
population subdivision using moments, without the need to recode python scripts
(the only thing that needs to be hard-coded in the scripts are settings for 
mutation rate and generation time for your species, see moments_scripts_README.txt)

Link to original paper: http://www.genetics.org/content/early/2017/05/08/genetics.117.200493
Moments manual (may change with updates): https://bitbucket.org/simongravel/moments/raw/efc4da3047226e3662dd43b525e41c85b93e90fd/doc/manual/manual.pdf

#------------
# installing moments

# visit https://bitbucket.org/simongravel/moments

# On a Mac: in Enthought Canopy terminal, say
git clone https://bitbucket.org/simongravel/moments.git
cd moments
sudo python setup.py install 

# installing locally on a Linux cluster:
git clone https://bitbucket.org/simongravel/moments.git
cd moments
python setup.py build_ext --inplace

#------------------------
# PREPARING DATA

# unpacking vcf file to play with
tar vxf coral.tgz

# converting vcf into dadi/moments format
# pops.txt is a two-column table giving correspondence 
# between samples (in vcf file) and populations
perl vcf2dadi.pl coral.vcf pops.txt 

# calculating optimal projections: number of AFS bins giving the most segregating sites
# (aim for slightly lower than 2x number of individuals)
# arguments: datafile, population, lowest projection to try, highest projection to try
# the last pair of numbers printed to STDOUT is the optimal projection number 
# and the number of segregating sites under this projection
python projections.py coral_dadi.data K 30 40  
python projections.py coral_dadi.data W 30 40  

#-----------------------
#  DATA EXPLORATION

# plotting 1d AFS for K population
# arguments: datafile, population, projections
python 1dAFS.py coral_dadi.data K 36

# plotting 2d AFS for W and K populations
# arguments: datafile, pop1, pop2, projections for pop1, projections for pop2
python 2dAFS.py coral_dadi.data W K 32 36

# plotting folded 2dAFS
python 2dAFS_fold.py coral_dadi.data W K 32 36

#==================
#  FITTING MODELS

# main models (S2M, IM2 and their folded or "hot" versions) must be run
# multiple times on the same dataset to find parameters producing 
# best likelihood. Initial parameter setting can be the same; they will be perturbed
# when starting the model (more strongly perturbed in "hot" version)
# Outputs on main model scripts:
#	- resulting parameter values and their uncertainty SDs based on bootstrap, printed to STDOUT 
#     (uncertainties are for all parameters and theta);
#   - pdf quad-plot of data - model - residuals
#		(examine them for correlated residuals in the AFS - 
#		indication that the model does not capture something)
#   - png of inferred demographic scenario 
#		to properly scale this graph, modify mutation rate and generation time within 
#		the model script, lines ~18-22. Currently set to:
#		# mutation rate per sequenced portion of genome per generation
#		mu=0.018
#		# generation time, in thousand years
#		gtime=0.005 

# Each run is assigned a random ID number to associate printed parameter values 
# with quad-plots and demographic cartoons

# If your reads are mapped to an outgroup genome, which allows identifying which SNP 
# state is ancestral, use basic (polarized) versions of models.
# if not (mapped to same-species genome or called genotypes de novo), use "fold" versions. 
# Example commands below assume polarized data, using basic models.
# The folded models are run in the same way but have one less argument 
# (percent of ancestral state mis-identification, 
# which is the last argument in basic models - remove it when running folded models)

#-----------------------
# split with asymmetrical migration (S2M) model

# arguments: 
#	datafile, pop1, pop2, projections for pop1, projections for pop2
#	pop1 size, pop2 size, 
#   time of split,
#	migration from pop2 to pop1 
# 		(number of individuals in pop1 that are immigrants from pop2 every generation),
#	migration from pop1 to pop2,
#	percent of misidentified ancestral states
python S2M.py coral_dadi.data W K 32 36 1 1 1 10 10 0.01

# "S2M" without split
# arguments: 
#	datafile, pop1, pop2, projections for pop1, projections for pop2
#	final pop size 
#   time of size change,
#	percent of misidentified ancestral states
python te2d_S2Mnull.py coral_dadi.data W K 32 36 1 1 0.01


# S2M without migration (can the AFS be explained by very recent split instead?)
# arguments: 
#	datafile, pop1, pop2, projections for pop1, projections for pop2
#	pop1 size, pop2 size, 
#   time of split,
#	percent of misidentified ancestral states
python S2M_noMig.py coral_dadi.data W K 32 36 1 1 1 0.01


# S2M with symmetrical migration
# arguments: 
#	datafile, pop1, pop2, projections for pop1, projections for pop2
#	pop1 size, pop2 size, 
#   time of split,
#   migration (same in both directions),
#	percent of misidentified ancestral states
python S2M_symmMig.py coral_dadi.data W K 32 36 1 1 1 10 0.01

# use AICweights.R to compare models by delta-AIC criterion

#-------------
# IM2 model - like S2M but with growth

# split with asymmetrical migration and growth (IM2) model
# (this one might take a long time!)
# arguments: 
#	datafile, pop1, pop2, projections for pop1, projections for pop2
#	initial pop1 size, initial pop2 size, 
#	final pop1 size, final pop2 size, 
#   time of split,
#	migration from pop2 to pop1,
#	migration from pop1 to pop2,
#	percent of misidentified ancestral states
python IM2.py coral_dadi.data W K 32 36 1 1 1 1 5 10 10 0.01

# IM2 without split
# arguments: 
#	datafile, pop1, pop2, projections for pop1, projections for pop2
#	initial pop size
#	final pop size 
#   time of growth start (= initial size change)
#	percent of misidentified ancestral states
python growth2d_IM2null.py coral_dadi.data W K 32 36 1 1 1 0.01

# IM2 without migration (can the AFS be explained by very recent split instead?)
# arguments: 
#	datafile, pop1, pop2, projections for pop1, projections for pop2
#	initial pop1 size, initial pop2 size, 
#	final pop1 size, final pop2 size, 
#   time of split,
#	percent of misidentified ancestral states
python IM2_noMig.py coral_dadi.data W K 32 36 1 1 1 1 5 0.01

# IM2 with symmetrical migration
# arguments: 
#	datafile, pop1, pop2, projections for pop1, projections for pop2
#	initial pop1 size, initial pop2 size, 
#	final pop1 size, final pop2 size, 
#   time of split,
#   migration (same in both directions),
#	percent of misidentified ancestral states
python IM2_symmMig.py coral_dadi.data W K 32 36 1 1 1 1 1 10 0.01

# use AICweights.R to compare models by delta-AIC criterion



