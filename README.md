# Multimodel demographic inference using *moments* 

[*Moments*](https://bitbucket.org/simongravel/moments/src/master/) is AFS analysis method superficially similar to *dadi* but operating on different math (ordinary differential equations rather than diffusion approximation) It is considerably faster than *dadi*, handles up to five populations simultaneously, has easy functions for bootstrapping and estimating parameter uncertainties, and plots cartoons of inferred demographic scenarios.

The problem with *moments* and *dadi* is that they are designed to evaluate a pre-specified demographic model; they cannot find the model structure that best fits the data. The idea of this repository is pretty simple: to fit all possible two-population models we can possible think of to the experimental 2dSFS and  use Akaike Information Criterion to select the best one.

See [GADMA](https://github.com/ctlab/GADMA) for the alternative solution to this problem. Compared to GADMA, we are far less elegant but somewhat more flexible (we can incorporate essentially any model, for example, involving selection and heterogeneous introgression rates across the genome). Our approach also lets the user evaluate how much better the winning model is compared to "null" alternatives (for example, models with no population split, with symmetrical migration, or with constant population sizes). Our disadvantage is a ridiculously huge number of model-fit runs that we have to perform; the good news is, all this can be done in parallel, and each single-model run takes no more than 2 hours.

## Stages of the method ##
- The first step is **model selection**, where we run all possible models on 10 bootstrapped SFS. We actually run each model on each bootstrap six times, to make sure the model converges to its best likelihood at least once. Then we use an `R` script `model_selection_10boots.R` to select the best-fitted instance (out of 6) for each model for each bootstrap, and compare the AIC scores for all models. The best model is the one with the *lowest median AIC score among bootstrap replicates*.  
- The second step is running the winning model on 100 bootstrapped SFS, to **evaluate parameter uncertainties**. Once again, we will have to do 6 model runs for each bootstrap. The parameter meanings and uncertainties are deciphered by the second R script that we have, `bestmodel_bootstrap.R`.

## Model selection ##
Let's assume we have a 2dSFS formatted for *moments* or *dadi*, which is no more than a string of numbers with a header line giving two numbers, the dimensions of the spectrum ( 2 x N + 1 for each of the two populations, where N is the number of sampled diploids). See **Appendix** for instructions how to obtain bootstrapped 2dSFS from [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD)

See **moments_scripts_README.txt** for instructions.
## Appendix ## 

[Link to original *Moments* paper]( http://www.genetics.org/content/early/2017/05/08/genetics.117.200493)

[*Moments* manual](https://bitbucket.org/simongravel/moments/raw/efc4da3047226e3662dd43b525e41c85b93e90fd/doc/manual/manual.pdf) (link may change with updates)
