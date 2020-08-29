# Multimodel demographic inference using *moments* 

[*Moments*](https://bitbucket.org/simongravel/moments/src/master/) is the site frequency spectrum (SFS) analysis method superficially similar to *dadi* but operating on different math (ordinary differential equations rather than diffusion approximation). It is considerably faster than *dadi* (although somewhat less accurate), handles up to five populations simultaneously, and plots cartoons of inferred demographic scenarios.

The problem with *moments* and *dadi* is, while they can evaluate the fit of a pre-specified demographic model, they are not designed to search for the general model structure that best fits the data (i.e., was there population splitting or not, with or without migration, with or without additional growth periods). The idea of this repository to solving this problem is pretty simple: to fit all two-population models we can possibly think of to the experimental 2dSFS and use Akaike Information Criterion to select the best one.

See [GADMA](https://github.com/ctlab/GADMA) for the alternative solution to this problem. Compared to GADMA, we are far less elegant but somewhat more flexible (we can incorporate essentially any model, for example, involving selection and heterogeneous introgression rates across the genome). Our approach also lets the user evaluate how much better the winning model is compared to "null" alternatives (for example, models with no population split, with symmetrical migration, or with constant population sizes). Our disadvantage is a ridiculously huge number of model-fit runs that we have to perform; the good news is, all this can be done in parallel, and each single-model run takes no more than 2 hours.

## Installation ##
First of all, install *moments*. The example below would clone it into root directory and install it for a specific user.
```bash
cd
git clone https://bitbucket.org/simongravel/moments.git 
cd moments
python setup.py build_ext --inplace
cd
```
Then, clone this repository and copy all the `*.py` files from `~/AFS-analysis-with-moments/multimodel_inference/nodadi.python2/` or from `/AFS-analysis-with-moments/multimodel_inference/nodadi.python3/` (depending on your `python` version) to where you keep your executables (for example, `~/bin`). 
> NOTE: all code examples here assume the repository is cloned in the home directory, `~/`. If you cloned it elsewhere, make sure to replace `~/` in all examples with the actual path.

## Overview of the method ##
- The first step is **model selection**, where we run all possible models on 10 bootstrapped SFS. We actually run each model on each bootstrap six times, to make sure the model converges to its best likelihood at least once. Then we use an `R` script `model_selection_10boots.R` to select the best-fitted instance (out of 6) for each model for each bootstrap, and compare the AIC scores for all models. The best model is the one with the *lowest median AIC score among bootstrap replicates*.  
- The second step is running the winning model on 100 bootstrapped SFS, to **evaluate parameter uncertainties**. Once again, we will have to do 6 model runs for each bootstrap. The parameter meanings and uncertainties are deciphered by the second R script that we have, `bestmodel_bootstrap.R`.

## Model selection ##
Let's assume we have ten bootstrapped 2dSFS formatted for *moments* or *dadi*. Such file is nothing more than a line of numbers with a header line giving the dimensions of the spectrum ( 2 x N + 1 for each of the two populations, where N is the number of sampled diploids). See **Appendix** for instructions how to obtain bootstrapped 2dSFS from [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD).
Bootstrapped sfs shoudl be named like `p12_1.sfs`, `p12_2.sfs` etc. where `p12` is the name of population contrast.

```bash
cd [where your boostrapped SFS files are]
cp ~/AFS-analysis-with-moments/multimodel_inference/allmodels_unfolded allmodels
# if your SFS needs to be folded, use this line instead:
# cp ~/AFS-analysis-with-moments/multimodel_inference/allmodels_folded allmodels
NREPS=6
>mods
for i in `seq 1 $NREPS`;do 
cat allmodels >>mods;
done

CONTRAST=p12 # name of population comparison, should be matching the leading part of the bootstapped SFS names
ARGS="p1 p2 16 16 0.02 0.005" # pop1, pop2, projection for pop1, projection for pop2, mutation rate (per genotyped portion of the genome per generation), generation time in thousands of years. Population names can be anything. For ANGSD-derived SFS, projections should be 0.8*2N for each population (ronded to integer); in the case shown here, each population was represented by 10 individuals.

>modsel
for B in `seq 1 10`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS" >>args;
done;
paste mods args -d " " >>modsel;
done
```

See **moments_scripts_README.txt** for instructions.
## Appendix ## 

[Link to original *Moments* paper]( http://www.genetics.org/content/early/2017/05/08/genetics.117.200493)

[*Moments* manual](https://bitbucket.org/simongravel/moments/raw/efc4da3047226e3662dd43b525e41c85b93e90fd/doc/manual/manual.pdf) (link may change with updates)
