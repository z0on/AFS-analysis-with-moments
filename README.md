# Multimodel demographic inference using *moments* 

[*Moments*](https://bitbucket.org/simongravel/moments/src/master/) is the site frequency spectrum (SFS) analysis method superficially similar to *dadi* but operating on different math (ordinary differential equations rather than diffusion approximation). It is considerably faster than *dadi* (although somewhat less accurate), handles up to five populations simultaneously, and plots cartoons of inferred demographic scenarios.

The problem with *moments* and *dadi* is, they can evaluate the fit of a pre-specified demographic model but are **not designed to search for the general model structure** that best fits the data (i.e., was there population split or not, is there any migration and if yes, is it symmetrical or not, were there additional growth periods before or after spit, etc). The idea of this repository to solving this problem is pretty simple: to fit all the two-population models we can possibly think of to our experimental 2dSFS and use Akaike Information Criterion to select the best one.

See [GADMA](https://github.com/ctlab/GADMA) for the alternative solution to this problem. Compared to GADMA, we are far less elegant but somewhat more flexible (we can incorporate essentially any model, for example, involving selection and heterogeneous introgression rates across the genome). Our approach also lets the user evaluate how much better the winning model is compared to certain "null" alternatives (for example, models with no population split, with symmetrical migration, or with constant population sizes), which provides statistical evidence for general aspects of the model structure. Our disadvantage is a ridiculously huge number of model-fit runs that we have to perform; the good news is, all this can be done in parallel, and each single-model run takes no more than 2 hours.

## Installation ##
First of all, install *moments*. The example below would clone it into root directory and install it for a specific user.
```bash
cd
git clone https://bitbucket.org/simongravel/moments.git 
cd moments
python setup.py build_ext --inplace
cd
```
Then, clone this repository and copy all the `*.py` files from `~/AFS-analysis-with-moments/multimodel_inference/py2/` or from `/AFS-analysis-with-moments/multimodel_inference/py3/` (depending on your `python` version) to where you keep your executables (for example, `~/bin`). 
> NOTE: all code examples here assume the repository is cloned in the home directory, `~/`. If you cloned it elsewhere, make sure to replace `~/` in all examples with the actual path.

## Overview of the method ##
- The first step is **model selection**, where we run all possible models on 10 bootstrapped SFS. We actually run each model on each bootstrap six times (six random restarts), to make sure the model converges to its best likelihood at least once. Then we use an `R` script `model_selection_10boots.R` to select the best-fitted instance (out of 6) for each model for each bootstrap, and compare the AIC scores for all models. The best model is the one with the *lowest median AIC score among bootstrap replicates*.  
- The second step is running the winning model on 100 bootstrapped SFS, to **evaluate parameter uncertainties**. Once again, we will have to do 6 random restarts for each bootstrap. The parameter meanings and uncertainties are deciphered by the second R script that we have, `bestmodel_bootstrap.R`.

## Model selection ##
Let's assume we have ten bootstrapped 2dSFS formatted for *moments* or *dadi* (See **Appendix** for instructions how to obtain bootstrapped 2dSFS from [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD)).Such file is nothing more than a line of numbers with a header line giving the dimensions of the spectrum ( 2 x N + 1 for each of the two populations, where N is the number of sampled diploids). 
Bootstrapped sfs shoudl be named like `p12_1.sfs`, `p12_2.sfs` etc. where `p12` is the name of population contrast.

```bash
cd [where your boostrapped SFS files are]
cp ~/AFS-analysis-with-moments/multimodel_inference/allmodels_unfolded allmodels
# if your SFS needs to be folded, use this line instead:
# cp ~/AFS-analysis-with-moments/multimodel_inference/allmodels_folded allmodels
NREPS=6 # number of random restarts per model per bootstrap rep
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
>${CONTRAST}.stdout
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS >>${CONTRAST}.stdout" >>args;
done;
paste mods args -d " " >>modsel;
done
```
Run all commands in `modsel` file. This is the most computaitonally intensive thing I have ever done - there are 6 x 108 x 10 model runs, requiring 2 hours each. Best run these on an HPC cluster, in parallel! Note that all the screen output is collected in a file, `p12.stdout` in this case.

Then, to extract results:
```bash
CONTRAST=p12 # don't forget to regenerate $CONTRAST variable if you had to re-login!
grep RESULT ${CONTRAST}.stdout -A 4 | grep -E "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep RESULT >${CONTRAST}.res
```
Lastly, run `modelSelection.R` on the file `${CONTRAST}.res`:
```bash
Rscript modelSelection.R infile=${CONTRAST}.res
```
Two plots will be generated. The first one is the boxplot of best AIC scores for each model for all bootstrap replicates:
![all boxplots](all_boxplots.png)

And the second one is the plot of just the AIC medians for the top 10 models:

![top10](top10_medians.png)

The script also outputs the text file named `[infile].[modelname]`, **where `[modelname]` is the name of the winning model**. This file contains the fitted parameter values for the winning model, which will be used at the next stage as "guiding values" for random restarts.

## Bootstrapping the winning model ## 
Assuming we have 100 boostrapped SFS (See **Appendix** for instructions how to obtain bootstrapped 2dSFS from [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD)), we are now going to run just the winning model on all of them. We will do 6 random restarts per bootstrap replicate, using parameters output by `modelSelection.R` as "guides": these values will be randomly perturbed (changed on average 2-fold up or down) at the start of each model run.

```bash
CONTRAST=p12 
WINNER="sc3ielsm1" # change this to your winning model, according to stage 1
ARGS="p1 p2 16 16 0.02 0.005 ${CONTRAST}.res.${WINNER}"

NREPS=6 # number of random restarts per bootstrap rep
>mods
for i in `seq 1 $NREPS`;do 
echo ${WINNER}.py >>mods;
done

>winboots
for B in `seq 1 100`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
>${CONTRAST}.boots
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS >>${CONTRAST}.boots" >>args;
done;
paste mods args -d " " >>winboots;
done
```
Run all commands in `winboots`, collecting all the text output in a file `p12.boots`. Then extract results from that file and run `bestmodel_bootstrap.R` on them:
```bash
grep RESULT ${CONTRAST}.boots -A 4 | grep -E "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep RESULT >${CONTRAST}.boots.res
Rscript bestmodel_bootstrap.R infile=${CONTRAST}.boots.res
```
>Note: Additonal options to `bestmodel_bootstrap.R` are:
- `topq`: top quantile cutoff. Only boostrap runs in this top quantile will be summarized. Default 0.25.
- `path2models`: path to the subdir `multimodel_inference` within this cloned repository. Default `~/AFS-analysis-with-moments/multimodel_inference/`.

This will generate the histogram of likelihoods and red line for top=quantile cutoff:
![boots histogram](boothist.png)
...and, finally, boxplots for parameter estimates:
![params](boot_params.png)

The script also saves an RData bundle containing the summary dataframe (medians, 25% quantile, 75% quantile for all parameters) and the big dataframe containing all summarized bootstrap data.



## Appendix ## 

[Link to original *Moments* paper]( http://www.genetics.org/content/early/2017/05/08/genetics.117.200493)

[*Moments* manual](https://bitbucket.org/simongravel/moments/raw/efc4da3047226e3662dd43b525e41c85b93e90fd/doc/manual/manual.pdf) (link may change with updates)
