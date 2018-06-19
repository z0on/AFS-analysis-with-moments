# INSTALLATIONS

------- Moments: 

cd
git clone https://bitbucket.org/simongravel/moments.git 
cd moments
python setup.py build_ext --inplace

# add this to .bashrc, section 2:
  export PYTHONPATH=$PYTHONPATH:$HOME/moments
# re-login

cds
cd RAD

------- ANGSD (if this does not work, look for changes in installation procedure on ANGSD github site) : 

# install xz first from https://tukaani.org/xz/
cd
wget https://tukaani.org/xz/xz-5.2.3.tar.gz --no-check-certificate
tar vxf xz-5.2.3.tar.gz 
cd xz-5.2.3/
./configure --prefix=$HOME/xz-5.2.3/
make
make install

# edit .bashrc:
nano .bashrc
   export LD_LIBRARY_PATH=$HOME/xz-5.2.3/lib:$LD_LIBRARY_PATH
   export LIBRARY_PATH=$HOME/xz-5.2.3/lib:$LIBRARY_PATH
   export C_INCLUDE_PATH=$HOME/xz-5.2.3/include:$C_INCLUDE_PATH
logout
# re-login

# now, install htslib:
cd
git clone https://github.com/samtools/htslib.git
cd htslib
make CFLAGS=" -g -Wall -O2 -D_GNU_SOURCE -I$HOME/xz-5.2.3/include"

cd
git clone https://github.com/ANGSD/angsd.git 
cd angsd
make HTSSRC=../htslib

# now adding ANGSD to $PATH
cd
nano .bashrc
# section 2:
   export PATH=$HOME/angsd:$PATH
   export PATH=$HOME/angsd/misc:$PATH
# save (Ctl-O, Ctl-X)


#==========================
# ANDSD => SFS for demographic analysis

# make separate files listing bams for each population (NB: without clones or replicates!)
# Let's assume we have two populations, pop0 and pop1, 20 individuals each, with corresponding bams listed in pop0.bams and pop1.bams

# filtering sites to work on - use only filters that do not distort allele frequency
# set minInd to 0.9 x (pop0 + pop1 sizes) 
cat pop0.bams pop1.bams > all.bams
FILTERS="-minMapQ 30 -minQ 35 -minInd 36 -doHWE 1 -sb_pval 1e-2 -hetbias_pval 1e-2 -skipTriallelic 1"
DOS="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -dogeno 3 -doPost 2"
angsd -b all.bams -GL 1 $FILTERS $DOS -P 1 -out sfsSites

# extracting and indexing list of sites to make SFS from 
# filtering out sites where heterozygote counts comprise more than 50% of all counts (likely lumped paralogs)
zcat sfsSites.snpStat.gz | awk '($3+$4+$5+$6)>0' | awk '($12+$13+$14+$15)/($3+$4+$5+$6)<0.5' | cut -f 1,2  >sites2do 
angsd sites index sites2do

# estimating site frequency likelihoods for each population 
export GENOME_REF=mygenome.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF""
# In the following lines, set minInd to 80-90% of 2x pop sample size; i.e., if sample size is 20 set to 2*20*0.9: 36)
angsd -sites sites2do -b pop0.bams -GL 1 -P 1 $TODO -out pop0
angsd -sites sites2do -b pop1.bams -GL 1 -P 1 $TODO -out pop1

# generating per-population SFS
realSFS pop0.saf.idx >pop0.sfs
realSFS pop1.saf.idx >pop1.sfs

# generating dadi-like posterior counts based on sfs priors
realSFS dadi pop0.saf.idx pop1.saf.idx -sfs pop0.sfs -sfs pop1.sfs -ref $GENOME_REF -anc $GENOME_REF >dadiout

# converting to dadi-snp format understood by dadi an Moments:
# (numbers after the input file name are numbers of individuals sampled per population)
realsfs2dadi.pl dadiout 20 20 >2pops_dadi.data

#=====================
# 2d AFS analysis using Moments

# get Misha's Moments scripts collection
git clone https://github.com/z0on/AFS-analysis-with-moments.git
move all files from AFS-analysis-with-moments to one of your $PATH locations, or set $PATH to include the AFS-analysis-with-moments directory

# print folded 2d SFS - for denovo or when mapping to genome of the studied species
# (change numbers to 2 x 0.9 x number of samples for in each pop):
2dAFS_fold.py 2pops_dadi.data pop0 pop1 36 36

# print unfolded 2d SFS - if mapping to genome of sister species
# (change numbers to 2 x 0.9 x number of samples for in each pop):
2dAFS.py 2pops_dadi.data pop0 pop1 36 36

# ------ multimodel inference: fit a diversity of 2-population models, then select the best one based on AIC.
# there are models with a period of exponential growth ("IM" models), models with one, two or three different size and/or migration rate epochs ("SC" models, including models with no migration in some epochs), models with symmetrical migration ("sm" models, in other cases migration is asymmetrical), models with two types of genomic loci ("genomic islands") introgressing at different rates ("i" models), and some fun combinations thereof.
# differences between models are summarized in excel table moments_multimodels.xls

# read about multimodel inference here: 
# https://pdfs.semanticscholar.org/a696/9a3b5720162eaa75deec3a607a9746dae95e.pdf

# this HAS to be parallelized - we need to fit ~100 models 10 times to make sure each model converges at its best fit at least once.

# input line: the last four numbers are:
# - projections (2 x 0.9 x number of samples) for in each pop;
# - mutation rate per gamete per generation
# - generation time, in thousand years
IN="2pops_dadi.data pop0 pop1 36 36 0.02 0.005"

# after setting the $IN variable, execute all commands listed in the text file "allmodels_unfolded" (if your alleles are polarized into ancestral and derived, for example by mapping to a sister species genome) or "allmodels_folded" 

# collecting results while fixing broken lines
cat *.mom | perl -pe 's/RESULT(.+)(\d)\n/RESULT$1$2/' |perl -pe 's/RESULT(.+)(\d)\n/RESULT$1$2/' |perl -pe 's/RESULT(.+)(\d)\n/RESULT$1$2/' | perl -pe 's/RESULT(.+)([\d\s])\n/RESULT$1$2/' | grep RESULT > mmods.res

# extracting likelihoods and parameter numbers for AIC:
cut -f 2,3,4,5 -d " " mmods.res >likes

# use R script deltaAIC_multimodels.R to find best-fitting model.
# using model ID number that the R script will list: 
# - examine the model's graphic output (*.pdf of actual and modeled SFS, and *.png of the model graph)
# - grep fitted model parameters and their SDs from *.mom files

# the order of parameters are listed in files unfolded_params and folded_params. Typically pop size parameters are first, then times, then migration rates, then the fraction of genomic "islands" (in "i"  models), then percentage of misidentified ancestral states (in unfolded models).

