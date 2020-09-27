
if (length(commandArgs(trailingOnly=TRUE))<1) {
	options(warning.length=8000)
	stop("
	
Writes a list of commands to perform SFS model selection. Uses nreps random restarts fo each model on 10 bootstrapped SFS.
Writes to the file [pop.contrast].modsel.runs

Assume we have ten bootstrapped 2dSFS formatted for moments or dadi. Such file is just a line of numbers with a header line 
giving the dimensions of the spectrum ( 2 x N + 1 for each of the two populations, where N is the number of sampled diploids). 
Bootstrapped SFS files should be named like p12_1.sfs, p12_2.sfs, etc. where p12 is the name of population contrast.

contrast=\"p12\"                   name of the contrast: leading part of the bootstrapped sfs file names, for example, p12_1.sfs, p12_2.sfs etc.

nreps=6                            number of random restarts for each model for each bootstrap replicate 
nboots=10                          number of bootstrap replicates to use

args=[list of arguments]           names of pop1 and pop2, projection for pop1, projection for pop2, mutation rate (per genotyped portion of the 
                                   genome per generation), generation time in thousands of years. Population names can be anything. 
                                   For ANGSD-derived SFS, projections should be 0.8*2N for each population (ronded to integer); 
                                   in the case corresponding to the default setting each population was represented by 10 individuals.
                                   Example: args=\"p1 p2 16 16 0.02 0.005\"

folded=FALSE                       whether the analysis was using folded SFS
				                
path2models=\"~/AFS-analysis-with-moments/multimodel_inference/\"      path to the cloned repository

")
}

ctrst=grep("contrast=",commandArgs())
if (length(ctrst)==0) { stop ("specify contrast name\nRun script without arguments for details\n") }
contrast =sub("contrast=","", commandArgs()[ctrst])

args=grep("args=",commandArgs())
if (length(args)==0) { stop ("specify arguments for SFS model runs\nRun script without options to see details\n") }
args =sub("args=","", commandArgs()[args])

path2models =grep("path2models=",commandArgs())
if(length(path2models)>0) { path2models=sub("path2models=","", commandArgs()[path2models]) } else { path2models="~/AFS-analysis-with-moments/multimodel_inference/" }

nr=grep("nreps=",commandArgs())
if(length(nr)>0) { nreps=as.numeric(sub("nreps=","", commandArgs()[nr])) } else { nreps=6 }

nboots =grep("nboots=",commandArgs())
if(length(nboots)>0) { nreps=as.numeric(sub("nboots=","", commandArgs()[nboots])) } else { nboots=10 }

if(length(grep("folded=T",commandArgs()))>0) { folded=TRUE } else { folded=FALSE }

# contrast="p12"
# path2models="~/AFS-analysis-with-moments/multimodel_inference/"
# folded=T
# args="p1 p2 16 16 0.02 0.005"
# nreps=3
# nboots=10

if (folded) { 
	mods=scan(paste(path2models,"allmodels_folded",sep=""),what="character") 
	} else {
	mods=scan(paste(path2models,"allmodels_unfolded",sep=""),what="character") 		
	}

args2=c()
for (b in 1:nboots) {
	for (n in 1:nreps) {
		for (m in mods) {
			bname=paste(contrast,"_",b,".sfs",sep="")
			args2=c(args2,paste("sleep ",n," && ",m,".py ",bname," ",args," >>",contrast,".modsel",sep=""))
		}
	}
}

writeLines(args2,paste(contrast,".modsel.runs",sep=""))


