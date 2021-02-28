# read about multimodel inference here: 
# https://pdfs.semanticscholar.org/a696/9a3b5720162eaa75deec3a607a9746dae95e.pdf

if (length(commandArgs(trailingOnly=TRUE))<1) {
	options(warning.length=8000)
	stop("
	
Summarizes bootstrap replicates for the best model.

bootRes=[filename]    results of bootstraps (likely named as [pop.contrast].winboots)
                     
folded=FALSE    whether the analysis was using folded SFS
				
topq=0.9        top quantile to summarize. For example 0.9 means that only the best-likelihood 90% 
                of bootstrap replicates will be used to summarize paramter values.
                
path2models=\"~/AFS-analysis-with-moments/work/\"   path to the cloned repository

")
}

# ----------- reading input parameters

infl=grep("bootRes=",commandArgs())
if (length(infl)==0) { stop ("specify input file (infile=filename)\nRun script without arguments to see all options\n") }
bootRes=sub("bootRes=","", commandArgs()[infl])

topq =grep("topq=",commandArgs())
if(length(topq)>0) { topq=as.numeric(sub("topq=","", commandArgs()[topq])) } else { topq=0.9 }

path2models =grep("path2models=",commandArgs())
if(length(path2models)>0) { path2models=sub("path2models=","", commandArgs()[path2models]) } else { path2models="~/AFS-analysis-with-moments/work/" }

if(length(grep("folded=T",commandArgs()))>0) { folded=TRUE } else { folded=FALSE }

# ----------- reading data

require(ggplot2)

'
setwd("/home/misha/Dropbox/yap_hetcheck_2020")
bootRes="ny.winboots"
folded=F
topq=0.9
path2models="/home/misha/Dropbox/Documents/perl_bin/moments_scripts/multimodel_inference/py3_v2/GA/"

'

system(paste("grep RESULT ", bootRes," -A 4 | grep -v Launcher | grep -E \"[0-9]|\\]\" | perl -pe 's/^100.+\\.o\\d+\\S//' | perl -pe 's/\\n//' | perl -pe 's/[\\[\\]]//g' | perl -pe 's/RESULT/\\nRESULT/g' | grep RESULT >", bootRes,".res",sep=""))

infile=paste(bootRes,".res",sep="")

npl=read.table(infile)
pdf(paste(bootRes,"_plots.pdf",sep=""),height=3, width=8)

#------ retrieving parameter names, applying them

if (folded) { 
  pa=readLines(paste(path2models,"folded_params",sep=""))
} else {
  pa=readLines(paste(path2models,"unfolded_params",sep=""))
}
wmod=as.character(npl[1,2])
npl=npl[,-c(1:2)]
#params=c(strsplit(gsub("[ \t]","",pa[grep(paste0(wmod,"_ga.py"),pa)]),split="[:,]")[[1]][-1],"theta")
params=c(strsplit(pa[grep(wmod,pa)],split="[ ,]")[[1]][-1],"theta")
names(npl)=c("id","np","ll","boot","p1","p2",params)
npl$boot=factor(npl$boot)

#------ finding best likelihood for each bootstrap rep
maxlike=c()
for (b in 1:length(levels(npl$boot))) {
	bb=levels(npl$boot)[b]
	sub=subset(npl,boot==bb)
	maxlike=data.frame(rbind(maxlike,sub[sub$ll==max(sub$ll),]))
}
#head(maxlike)

# ---- leaving only topq % of bootstraps

hist(maxlike$ll,breaks=50,main="bootstrap likes")
abline(v=quantile(maxlike$ll,(1-topq)),col="red")
#abline(v=quantile(maxlike$ll,topq+(1-topq)/2),col="red")
maxlike=maxlike[maxlike$ll>quantile(maxlike$ll,(1-topq)),]
#hist(maxlike$ll,breaks=20)

# ---- log-transforming migration rates

migrations=grep("^m",names(maxlike))
maxlike[,migrations]=apply(maxlike[,migrations],2,log,base=10)

# ----- "folding" genomic islands such that island migration is lower on average
# 
# if(length(grep("i",wmod))>0) {
# 
# 	mxl=maxlike
# 	migrations=grep("^m",names(maxlike))
# 	migrations.i=grep("^m.*i",names(maxlike))
# 	migrations=migrations[!(migrations %in% migrations.i)]
# 
# 	mm=apply(mxl[,migrations],1,mean)
# 	mi=apply(mxl[,migrations.i],1,mean)
# 
# 	maxlike[mm<mi,migrations]=mxl[mm<mi,migrations.i]
# 	maxlike[mm<mi,migrations.i]=mxl[mm<mi,migrations]
# 	maxlike[mm<mi,"P"]=1-maxlike[mm<mi,"P"]
# }
# 

# ---- making a long table, defining parameter types

ms=stack(maxlike[,7:ncol(maxlike)])
names(ms)[2]="parameter"
ms$type="misc"
ms$type[grep("^nu",ms$parameter)]="nu"
ms$type[grep("^T",ms$parameter)]="T"
ms$type[grep("^m",ms$parameter)]="log10.migr"
ms$type[grep("^P|^F",ms$parameter)]="prop.specials"
ms$type[ms$parameter=="theta"]="theta"
ms$type[ms$parameter=="p_misid"]="p_misid"
ms$type=factor(ms$type,levels=c("nu","T","log10.migr","prop.specials","theta","p_misid"))

# ---- plotting, saving results

pp=ggplot(ms,aes(parameter,values))+geom_boxplot()+theme(axis.text.x = element_text(angle = 45,hjust=1))+facet_wrap(~type,scale="free",nrow=1)
plot(pp)
medians=apply(maxlike[,7:ncol(maxlike)],2,median)
q25=apply(maxlike[,7:ncol(maxlike)],2,quantile,prob=0.25)
q75=apply(maxlike[,7:ncol(maxlike)],2,quantile,prob=0.75)
mres=data.frame(cbind(medians, q25,q75))
save(mres,maxlike,file=paste(bootRes,"_bootres.RData",sep=""))
print(medians)
dev.off()

# ---- finding represenative run (most similar to medians)

idpara=t(maxlike[,7:ncol(maxlike)])
cors=c()
for(i in 1:ncol(idpara)){
	cors=c(cors,cor(medians,idpara[,i]))
}
bestrun=maxlike$id[which(cors==max(cors))[1]]
system(paste("cp *_",bestrun,"_*pdf ", bootRes,"_representativeModel.pdf",sep=""))
system(paste("cp *_",bestrun,".png ", bootRes,"_representativeModel.png",sep=""))
message("representative run: ",bestrun)
