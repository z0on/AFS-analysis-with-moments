# read about multimodel inference here: 
# https://pdfs.semanticscholar.org/a696/9a3b5720162eaa75deec3a607a9746dae95e.pdf

if (length(commandArgs(trailingOnly=TRUE))<1) {
	options(warning.length=8000)
	stop("
	
Ranks models based on AIC. Models are assumed to have been run on 10+ boostrapped replicates, each model with 6+ random starts.

infile=[filename]    results of model runs. Obtained by summarizing the STDOUT output of model fitting
                     by running 
                     
grep RESULT manymodels.out -A 4 | grep -E \"[0-9]|\\]\" | perl -pe 's/^100.+\\.o\\d+\\S//' | perl -pe 's/\\n//' | perl -pe 's/[\\[\\]]//g' | perl -pe 's/RESULT/\\nRESULT/g' | grep RESULT >manymodels.res

topq=0.25        top quantile to summarize. 0.25 means that only the best-likelihood 25% 
                of bootstrap replicates will be used to summarize paramter values.

")
}

infl=grep("infile=",commandArgs())
if (length(infl)==0) { stop ("specify input file (infile=filename)\nRun script without arguments to see all options\n") }
infile=sub("infile=","", commandArgs()[infl])

#infile="c14.res"

library(ggplot2)

system(paste("cut -f 2,3,4,5,6 -d ' ' ",infile," > ",infile,".likes",sep=""))
npl=read.table(paste(infile,".likes",sep=""))
system(paste("rm ",infile,".likes",sep=""))

names(npl)=c("model","id","npara","ll","boot")
#head(npl)
#npl=npl[grep(infile,npl$boot),]
aics=list()
for (b in 1:length(unique(npl$boot))) {
	bb=levels(npl$boot)[b]
	nplb=subset(npl,boot==bb)
	maxlike=c();nmod=c()
	for (m in unique(nplb$model)) {
		sub=subset(nplb,model==m)
		nmod=c(nmod,nrow(sub))
		maxlike=data.frame(rbind(maxlike,sub[sub$ll==max(sub$ll),]))
	}
	npara=maxlike$npara
	likes=maxlike$ll
	aic=2*npara-2*likes
	aicc=data.frame(cbind(model=as.character(unique(nplb$model))))
	aicc$aic=aic
	aicc$nmod=nmod
	aicc$boot=bb
	aics[[b]]=aicc
}
awt=data.frame(do.call(rbind,aics))
models=unique(awt$model)
med=c()
for (m in models) {
	ss=subset(awt,model==m)
	med=c(med,median(ss$aic))
}
modmed=data.frame(cbind(mod=as.character(models)))
modmed$med=med
modmed=modmed[order(med),]
modmed$mod=factor(modmed$mod,levels=modmed$mod)
awt$model=factor(awt$model,levels=modmed$mod)

pdf(width=2.2,height=2,file=paste(infile,"_modsel_top10medians.pdf",sep=""))
pp=ggplot(modmed[1:10,],aes(mod, med))+geom_point()+theme(axis.text.x = element_text(angle = 45,hjust=1))
#ggplot(modmed,aes(mod, med))+geom_point()+theme(axis.text.x = element_text(angle = 45,hjust=1))
plot(pp)
dev.off()
pdf(width=15,height=5,file=paste(infile,"_modsel_allBoxplots.pdf",sep=""))
pp=ggplot(awt,aes(model,aic))+geom_boxplot()+scale_y_log10()+theme(axis.text.x = element_text(angle = 45,hjust=1))
plot(pp)
dev.off()

# ----- extracting name and parameters of the winning model

winner=as.character(modmed[1,1])

system(paste("grep ",winner," ", infile," > ",infile,".winmod",sep=""))
npl0=read.table(paste(infile,".winmod",sep=""))
system(paste("rm ",infile,".winmod",sep=""))

npl0=npl0[which(npl0$V5==max(npl0$V5)),]
params=as.vector(npl0[1,c(9:(ncol(npl0)-1))])
write.table(params,file=paste(infile,".",winner,sep=""),quote=F,col.names=F,row.names=F)

print(as.character(modmed[1,1]),quote=F)

