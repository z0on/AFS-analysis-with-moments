# read about multimodel inference here: 
# https://pdfs.semanticscholar.org/a696/9a3b5720162eaa75deec3a607a9746dae95e.pdf

npl=read.table("c34wins.res")

#------ retrieving parameter names, applying them

pa=readLines("~/AFS-analysis-with-moments/multimodel_inference/unfolded_params")
wmod=sub("sc","SC",npl[1,2])
npl=npl[,-c(1:2)]
params=c(strsplit(gsub("[ \t]","",pa[grep(wmod,pa)]),split="[:,]")[[1]][-1],"theta")
names(npl)=c("id","np","ll","boot","p1","p2",params)
head(npl)

#------ finding best likelihood for each bootstrap rep
maxlike=c()
for (b in 1:length(levels(npl$boot))) {
	bb=levels(npl$boot)[b]
	sub=subset(npl,boot==bb)
	maxlike=data.frame(rbind(maxlike,sub[sub$ll==max(sub$ll),]))
}
head(maxlike)

# ---- removing bottom 25% of bootstraps

hist(maxlike$ll,breaks=20)
maxlike=maxlike[maxlike$ll>quantile(maxlike$ll,0.75),]
hist(maxlike$ll,breaks=20)

# ----- "folding" genomic islands

fold=rep(FALSE,nrow(maxlike))
fold[maxlike$P>0.5]=TRUE

mxl=maxlike
migrations=grep("^m",names(maxlike))
migrations.i=grep("^m.*i",names(maxlike))
migrations= migrations[!(migrations %in% migrations.i)]

maxlike[fold,migrations]=mxl[fold,migrations.i]
maxlike[fold,migrations.i]=mxl[fold,migrations]
maxlike[fold,"P"]=1-maxlike[fold,"P"]

# ---- log-transforming migration rates

migrations=grep("^m",names(maxlike))
maxlike[,migrations]=apply(maxlike[,migrations],2,log,base=10)

# ---- makig a long table, defining parameter types

ms=stack(maxlike[,7:ncol(maxlike)])
names(ms)[2]="parameter"
ms$type="misc"
ms$type[grep("^nu",ms$parameter)]="nu"
ms$type[grep("^T",ms$parameter)]="T"
ms$type[grep("^m",ms$parameter)]="log10.migr"
ms$type[ms$parameter=="P"]="prop.islands"
ms$type[ms$parameter=="theta"]="theta"
ms$type[ms$parameter=="p_misid"]="p_misid"
ms$type=factor(ms$type,levels=c("nu","T","log10.migr","prop.islands","theta","p_misid"))

# ---- plotting
library(ggplot2)
ggplot(ms,aes(parameter,values))+geom_boxplot()+theme(axis.text.x = element_text(angle = 45,hjust=1))+facet_wrap(~type,scale="free",nrow=1)
