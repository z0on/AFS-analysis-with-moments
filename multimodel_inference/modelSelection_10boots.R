# read about multimodel inference here: 
# https://pdfs.semanticscholar.org/a696/9a3b5720162eaa75deec3a607a9746dae95e.pdf
library(ggplot2)

setwd("~/Dropbox/keysMcav/mcav_aug2020/")
weight.cutoff=0 # lowest Akaike weight to report a model

ins=c("c12","c13","c14","c23","c24","c34")
inn="c13"
for (inn in ins) {
	infile=paste(inn,".likes",sep="")
	npl=read.table(infile)
	names(npl)=c("model","id","npara","ll","boot")
	head(npl)
	npl=npl[grep(inn,npl$boot),]
	# only taking non-island models
	# npl=npl[grep("i",npl$model,invert=T),]
	
	#subset(npl,model=="sc3m" & boot=="c12_8.sfs")
	
	#npl$model=as.factor(npl$model)
	
	# finding minimum likelihood for each model
	aics=list()
	for (b in 1:length(unique(npl$boot))) {
		bb=levels(npl$boot)[b]
		nplb=subset(npl,boot==bb)
		maxlike=c();nmod=c();m="sc3m"
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
	
	pdf(width=2.2,height=2,file=paste(inn,"_modsel_top10medians.pdf",sep=""))
	pp=ggplot(modmed[1:10,],aes(mod, med))+geom_point()+theme(axis.text.x = element_text(angle = 45,hjust=1))
	#ggplot(modmed,aes(mod, med))+geom_point()+theme(axis.text.x = element_text(angle = 45,hjust=1))
	plot(pp)
	dev.off()
	pdf(width=15,height=5,file=paste(inn,"_modsel_allBoxplots.pdf",sep=""))
	pp=ggplot(awt,aes(model,aic))+geom_boxplot()+scale_y_log10()+theme(axis.text.x = element_text(angle = 45,hjust=1))
	plot(pp)
	dev.off()
}
# awt20=subset(awt,model %in% modmed[1:20,"mod"])
# ggplot(awt20,aes(model,aic))+geom_boxplot()+scale_y_log10()+theme(axis.text.x = element_text(angle = 45,hjust=1))
