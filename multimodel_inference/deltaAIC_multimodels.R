# read about multimodel inference here: 
# https://pdfs.semanticscholar.org/a696/9a3b5720162eaa75deec3a607a9746dae95e.pdf

dir <- "set1/mcav/analysis_sfs/moments/c13/"
weight.cutoff=1e-2 # lowest Akaike weight to report a model
npl=read.table(paste0(dir,"likes"))

names(npl)=c("model","id","npara","ll")
#npl$model=as.factor(npl$model)

# finding minimum likelihood for each model
maxlike=c()
for (m in unique(npl$model)) {
  sub=subset(npl,model==m)
  maxlike=data.frame(rbind(maxlike,sub[sub$ll==max(sub$ll),]))
}
npara=maxlike$npara
likes=maxlike$ll

aic=2*npara-2*likes
delta.aic= aic-min(aic)
delta.aic[delta.aic>100]=100

# weights of evidence
exps=exp(-delta.aic/2)
wts=exps/sum(exps)
wts
maxlike$wts=wts
maxlike=maxlike[order(maxlike$wts,decreasing=T),]
maxlike=subset(maxlike,wts>weight.cutoff)
maxlike$model=factor(maxlike$model,levels=maxlike$model[order(maxlike$wts,decreasing=T)])
plot(wts~model,maxlike,las=2,log="y",xlab="",ylab="evidence weight")
maxlike

# report parameter estimates for top models
maxlike$model <- as.character(maxlike$model)

library(stringr)
mod.results <- apply(read.table(paste0(dir,"mmods.res2"), sep = ":", as.is = T), 1, function(x) gsub("\\[|\\]|\\[ | \\]", "", x))
mod.params.v <- mod.results[grep(paste(with(maxlike, paste(model, id, sep = " ")), collapse = "|"), mod.results)]
mod.params.v <- gsub("\\s+", " ", str_trim(mod.params.v))
mod.params.l <- sapply(as.list(mod.params.v), function(x) strsplit(x, split = " "))

param.ids <- read.table("set1/mcav/analysis_sfs/moments/folded_params", sep = ":", as.is = T)

params.out <- list()
for(i in 1:length(mod.params.l)){
  model.id <- mod.params.l[[i]][2]
  mod.index <- grep(paste0(mod.params.l[[i]][2],".py"), param.ids[,1], ignore.case = T)
  all.params <- unlist(strsplit(param.ids[mod.index,2], split = ","))
  nparams <- length(all.params)
  out.df <- data.frame(Parameter = all.params, Value = mod.params.l[[i]][seq(nparams)+8], Uncertainty = mod.params.l[[i]][seq(nparams)+(nparams+9)])
  out.df$Parameter <- gsub(" ", "", as.character(out.df$Parameter))
  params.out[[model.id]] <- out.df
  write.table(out.df, paste0(dir,model.id, "_params.txt"), quote = F, sep = "\t", row.names = F)
}
params.out


