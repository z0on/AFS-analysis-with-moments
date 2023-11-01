#!/usr/bin/env Rscript
if (length(commandArgs(trailingOnly=TRUE))<1) {
  options(warning.length=8000)
  stop("

Takes .mafs.gz output from angsd, produces NJACK replicates of sites lists
sampling PROP proportion of contigs (listed as chromo in the .mafs table) 
wihtout replacement (for making SFS replicates)

arguments:

mafs=[filename]   mafs.gz output from angsd (use `-doMajorMinor 1 -doMaf 1` to make that)
NJACK=[integer]   how many replicates to make
PROP=0.9          proportion of contigs to sample


Output:           Sites files named [mafs lead name]_1.sites, [mafs lead name]_2.sites etc.

Mikhail Matz, matz@utexas.edu, October 2023

")
}

# require(dplyr)
# require(R.utils)
# require(data.table)
#

mff =grep("mafs=",commandArgs())
if (length(mff)==0) { stop ("specify mafs.gz file (mafs=filename)\nRun script without arguments to see all options\n") }
mff=sub("mafs=","", commandArgs()[mff])

nj =grep("NJACK=",commandArgs())
if (length(nj)==0) { nj=10 } else { nj =as.numeric(sub("NJACK=","", commandArgs()[nj])) }

pp=grep("PROP=",commandArgs())
if (length(pp)==0) { pp=0.9 } else { pp =as.numeric(sub("PROP=","", commandArgs()[pp])) }

message("making ", nj," replicates\n");
message("sampling ",pp," fraction of all contigs each rep\n")

#mff="~/Dropbox/D.mafs.gz";pp=0.9;nj=2
mafs=read.table(mff,header=TRUE)[,1:2]
cc=unique(mafs$chromo)

for (r in 1:nj){
  message("replicate ",r,"\n")
  boot=sample(cc,size=round(length(cc)*pp))
  bt=mafs[mafs$chromo %in% boot,]
  outname=paste(sub(".mafs.gz","",mff),"_",r,".sites",sep="")
  write.table(bt,file=outname,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
}

