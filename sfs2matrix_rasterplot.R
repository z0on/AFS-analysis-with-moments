# converts one-liner 2dSFS to matrix that can be plotted
sfs2matrix=function(sfs,n1,n2,zero.ends=TRUE) {
  dd=matrix(ncol=2*n1+1,nrow=2*n2+1,sfs)
  if(zero.ends==TRUE) { dd[1,1]=dd[2*n2+1,2*n1+1]=0 }
  return(apply(dd,2,rev))
}

sfs=scan("~/Dropbox/porites_adultJuv_oct2023/aj4.sfs",skip=1)
sfsm=sfs2matrix(sfs,6,8)
library(raster)
plot(raster(log(sfsm+1,10)))

