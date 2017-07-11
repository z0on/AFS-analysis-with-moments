# Read about delta-AIC here:
# http://theses.ulaval.ca/archimede/fichiers/21842/apa.html

# can be used with more than 2 models in the same time

# model 1 :
npar1=7
like1=-2000

# model 2 :
npar2=5
like2=-2004

# model 3 :
npar3=4
like3=-2002

# and so on for each model

aic1=2*npar1-2*like1 
aic2=2*npar2-2*like2
aic3=2*npar3-2*like3
# and so on for each model

aic=c(aic1,aic2,aic3) # add model if needed

aic.weights=exp(-(aic-min(aic)))/sum(exp(-(aic-min(aic))))

# Akaike weights: strength of evidence for each model
aic.weights