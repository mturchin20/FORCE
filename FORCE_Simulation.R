### Clear Console ###
cat("\014")

### Clear Environment ###
rm(list = ls(all = TRUE))

### Load in the R libraries ###
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
#library(Rcpp, lib="/users/mturchin/data/mturchin/Software/Module_R-3.4.0")
#library(RcppArmadillo, lib="/users/mturchin/data/mturchin/Software/Module_R-3.4.0")
library(RcppParallel)
library(CompQuadForm)
library(Matrix)
library(MASS)
library(truncnorm)

### Load in the MAPIT Functions ###
sourceCpp("/users/mturchin/LabMisc/RamachandranLab/FORCE/FORCE.cpp") 

### Set the OpenMP Envrionment ###
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")

######################################################################################
######################################################################################
######################################################################################

ind = 2e3; nsnp = 1e4;
maf <- 0.05 + 0.45*runif(nsnp)
X   <- (runif(ind*nsnp) < maf) + (runif(ind*nsnp) < maf)
X   <- matrix(as.double(X),ind,nsnp,byrow = TRUE)
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)

pve = 0.6; ncausal= 100 #Set 1 of causal SNPs

## Select Causal SNPs ###
snp.ids = 1:nsnp
pthwy1 = sample(snp.ids, ncausal/2, replace=F)
pthwy2 = sample(snp.ids[-pthwy1], ncausal/2, replace=F)

### Simulate the Additive Effects ###
Xmarginal = X[,c(pthwy1,pthwy2)]
beta=rnorm(dim(Xmarginal)[2])
y_marginal=c(Xmarginal%*%beta)
beta=beta*sqrt(pve/var(y_marginal))
y_marginal=Xmarginal%*%beta

# Simulate the (Environmental) Error/Noise
y_err=rnorm(ind)
y_err=y_err*sqrt((1-pve)/var(y_err))

y=y_marginal+y_err

### Check dimensions and add SNP names ###
dim(X); dim(y)

######################################################################################
######################################################################################
######################################################################################

### Determine a Preset List of Pathways (i.e. Indices) ###
regions = list()
for(k in 1:98){
  regions[[k]] = sort(sample(snp.ids[-c(pthwy1,pthwy2)], ncausal, replace=F))
}
regions[[99]] = sort(pthwy1); regions[[100]] = sort(pthwy2)

######################################################################################
######################################################################################
######################################################################################

### Set the number of cores ###
cores = detectCores()

### Run MAPIT: Full Version with Satterthwaite Method###
ptm <- proc.time() #Start clock
vc.mod = FORCE(t(X),y,regions,cores = cores)
proc.time() - ptm #Stop clock

### Apply Davies Exact Method ###
vc.ts = vc.mod$Est
names(vc.ts) = paste("Pathway",c(1:length(regions)),sep="_")

pvals = c()
for(i in 1:length(vc.ts)){
  lambda = sort(vc.mod$Eigenvalues[,i],decreasing = T)
  Davies_Method = davies(vc.mod$Est[i], lambda = lambda, acc=1e-8)
  pvals[i] = 2*min(1-Davies_Method$Qq,Davies_Method$Qq)
  names(pvals)[i] = names(vc.ts[i])
}
summary(pvals)
