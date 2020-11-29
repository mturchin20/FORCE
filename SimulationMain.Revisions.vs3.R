#library("MAPITR")
library("devtools")
devtools::load_all("/users/mturchin/LabMisc/RamachandranLab/MAPITR")

args <- commandArgs()
#print(args)
X.File <- args[6]
Genes.File <- args[7]
Covars.File <- args[8]
Genes.Analysis.File <- args[9]
Output1.File <- args[10]
seed.value <- as.numeric(as.character(args[11]))
n.datasets <- as.numeric(as.character(args[12]))
pve <- as.numeric(as.character(args[13]))
rho <- as.numeric(as.character(args[14]))
pc.var <- as.numeric(as.character(args[15]))
#ngenes <- as.numeric(as.character(args[16]))
#nsnps <- as.numeric(as.character(args[16]))
ncausaltotal <- as.numeric(as.character(args[16]))
ncausal1 <- as.numeric(as.character(args[17]))
ncausal2 <- as.numeric(as.character(args[18]))
ncausal3 <- as.numeric(as.character(args[19]))
#ncausal3 = 1-(ncausal1+ncausal2) #Remaining SNPs to be allocated to G3
#nmask <- as.numeric(as.character(args[20]))

set.seed(seed.value)

X <- as.matrix(read.table(X.File, header=T));
Genes <- read.table(Genes.File, header=F);
Covars <- read.table(Covars.File, header=T);
Genes.Analysis <- read.table(Genes.Analysis.File, header=F);
PCs <- as.matrix(Covars[,(ncol(Covars)-9):ncol(Covars)]);

#Genes <- Genes[1:20,]
#Genes.Analysis <- Genes.Analysis[1:20,]

Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)
ind = nrow(X); nsnp = ncol(X)

### Run the Analysis ###
for(i in 1:n.datasets) {
  
  #Select Causal Pathways

  genes.ids = 1:nrow(Genes)
  s1.genes.ids <- sample(genes.ids, ncausal1, replace=F);
  s2.genes.ids <- sample(genes.ids[-s1.genes.ids], ncausal2, replace=F);
  s3.genes.ids <- sample(genes.ids[c(-s1.genes.ids,-s2.genes.ids)], ncausal3, replace=F);
  genes.pulled <- rbind(cbind(as.character(Genes[s1.genes.ids,1]), rep("Epi1", length(s1.genes.ids))), cbind(as.character(Genes[s2.genes.ids,1]), rep("Epi2", length(s2.genes.ids))), cbind(as.character(Genes[s3.genes.ids,1]), rep("Add", length(s3.genes.ids))));
  s1 <- c(); s1.sizes <- c(); s2 <- c(); s2.sizes <- c(); s3 <- c(); s3.sizes <- c();
  for (j in 1:length(s1.genes.ids)) {
	gene.temp <- Genes[s1.genes.ids[j],];
	gene.temp.SNPs <- unlist(strsplit(as.character(gene.temp[1,2]), ","));
  	gene.temp.SNPs.pulled <- sample(gene.temp.SNPs, round(ncausaltotal * length(gene.temp.SNPs)), replace=F);
	s1 <- c(s1, gene.temp.SNPs.pulled); s1.sizes <- c(s1.sizes, length(gene.temp.SNPs.pulled));
  }
  for (j in 1:length(s2.genes.ids)) {
	gene.temp <- Genes[s2.genes.ids[j],];
	gene.temp.SNPs <- unlist(strsplit(as.character(gene.temp[1,2]), ","));
  	gene.temp.SNPs.pulled <- sample(gene.temp.SNPs, round(ncausaltotal * length(gene.temp.SNPs)), replace=F);
	s2 <- c(s2, gene.temp.SNPs.pulled); s2.sizes <- c(s2.sizes, length(gene.temp.SNPs.pulled));
  }
  for (j in 1:length(s3.genes.ids)) {
	gene.temp <- Genes[s3.genes.ids[j],];
	gene.temp.SNPs <- unlist(strsplit(as.character(gene.temp[1,2]), ","));
  	gene.temp.SNPs.pulled <- sample(gene.temp.SNPs, round(ncausaltotal * length(gene.temp.SNPs)), replace=F);
	s3 <- c(s3, gene.temp.SNPs.pulled); s3.sizes <- c(s3.sizes, length(gene.temp.SNPs.pulled));
  }
  print(c(ncausaltotal,ncausal1,ncausal2,ncausal3,length(s1.genes.ids),length(s2.genes.ids),length(s3.genes.ids),length(s1),length(s2),length(s3)));
  print(s1.sizes); print(s2.sizes); print(s3.sizes);
  
  ### Simulate the Additive Effects ###
  SNPs.additive = c(s1,s2,s3);
#  print(dim(X))
#  print(SNPs.additive)
  Xmarginal = X[,SNPs.additive]
  beta=rnorm(dim(Xmarginal)[2])
  y_marginal=c(Xmarginal%*%beta)
  beta=beta*sqrt(pve*rho/var(y_marginal))
  y_marginal=Xmarginal%*%beta
  
  ### Simulate Pairwise Interaction matrix ###
  y_epi <- 0;
  if ((ncausal1 > 0) && (rho < 1)) { 
	  Xepi = c(); b = c();
	  for(j in 1:length(s1)){
	      Xepi = cbind(Xepi,X[,s1[j]]*X[,s2]) 
	  }
	  
	  ### Simulate the Pairwise Effects ###
	  beta=rnorm(dim(Xepi)[2])
	  y_epi=c(Xepi%*%beta)
	  beta=beta*sqrt(pve*(1-rho)/var(y_epi))
	  y_epi=Xepi%*%beta
   }

  ### Define the effects of the PCs, if included in simulation ###
  y_pcs <- 0;
  if (pc.var > 0) { 
	  beta=rnorm(dim(PCs)[2])
	  y_pcs=c(PCs%*%beta)
	  beta=beta*sqrt(pc.var/var(y_pcs))
	  y_pcs=PCs%*%beta
  }

  ### Simulate the (Environmental) Error/Noise ###
  y_err=rnorm(ind)
  y_err=y_err*sqrt((1-pve-pc.var)/var(y_err))
  
  ### Simulate the Total Phenotypes ###
  y=y_marginal+y_epi+y_pcs+y_err
  y=(y-mean(y))/sd(y)
  
  ### Check dimensions ###
  print(seed.value); print(dim(X)); print(length(y));
  
  ######################################################################################
  ######################################################################################
  ######################################################################################

#  if (nmask > 0) {
#
#  remove G3 SNPs from X
#
#  }

  ptm <- proc.time() #Start clock
#  MAPITR_Output <- MAPITR(X,y,Genes.Analysis,Covariates=PCs) 
  MAPITR_Output <- MAPITR(X,y,Genes.Analysis) 
  print(proc.time() - ptm) #Stop clock
#  print(head(MAPITR_Output$Results))

  write.table(y, paste(Output1.File, ".Simulation.Pheno.txt", sep=""), quote=FALSE, col.name=FALSE, row.name=FALSE);
  write.table(genes.pulled, paste(Output1.File, ".Simulation.nGenes.txt", sep=""), quote=FALSE, col.name=FALSE, row.name=FALSE);
  write.table(s1, paste(Output1.File, ".Simulation.nCausal1.txt", sep=""), quote=FALSE, col.name=FALSE, row.name=FALSE);
  write.table(s2, paste(Output1.File, ".Simulation.nCausal2.txt", sep=""), quote=FALSE, col.name=FALSE, row.name=FALSE);
  write.table(MAPITR_Output$Results, paste(Output1.File, ".Results.Output.txt", sep=""), quote=FALSE, col.name=TRUE, row.name=FALSE);

}
