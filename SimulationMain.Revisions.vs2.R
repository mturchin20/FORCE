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
nsnps <- as.numeric(as.character(args[16]))
ncausal1 <- as.numeric(as.character(args[17]))
ncausal2 <- as.numeric(as.character(args[18]))
#ncausal3 <- args[19]
ncausal3 = 1-(ncausal1+ncausal2) #Remaining SNPs to be allocated to G3
#nmask <- as.numeric(as.character(args[20]))

set.seed(seed.value)

X <- as.matrix(read.table(X.File, header=T));
Genes <- read.table(Genes.File, header=F);
Covars <- read.table(Covars.File, header=T);
Genes.Analysis <- read.table(Genes.Analysis.File, header=F);
PCs <- as.matrix(Covars[,(ncol(Covars)-9):ncol(Covars)]);

Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)
ind = nrow(X); nsnp = ncol(X)

### Run the Analysis ###
for(i in 1:n.datasets) {
  
  #Select Causal Pathways
  gene.ids = 1:nrow(Genes)
#  genes.pulled.ids.temp <- sample(gene.ids, 1, replace=F);
  gene.ids.sampled <- sample(gene.ids);
  snp.total <- 0;
  genes.pulled.ids <- c();
  genes.pulled.SNPs <- c();
  snp.total.running <- c(0);
  numloop <- 1;
  while ((snp.total <= 1000) || (numloop >= length(gene.ids))) {
  	gene.temp <- Genes[gene.ids.sampled[numloop],];
  	genes.pulled.ids <- c(genes.pulled.ids, gene.ids.sampled[numloop]);
	genes.pulled.SNPs <- c(genes.pulled.SNPs, unlist(strsplit(as.character(gene.temp[1,2]), ",")));
  	snp.total <- length(genes.pulled.SNPs);
	snp.total.running <- c(snp.total.running, snp.total);
	numloop <- numloop + 1;
  }
  genes.pulled <- Genes[genes.pulled.ids,];
  print(c(numloop, snp.total, length(genes.pulled.ids), snp.total.running[length(snp.total.running)-1]));

#  for (j in 1:nrow(genes.pulled)) {
#	genes.pulled.SNPs <- c(genes.pulled.SNPs, unlist(strsplit(as.character(genes.pulled[j,2]), ",")));
#  }
  genes.pulled.SNPs.uniq <- unique(genes.pulled.SNPs);
  genes.pulled.SNPs.uniq.ids <- 1:length(genes.pulled.SNPs.uniq)
  s1.ids=sample(genes.pulled.SNPs.uniq.ids, round(ncausal1*length(genes.pulled.SNPs.uniq)), replace=F)
  s2.ids=sample(genes.pulled.SNPs.uniq.ids[-s1.ids], round(ncausal2*length(genes.pulled.SNPs.uniq)), replace=F)
  s3.ids.size <- round(ncausal3*length(genes.pulled.SNPs.uniq)); if (s3.ids.size > length(genes.pulled.SNPs.uniq.ids[-c(s1.ids,s2.ids)])) { s3.ids.size <- s3.ids.size - 1; }
  s3.ids=sample(genes.pulled.SNPs.uniq.ids[-c(s1.ids,s2.ids)], s3.ids.size, replace=F)
  s1 <- genes.pulled.SNPs.uniq[s1.ids]
  s2 <- genes.pulled.SNPs.uniq[s2.ids]
  s3 <- genes.pulled.SNPs.uniq[s3.ids]

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
  if (ncausal1 > 0) { 
	  Xepi = c(); b = c();
	  for(j in 1:ncausal1){
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
  print(dim(X)); print(length(y));
  
  ######################################################################################
  ######################################################################################
  ######################################################################################

#  if (nmask > 0) {
#
#  remove G3 SNPs from X
#
#  }

  ptm <- proc.time() #Start clock
  MAPITR_Output <- MAPITR(X,y,Genes.Analysis[1:10,],Covariates=PCs) 
  print(proc.time() - ptm) #Stop clock
#  print(head(MAPITR_Output$Results))

  write.table(y, paste(Output1.File, ".Simulation.Pheno.txt", sep=""), quote=FALSE, col.name=FALSE, row.name=FALSE);
  write.table(genes.pulled, paste(Output1.File, ".Simulation.nGenes.txt", sep=""), quote=FALSE, col.name=FALSE, row.name=FALSE);
  write.table(s1, paste(Output1.File, ".Simulation.nCausal1.txt", sep=""), quote=FALSE, col.name=FALSE, row.name=FALSE);
  write.table(s2, paste(Output1.File, ".Simulation.nCausal2.txt", sep=""), quote=FALSE, col.name=FALSE, row.name=FALSE);
  write.table(MAPITR_Output$Results, paste(Output1.File, ".Results.Output.txt", sep=""), quote=FALSE, col.name=TRUE, row.name=FALSE);

}
