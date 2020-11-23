#library("MAPITR")
library("devtools")
devtools::load_all("/users/mturchin/LabMisc/RamachandranLab/MAPITR")

args <- commandArgs()
print(args)
X.File <- args[6]
Genes.File <- args[7]
Covars.File <- args[8]
Output1.File <- args[9]
seed.value <- as.numeric(as.character(args[10]))
n.datasets <- as.numeric(as.character(args[11]))
pve <- as.numeric(as.character(args[12]))
rho <- as.numeric(as.character(args[13]))
pc.var <- as.numeric(as.character(args[14]))
ngenes <- as.numeric(as.character(args[15]))
ncausal1 <- as.numeric(as.character(args[16]))
ncausal2 <- as.numeric(as.character(args[17]))
#ncausal3 <- args[18]

set.seed(seed.value)

#X <- as.matrix(read.table("/users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v3.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.ForSimulations.chr16.raw.edit.Rheaders.gz", header=T));
#X <- as.matrix(read.table("/users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran40000/Imputation/mturchin20/ukb_chrAll_v3.British.Ran40000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.ForSimulations.chr16.raw.edit.Rheaders.gz", header=T));
#Genes <- read.table("/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/Simulations/Data/ukb_chrAll_v3.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus20kb.SimFormat.Chr16.Rformat.txt", header=F);
#Genes <- read.table("/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/Simulations/Data/ukb_chrAll_v3.British.Ran40000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus20kb.SimFormat.Chr16.Rformat.txt", header=F);
#Covars <- read.table("/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v3.African.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.wAC.txt", header=T);
#Covars <- read.table("/users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran40000/mturchin20/ukb_chrAll_v3.British.Ran40000.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.wAC.txt", header=T);
#PCs <- Covars[,(ncol(Covars)-9):ncol(Covars)];
X <- as.matrix(read.table(X.File, header=T));
Genes <- read.table(Genes.File, header=F);
Covars <- read.table(Covars.File, header=T);
PCs <- as.matrix(Covars[,(ncol(Covars)-9):ncol(Covars)]);

Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)
ind = nrow(X); nsnp = ncol(X)

### Define the Simulation Parameters ###
n.datasets = 1 #Total Number of Simulations 
pve = 0.6; #Heritability of the trait
rho = 0.8; #Proportion of the heritability caused by additive effects {0.8, 0.5}
pc.var = 0.1

### Set Up Causal Pathways in Three Groups
ngenes = 20; #number of genes needed to get around 1000 SNPs pulled for G1, G2, and G3
ncausal1 = .05; ncausal2 = .2 #Percent of SNPs needed for G1 and G2
ncausal3 = 1-(ncausal1+ncausal2) #Remaining SNPs to be allocated to G3

### Create a list to save the final Results ###
pval_mat = matrix(nrow = nrow(Genes),ncol = n.datasets); pval_mat <- cbind(Genes[,1], pval_mat);
genes_chosen = matrix(nrow = ngenes, ncol = n.datasets); 
G1_snps = matrix(nrow = ncausal1,ncol = n.datasets)
G2_snps = matrix(nrow = ncausal2,ncol = n.datasets)

### Run the Analysis ###
for(i in 1:n.datasets) {
  
  #Select Causal Pathways
  gene.ids = 1:nrow(Genes)
  genes.pulled.ids <- sample(gene.ids, ngenes, replace=F);
  genes.pulled <- Genes[genes.pulled.ids,];

  genes.pulled.SNPs <- c();
  for (j in 1:nrow(genes.pulled)) {
	genes.pulled.SNPs <- c(genes.pulled.SNPs, unlist(strsplit(as.character(genes.pulled[j,2]), ",")));
  }
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
  Xepi = c(); b = c()
  for(j in 1:ncausal1){
      Xepi = cbind(Xepi,X[,s1[j]]*X[,s2]) 
  }
  
  ### Simulate the Pairwise Effects ###
  beta=rnorm(dim(Xepi)[2])
  y_epi=c(Xepi%*%beta)
  beta=beta*sqrt(pve*(1-rho)/var(y_epi))
  y_epi=Xepi%*%beta
  
  ### Define the effects of the PCs ###
  beta=rnorm(dim(PCs)[2])
  y_pcs=c(PCs%*%beta)
  beta=beta*sqrt(pc.var/var(y_pcs))
  y_pcs=PCs%*%beta
  
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
 
  vc.mod = InterPath(t(X),y,regions,cores = cores)
  
  MAPITR_Output <- MAPITR(X, y, Genes) 
#  MAPITR_SimData_Genotypes, MAPITR_SimData_Phenotype, MAPITR_SimData_Pathways)
  
  ### Save Results ###
#  pval_mat[,j] = pvals
#  G1_snps[,j] = Pthwys_1
#  G2_snps[,j] = Pthwys_2

}

#Save final Results
#Final = list(pval_mat,G1_snps,G2_snps)

#### Save the Results ###
#file = "/home/lcrawfo1/Results/InterPath/PopStruct_Sim_rho8_S1.RData"
#save(Final, file = file)

