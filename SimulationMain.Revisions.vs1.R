library(MAPITR)

set.seed(11151990)

args <- commandArgs()
print(args)
rho <- args[1]
pc.var = 0.10 
#Variance Explained by PCs (Population Structure)

#X <- read.table("/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/Simulations/Data/ukb_chrAll_v3.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Chr16.txt.gz", header=T);
X <- read.table("/users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v3.African}QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.ForSimulations.chr16.raw.edit.Rheaders.gz", header=T);
##X.Headers <- read.table("/users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v3.African}QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.ForSimulations.chr16.raw.edit.SNPIDs.Rformat.txt", header=F);
Genes <- read.table("/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/Simulations/Data/ukb_chrAll_v3.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus20kb.SimFormat.Chr16.Rformat.txt", header=F);
#PCs <- read.table("", header=F);

#maf=apply(X, 2, mean)/2
#X = X[,maf>0.01]; dim(X)
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)

ind = nrow(X); nsnp = ncol(X)

### Define the Simulation Parameters ###
n.datasets = 100 #Total Number of Simulations 
pve = 0.6; #Heritability of the trait
rho = 0.8; #Proportion of the heritability caused by additive effects {0.8, 0.5}

### Set Up Causal Pathways in Three Groups
ngenes = 20; #number of genes needed to get around 1000 SNPs pulled for G1, G2, and G3
ncausal1 = .05; ncausal2 = .2 #Percent of SNPs needed for G1 and G2
ncausal3 = 1-(ncausal1+ncausal2) #Remaining SNPs to be allocated to G3

### Create a list to save the final Results ###
pval_mat = matrix(nrow = npthwy,ncol = n.datasets); rownames(pval_mat) = names(gene_list)
G1_snps = matrix(nrow = ncausal1,ncol = n.datasets)
G2_snps = matrix(nrow = ncausal2,ncol = n.datasets)

### Run the Analysis ###
for(i in 1:n.datasets){
  
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
  s3.ids.size <- round(ncausal3*length(genes.pulled.SNPs.uniq)); if (s3.size > length(genes.pulled.SNPs.uniq.ids[-c(s1.ids,s2.ids)])) { s3.size <- s3.size - 1; }
  s3.ids=sample(genes.pulled.SNPs.uniq.ids[-c(s1.ids,s2.ids)], s3.ids.size, replace=F)
  s1 <- genes.pulled.SNPs.uniq[s1.ids]
  s2 <- genes.pulled.SNPs.uniq[s2.ids]
  s3 <- genes.pulled.SNPs.uniq[s3.ids]

  ### Simulate the Additive Effects ###
  SNPs.additive = c(s1,s2,s3);
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
  }
  
  ### Simulate the Pairwise Effects ###
  beta=rnorm(dim(Xepi)[2])
  y_epi=c(Xepi%*%beta)
  beta=beta*sqrt(pve*(1-rho)/var(y_epi))
  y_epi=Xepi%*%beta
  
  ### Define the effects of the PCs ###
  beta=rnorm(dim(PCs)[2])
  y_pcs=c(PCs%*%beta)
  beta=beta*sqrt(pc.var/var(z_pcs))
  y_pcs=PCs%*%beta
  
  ### Simulate the (Environmental) Error/Noise ###
  y_err=rnorm(ind)
  y_err=y_err*sqrt((1-pve-pc.var)/var(y_err))
  
  ### Simulate the Total Phenotypes ###
  y=y_marginal+y_epi+y_pcs+y_err
  y=(y-mean(y))/sd(y)
  
  ### Check dimensions ###
  dim(X); length(y)
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Rewrite the Pathway List to Give Column IDs Instead of SNP Names ###
  regions = matrix(nrow = nsnp, ncol = length(gene_list))
  rownames(regions) = colnames(X)
  colnames(regions) = names(gene_list)
  for(k in 1:length(gene_list)){
      regions[which(colnames(X)%in%gene_list[[k]]),k] = 1
  }
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Set the number of cores ###
  cores = detectCores()
  
  ### Run InterPath ###
  ptm <- proc.time() #Start clock
  vc.mod = InterPath(t(X),y,regions,cores = cores)
  proc.time() - ptm #Stop clock
  
  ### Apply Davies Exact Method ###
  vc.ts = vc.mod$Est
  names(vc.ts) = colnames(regions)
  
  pvals = c()
  for(i in 1:length(vc.ts)){
      lambda = sort(vc.mod$Eigenvalues[,i],decreasing = T)
      Davies_Method = davies(vc.mod$Est[i], lambda = lambda, acc=1e-8)
      pvals[i] = 2*min(1-Davies_Method$Qq,Davies_Method$Qq)
      names(pvals)[i] = names(vc.ts[i])
  }
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Find power for the first group of SNPs ###
  Pthwys_1 = colnames(regions)[s1]
  
  ### Find power for the second group of SNPs ###
  Pthwys_2 = colnames(regions)[s2]
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Save Results ###
  pval_mat[,j] = pvals
  G1_snps[,j] = Pthwys_1
  G2_snps[,j] = Pthwys_2
  
  ### Report Status ###
  cat("Completed Dataset", j, "\n", sep = " ")
}

#Save final Results
Final = list(pval_mat,G1_snps,G2_snps)

### Save the Results ###
file = "/home/lcrawfo1/Results/InterPath/PopStruct_Sim_rho8_S1.RData"
save(Final, file = file)

