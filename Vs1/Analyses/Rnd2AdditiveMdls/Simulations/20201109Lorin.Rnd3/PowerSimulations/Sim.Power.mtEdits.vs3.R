### Clear Console ###
#cat("\014")

### Clear Environment ###
#rm(list = ls(all = TRUE))

### Load in the R libraries ###
#library(doParallel)
#library(Rcpp)
#library(RcppArmadillo)
#library(RcppParallel)
#library(CompQuadForm)
#library(Matrix)
#library(MASS)
#library(truncnorm)
#library(data.table)
#library("devtools")
#devtools::load_all("/users/mturchin/LabMisc/RamachandranLab/MAPITR")
library("MAPITR")
#devtools::load_all("/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/Simulations/20201109Lorin.Rnd3/TestMAPITR/MAPITR")

### Load in the MAPIT Functions ###
#sourceCpp("~/Desktop/Simulations/Code/InterPath.cpp")
##sourceCpp("/gpfs/home/mturchin/LabMisc/RamachandranLab/MAPITR/src/MAPITR.cpp")
#sourceCpp("/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Vs2.download2.20210207.mtEdits.cpp")

### Set the OpenMP Envrionment ###
#Sys.setenv("PKG_CXXFLAGS"="-fopenmp")

######################################################################################
######################################################################################
######################################################################################

args <- commandArgs()
#print(args)
X.File <- args[6]
Genes.File <- args[7]
Covars.File <- args[8]
##old genes file was here <- args[9]
Output1.File <- args[10]
seed.value <- as.numeric(as.character(args[11]))
rounds.start <- as.numeric(as.character(args[12]))
pve <- as.numeric(as.character(args[13]))
rho <- as.numeric(as.character(args[14]))
pc.var <- as.numeric(as.character(args[15]))
ncausaltot <- as.numeric(as.character(args[16]))
ncausal1 <- as.numeric(as.character(args[17]))
ncausal2 <- as.numeric(as.character(args[18]))
ncausal3 <- as.numeric(as.character(args[19]))

Covars <- read.table(Covars.File, header=T);
Covars.PCs <- Covars[,(ncol(Covars)-9):ncol(Covars)];
print(head(Covars.PCs))

#6	7		8		9	10		11	12	13	14	15	16		17		18	19
#$X_File1 $Genes_File1 $Covars_File1 $Genes_File2 $Output1_File1 $Seed1 $Datasets1 $PVE1 $Rho1 $PCs_var1 $ncausaltotal1 $nCausal1a $nCausal2a $nCausal3a"

##X <- as.matrix(read.table(X.File, header=T));
##Genes <- read.table(Genes.File, header=F);

### Set the Seed for the analysis ###
set.seed(seed.value)

### Load in the Gene IDs ###
#don't think this is actually used anywhere?
#load("~/Desktop/Simulations/Data/gene_ids.RData") 

### Load in the Gene-SNP Mapping List ###
#load("~/Desktop/Simulations/Data/gene_snp_list.RData")
##gene_snp_list_full <- read.table("/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/Simulations/Data/ukb_chrAll_v3.British.Ran4000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.wRowPos.Regions.ExonicPlus20kb.SimFormat.Chr16.perSNPs.noGeneDups.GenesFormat.Rformat.txt", header=F);
gene_snp_list_full <- read.table(Genes.File, header=F);
gene_snp_list <- list()
for (i in 1:nrow(gene_snp_list_full)) {
	gene_snp_list[[as.character(gene_snp_list_full[i,1])]] <- strsplit(as.character(gene_snp_list_full[i,2]), ",")
}

gene_list = list()
for(i in 1:length(gene_snp_list)){
  x = unlist(gene_snp_list[[i]])
  gene_list[[i]] = x[!is.na(x)]
  names(gene_list)[i] = names(gene_snp_list)[i]
}

rm(gene_snp_list)

######################################################################################
######################################################################################
######################################################################################

### Load in the Geno Type Data ###
#load("/home/lcrawfo1/Data/UK_Biobank/chromosome16_snps.RData")
##chromosome16_snps <- as.matrix(read.table("/users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran4000/Imputation/mturchin20/ukb_chrAll_v3.British.Ran4000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.ForSimulations.chr16.raw.edit.Rheaders.gz", header=T));
chromosome16_snps <- as.matrix(read.table(X.File, header=T));
X = chromosome16_snps; rm(chromosome16_snps)
colnames(X) = sapply(colnames(X),function(x) unlist(strsplit(x,split = "_"))[1])
maf=apply(X, 2, mean)/2
X = X[,maf>0.01]; dim(X)
unique.snps = apply(X,2,function(x) length(unique(x)))
X = X[,unique.snps>1]; dim(X)
duplicated.snps = duplicated(t(X))
X = X[,duplicated.snps==FALSE]; dim(X)
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)

ind = nrow(X); nsnp = ncol(X)

######################################################################################
######################################################################################
######################################################################################

### Make Sure that Gene List and Genotype SNPs Match Up ###
for(i in 1:length(gene_list)){
  x = gene_list[[i]]
  gene_list[[i]] = x[x%in%colnames(X)] 
}

### Remove Duplicate SNPs ### 
repeated_snps = unlist(gene_list)[duplicated(unlist(gene_list))]
for(i in 1:length(gene_list)){
  x = gene_list[[i]]
  gene_list[[i]] = x[which(x%in%repeated_snps==FALSE)]  
}

rm(repeated_snps)

### Remove Any Gene with One or Fewer SNPs ###
gene_list = gene_list[unlist(lapply(gene_list,function(x){length(x)}))>=2]

npthwy = length(gene_list) #Final Number of Genes/Pathways 

######################################################################################
######################################################################################
######################################################################################

### Define the Simulation Parameters ###
##n.datasets = 2 #Total Number of Simulations 
##pve = 0.6; #Heritability of the trait
##rho = 0.5; #Proportion of the heritability caused by additive effects {0.8, 0.5}

### Set Up Causal Pathways in Three Groups
##ncausal1 = 5; ncausal2 = 50 #Pathways in each group
##ncausal3 = 150-(ncausal1+ncausal2) #Total Pathways

### Create a list to save the final Results ###
pval_mat = matrix(nrow = npthwy,ncol = 10); rownames(pval_mat) = names(gene_list)
G1_genes = matrix(nrow = ncausal1,ncol = 10)
G2_genes = matrix(nrow = ncausal2,ncol = 10)
##G1_snps = matrix(nrow = ncausal1,ncol = 10)
##G2_snps = matrix(nrow = ncausal2,ncol = 10)
##Add_snps = matrix(nrow = ncausal3,ncol = 10)
stats <- c()

#do a stderr print out here

Counter1 <- 1;
### Run the Analysis ###
for(j in rounds.start:(rounds.start+9)) {
  
  #Select Causal Pathways
  pthwy.ids = 1:npthwy
  s1=sample(pthwy.ids, ncausal1, replace=F)
  s2=sample(pthwy.ids[-s1], ncausal2, replace=F)
  s3=sample(pthwy.ids[-c(s1,s2)], ncausal3, replace=F)
  
  ### Simulate the Additive Effects ###
  snps = unlist(gene_list[c(s1,s2,s3)])
  Xmarginal = X[,snps]
  beta=rnorm(dim(Xmarginal)[2])
  y_marginal=c(Xmarginal%*%beta)
  beta=beta*sqrt(pve*rho/var(y_marginal))
  y_marginal=Xmarginal%*%beta
  
  ### Simulate Pairwise Interaction matrix ###
  Xepi = c(); b = c()
  ##causal_pthwys = gene_list[c(s1,s2,s3)]
  snps.s1.use.total <- c();
  snps.s2.use.total <- c();
  for (k in 1:ncausal2) { 
  	snps.s2.full.temp = unlist(gene_list[s2[k]])
	snps.s2.use.total <- c(snps.s2.use.total, sample(snps.s2.full.temp, ceiling(length(snps.s2.full.temp) * ncausaltot), replace=F))
  }
  print(head(snps.s2.use.total))
  print("yaya1")
  for(i in 1:ncausal1){
    snps.s1.full = unlist(gene_list[s1[i]])
    snps.s1.use = sample(snps.s1.full, ceiling(length(snps.s1.full) * ncausaltot), replace=F)
    print(head(snps.s1.full)); print(head(snps.s1.use));
    print("yaya2")
    for(k in 1:length(snps.s1.use)){
      Xepi = cbind(Xepi,X[,snps.s1.use[k]]*X[,snps.s2.use.total]) 
    }
    snps.s1.use.total <- c(snps.s1.use.total, snps.s1.use)
  }
  
  ### Simulate the Pairwise Effects ###
  beta=rnorm(dim(Xepi)[2])
  y_epi=c(Xepi%*%beta)
  beta=beta*sqrt(pve*(1-rho)/var(y_epi))
  y_epi=Xepi%*%beta
  
  ### Simulate the (Environmental) Error/Noise ###
  y_err=rnorm(ind)
  y_err=y_err*sqrt((1-pve)/var(y_err))
  
  ### Simulate the Total Phenotypes ###
  y=y_marginal+y_epi+y_err
  y=(y-mean(y))/sd(y)
  
  ### Check dimensions ###
  dim(X); dim(y)
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Rewrite the Pathway List to Give Column IDs Instead of SNP Names ###
#  regions = matrix(nrow = nsnp, ncol = length(gene_list))
#  rownames(regions) = colnames(X)
#  colnames(regions) = names(gene_list)
  regions <- c();
  X.seq <- seq(1,ncol(X),by=1);
  for(i in 1:length(gene_list)){
#      regions[which(colnames(X)%in%gene_list[[i]]),i] = 1
  	regions <- rbind(regions, c(names(gene_list)[i], paste(X.seq[which(colnames(X)%in%gene_list[[i]])], collapse=",")));
  }
  print(gene_list[[1]])
  print(regions[1,])
  print(gene_list[[2]])
  print(regions[2,])
##  regions.vs2 <- regions[,130:180]

  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Set the number of cores ###
#  cores = detectCores()

	y.old <- y
	y.lm <- lm(y ~ as.matrix(Covars.PCs) - 1)
	y.new <- residuals(y.lm)
	y <- y.new

#	print(head(y.old)); print(head(y.new)); print(head(y)); print(head(as.matrix(Covars.PCs))); print(y.lm); print(summary(y.lm));

  ### Run InterPath ###
  ptm <- proc.time() #Start clock
##  vc.mod = InterPath(t(X),y,regions.vs2,cores = cores)
#  vc.mod = MAPITR(X,y,regions,OpenMP=TRUE)
  vc.mod = MAPITR(X,y,regions,Covariates=as.matrix(Covars.PCs),OpenMP=TRUE)
  proc.time() - ptm #Stop clock
  
#  ### Apply Davies Exact Method ###
#  vc.ts = vc.mod$Est
#  names(vc.ts) = colnames(regions.vs2)
#  
#  pvals = c()
#  for(i in 1:length(vc.ts)){
#      lambda = sort(vc.mod$Eigenvalues[,i],decreasing = T)
#      Davies_Method = davies(vc.mod$Est[i], lambda = lambda, acc=1e-8)
#      pvals[i] = 2*min(1-Davies_Method$Qq,Davies_Method$Qq)
#      names(pvals)[i] = names(vc.ts[i])
#  }
  pvals <- vc.mod$Results[,2]
  names(pvals) <- vc.mod$Results[,1]

  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Find power for the first group of SNPs ###
#  Pthwys_1 = colnames(regions)[s1]
  Pthwys_1 = regions[s1,1]
  
  ### Find power for the second group of SNPs ###
#  Pthwys_2 = colnames(regions)[s2]
  Pthwys_2 = regions[s2,1]
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  snps = unlist(gene_list[c(s1,s2,s3)])
  
  ### Save Results ###
  pval_mat[,Counter1] = pvals
  G1_genes[,Counter1] = Pthwys_1
  G2_genes[,Counter1] = Pthwys_2
##  G1_snps[,Counter1] = lapply(gene_list[s1], sum)
##  G2_snps[,Counter1] = lapply(gene_list[s2], sum)
##  Add_snps[,Counter1] = lapply(gene_list[s3], sum)
  stats <- rbind(stats, c(length(unlist(gene_list[s1])), length(unlist(gene_list[s2])), length(unlist(gene_list[s3]))));
  stats <- rbind(stats, c(length(snps.s1.use.total), length(snps.s2.use.total), "NA"))

  ### Report Status ###
  cat("Completed Dataset", j, "\n", sep = " ")

  Counter1 <- Counter1 + 1;
}

  colnames(pval_mat) <- rounds.start:(rounds.start+9)	
  colnames(G1_genes) <- rounds.start:(rounds.start+9)	
  colnames(G2_genes) <- rounds.start:(rounds.start+9)	

  write.table(pval_mat, paste(Output1.File, ".Results.pValues.txt", sep=""), quote=FALSE, col.name=TRUE, row.name=TRUE);
  write.table(G1_genes, paste(Output1.File, ".Results.nCausal1.txt", sep=""), quote=FALSE, col.name=TRUE, row.name=FALSE);
  write.table(G2_genes, paste(Output1.File, ".Results.nCausal2.txt", sep=""), quote=FALSE, col.name=TRUE, row.name=FALSE);
  write.table(stats, paste(Output1.File, ".Results.statistics.txt", sep=""), quote=FALSE, col.name=FALSE, row.name=FALSE);

##Save final Results
#Final = list(pval_mat,G1_genes,G2_genes)

#### Save the Results ###
#file = "~/Desktop/Simulations/Results/Sim_rho5_S4.RData"
#save(Final, file = file)
  
