#!/bin/sh

###20180605 -- FORCE

#20180605
#TOC

 - 20180605: Vs1


#20180605
#Conda environment setup/details/etc
#See `/users/mturchin/RamachandranLab.CCV_General.Workbook.vs1.sh` for initial setup 

conda create -n FORCE
source activate FORCE
conda update conda
conda install R perl java-jdk
conda install git
conda install plink bedtools vcftools vcftools bwa samtools picard gatk imagemagick gnuplot eigensoft tabix
conda install r-base r-devtools r-knitr r-testthat r-cairo r-ashr r-rcolorbrewer r-essentials r-extrafont fonts-anaconda
conda install eigen boost gcc
##20180311 (From MultiEthnicGWAS)
##conda install flashpca -- failed
#conda install eigen
##conda install spectra -- installed a Python package called spectra, not what was intended
#conda install boost
##conda install libgomp -- failed
#conda install gcc
conda install r-doParallel r-Rcpp r-RcppArmadillo r-RcppParallel r-CompQuadForm r-Matrix r-MASS r-truncnorm
#20180611
conda install armadillo fortran
conda install r-data.table r-bigmemory
conda install julia
##conda install r-coop -- failed, install via 'install.packages("coop")'
##conda install r-rbenchmark -- failed, install via 'install.packages("rbenchmark")'
##conda install r-cpgen -- failed, install via 'install.packages("cpgen")'
conda install docker


#20180607
#Vs1

mkdir /users/mturchin/LabMisc/RamachandranLab/FORCE
mkdir /users/mturchin/LabMisc/RamachandranLab/FORCE/Vs1
mkdir /users/mturchin/LabMisc/RamachandranLab/FORCE/Vs1/Analyses
cd /users/mturchin/LabMisc/RamachandranLab/FORCE; git clone https://github.com/mturchin20/Ramachandran_FORCE
mkdir /users/mturchin/data/mturchin/Broad
mkdir /users/mturchin/data/mturchin/Broad/MSigDB
#NOTE -- need to log-in with an e-mail address and then DL directly from webpage, so cannot really do wget setup (I think anyways); so files were downloaded to MacBook Air and then scp'd over
#cd /users/mturchin/data/mturchin/Broad/MSigDB
#wget 
#From MackBook Air
#scp -p /Users/mturchin20/Downloads/*gmt mturchin@ssh.ccv.brown.edu:/users/mturchin/data/mturchin/Broad/MSigDB/.

#NOTE -- copy and pasted 'FORCE_Simulation.R' and 'FORCE.cpp' from Slack via Lorin

cp -p /users/mturchin/LabMisc/RamachandranLab/FORCE/FORCE.cpp /users/mturchin/LabMisc/RamachandranLab/FORCE/Vs1/FORCE.mtEdits.vs1.cpp
cp -p /users/mturchin/LabMisc/RamachandranLab/FORCE/FORCE_Simulation.R /users/mturchin/LabMisc/RamachandranLab/FORCE/Vs1/FORCE.mtEdits.vs1.R
cp -p /users/mturchin/LabMisc/RamachandranLab/FORCE/FORCE_Simulation.R /users/mturchin/LabMisc/RamachandranLab/FORCE/Vs1/FORCE.Source.vs1.R

#See /users/mturchin/.bash_profile for how eventually got the below g++ call to work/function to get the stand-alone armadillo c++ code going
g++ Test1.cpp -o Test1.out  -O2 /gpfs/data/sramacha/mturchin/miniconda2/envs/FORCE/lib/libgfortran.so.4 /gpfs/data/sramacha/mturchin/miniconda2/envs/FORCE/lib/libarmadillo.so.8 -larmadillo -lopenblas -larpack -lgfortran -L/gpfs/data/sramacha/mturchin/miniconda2/envs/FORCE/lib

mkdir /users/mturchin/data/mturchin/20180611_annovar
mkdir /users/mturchin/data/mturchin/20180611_annovar/humandb
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb -webfrom annovar refGene /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg19
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb genomicSuperDups /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg19
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb -webfrom annovar snp138 /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg19
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb -webfrom annovar avsift /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg19
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb -webfrom annovar ljb26_all /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg19
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb -webfrom annovar esp6500_all /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg19
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb -webfrom annovar 1000g2012apr /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg19
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb phastConsElements46way /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg19
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb -webfrom annovar dbnsfp33a /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg19
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb tfbsConsSites /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg19
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb -webfrom annovar refGene /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg18
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb genomicSuperDups /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg18
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb -webfrom annovar snp138 /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg18
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb -webfrom annovar avsift /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg18
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb -webfrom annovar ljb26_all /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg18
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb -webfrom annovar esp6500_all /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg18
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb -webfrom annovar 1000g2012apr /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg18
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb phastConsElements44way /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg18
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb -webfrom annovar dbnsfp33a /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg18
perl /users/mturchin/Software/20180611_annovar/annotate_variation.pl --downdb tfbsConsSites /users/mturchin/data/mturchin/20180611_annovar/humandb/ -build hg18

#   3396 African
# 442688 British (Ran4000)
#   4519 Caribbean
#   1574 Chinese
#   5951 Indian
###  13213 Irish
#   1837 Pakistani
#   1662 Prefer_not_to_answer

UKBioBankPops=`echo "African;African;Afr British;British;Brit British;British.Ran4000;Brit4k Caribbean;Caribbean;Carib Chinese;Chinese;Chi Indian;Indian;Indn Irish;Irish;Irish Pakistani;Pakistani;Pkstn"`;

// [[Rcpp::export]]
arma::mat GetLinearKernel(arma::mat X){
    double p = X.n_rows;
    return X.t()*X/p;
}

ind = 2e3; nsnp = 1e4;
maf <- 0.05 + 0.45*runif(nsnp)
X   <- (runif(ind*nsnp) < maf) + (runif(ind*nsnp) < maf)
X   <- matrix(as.double(X),ind,nsnp,byrow = TRUE)
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)

for (i in c(1e4, 1e5, 1e6)) {
	print(i);
	ind = 2e3; nsnp = i;
	maf <- 0.05 + 0.45*runif(nsnp);
	X   <- (runif(ind*nsnp) < maf) + (runif(ind*nsnp) < maf);
	X   <- matrix(as.double(X),ind,nsnp,byrow = TRUE);
	ptm <- proc.time();
	X.cov <- X %*% t(X);
	print(proc.time() - ptm);
}
set.seed(123); for (i in c(1e4, 1e5, ncol(Data3))) {
	print(i);
	Data3.sub <- Data3[,sample(ncol(Data3),i)];
	ptm <- proc.time();
	Data3.sub.cov <- as.matrix(Data3.sub) %*% t(as.matrix(Data3.sub));
	print(proc.time() - ptm);
}
sourceCpp("/users/mturchin/LabMisc/RamachandranLab/FORCE/FORCE.cpp"); set.seed(123); for (i in c(1e4, 1e5, ncol(Data3))) {
	print(i);
	Data3.sub <- Data3[,sample(ncol(Data3),i)];
	ptm <- proc.time();
	Data3.sub.cov.arma <- GetLinearKernel(t(as.matrix(Data3.sub))); 
	print(proc.time() - ptm);
}
library("coop"); set.seed(123); for (i in c(1e4, 1e5, ncol(Data3))) {
	print(i);
	Data3.sub <- Data3[,sample(ncol(Data3),i)];
	ptm <- proc.time();
	Data3.sub.cov.coop <- covar(t(Data3.sub)); 
	print(proc.time() - ptm);
}
library("cpgen"); set.seed(123); for (i in c(1e4, 1e5, 2e5, 3e5, ncol(Data3))) {
	print(i);
	Data3.sub <- as.matrix(Data3[,sample(ncol(Data3),i)]);
	ptm <- proc.time();
	Data3.sub.cov.coop <- ccov(t(Data3.sub)); 
	print(proc.time() - ptm);
}
set.seed(123); for (i in c(1e4, 1e5, 2e5, 3e5, ncol(Data3))) {
	print(i);
	Data3.sub <- as.matrix(Data3[,sample(ncol(Data3),i)]);
	ptm <- proc.time();
	Data3.sub.cov.coop <- 1/nrow(Data3.sub) * tcrossprod(as.matrix(Data3.sub)) 
	print(proc.time() - ptm);
}
system.time({ x <- replicate(5e3, rnorm(5e3)); tcrossprod(x) })
rbenchmark::benchmark(covar(t(X)), (X %*% t(X)) / nrow(X), replications=10, columns=c("test", "replications", "elapsed", "relative"))
rbenchmark::benchmark(ccov(t(X)), (X %*% t(X)) / nrow(X), replications=10, columns=c("test", "replications", "elapsed", "relative"))

join <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim | awk '{ print $2 "\t" $5 "\t" $6 }' | sort -k 1,1 | uniq) <(zcat /users/mturchin/data/mturchin/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz | awk '{ print $3 "\t" $4 "\t" $5 }' | sort -k 1,1 | uniq) | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.CompHRC.txt.gz
zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.CompHRC.txt.gz | awk '{ if ((($2 != $4) && ($3 != $5)) && (($2 != $5) && ($3 != $4))) { print $1 } }' > /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.CompHRC.flipSNPs.rsIDs 

mkdir /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation
mkdir /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20

for i in {1..22}; do
	echo $i;

	plink --bfile /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop --chr $i --recode vcf --out /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chr${i}_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop
	/users/mturchin/Software/vcftools_0.1.13/bin/vcf-sort /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chr${i}_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.vcf | bgzip -c > /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chr${i}_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz
	rm /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chr${i}_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.vcf /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chr${i}_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.log
done	

#MacBook Air
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chr*_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/ 

cd /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20
wget
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_1.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_2.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_3.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_4.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_5.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_6.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_7.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_8.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_9.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_10.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_11.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_12.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_13.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_14.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_15.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_16.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_17.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_18.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_19.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_20.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_21.zip -p 
7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_22.zip -p 

wget https://imputationserver.sph.umich.edu/share/results/bb235fb15db46210e9321baa1354e409/chr_1.zip https://imputationserver.sph.umich.edu/share/results/4c64843d6103fe5b810f39ad8811c48a/chr_1.log
7za x /users/mturchin/data/POPRES/NHGRI/POPRES/phs000145v2/p2/mturchin20/UKBHeightRspnd/Imputation/UK/chr_1.zip -p'VvnB1J2k(KhvOk'

#From http://csg.sph.umich.edu/abecasis/mach/tour/imputation.html
#r2 cutoff .3
for i in {1..22}; do
	echo $i;
	
	ln -s /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr${i}.dose.vcf.gz  /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chr${i}_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.vcf.gz
	ln -s /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr${i}.dose.vcf.gz.tbi  /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chr${i}_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.vcf.gz.tbi
	ln -s /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr${i}.info.gz  /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chr${i}_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.info.gz

done

zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chr*_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.info.gz | awk '{ if ($7 > .3) { 
print $2 } }' | grep -v SNP | sort | uniq -u | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.info.r2gt3.noDups.rsIDs.gz
cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim | awk '{ print $2 }' > /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim.rsIDs
join <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.info.r2gt3.noDups.rsIDs.gz | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim.rsIDs | sort -k 1,1) > /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim.ImptHRC.info.r2gt3.noDups.rsIDs 

zcat /users/mturchin/data/POPRES/NHGRI/POPRES/phs000145v2/p2/mturchin20/UKBHeightRspnd/Imputation/UK/POPRES_Genotypes_QC2_v2.United_Kingdom.chr*.flip.hg19.ImptHRC.info.gz | awk '{ if ($7 > .3) { print $1 } } ' | sed 's/:/_/g' | grep -v SNP | sort | uniq -u | gzip > /users/mturchin/data/POPRES/NHGRI/POPRES/phs000145v2/p2/mturchin20/UKBHeightRspnd/Imputation/UK/POPRES_Genotypes_QC2_v2.United_Kingdom.chrAll.flip.hg19.ImptHRC.info.r2gt3.noDups.ChrBP.gz

for i in {1..22}; do
        echo $i
	sbatch -t 72:00:00 --mem 10g --qos=Normal -o blah -e blah --comment "ukbPpRs1kGMatch $i" <(echo -e '#!/bin/sh';
	echo -e "\nvcftools --gzvcf /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chr${i}_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.vcf.gz --plink --snps /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim.ImptHRC.info.r2gt3.noDups.rsIDs --out /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chr${i}_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp"; \
	echo -e "\nplink --file /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chr${i}_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp --make-bed --out /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chr${i}_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp"; \
	echo -e "\nrm /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chr${i}_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.ped /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chr${i}_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.map";)	
#	rm /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chr${i}_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp* 
 
	echo -e "\nvcftools --gzvcf /users/mturchin/data/POPRES/NHGRI/POPRES/phs000145v2/p2/mturchin20/UKBHeightRspnd/Imputation/UK/POPRES_Genotypes_QC2_v2.United_Kingdom.chr${i}.flip.hg19.ImptHRC.dose.vcf.gz --plink --snps /users/mturchin/data/POPRES/NHGRI/POPRES/phs000145v2/p2/mturchin20/UKBHeightRspnd/Imputation/POPRES_Genotypes_QC2_v2.UKItalyPortugal.chrAll.flip.hg19.ImptHRC.info.r2gt3.ukb1kGMatch.ChrSemiBP --out /users/mturchin/data/POPRES/NHGRI/POPRES/phs000145v2/p2/mturchin20/UKBHeightRspnd/Imputation/UK/POPRES_Genotypes_QC2_v2.United_Kingdom.chr${i}.flip.hg19.ImptHRC.dose.plinkTemp"; \ 
	echo -e "\nplink --file /users/mturchin/data/POPRES/NHGRI/POPRES/phs000145v2/p2/mturchin20/UKBHeightRspnd/Imputation/UK/POPRES_Genotypes_QC2_v2.United_Kingdom.chr${i}.flip.hg19.ImptHRC.dose.plinkTemp --make-bed --out /users/mturchin/data/POPRES/NHGRI/POPRES/phs000145v2/p2/mturchin20/UKBHeightRspnd/Imputation/UK/POPRES_Genotypes_QC2_v2.United_Kingdom.chr${i}.flip.hg19.ImptHRC.dose.plinkTemp"; \ 
	echo -e "\nrm /users/mturchin/data/POPRES/NHGRI/POPRES/phs000145v2/p2/mturchin20/UKBHeightRspnd/Imputation/UK/POPRES_Genotypes_QC2_v2.United_Kingdom.chr${i}.flip.hg19.ImptHRC.dose.plinkTemp.ped /users/mturchin/data/POPRES/NHGRI/POPRES/phs000145v2/p2/mturchin20/UKBHeightRspnd/Imputation/UK/POPRES_Genotypes_QC2_v2.United_Kingdom.chr${i}.flip.hg19.ImptHRC.dose.plinkTemp.map";)
        rm /users/mturchin/data/POPRES/NHGRI/POPRES/phs000145v2/p2/mturchin20/UKBHeightRspnd/Imputation/UK/POPRES_Genotypes_QC2_v2.United_Kingdom.chr${i}.flip.hg19.ImptHRC.dose.plinkTemp* /users/mturchin/data/POPRES/NHGRI/POPRES/phs000145v2/p2/mturchin20/UKBHeightRspnd/Imputation/Italy/POPRES_Genotypes_QC2_v2.Italy.chr${i}.flip.hg19.ImptHRC.dose.plinkTemp* /users/mturchin/data/POPRES/NHGRI/POPRES/phs000145v2/p2/mturchin20/UKBHeightRspnd/Imputation/Portugal/POPRES_Genotypes_QC2_v2.Portugal.chr${i}.flip.hg19.ImptHRC.dose.plinkTemp*
done












for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep British | grep -v Ran`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`

        echo $ancestry1 $ancestry2 /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.MergeList.Vs2.txt

        cat /dev/null > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.MergeList.Vs2.txt

        for chr in {2..22}; do
                echo "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${chr}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bed /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${chr}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${chr}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.fam" >> /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.MergeList.Vs2.txt
        done

done

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep British | grep -v Ran`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`

        echo $ancestry1 $ancestry2 /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.MergeList.Vs2.txt

        sbatch -t 24:00:00 --mem 100g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.MergeList.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.MergeList.error --comment "$ancestry1 $ancestry2" <(echo -e '#!/bin/sh'; echo -e "\nplink --bfile /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr1_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop --merge-list /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.MergeList.Vs2.txt --recodeAD --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop")

done

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`

        for i in {1..22}; do
                echo $ancestry1 $ancestry2 $i

                sbatch -t 72:00:00 --mem 20g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${i}_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.clump.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${i}_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.clump.error --comment "$ancestry1 $ancestry2 $i" <(echo -e '#!/bin/sh'; \
                echo -e "\nplink --bfile /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop --indep-pairwise 1000 50 .1 --exclude range /users/mturchin/Software/flashpca/exclusion_regions_hg19.txt --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop"; \
		echo -e "\nplink --bfile /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop --extract /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.prune.in --make-bed --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned") 

        done
done

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`

	rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr*_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.prune.* /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr*_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.log /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr*_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.log /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr*_v2.African.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.clump.output /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr*_v2.African.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.clump.error
        for chr in {2..22}; do
                echo "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${chr}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.bed /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${chr}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.bim /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${chr}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.fam" 
        done > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.MergeList.Vs2.txt
done

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`

        echo $ancestry1 $ancestry2 /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.MergeList.Vs2.txt

        sbatch -t 24:00:00 --mem 100g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.MergeList.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.MergeList.error --comment "$ancestry1 $ancestry2" <(echo -e '#!/bin/sh'; echo -e "\nplink --bfile /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr1_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned --merge-list /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.MergeList.Vs2.txt --recodeAD --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned")

done

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`

        echo $ancestry1 $ancestry2 /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.MergeList.Vs2.txt

	gzip -f /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.raw
	rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.bed /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.bim /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.fam /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.log

done

#From https://cran.r-project.org/web/packages/coop/vignettes/coop.pdf
library("data.table"); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane \'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);\'', header=T);
library("data.table"); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.raw.gz | perl -lane \'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);\'', header=T);

ptm <- proc.time();
Data3.cov2 <- 1/nrow(Data3) * tcrossprod(scale(as.matrix(Data3), TRUE, FALSE))
print(proc.time() - ptm);





for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep British | grep Ran4000`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
        ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
		
	echo $pheno1 $ancestry1 $ancestry2 $ancestryr3

	if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses ]; then
		mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses
	fi
	if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/EpiPath ]; then
		mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/EpiPath
	fi

done





for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep British | grep Ran4000`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
        ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
		
	echo $pheno1 $ancestry1 $ancestry2 $ancestryr3

	if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses ]; then
		mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses
	fi
	if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/EpiPath ]; then
		mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/EpiPath
	fi

done



~~~
#20180615
(FORCE) [  mturchin@node621  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.CompHRC.flipSNPs.rsIDs | wc
      0       0       0



~~~







