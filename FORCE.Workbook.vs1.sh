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
##conda install docker -- failed


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

#20180618 NOTE -- overall impression from the first round of these explorations below: a) base R BLAS is not good (as is well known by this poitn) b) load R from Oscar since conda R currently does not implement any of the better BLAS libraries (eg OpenBLAS or MKL) c) tcrossprod() outperforms GetLinearKernal() (apparently), and in generally appears to be best option d) doing the for loop thing with tcrossprod() though leads to a seg fault by the second loop, not sure why e) ccov and covar don't seem to either really matter or make much of a difference (at least top-level enough that once I found out tcrossprod() and the Oscar R combination worked, I stuck with that; I don't think I checked whether ccov/covar did better on the Oscar R module, so that may make a difference tbh) f) covar needs full data (no missing genotypes) I think whereas tcrossprod() seems to handle NAs somehow/someway g) GetLinearKernel() and tcrossprod() do not handle NAs, eg by having them it just makes every resulting new matrix entry NA; need to use imputed data or data removed of any missing genotypes h) the amount of memory for getting the tcrossprod() result on the African raw dataset is about ~18-19gb 
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
sourceCpp("/users/mturchin/LabMisc/RamachandranLab/FORCE/FORCE.cpp"); set.seed(123); for (i in c(1e4, 1e5, 2e5, 3e5, ncol(Data3))) {
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
	Data3.sub.cov.coop <- 1/nrow(Data3.sub) * tcrossprod(scale(as.matrix(Data3.sub), TRUE, FALSE)) 
	print(proc.time() - ptm);
}
system.time({ x <- replicate(5e3, rnorm(5e3)); tcrossprod(x) })
rbenchmark::benchmark(covar(t(X)), (X %*% t(X)) / nrow(X), replications=10, columns=c("test", "replications", "elapsed", "relative"))
rbenchmark::benchmark(ccov(t(X)), (X %*% t(X)) / nrow(X), replications=10, columns=c("test", "replications", "elapsed", "relative"))

#From https://cran.r-project.org/web/packages/coop/vignettes/coop.pdf
library("data.table"); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane \'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);\'', header=T);
#library("data.table"); Data4 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.pruned.raw.gz | perl -lane \'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);\'', header=T);
#library("Rcpp"); library("RcppArmadillo"); sourceCpp("/users/mturchin/LabMisc/RamachandranLab/FORCE/FORCE.cpp"); 

ptm <- proc.time();
Data3.cov2 <- 1/nrow(Data3) * tcrossprod(scale(as.matrix(Data3), FALSE, FALSE))
print(proc.time() - ptm);
rm(Data3.cov2);
ptm <- proc.time();
Data4.cov2 <- 1/nrow(Data4) * tcrossprod(scale(as.matrix(Data4), FALSE, FALSE))
print(proc.time() - ptm);
rm(Data4.cov2);

ptm <- proc.time();
Data3.cov3 <- GetLinearKernel(t(as.matrix(Data3)));
print(proc.time() - ptm);
rm(Data3.cov2);
ptm <- proc.time();
Data4.cov3 <- GetLinearKernel(t(as.matrix(Data4)));
print(proc.time() - ptm);
rm(Data4.cov2);

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

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep British | grep -v Ran4000`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`

        echo $ancestry1 $ancestry2 
	
	if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation ]; then
		mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation
	fi
	if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20 ]; then
		mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20
	fi

        for i in {1..22}; do
        	plink --bfile /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop --recode vcf --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop
		/users/mturchin/Software/vcftools_0.1.13/bin/vcf-sort /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.vcf | bgzip -c > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz
		rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.vcf /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.log
	done

done

#MacBook Air
#mkdir /Volumes/NO\ NAME/African /Volumes/NO\ NAME/British /Volumes/NO\ NAME/British.Ran4000 /Volumes/NO\ NAME/Caribbean; mkdir /Volumes/NO\ NAME/Chinese /Volumes/NO\ NAME/Indian /Volumes/NO\ NAME/Irish /Volumes/NO\ NAME/Pakistani 
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chr*_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/African 
scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/British/British/mturchin20/ukb_chr*_v2.British.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/British 
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran4000/mturchin20/ukb_chr*_v2.British.Ran4000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/British.Ran4000
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/mturchin20/ukb_chr*_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/Caribbean 
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/Chinese/Chinese/mturchin20/ukb_chr*_v2.Chinese.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/Chinese 
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/Indian/Indian/mturchin20/ukb_chr*_v2.Indian.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/Indian 
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/Irish/Irish/mturchin20/ukb_chr*_v2.Irish.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/Irish 
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/Pakistani/Pakistani/mturchin20/ukb_chr*_v2.Pakistani.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/Pakistani 

cd /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20
mkdir /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/Part1; cd /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/Part1; wget https://imputationserver.sph.umich.edu/share/results/56d545ecfadae67f277a19213300b028/qcreport.html https://imputationserver.sph.umich.edu/share/results/a2a92469a201d8b506f6eb8bac463927/chr_1.zip https://imputationserver.sph.umich.edu/share/results/d22192341e5b80d9d3de4f43a468e045/chr_2.zip https://imputationserver.sph.umich.edu/share/results/fd562d2058434b03ce4e6e3bd2cbfc39/chr_3.zip https://imputationserver.sph.umich.edu/share/results/3d7be1309bff37a51e95e82e9052f851/chr_4.zip https://imputationserver.sph.umich.edu/share/results/1da2b3deab2c57f1119933744786a40d/chr_5.zip https://imputationserver.sph.umich.edu/share/results/750e64529a13589d69c6e5a9c61246d/statistics.txt https://imputationserver.sph.umich.edu/share/results/280912b3f0d2c522122e6b59405e5292/chr_1.log https://imputationserver.sph.umich.edu/share/results/6072bebce7e8d6126d7dbd26ebd30a7c/chr_2.log https://imputationserver.sph.umich.edu/share/results/11a8a3ef083ef05b4948e381de00a254/chr_3.log https://imputationserver.sph.umich.edu/share/results/234f826498e5259802adc3dcd1c917ac/chr_4.log https://imputationserver.sph.umich.edu/share/results/8408aa6a10e95597e3b73b4bec522738/chr_5.log
mkdir /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/Part2; cd /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/Part2; wget https://imputationserver.sph.umich.edu/share/results/5cf2b405a9b41f8c83f3ab109b22037f/qcreport.html https://imputationserver.sph.umich.edu/share/results/8f313c0abfbd6c95ea458a4a0b59f31/chr_10.zip https://imputationserver.sph.umich.edu/share/results/cceff332dc42f913b66199d00d3c6aa5/chr_11.zip https://imputationserver.sph.umich.edu/share/results/ee3af88208dafd3e2520ab9ea5a3b5ec/chr_12.zip https://imputationserver.sph.umich.edu/share/results/c13e4f1c337d4a58bcfad95ffe6fef7f/chr_6.zip https://imputationserver.sph.umich.edu/share/results/5fc5760498f5920b705d83c3576a4b19/chr_7.zip https://imputationserver.sph.umich.edu/share/results/195043be82e442dbed73e19f7b59b8df/chr_8.zip https://imputationserver.sph.umich.edu/share/results/889f456c90e1c5b78833c47acf540534/chr_9.zip https://imputationserver.sph.umich.edu/share/results/b1344d2092d9d918d70a2f4cbcc5bdcc/statistics.txt https://imputationserver.sph.umich.edu/share/results/ab867ba77fd62a3b1c9e83b3daaf14da/chr_10.log https://imputationserver.sph.umich.edu/share/results/992d1fd9b515fa5d491d9f912caedf91/chr_11.log https://imputationserver.sph.umich.edu/share/results/671732668b1d779c495087f8ac54d76d/chr_12.log https://imputationserver.sph.umich.edu/share/results/46c0e94d1de168e495a4798f12d2ff83/chr_6.log https://imputationserver.sph.umich.edu/share/results/3bc3abf96c478d5bacf3d7295d55db14/chr_7.log https://imputationserver.sph.umich.edu/share/results/7583d120cdb5a4286ad96278793cca90/chr_8.log https://imputationserver.sph.umich.edu/share/results/9ce5c8948e7bc075f8c37dad5fc7fc9e/chr_9.log
mkdir /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/Part3; cd /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/Part3; wget https://imputationserver.sph.umich.edu/share/results/3f3b00495ab6bde28a35c6c7e65eb455/qcreport.html https://imputationserver.sph.umich.edu/share/results/36fd7ae63aefe201b29de33b380fab31/chr_13.zip https://imputationserver.sph.umich.edu/share/results/13c53d8fae0e47de98815673fd173650/chr_14.zip https://imputationserver.sph.umich.edu/share/results/2b0585da89059cf6dce1697c1b6366aa/chr_15.zip https://imputationserver.sph.umich.edu/share/results/6777162f9eaa2c59bd3254f29bcc9050/chr_16.zip https://imputationserver.sph.umich.edu/share/results/23fe40a43c9e74959059b01885543f19/chr_17.zip https://imputationserver.sph.umich.edu/share/results/74384e9b46fcd5f88bbcc5cdf40710fc/chr_18.zip https://imputationserver.sph.umich.edu/share/results/dcd89dbf47199ab478637446a40e1833/chr_19.zip https://imputationserver.sph.umich.edu/share/results/7f364000193dcd6eee14f0552c6a01a6/chr_20.zip https://imputationserver.sph.umich.edu/share/results/f74aa2ea6480625bf25329e30cea34e2/chr_21.zip https://imputationserver.sph.umich.edu/share/results/a0d6ab789b2d5d7f8c416e58fff23c30/chr_22.zip https://imputationserver.sph.umich.edu/share/results/f20b824f76d4fa922bb94077690d6664/statistics.txt https://imputationserver.sph.umich.edu/share/results/6f714f1bea5629d7398c8e7847af0b49/chr_13.log https://imputationserver.sph.umich.edu/share/results/d577d17be363d87a14177b1086c4c60d/chr_14.log https://imputationserver.sph.umich.edu/share/results/898fac76a61cbffcc644bcff9d786762/chr_15.log https://imputationserver.sph.umich.edu/share/results/601749419847724c533ad8e0317caeb6/chr_16.log https://imputationserver.sph.umich.edu/share/results/752ad5fe36038183c631d047ca961008/chr_17.log https://imputationserver.sph.umich.edu/share/results/e5f020bc64c1fbe6d15f3e6a086c115d/chr_18.log https://imputationserver.sph.umich.edu/share/results/2baab818142d280d7e4b9ea258da0f7/chr_19.log https://imputationserver.sph.umich.edu/share/results/ff6bfa773a2fc7f2601867d6d445a205/chr_20.log https://imputationserver.sph.umich.edu/share/results/c4a3f75bfaf8c0a9078c26e2483a4cf9/chr_21.log https://imputationserver.sph.umich.edu/share/results/26765d4cbd03c8faa049d1d1af5a289f/chr_22.log

mv /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/Part1/*zip /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20
mv /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/Part2/*zip /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20
mv /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/Part3/*zip /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20
for i in {1..5}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_${i}.zip -p'OU6mc0gyOpFOTf'
done
for i in {6..12}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_${i}.zip -p'Nhs2xRH6}Kig)X'
done
for i in {13..22}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/chr_${i}.zip -p'ECnPNblYi3Xoq0'
done







cd /users/mturchin/data/ukbiobank_jun17/subsets/British/British/Imputation/mturchin20
wget 
for i in {1..22}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/British/British/Imputation/mturchin20/chr_${i}.zip -p'
done
cd /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran4000/Imputation/mturchin20
wget
for i in {1..22}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran4000/Imputation/mturchin20/chr_${i}.zip -p'
done
cd /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20
wget
for i in {1..22}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/chr_${i}.zip -p'
done
cd /users/mturchin/data/ukbiobank_jun17/subsets/Chinese/Chinese/Imputation/mturchin20
wget
for i in {1..22}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/Chinese/Chinese/Imputation/mturchin20/chr_${i}.zip -p'
done
cd /users/mturchin/data/ukbiobank_jun17/subsets/Irish/Irish/Imputation/mturchin20
wget
for i in {1..22}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/Irish/Irish/Imputation/mturchin20/chr_${i}.zip -p'
done
cd /users/mturchin/data/ukbiobank_jun17/subsets/Indian/Indian/Imputation/mturchin20
wget
for i in {1..22}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/Indian/Indian/Imputation/mturchin20/chr_${i}.zip -p'
done
cd /users/mturchin/data/ukbiobank_jun17/subsets/Pakistani/Pakistani/Imputation/mturchin20
wget
for i in {1..22}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/Pakistani/Pakistani/Imputation/mturchin20/chr_${i}.zip -p'
done












#From http://csg.sph.umich.edu/abecasis/mach/tour/imputation.html
#r2 cutoff .3
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
        ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
		
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3

	for i in {1..22}; do
		echo $i;

		ln -s /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/chr${i}.dose.vcf.gz /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.vcf.gz	
		ln -s /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/chr${i}.dose.vcf.gz.tbi /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.vcf.gz.tbi	
		ln -s /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/chr${i}.info.gz /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.info.gz	

	done
done

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
        ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
		
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3

	zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr*_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.info.gz | awk '{ if ($7 > .3) { print $1 } }' | sed 's/:/_/g' | grep -v SNP | sort | uniq -u | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.info.r2gt3.noDups.ChrBPs.gz
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr*_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim | awk '{ print $1 "_" $4 }' | sort | uniq > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim.ChrBPs
	join <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.info.r2gt3.noDups.ChrBPs.gz | sort) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim.ChrBPs | sort) | sed 's/_/:/g' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim.ImptHRC.info.r2gt3.noDups.ChrBPs

done

#From https://www.biostars.org/p/46060/ & https://sourceforge.net/p/vcftools/mailman/message/29115811/
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
        ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
		
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3

	for i in {1..22}; do
		echo $i
		sbatch -t 72:00:00 --mem 10g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.error --comment "$ancestry1 $ancestry2 $i" <(echo -e '#!/bin/sh';
		echo -e "\nvcftools --gzvcf /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.vcf.gz --plink-tped --snps /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim.ImptHRC.info.r2gt3.noDups.ChrBPs --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp"; \
		echo -e "\nplink --tfile /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp --make-bed --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp"; \
		echo -e "\nrm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.tped /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.tfam";)	
#		rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp* 
#		echo -e "\nplink --tfile /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp --exclude /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.missnp --make-bed --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp"; \
	done
done 

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`

        echo $ancestry1 $ancestry2 /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.MergeList.Vs2.txt

        for chr in {2..22}; do
                echo "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${chr}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.bed /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${chr}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.bim /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${chr}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.fam" 
        done > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.MergeList.Vs2.txt

done

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`

        echo $ancestry1 $ancestry2 /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.MergeList.Vs2.txt

        sbatch -t 24:00:00 --mem 100g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.MergeList.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.MergeList.error --comment "$ancestry1 $ancestry2" <(echo -e '#!/bin/sh'; echo -e "\nplink --bfile /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr1_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp --merge-list /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.MergeList.Vs2.txt --recodeAD --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose")

done

~for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
~        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
~        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
~
~        for i in {1..22}; do
~                echo $ancestry1 $ancestry2 $i
~
~                sbatch -t 72:00:00 --mem 20g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.clump.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.clump.error --comment "$ancestry1 $ancestry2 $i" <(echo -e '#!/bin/sh'; \
~                echo -e "\nplink --bfile /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose --indep-pairwise 1000 50 .1 --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose"; \
~		echo -e "\nplink --bfile /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose --extract /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.prune.in --make-bed --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned") 
~
~        done
~done
~
~for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
~        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
~        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
~
~	rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr*_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.prune.* /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr*_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.log /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr*_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned.log /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr*_v2.African.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.clump.output /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr*_v2.African.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.clump.error
~        for chr in {2..22}; do
~                echo "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${chr}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned.bed /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${chr}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned.bim /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${chr}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned.fam" 
~        done > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned.MergeList.Vs2.txt
~done
~
~for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
~        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
~        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
~
~        echo $ancestry1 $ancestry2 /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned.MergeList.Vs2.txt
~
~        sbatch -t 24:00:00 --mem 100g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned.MergeList.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned.MergeList.error --comment "$ancestry1 $ancestry2" <(echo -e '#!/bin/sh'; echo -e "\nplink --bfile /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr1_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned --merge-list /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned.MergeList.Vs2.txt --recodeAD --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned")
~
~done






for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`

        echo $ancestry1 $ancestry2 /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned.MergeList.Vs2.txt

	gzip -f /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw
#	gzip -f /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned.raw
	rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned.bed /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned.bim /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned.fam /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned.log

done








##zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' | head -n 10
##join <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' | sort) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim |  
##cat file1 | R -q -e "Data1 <- read.table(file('stdin'), header=F); Data2 <- Data1[,1]; Data3 <- Data1[,2:ncol(Data1)]; colnames(Data2) <- c(\"SNPID\"); colnames(Data3) <- c(\"SNPID\", ....\"POS\"?....); Data4 <- merge(Data2, Data3, by=\"SNPID\"); write.table(Data4, quote=FALSE, row.name=FALSE, col.name=FALSE);" | grep -v \> | 
zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.gz | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' > /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs
cmp <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs | sed 's/_/ /g' | awk '{ print $1 }') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | awk '{ print $2 }') | wc
paste <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs | sed 's/_/ /g' | awk '{ print $1 }') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | awk '{ print $2 }') | awk '{ if ($1 == $2) { print $0 } } ' | wc

#/users/mturchin/data/mturchin/Broad/MSigDB
#/users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt
#perl /mnt/gluster/home/mturchin20/Software/20170208DL_20160204Vs_annovar/table_annovar.pl /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs2/bmass.GlobalLipids2013.Vs2.NewSNPs.SNPs.vs1.AnnovarFormat /mnt/gluster/home/mturchin20/Software/20170208DL_20160204Vs_annovar/humandb/ -buildver hg18 -out /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs2/bmass.GlobalLipids2013.Vs2.NewSNPs.SNPs.vs1.AnnovarFormat.TableAnnovar -remove -protocol refGene,phastConsElements44way,tfbsConsSites,gwasCatalog,snp129,snp138,1000g2014oct_all,esp6500siv2_all,exac03,dbnsfp33a -operation g,r,r,r,f,f,f,f,f,f -nastring NA

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
        ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
		
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3

	if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses ]; then
		mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses
	fi
	if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/FORCE ]; then
		mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/FORCE
	fi

	zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.gz | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs
	paste <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs | sed 's/_/ /g' | awk '{ print $1 }') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | awk '{ print $2 }') | awk '{ if ($1 == $2) { print $0 } } ' | wc
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs | wc
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | wc
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs | awk '{ print $1 "\t" NR }' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDsRowPos

done

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
        ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
		
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3

	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | awk '{ print $1 "\t" $4 "\t" $4 "\t" $5 "\t" $6 }' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat 
	perl /users/mturchin/Software/20180611_annovar/table_annovar.pl /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat /users/mturchin/data/mturchin/20180611_annovar/humandb/ -buildver hg19 -out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar -remove -protocol refGene -operation g -nastring NA
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | perl -F\\t -lane 'if ($. == 1) { %hash1; } if ($hash1{$F[6]}) { push(@{$hash1{$F[6]}}, $F[0] . ":" . $F[1] . "\t" . $F[5] . "," . $F[7] . "," . $F[8]); } else { $hash1{$F[6]} = [($F[0] . ":" . $F[1] . "\t" . $F[5] . "," . $F[7] . "," . $F[8])]; } if (eof()) { foreach my $gene1 (keys %hash1) { foreach my $snp1 (@{$hash1{$gene1}}) { print $gene1, "\t", $snp1; } } };' | sed 's/ /_/g ' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.txt
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.txt | perl -lane 'my @genes1 = split(/;/, $F[0]); my @info1 = split(/,/, $F[2]); my @dists1 = split(/;/, $info1[1]); if (($F[0] =~ m/;/) && ($info1[0] eq "intergenic")) { print $genes1[0], "\t", $F[1], "\t", $info1[0] . "_upstream", ",", $dists1[0], ",", $info1[2]; print $genes1[1], "\t", $F[1], "\t", $info1[0] . "_downstream", ",", $dists1[1], ",", $info1[2]; } elsif ($F[0] =~ m/;/) { print $genes1[0], "\t", $F[1], "\t", $info1[0], ",", $info1[1], ",", $info1[2]; print $genes1[1], "\t", $F[1], "\t", $info1[0], ",", $info[1], ",", $info1[2]; } else { print join("\t", @F); }' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.txt
	join <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDsRowPos | sed 's/_/ /g' | awk '{ print $1 "\t" $3 }' | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.txt | awk '{ print $2 "\t" $1 "\t" $3 }' | sort -k 1,1) | awk '{ print $3 "\t" $1 "\t" $4 "\t" $2 }' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt 

done

#nonsynonymous, exonic, exonic + intronic + UTR, exonic + intronic + UTR + upstream/downstream + 20kb away
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`; ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`; echo $pheno1 $ancestry1 $ancestry2 $ancestry3;

	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep nonsynonymous | sort -k 1,1 | perl -lane 'if ($. == 1) { @gene1; push(@gene1, $F[0]); push(@gene1, $F[3]); } else { if ($F[0] ne $gene1[0]) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); @gene1 = (); push(@gene1, $F[0]); push(@gene1, $F[3]); } else { push(@gene1, $F[3]); } if (eof()) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); } };' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/FORCE/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.NonSyn.txt
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep exonic | sort -k 1,1 | perl -lane 'if ($. == 1) { @gene1; push(@gene1, $F[0]); push(@gene1, $F[3]); } else { if ($F[0] ne $gene1[0]) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); @gene1 = (); push(@gene1, $F[0]); push(@gene1, $F[3]); } else { push(@gene1, $F[3]); } if (eof()) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); } };' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/FORCE/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.Exonic.txt
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep -E 'exonic|intronic|UTR' | sort -k 1,1 | perl -lane 'if ($. == 1) { @gene1; push(@gene1, $F[0]); push(@gene1, $F[3]); } else { if ($F[0] ne $gene1[0]) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); @gene1 = (); push(@gene1, $F[0]); push(@gene1, $F[3]); } else { push(@gene1, $F[3]); } if (eof()) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); } };' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/FORCE/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus.txt
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep -E 'exonic|intronic|UTR|upstream|downstream' | sort -k 1,1 | perl -lane 'my @info1 = split(/,/, $F[2]); if ($info1[0] =~ m/intergenic/) { my @dists1 = split(/=/, $info1[1]); if ($dists1[1] <= 20000) { print join("\t", @F); } } else { print join("\t", @F); }' | perl -lane 'if ($. == 1) { @gene1; push(@gene1, $F[0]); push(@gene1, $F[3]); } else { if ($F[0] ne $gene1[0]) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); @gene1 = (); push(@gene1, $F[0]); push(@gene1, $F[3]); } else { push(@gene1, $F[3]); } if (eof()) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); } };' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/FORCE/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus20kb.txt
	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/FORCE/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.NonSyn.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(",", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/FORCE/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.NonSyn.c2.txt
	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/FORCE/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.Exonic.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(",", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/FORCE/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.Exonic.c2.txt
	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/FORCE/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(",", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/FORCE/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus.c2.txt
	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/FORCE/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus20kb.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(",", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/FORCE/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus20kb.c2.txt
done

#	for k in `cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | perl -F\\\t -lane 'print $F[6];' | sort | uniq`; do
#		cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | perl -F\\t -slane 'if ($F[6] eq $gene1) { print $gene1, "\t", $F[0], ":", $F[1]; }' -- -gene1=$k
#	done > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.txt

#From https://blog.rstudio.com/2016/03/29/feather/, https://blog.dominodatalab.com/the-r-data-i-o-shootout/, https://stackoverflow.com/questions/1727772/quickly-reading-very-large-tables-as-dataframes
/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/FORCE/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.NonSyn.c2.txt

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
        ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
		
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3

	R -q -e "library(\"data.table\"); library(\"feather\"); \
	Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane \'if (\$. == 1) { @vals1; for (my \$i = 6; \$i <= \$#F; \$i++) { if (\$F[\$i] =~ m/HET/) { \$PH = 1 } else { push(@vals1, \$i); } } } print join(\"\t\", @F[@vals1]);\'', header=T); \
	Data3.cov <- 1/nrow(Data3) * tcrossprod(as.matrix(Data3)); \
	write_feather(Data3, \"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/FORCE/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.edit.fthr\"); \
	write_feather(Data3.cov, \"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/FORCE/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.edit.cov.fthr\");"
done 

Data3.cov2 <- 1/nrow(Data3) * tcrossprod(scale(as.matrix(Data3), FALSE, FALSE))

library("data.table"); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane \'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);\'', header=T);

done

R -q -e "library(\"data.table\"); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz'); Data3[1:5,1:5]"
R -q -e "library(\"data.table\"); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz'); Data3[1:5,1:5]"








~~~
#20180615
(FORCE) [  mturchin@node621  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.CompHRC.flipSNPs.rsIDs | wc
      0       0       0
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' | wc        
 379088  379088 4722679
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' | head -n 10
rs116720794_T
rs3131972_G
rs12184325_T
rs3131962_G
rs114525117_A
rs3115850_C
rs115991721_G
rs12562034_A
rs116390263_T
rs4040617_G
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | head -n 10
1       1:729632        0       729632  T       C
1       1:752721        0       752721  G       A
1       1:754105        0       754105  T       C
1       1:756604        0       756604  G       A
1       1:759036        0       759036  A       G
1       1:761147        0       761147  C       T
1       1:767096        0       767096  G       A
1       1:768448        0       768448  A       G
1       1:772927        0       772927  T       C
1       1:779322        0       779322  G       A
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' | head -n 10      1:729632_T                                                                                                                                                                                                                      1:752721_G                                                                                                                                                                                                                      1:754105_T                                                                                                                                                                                                                      1:756604_G                                                                                                                                                                                                                      1:759036_A                                                                                                                                                                                                                      1:761147_C
1:767096_G
1:768448_A
1:772927_T
1:779322_G
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' | tail -n 10
22:51161019_T
22:51161093_T
22:51161620_T
22:51162059_A
22:51165664_G
22:51171497_A
22:51175626_G
22:51193629_G
22:51217954_A
22:51224208_A
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | tail -n 10
22      22:51161019     0       51161019        T       C
22      22:51161093     0       51161093        T       C
22      22:51161620     0       51161620        T       C
22      22:51162059     0       51162059        A       G
22      22:51165664     0       51165664        G       A
22      22:51171497     0       51171497        A       G
22      22:51175626     0       51175626        G       A
22      22:51193629     0       51193629        G       A
22      22:51217954     0       51217954        A       G
22      22:51224208     0       51224208        A       G
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' > /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cmp /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | awk '{ print $2 }')
/users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs /dev/fd/63 differ: char 9, line 1
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$paste /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | awk '{ print $2 }') | head -n 10
1:729632_T      1:729632
1:752721_G      1:752721
1:754105_T      1:754105
1:756604_G      1:756604
1:759036_A      1:759036
1:761147_C      1:761147
1:767096_G      1:767096
1:768448_A      1:768448
1:772927_T      1:772927
1:779322_G      1:779322
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cmp <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs | sed 's/_/ /g' | awk '{ print $1 }') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | awk '{ print $2 }') | wc
      0       0       0
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$paste <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs | sed 's/_/ /g' | awk '{ print $1 }') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | awk '{ print $2 }') | awk '{ if ($1 == $2) { print $0 } } ' | wc
 373401  746802 8690114
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs | wc
 373401  373401 5091859
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | wc
 373401 2240406 10930520
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | perl -F\\t -lane 'print $#F;' | sort | uniq -c             
 373402 9
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | head -n 10
KEGG_GLYCOLYSIS_GLUCONEOGENESIS http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_GLYCOLYSIS_GLUCONEOGENESIS ACSS2   GCK     PGK2    PGK1    PDHB    PDHA1   PDHA2   PGM2    TPI1    ACSS1   FBP1    ADH1B   HK2     ADH1C  HK1      HK3     ADH4    PGAM2   ADH5    PGAM1   ADH1A   ALDOC   ALDH7A1 LDHAL6B PKLR    LDHAL6A ENO1    PKM2    PFKP    BPGM    PCK2    PCK1    ALDH1B1 ALDH2   ALDH3A1 AKR1A1  FBP2    PFKM    PFKL    LDHC    GAPDH   ENO3   ENO2     PGAM4   ADH7    ADH6    LDHB    ALDH1A3 ALDH3B1 ALDH3B2 ALDH9A1 ALDH3A2 GALM    ALDOA   DLD     DLAT    ALDOB   G6PC2   LDHA    G6PC    PGM1    GPI
KEGG_CITRATE_CYCLE_TCA_CYCLE    http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_CITRATE_CYCLE_TCA_CYCLE    IDH3B   DLST    PCK2    CS      PDHB    PCK1    PDHA1   LOC642502       PDHA2   LOC283398       FH      SDHD   OGDH     SDHB    IDH3A   SDHC    IDH2    IDH1    ACO1    ACLY    MDH2    DLD     MDH1    DLAT    OGDHL   PC      SDHA    SUCLG1  SUCLA2  SUCLG2  IDH3G   ACO2
KEGG_PENTOSE_PHOSPHATE_PATHWAY  http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_PENTOSE_PHOSPHATE_PATHWAY  RPE     RPIA    PGM2    PGLS    PRPS2   FBP2    PFKM    PFKL    TALDO1  TKT     FBP1    TKTL2   PGD     RBKS   ALDOA    ALDOC   ALDOB   H6PD    LOC729020       PRPS1L1 PRPS1   DERA    G6PD    PGM1    TKTL1   PFKP    GPI
KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS   http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS   UGT1A10 UGT1A8  RPE     UGT1A7  UGT1A6  UGT2B28 UGT1A5  CRYL1   UGDH    UGT2A1 GUSB     UGT1A9  DCXR    LOC729020       DHDH    UGT2B11 UGP2    XYLB    UGT2B10 AKR1B1  UGT2B7  UGT2B4  UGT2A3  UGT1A4  UGT2B17 UGT1A1  UGT1A3  UGT2B15
KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM    http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM    MPI     PMM2    PMM1    FBP2    PFKM    GMDS    PFKFB4  PFKL    MTMR6   TPI1    PHPT1   PFKFB3 FUK      PFKFB2  MTMR1   PFKFB1  AKR1B10 FPGT    KHK     FBP1    MTMR2   HK2     HK3     HK1     ALDOA   ALDOC   ALDOB   MTMR7   TSTA3   AKR1B1  SORD    GMPPA   PFKP    GMPPB
KEGG_GALACTOSE_METABOLISM       http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_GALACTOSE_METABOLISM       GCK     GALK1   GLB1    GALE    B4GALT1 PGM2    LALBA   PFKM    PFKL    MGAM    HK2     HK1     HK3     GALT   G6PC2    GLA     GANC    LCT     GALK2   G6PC    UGP2    PGM1    AKR1B1  B4GALT2 GAA     PFKP
KEGG_ASCORBATE_AND_ALDARATE_METABOLISM  http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_ASCORBATE_AND_ALDARATE_METABOLISM  UGT1A10 UGT1A8  UGT1A7  UGT1A6  ALDH1B1 UGT2B28 ALDH2   UGT1A5  MIOX    UGDH    UGT2A1  ALDH9A1ALDH3A2  UGT1A9  ALDH7A1 UGT2B11 UGT2B10 UGT2B7  UGT2B4  UGT2A3  UGT1A4  UGT1A1  UGT2B17 UGT1A3  UGT2B15
KEGG_FATTY_ACID_METABOLISM      http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_FATTY_ACID_METABOLISM      CPT1A   CPT1C   ACADS   ALDH1B1 ACADSB  ACADL   ALDH2   ACADM   CYP4A11 ACAT2   ACADVL  ACAT1   ACAA2   HADH   HADHB    HADHA   CYP4A22 ADH7    ADH6    ACSL6   ADH1B   ADH1C   ADH4    ECHS1   ADH5    ALDH9A1 ALDH3A2 ACSL5   ADH1A   EHHADH  GCDH    ALDH7A1 ACOX3   ACSL1   ACAA1   CPT2    CPT1B   ACOX1   ECI2    ECI1    ACSL3   ACSL4
KEGG_STEROID_BIOSYNTHESIS       http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_STEROID_BIOSYNTHESIS       SOAT1   LSS     SQLE    EBP     CYP51A1 DHCR7   CYP27B1 DHCR24  HSD17B7 MSMO1   FDFT1   SC5DL   LIPA    CEL    TM7SF2   NSDHL   SOAT2
KEGG_PRIMARY_BILE_ACID_BIOSYNTHESIS     http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_PRIMARY_BILE_ACID_BIOSYNTHESIS     CYP46A1 SLC27A5 BAAT    CYP7B1  AKR1C4  HSD17B4 SCP2    AKR1D1  ACOX2   HSD3B7  CYP27A1 AMACR  CYP7A1   CYP8B1  CYP39A1 CH25H
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -ane 'print join("\n", @F[2..$#F]), "\n";' | wc
 444687  444687 2767870
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -ane 'print join("\n", @F[2..$#F]), "\n";' | sort | uniq | wc
  21095   21095  144468
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | perl -F\\t -lane 'print $F[6];' | sort | uniq | wc
  33331   33331  341748
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.IntergenicSplit.wRowPos.txt | sed 's/,/ /g' | awk '{ print $3 }' | sort | uniq -c
   4480 UTR3
     45 UTR3_downstream
     45 UTR3_upstream
   1024 UTR5
      3 UTR5;UTR3_downstream
      3 UTR5;UTR3_upstream
     20 UTR5_downstream
     20 UTR5_upstream
   3061 downstream
    101 downstream_downstream
    101 downstream_upstream
  23451 exonic
     15 exonic;splicing_downstream
     15 exonic;splicing_upstream
    140 exonic_downstream
    140 exonic_upstream
 180350 intergenic_downstream
 180350 intergenic_upstream
 132282 intronic
   1288 intronic_downstream
   1288 intronic_upstream
   2288 ncRNA_exonic
      1 ncRNA_exonic;splicing_downstream
      1 ncRNA_exonic;splicing_upstream
     12 ncRNA_exonic_downstream
     12 ncRNA_exonic_upstream
  21109 ncRNA_intronic
    490 ncRNA_intronic_downstream
    490 ncRNA_intronic_upstream
     15 ncRNA_splicing
    109 splicing
   2885 upstream
    129 upstream;downstream_downstream
    129 upstream;downstream_upstream
    103 upstream_downstream
    103 upstream_upstream
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.txt | wc                                                 
 556099 1663603 26951222
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.txt | sed 's/,/ /g' | awk '{ print $3 }' | sort | uniq -c
      1 Func.refGene
   4384 NA
   4480 UTR3
   1024 UTR5
   3061 downstream
  23451 exonic
 180350 intergenic_downstream
 180350 intergenic_upstream
 132282 intronic
   2288 ncRNA_exonic
  21109 ncRNA_intronic
     15 ncRNA_splicing
     42 nonsynonymous_SNV
    109 splicing
    266 synonymous_SNV
      2 unknown
   2885 upstream
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep nonsynonymous | sort -k 1,1 | awk '{ print $1 }' | sort | uniq -c | awk '{ if ($1 > 1) { print $0 } } ' | wc
    503    1006    7288
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep nonsynonymous | sort -k 1,1 | awk '{ print $1 }' | sort | uniq -c | awk '{ if ($1 > 1) { print $0 } } ' | R -q -e "Data1 <- read.table(file('stdin'), header=F); table(Data1[,1]);"
> Data1 <- read.table(file('stdin'), header=F); table(Data1[,1]);

  2   3   4   5   6   7   9  10  11  17 
344  97  33  10   8   5   2   2   1   1 
> 
> 
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep nonsynonymous | sort -k 1,1 | head -n 10
A1BG    19:58864479     exonic,NA,nonsynonymous_SNV     351099
A2ML1   12:9020489      exonic,NA,nonsynonymous_SNV     248900
A4GALT  22:43089849     exonic,NA,nonsynonymous_SNV     371204
A4GNT   3:137843476     exonic,NA,nonsynonymous_SNV     73938
AADAC   3:151545601     exonic,NA,nonsynonymous_SNV     75591
AADACL2 3:151463421     exonic,NA,nonsynonymous_SNV     75582
AARS2   6:44275011      exonic,NA,nonsynonymous_SNV     140933
ABCA1   9:107586753     exonic,NA,nonsynonymous_SNV     204952
ABCA1   9:107620867     exonic,NA,nonsynonymous_SNV     204966
ABCA10  17:67178316     exonic,NA,nonsynonymous_SNV     324862
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep -E 'exonic|intronic|UTR|upstream|downstream' | sort -k 1,1 | perl -lane 'my @info1 = split(/,/, $F[2]); if ($info1[0] =~ m/intergenic/) { my @dists1 = split(/=/, $info1[1]); if ($dists1[1] <= 20000) { print join("\t", @F); } } else { print join("\t", @F); }' | sed 's/,/ /g' | sed 's/=/ /'g | grep dist | awk '{ print $5 }' | sort | uniq -c | sort -g -k 2,2 | head -n 10
   2394 NONE
     13 1
      8 2
      7 3
      7 4
     10 5
      8 6
      6 7
      4 8
     10 9
(FORCE) [  mturchin@login002  ~/LabMisc/RamachandranLab/FORCE]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep -E 'exonic|intronic|UTR|upstream|downstream' | sort -k 1,1 | perl -lane 'my @info1 = split(/,/, $F[2]); if ($info1[0] =~ m/intergenic/) { my @dists1 = split(/=/, $info1[1]); if ($dists1[1] <= 20000) { print join("\t", @F); } } else { print join("\t", @F); }' | sed 's/,/ /g' | sed 's/=/ /'g | grep dist | awk '{ print $5 }' | sort | uniq -c | sort -rg -k 2,2 | head -n 10
      3 20000
      5 19999
      4 19998
      6 19997
      1 19996
      6 19995
      6 19994
      3 19993
      6 19992
      3 19991



~~~














