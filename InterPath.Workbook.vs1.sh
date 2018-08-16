#!/bin/sh

###20180605 -- InterPath

#20180605
#TOC

 - 20180605: Vs1

#20180813
#See https://stackoverflow.com/questions/4089430/how-can-i-determine-the-url-that-a-local-git-repository-was-originally-cloned-fr, https://stackoverflow.com/questions/42830557/git-remote-add-origin-vs-remote-set-url-origin/42830632, https://conda.io/docs/user-guide/tasks/manage-environments.html#cloning-an-environment
#NOTE -- Made switch from `InterPath` to `InterPath`, via ':.,$ s/InterPath/InterPath/g', and manually follow-up with directories and such; made changes to git repo as well


#20180605
#Conda environment setup/details/etc
#See `/users/mturchin/RamachandranLab.CCV_General.Workbook.vs1.sh` for initial setup 

conda create -n InterPath
source activate InterPath
#conda update conda
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
#20180611 #20180814 NOTE -- not sure if I actually installed fortran this way, since I tried to reinstall using these same commands and conda said fortran was not available? Briefly looking online suggests that installing 'gcc' also installs 'gfortran', which is what you want?
conda install armadillo fortran
conda install r-data.table r-bigmemory
conda install julia
##conda install r-coop -- failed, install via 'install.packages("coop")'
##conda install r-rbenchmark -- failed, install via 'install.packages("rbenchmark")'
##conda install r-cpgen -- failed, install via 'install.packages("cpgen")'
##conda install docker -- failed
#20180815 NOTE -- re-install for movement into 'InterPath' from 'FORCE'
###conda install R perl java-jdk
###conda install git plink bedtools vcftools vcftools bwa samtools picard gatk imagemagick gnuplot eigensoft tabix
###conda install r-base r-devtools 
###conda install r-knitr r-testthat r-cairo r-ashr 
###conda install r-rcolorbrewer r-essentials r-extrafont fonts-anaconda 
###conda install eigen boost gcc 
###conda install r-doParallel r-Rcpp r-RcppArmadillo r-RcppParallel r-CompQuadForm r-Matrix r-MASS r-truncnorm 
###conda install armadillo r-data.table r-bigmemory julia



#20180607
#Vs1

mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath
mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1
mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses
cd /users/mturchin/LabMisc/RamachandranLab/InterPath; git clone https://github.com/mturchin20/Ramachandran_InterPath
mkdir /users/mturchin/data/mturchin/Broad
mkdir /users/mturchin/data/mturchin/Broad/MSigDB
#NOTE -- need to log-in with an e-mail address and then DL directly from webpage, so cannot really do wget setup (I think anyways); so files were downloaded to MacBook Air and then scp'd over
#cd /users/mturchin/data/mturchin/Broad/MSigDB
#wget 
#From MackBook Air
#scp -p /Users/mturchin20/Downloads/*gmt mturchin@ssh.ccv.brown.edu:/users/mturchin/data/mturchin/Broad/MSigDB/.

#NOTE -- copy and pasted 'InterPath_Simulation.R' and 'InterPath.cpp' from Slack via Lorin

cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/InterPath.cpp /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.vs1.cpp
cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/InterPath_Simulation.R /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.vs1.R
cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/InterPath_Simulation.R /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Source.vs1.R

#See /users/mturchin/.bash_profile for how eventually got the below g++ call to work/function to get the stand-alone armadillo c++ code going
g++ Test1.cpp -o Test1.out  -O2 /gpfs/data/sramacha/mturchin/miniconda2/envs/InterPath/lib/libgfortran.so.4 /gpfs/data/sramacha/mturchin/miniconda2/envs/InterPath/lib/libarmadillo.so.8 -larmadillo -lopenblas -larpack -lgfortran -L/gpfs/data/sramacha/mturchin/miniconda2/envs/InterPath/lib

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
UKBioBankPops=`echo "African;African;Afr British;British.Ran4000;Brit4k Caribbean;Caribbean;Carib Chinese;Chinese;Chi Indian;Indian;Indn Irish;Irish;Irish Pakistani;Pakistani;Pkstn"`;

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
sourceCpp("/users/mturchin/LabMisc/RamachandranLab/InterPath/InterPath.cpp"); set.seed(123); for (i in c(1e4, 1e5, 2e5, 3e5, ncol(Data3))) {
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
#library("Rcpp"); library("RcppArmadillo"); sourceCpp("/users/mturchin/LabMisc/RamachandranLab/InterPath/InterPath.cpp"); 

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
wget https://imputationserver.sph.umich.edu/share/results/18af195da63a5ddac5bed2121ff91d23/qcreport.html https://imputationserver.sph.umich.edu/share/results/e4ebbc8805f1d1ca5cc907a6f427ab1d/chr_1.zip https://imputationserver.sph.umich.edu/share/results/d577ae0c30e12274d8d81018da1e742a/chr_10.zip https://imputationserver.sph.umich.edu/share/results/29335b7863323cc46487171a1733bcf1/chr_11.zip https://imputationserver.sph.umich.edu/share/results/9d6792a074b0069da818cb5ed42ff070/chr_12.zip https://imputationserver.sph.umich.edu/share/results/423ac9a702e190769fac8f1d8e46f714/chr_13.zip https://imputationserver.sph.umich.edu/share/results/a35d55128ce4f1504419cc25f322a987/chr_14.zip https://imputationserver.sph.umich.edu/share/results/7e57c480c1839fbaaf9db269580296b4/chr_15.zip https://imputationserver.sph.umich.edu/share/results/cc4988a76009b9152e63f1400d442565/chr_16.zip https://imputationserver.sph.umich.edu/share/results/907335ebc2aa249277d07d6f95749226/chr_17.zip https://imputationserver.sph.umich.edu/share/results/7fa773e847ab8af18f80bd84ba020007/chr_18.zip https://imputationserver.sph.umich.edu/share/results/1059f3a187b3c173d4f3c9405b3fb12b/chr_19.zip https://imputationserver.sph.umich.edu/share/results/d8b7cc1ad44aa127254737e243c48141/chr_2.zip https://imputationserver.sph.umich.edu/share/results/6fc5bb22eceb71687f0686e9f82dbe80/chr_20.zip https://imputationserver.sph.umich.edu/share/results/3699ad341f942e71631179ea4587948f/chr_21.zip https://imputationserver.sph.umich.edu/share/results/e8481e79f7a02ce8d063f88cbf454550/chr_22.zip https://imputationserver.sph.umich.edu/share/results/5f30403ceccfee49529e31ba1e923048/chr_3.zip https://imputationserver.sph.umich.edu/share/results/9a63534e084998d6cd2cac5ab5f0d3ac/chr_4.zip https://imputationserver.sph.umich.edu/share/results/d370075cc842a7ea528adabeacfe49ad/chr_5.zip https://imputationserver.sph.umich.edu/share/results/e02052f79b263a86f60a0c15f8a89893/chr_6.zip https://imputationserver.sph.umich.edu/share/results/b51808754155fa48133ba2e2087bca16/chr_7.zip https://imputationserver.sph.umich.edu/share/results/9d0c8422414a9934934df5509dabf442/chr_8.zip https://imputationserver.sph.umich.edu/share/results/6e566d76a1d2a36196c4423fce65fa35/chr_9.zip https://imputationserver.sph.umich.edu/share/results/80899845cb92ad45aba61420e2493210/statistics.txt https://imputationserver.sph.umich.edu/share/results/7cf83dcefd2d3b274cc56c78c9aa3556/chr_1.log https://imputationserver.sph.umich.edu/share/results/eae6ec69f4afed5ee633079de3e321c5/chr_10.log https://imputationserver.sph.umich.edu/share/results/8d7620b974c07c09863f5add9332a6ed/chr_11.log https://imputationserver.sph.umich.edu/share/results/6c99680acfce437df804628a0c0880ba/chr_12.log https://imputationserver.sph.umich.edu/share/results/e9b5253bc81db919ca8fa8f785091d05/chr_13.log https://imputationserver.sph.umich.edu/share/results/156d6f626acb1645697fdda999fe2ac/chr_14.log https://imputationserver.sph.umich.edu/share/results/28102eea804e3ee12ffdf094b14bbc31/chr_15.log https://imputationserver.sph.umich.edu/share/results/c3a9c675502248bfc355bfb567b70739/chr_16.log https://imputationserver.sph.umich.edu/share/results/255c44bdc1cf21aac6e2e48acfdd054d/chr_17.log https://imputationserver.sph.umich.edu/share/results/9f6c6ec7d3d106e909d45249b4fa38cd/chr_18.log https://imputationserver.sph.umich.edu/share/results/facaa2edfa7e1a99b54d41a8f4ca770e/chr_19.log https://imputationserver.sph.umich.edu/share/results/4f917056fa7f6e36d7a277720a665f9/chr_2.log https://imputationserver.sph.umich.edu/share/results/ad225cccb2986df6e3d0f8f5b7b359cf/chr_20.log https://imputationserver.sph.umich.edu/share/results/13ac54739be370d2bc3479b6150de7e7/chr_21.log https://imputationserver.sph.umich.edu/share/results/f272bc117cf56ef499948859d413db62/chr_22.log https://imputationserver.sph.umich.edu/share/results/c9f1baeb315f4d8106f88a130392f404/chr_3.log https://imputationserver.sph.umich.edu/share/results/b2dec836903d7f66d952918d4057a177/chr_4.log https://imputationserver.sph.umich.edu/share/results/e9bec888da06e2e96e8c19d9a535fcd5/chr_5.log https://imputationserver.sph.umich.edu/share/results/c69ff6b589196935ea49d20ef2ade91/chr_6.log https://imputationserver.sph.umich.edu/share/results/c5e09b62d45003ba297c0c0edc13f571/chr_7.log https://imputationserver.sph.umich.edu/share/results/a3734afc604653e8af4c1f755d66bcea/chr_8.log https://imputationserver.sph.umich.edu/share/results/51c5096eab99d970a2a828d8de3aca21/chr_9.log
for i in {1..22}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran4000/Imputation/mturchin20/chr_${i}.zip -p'c}iuH7gVEZ2Owh'
done
cd /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20
wget https://imputationserver.sph.umich.edu/share/results/a485b7b0b7052952a7c9621ff057a648/qcreport.html https://imputationserver.sph.umich.edu/share/results/6decd1a5d506724daf48ae4fd3d729c9/chr_1.zip https://imputationserver.sph.umich.edu/share/results/2fccb0963963791394a685fcd3f1dda1/chr_10.zip https://imputationserver.sph.umich.edu/share/results/2673843ed044f65268fa8f8ecde8f09d/chr_11.zip https://imputationserver.sph.umich.edu/share/results/6210d35edf635643310810d1c20c7d98/chr_12.zip https://imputationserver.sph.umich.edu/share/results/2dbdc1bcc22fc0caf05d4dcc047fc01e/chr_13.zip https://imputationserver.sph.umich.edu/share/results/3e7a895a2270baf2d7f94076bff4209a/chr_14.zip https://imputationserver.sph.umich.edu/share/results/62e524101365b92ad732b77d58dc9735/chr_15.zip https://imputationserver.sph.umich.edu/share/results/2967a835b8145401e1bcfe65f7a48348/chr_16.zip https://imputationserver.sph.umich.edu/share/results/4a77dc8e1535ce1a4c174919a6abea5f/chr_17.zip https://imputationserver.sph.umich.edu/share/results/bbe3e5183414820961604c230bc1f4a3/chr_18.zip https://imputationserver.sph.umich.edu/share/results/7ed32c5b53f5689d876ed8b41d1982cf/chr_19.zip https://imputationserver.sph.umich.edu/share/results/74f139cf2b941d8625ecbc7975b77c56/chr_2.zip https://imputationserver.sph.umich.edu/share/results/c027a3491de5d0bdb70a96d16a20c7eb/chr_20.zip https://imputationserver.sph.umich.edu/share/results/245444f21f2b3db1f1d3028286f6862d/chr_21.zip https://imputationserver.sph.umich.edu/share/results/edb6cde98d8343333517fa0bb7cd353/chr_22.zip https://imputationserver.sph.umich.edu/share/results/e6fd9a48817eaf83ba689d9397627112/chr_3.zip https://imputationserver.sph.umich.edu/share/results/5115bdfa022b4375c0aec20fb0253ecf/chr_4.zip https://imputationserver.sph.umich.edu/share/results/6683666699c0ccda33b0c63a55217e51/chr_5.zip https://imputationserver.sph.umich.edu/share/results/e97d27b89862e129c1ff514cea673711/chr_6.zip https://imputationserver.sph.umich.edu/share/results/12a39009117e3914a8af4cdf12926e7c/chr_7.zip https://imputationserver.sph.umich.edu/share/results/17f21694afb7de887b9bfa2f54a7d25a/chr_8.zip https://imputationserver.sph.umich.edu/share/results/437570d81c0cde6214cef1df89792b5a/chr_9.zip https://imputationserver.sph.umich.edu/share/results/e0953d319bc879ee1252a6637002f40d/statistics.txt https://imputationserver.sph.umich.edu/share/results/28258d8398858b1544f621b5428ed0ae/chr_1.log https://imputationserver.sph.umich.edu/share/results/f17c716bddb3c37c30e053e1e4068178/chr_10.log https://imputationserver.sph.umich.edu/share/results/d97bf7f986869a7fc28b604a3c3c0534/chr_11.log https://imputationserver.sph.umich.edu/share/results/ef089ae5540bd57699dd50d8146cbf13/chr_12.log https://imputationserver.sph.umich.edu/share/results/d31902cc5dfca42c84b23a4f1154efe4/chr_13.log https://imputationserver.sph.umich.edu/share/results/1a6e127938d46ed2e21634171805bb2a/chr_14.log https://imputationserver.sph.umich.edu/share/results/e72954d40cc648a32497ab76a3879c47/chr_15.log https://imputationserver.sph.umich.edu/share/results/91aea73db902c429c359bd394819474e/chr_16.log https://imputationserver.sph.umich.edu/share/results/9e01278d3a5f44dc2235260297b9f605/chr_17.log https://imputationserver.sph.umich.edu/share/results/571c491d461d8fdefaaf36cc746555eb/chr_18.log https://imputationserver.sph.umich.edu/share/results/b166c9e6a90fd1aea3a367b8fdb3d69/chr_19.log https://imputationserver.sph.umich.edu/share/results/2e0bd6e55afd4df7639f71e884034a5d/chr_2.log https://imputationserver.sph.umich.edu/share/results/df20d49017fc5c5f3f332bcbec121134/chr_20.log https://imputationserver.sph.umich.edu/share/results/e052113bd01c150dcee38ab6efd4d6ee/chr_21.log https://imputationserver.sph.umich.edu/share/results/4a91305b653af16ca3296a3c12ebfa34/chr_22.log https://imputationserver.sph.umich.edu/share/results/6e53f8918933312929fe3e6646b001a9/chr_3.log https://imputationserver.sph.umich.edu/share/results/f322d2bad3730aef4a0873fb1139f221/chr_4.log https://imputationserver.sph.umich.edu/share/results/d12437713af84c756c882dc19be46f67/chr_5.log https://imputationserver.sph.umich.edu/share/results/ce307a4ea18b7ce2d50c0d5f907f3b3d/chr_6.log https://imputationserver.sph.umich.edu/share/results/50425e45e17e250aff62d74de8fdd98/chr_7.log https://imputationserver.sph.umich.edu/share/results/522f6ef95a825e5a17e088a39a0a2a2e/chr_8.log https://imputationserver.sph.umich.edu/share/results/f408b06ec0e4aca0bfafbc4dac03742/chr_9.log
for i in {1..22}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/chr_${i}.zip -p'ReK7lYFndsE7'
done
cd /users/mturchin/data/ukbiobank_jun17/subsets/Chinese/Chinese/Imputation/mturchin20
wget https://imputationserver.sph.umich.edu/share/results/72fd72b0933b942c014cb84b4d3791da/qcreport.html https://imputationserver.sph.umich.edu/share/results/d34ff4d581219ab7c3a3ac8155886dbc/chr_1.zip https://imputationserver.sph.umich.edu/share/results/e960fe15cd643149e03c3e7948f9fed/chr_10.zip https://imputationserver.sph.umich.edu/share/results/1abb741440448737d11fcb26cfc4d4d6/chr_11.zip https://imputationserver.sph.umich.edu/share/results/be21465bf4b17f77da870679ebe709a0/chr_12.zip https://imputationserver.sph.umich.edu/share/results/2c5a33d41883fdca2e34f9b24c840107/chr_13.zip https://imputationserver.sph.umich.edu/share/results/4789cb094d0425b4043c4a5a0682e61e/chr_14.zip https://imputationserver.sph.umich.edu/share/results/cecc0c31b9d19f2db0a128e9b8d75287/chr_15.zip https://imputationserver.sph.umich.edu/share/results/419a04d7a603a61dc4382383f54325/chr_16.zip https://imputationserver.sph.umich.edu/share/results/bce04ca77fd2e79c75fe1f5840a1c77d/chr_17.zip https://imputationserver.sph.umich.edu/share/results/b2542cff81a0dfeb9523fdd2285799b6/chr_18.zip https://imputationserver.sph.umich.edu/share/results/b74dd2154a06572baa578c98815ff6e8/chr_19.zip https://imputationserver.sph.umich.edu/share/results/fa0ab64b9c178bad297f6d53fd50cac3/chr_2.zip https://imputationserver.sph.umich.edu/share/results/79f52165ef523e78bc5b22f349888077/chr_20.zip https://imputationserver.sph.umich.edu/share/results/60059ac3ac86a39d0b01f3fae5cd23e1/chr_21.zip https://imputationserver.sph.umich.edu/share/results/dbf4318da553da8354d2867699245d5c/chr_22.zip https://imputationserver.sph.umich.edu/share/results/f24465eb9040690614a13d375ea01e26/chr_3.zip https://imputationserver.sph.umich.edu/share/results/eda0b5faebd85404a480f13e58caf5b2/chr_4.zip https://imputationserver.sph.umich.edu/share/results/3c6cbcb7987fb0d592ba0ea8633830c/chr_5.zip https://imputationserver.sph.umich.edu/share/results/421541226f2f9cdff35a7486c02da4f2/chr_6.zip https://imputationserver.sph.umich.edu/share/results/499bd22ae71773085fa4008dc8e0a08e/chr_7.zip https://imputationserver.sph.umich.edu/share/results/7306a61b463828e5e2c4ed3071695320/chr_8.zip https://imputationserver.sph.umich.edu/share/results/9c92b491eb268b67f8837585b4374acc/chr_9.zip https://imputationserver.sph.umich.edu/share/results/b5908a724cb45a09b8d202e622060d08/statistics.txt https://imputationserver.sph.umich.edu/share/results/10756c72fc608357eda5d5a945d06ae3/chr_1.log https://imputationserver.sph.umich.edu/share/results/2cc9971ab5eb17be6aeddd5ead1eba1b/chr_10.log https://imputationserver.sph.umich.edu/share/results/5309d67a5a370422effd5da9b4445410/chr_11.log https://imputationserver.sph.umich.edu/share/results/a8072fa4d2b24e2375215d91c2fbb669/chr_12.log https://imputationserver.sph.umich.edu/share/results/12bb7d67ea34e88296b785e0bd37be9c/chr_13.log https://imputationserver.sph.umich.edu/share/results/4f0017e28ff6a07ae3413ec9e61ec6e9/chr_14.log https://imputationserver.sph.umich.edu/share/results/bbe8947650fe640d79eb171dc22700df/chr_15.log https://imputationserver.sph.umich.edu/share/results/caed589f216ff402123382a318939a4a/chr_16.log https://imputationserver.sph.umich.edu/share/results/88d6fba610f0d7e4d102e1e3fb8cce35/chr_17.log https://imputationserver.sph.umich.edu/share/results/5e75aa91f4b40c5b780cec2dc8126741/chr_18.log https://imputationserver.sph.umich.edu/share/results/99b44d7833cd2dd159fd9706887f6313/chr_19.log https://imputationserver.sph.umich.edu/share/results/ea14c6476e65a2584d549828561078ad/chr_2.log https://imputationserver.sph.umich.edu/share/results/828ce1048122096e6ce641ffba40235b/chr_20.log https://imputationserver.sph.umich.edu/share/results/b7c4f2cc743b0f3895c853cae3e03f0e/chr_21.log https://imputationserver.sph.umich.edu/share/results/7884ff42d714889317d7ffd297ce9da/chr_22.log https://imputationserver.sph.umich.edu/share/results/433bf557524905ee023083f5062abd75/chr_3.log https://imputationserver.sph.umich.edu/share/results/1cff446d7a2b75699eaf7f1a3c9c20c4/chr_4.log https://imputationserver.sph.umich.edu/share/results/ea7b48716f66e1afb4bc6411963941be/chr_5.log https://imputationserver.sph.umich.edu/share/results/3b9fdfd90410bee77b0ac5e7a790cbed/chr_6.log https://imputationserver.sph.umich.edu/share/results/73ee612f42f9a2e95fc8f24d0cbc07bd/chr_7.log https://imputationserver.sph.umich.edu/share/results/260b5bd8c020cbb2dea160879aee604c/chr_8.log https://imputationserver.sph.umich.edu/share/results/dde86bea56ec6f01fcb44981b05846bc/chr_9.log
for i in {1..22}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/Chinese/Chinese/Imputation/mturchin20/chr_${i}.zip -p'OhLHXUz7mhVfO3'
done
cd /users/mturchin/data/ukbiobank_jun17/subsets/Irish/Irish/Imputation/mturchin20
wget https://imputationserver.sph.umich.edu/share/results/d5df27d8eb41947a422d2715c5fb1a11/qcreport.html https://imputationserver.sph.umich.edu/share/results/b20b486baa0ed4b45bf8c522723509c5/chr_1.zip https://imputationserver.sph.umich.edu/share/results/572ef9960c8cb683ed490a627178dd8c/chr_10.zip https://imputationserver.sph.umich.edu/share/results/70dbf995ab04df8a1faa3ccdb9b92bb2/chr_11.zip https://imputationserver.sph.umich.edu/share/results/d66a6193ab41db26fda2069e44a46d34/chr_12.zip https://imputationserver.sph.umich.edu/share/results/a7ece25b02b2441df8d013f8d048a999/chr_13.zip https://imputationserver.sph.umich.edu/share/results/a0deca9a563dc4348a093201b1a19087/chr_14.zip https://imputationserver.sph.umich.edu/share/results/2cc537737eb1c5a331e2af151eae8963/chr_15.zip https://imputationserver.sph.umich.edu/share/results/b111157268f243d898ded4ba07e96e9f/chr_16.zip https://imputationserver.sph.umich.edu/share/results/51c65fb5c7fe077eff088327b3c150a0/chr_17.zip https://imputationserver.sph.umich.edu/share/results/2e1729bc745441ba93fb7dacc37c7bc5/chr_18.zip https://imputationserver.sph.umich.edu/share/results/144ad9824d68cb7786d8c0070ab08181/chr_19.zip https://imputationserver.sph.umich.edu/share/results/8e0dfb809dc323a1e515246ef0caf4e5/chr_2.zip https://imputationserver.sph.umich.edu/share/results/7445bf33f8fdbb21fdbffa4a57ebc9f7/chr_20.zip https://imputationserver.sph.umich.edu/share/results/3b454752a673208638b45be617382343/chr_21.zip https://imputationserver.sph.umich.edu/share/results/9633dacfa6d498a481be0a45b1854af6/chr_22.zip https://imputationserver.sph.umich.edu/share/results/4ed1f6546e6626242785da4956ce365a/chr_3.zip https://imputationserver.sph.umich.edu/share/results/26b28f462a51c4b4745369b9573958b4/chr_4.zip https://imputationserver.sph.umich.edu/share/results/bb7cd175ec841a45d6a38b8a5dd824ba/chr_5.zip https://imputationserver.sph.umich.edu/share/results/5444bed3443307f0ea74ffe1932c92c0/chr_6.zip https://imputationserver.sph.umich.edu/share/results/dbcf69353f4a1b09193c8d50a643af66/chr_7.zip https://imputationserver.sph.umich.edu/share/results/bcbdd1be75f9a4d78611cf60e4736354/chr_8.zip https://imputationserver.sph.umich.edu/share/results/ae5abcae9bdeed9db14f36aea68d8630/chr_9.zip https://imputationserver.sph.umich.edu/share/results/a790cae79120bc822b185e822a78cd45/statistics.txt https://imputationserver.sph.umich.edu/share/results/96975ce905d5803d08dd3c51f04c535e/chr_1.log https://imputationserver.sph.umich.edu/share/results/daa680f24483db1ccf8af04d8049771d/chr_10.log https://imputationserver.sph.umich.edu/share/results/154fb851dcda106ea27d690a20e05b48/chr_11.log https://imputationserver.sph.umich.edu/share/results/f36edf28ee94dde525f0b434bf2c1dd2/chr_12.log https://imputationserver.sph.umich.edu/share/results/8e2cb869aaa3f88adc36cb1f4fb6b511/chr_13.log https://imputationserver.sph.umich.edu/share/results/6de60bd26f6dbc39678e55a1a442f68b/chr_14.log https://imputationserver.sph.umich.edu/share/results/bd6cad65431891ec06f1ec0cfaebd2aa/chr_15.log https://imputationserver.sph.umich.edu/share/results/407c201dc9541cc4cd01d434168b92d6/chr_16.log https://imputationserver.sph.umich.edu/share/results/be85efb27216ea69f81e76389eec475d/chr_17.log https://imputationserver.sph.umich.edu/share/results/e8090ce9810a32f5f07e6565e3e90460/chr_18.log https://imputationserver.sph.umich.edu/share/results/61b0354d2b54a025cf1bc67d9617414b/chr_19.log https://imputationserver.sph.umich.edu/share/results/172ed05cfe2a6450331a693fd2236fb1/chr_2.log https://imputationserver.sph.umich.edu/share/results/2289dd89acfc4a7bd98c7c5622f8f1f6/chr_20.log https://imputationserver.sph.umich.edu/share/results/4bf20a70083807098e27898d098a8045/chr_21.log https://imputationserver.sph.umich.edu/share/results/c4791403ceb5b4f3f75478a412259c9e/chr_22.log https://imputationserver.sph.umich.edu/share/results/3aed3e52a1b08f86bc63da3a7c5ba615/chr_3.log https://imputationserver.sph.umich.edu/share/results/ee04b7a0f89a5a9e6070cbec2ee1ae91/chr_4.log https://imputationserver.sph.umich.edu/share/results/1e6787f53336b180acab0ac2c9b9bbf2/chr_5.log https://imputationserver.sph.umich.edu/share/results/8b5e9b096666533e1aae0e47598b0f50/chr_6.log https://imputationserver.sph.umich.edu/share/results/cf5f176fe09b549a7a611acfed3c5d38/chr_7.log https://imputationserver.sph.umich.edu/share/results/7b63091a404dfa1520890bc1da2bfdfa/chr_8.log https://imputationserver.sph.umich.edu/share/results/891b660fb644c92261113f84b253ee89/chr_9.log
for i in {1..22}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/Irish/Irish/Imputation/mturchin20/chr_${i}.zip -p'Ie0PdKwGPs2ui4'
done
cd /users/mturchin/data/ukbiobank_jun17/subsets/Indian/Indian/Imputation/mturchin20
wget https://imputationserver.sph.umich.edu/share/results/6a1bbb76e6cb43b790d87660db1543f8/qcreport.html https://imputationserver.sph.umich.edu/share/results/7971adf76bcc4b585602bddc72092cce/chr_1.zip https://imputationserver.sph.umich.edu/share/results/c2d8f85e9594e81ce338a3705f80fd8d/chr_10.zip https://imputationserver.sph.umich.edu/share/results/4555755805818507df6b55cf29bb09cd/chr_11.zip https://imputationserver.sph.umich.edu/share/results/2866929fa32ab51a67b011a59864fd0e/chr_12.zip https://imputationserver.sph.umich.edu/share/results/bfadf001b0801c82c421f6465caa633b/chr_13.zip https://imputationserver.sph.umich.edu/share/results/fff7da908616a000bc93c09646448b81/chr_14.zip https://imputationserver.sph.umich.edu/share/results/2c4f73b7cc4a49b169a940299bcfcf35/chr_15.zip https://imputationserver.sph.umich.edu/share/results/4b4ee48fd51035a2efb1d52d22b4ae99/chr_16.zip https://imputationserver.sph.umich.edu/share/results/b6fe0affba7c5d61701417ffb60d242/chr_17.zip https://imputationserver.sph.umich.edu/share/results/8ed8970caaff8ae2d1eb852707f1f0c4/chr_18.zip https://imputationserver.sph.umich.edu/share/results/75a76e1a90c68bacf89d9876b0092889/chr_19.zip https://imputationserver.sph.umich.edu/share/results/a620623334a33bad132eca7c8bb5ca72/chr_2.zip https://imputationserver.sph.umich.edu/share/results/6c868b4c9a58c0616e55b923f86365ef/chr_20.zip https://imputationserver.sph.umich.edu/share/results/1b38aeb03565d56a20b5f91c067d652a/chr_21.zip https://imputationserver.sph.umich.edu/share/results/6f39959029e7796b29402afba97f9055/chr_22.zip https://imputationserver.sph.umich.edu/share/results/4cee5b7e2e7af1f15c54f25b0244f861/chr_3.zip https://imputationserver.sph.umich.edu/share/results/f992ec599a39e0fbbc2725ee1d979077/chr_4.zip https://imputationserver.sph.umich.edu/share/results/e1092ef460bac48d373bec28f7aeedb7/chr_5.zip https://imputationserver.sph.umich.edu/share/results/814f945586d0237861e49e2e2d47ecc9/chr_6.zip https://imputationserver.sph.umich.edu/share/results/6493c285051dd38aa3220940f7f5f4df/chr_7.zip https://imputationserver.sph.umich.edu/share/results/a2c3a863199fd38be8d6bb8a831834c8/chr_8.zip https://imputationserver.sph.umich.edu/share/results/196e3d377482faf573c1b24d002429aa/chr_9.zip https://imputationserver.sph.umich.edu/share/results/8e17c38f822f66e2ce025eda58ea37fb/statistics.txt https://imputationserver.sph.umich.edu/share/results/2ca52e879ed00c3c26668cf4b46f07fa/chr_1.log https://imputationserver.sph.umich.edu/share/results/de6b69bc13320bc92c85f1fef656ad92/chr_10.log https://imputationserver.sph.umich.edu/share/results/93a2f725e36d5ccb55f7a87402b1020a/chr_11.log https://imputationserver.sph.umich.edu/share/results/dc788719a1c30b61fa4a1a64f1e700f2/chr_12.log https://imputationserver.sph.umich.edu/share/results/5f2c4bd16c50a9658b9f74bddd29293d/chr_13.log https://imputationserver.sph.umich.edu/share/results/3755fca0b0758f95ed016dcbbdc66242/chr_14.log https://imputationserver.sph.umich.edu/share/results/ff2922a4afeeb67740045881e11d874e/chr_15.log https://imputationserver.sph.umich.edu/share/results/9a3c90c07e6e31076bd572542cdb87bf/chr_16.log https://imputationserver.sph.umich.edu/share/results/95d9128b9ea90828f76cb72995fb2c50/chr_17.log https://imputationserver.sph.umich.edu/share/results/3b72f0bd381c3aadfc7a7b4a8b02ee7c/chr_18.log https://imputationserver.sph.umich.edu/share/results/121b01e348f912cbe5ad5e63f6f97539/chr_19.log https://imputationserver.sph.umich.edu/share/results/c65909e6103fc8bda7f74c4842b48c12/chr_2.log https://imputationserver.sph.umich.edu/share/results/7912fcb372b4bb6231627086fecf383c/chr_20.log https://imputationserver.sph.umich.edu/share/results/efc09b2c74541c308cd173ab08a79cfa/chr_21.log https://imputationserver.sph.umich.edu/share/results/f83e893920c952b264470f7631d50f0d/chr_22.log https://imputationserver.sph.umich.edu/share/results/6fd203b72b578b53743b97629eb2004d/chr_3.log https://imputationserver.sph.umich.edu/share/results/fb4727b559ebf5d48ccc39fa5f87e5e9/chr_4.log https://imputationserver.sph.umich.edu/share/results/33cf2cc23398df49afe374863baf4089/chr_5.log https://imputationserver.sph.umich.edu/share/results/92471d5695704d7b8c6cdee9fdba1b3f/chr_6.log https://imputationserver.sph.umich.edu/share/results/ba85d1b171aa0a0b3f83af5846c71b41/chr_7.log https://imputationserver.sph.umich.edu/share/results/faadd4484d3cc17e669216274652678a/chr_8.log https://imputationserver.sph.umich.edu/share/results/a8897f4740b8ee14369b35e374a0bd15/chr_9.log
for i in {1..22}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/Indian/Indian/Imputation/mturchin20/chr_${i}.zip -p'ZXFvfg9uUw4BoZ'
done
cd /users/mturchin/data/ukbiobank_jun17/subsets/Pakistani/Pakistani/Imputation/mturchin20 
wget https://imputationserver.sph.umich.edu/share/results/89b43b5e2690314033b3f0c6cd59a059/qcreport.html https://imputationserver.sph.umich.edu/share/results/be797334e15c635bd88242589eb26fe1/chr_1.zip https://imputationserver.sph.umich.edu/share/results/5d91caa3bcc16a9b38ad1dcab128bb9b/chr_10.zip https://imputationserver.sph.umich.edu/share/results/5dd297085806345d9af88bc2d194ba1c/chr_11.zip https://imputationserver.sph.umich.edu/share/results/84a5e2b9dd3d88c80427b41202f6df51/chr_12.zip https://imputationserver.sph.umich.edu/share/results/a199595f92db69551058f70963429001/chr_13.zip https://imputationserver.sph.umich.edu/share/results/2b2c82505e7b21e32729e3feb4b69347/chr_14.zip https://imputationserver.sph.umich.edu/share/results/98ae84299e667991b46451c1d8a54b35/chr_15.zip https://imputationserver.sph.umich.edu/share/results/ce602c06f2bbe2b0a760213cff91cb5e/chr_16.zip https://imputationserver.sph.umich.edu/share/results/b46491e54a500e3e2a0b77337cd1b5a6/chr_17.zip https://imputationserver.sph.umich.edu/share/results/8b0a1702e128789566219d4099cb57a8/chr_18.zip https://imputationserver.sph.umich.edu/share/results/96f36eef54aee0db8f60bb1aad86c5d8/chr_19.zip https://imputationserver.sph.umich.edu/share/results/1fce2ac575450eba4cd42549134177a0/chr_2.zip https://imputationserver.sph.umich.edu/share/results/b69ce6294470213b2a5f6db357e8ed35/chr_20.zip https://imputationserver.sph.umich.edu/share/results/cf0e18dfc413e0cde8cfc2ac020f74dd/chr_21.zip https://imputationserver.sph.umich.edu/share/results/d4d4824c70be4eb77bcc16e94bec12c1/chr_22.zip https://imputationserver.sph.umich.edu/share/results/153e4e9e5f5d3f0cb45531e56a15b81b/chr_3.zip https://imputationserver.sph.umich.edu/share/results/961f8b201a27d52383142317fcec197b/chr_4.zip https://imputationserver.sph.umich.edu/share/results/2e1f964fc3f65645ebd36adaf1b03362/chr_5.zip https://imputationserver.sph.umich.edu/share/results/585d1d39e30796a196981b0fb63af7ca/chr_6.zip https://imputationserver.sph.umich.edu/share/results/898a02419d5c54a2ef5812497b079e2b/chr_7.zip https://imputationserver.sph.umich.edu/share/results/efdfab33d8766adc6d90f803c3372d43/chr_8.zip https://imputationserver.sph.umich.edu/share/results/912f52fd4120c124d888b98f500d542/chr_9.zip https://imputationserver.sph.umich.edu/share/results/ef3c072d9b07558c0f97f58f83a5bde5/statistics.txt https://imputationserver.sph.umich.edu/share/results/5317b95204bdb7d6f7809a4a21befdbd/chr_1.log https://imputationserver.sph.umich.edu/share/results/eaf2dd544c7bd69f48e5bb895a3c121f/chr_10.log https://imputationserver.sph.umich.edu/share/results/36b12de97bdf419ec1a994ed9dacfcbe/chr_11.log https://imputationserver.sph.umich.edu/share/results/a10e9382946e1f8d987aa4cb548fc3af/chr_12.log https://imputationserver.sph.umich.edu/share/results/95d5b6c5521383515b5c0d4136b58828/chr_13.log https://imputationserver.sph.umich.edu/share/results/cb6138991310af861722b96530016da0/chr_14.log https://imputationserver.sph.umich.edu/share/results/35dbc3f17547c32c2149648f2820dfd4/chr_15.log https://imputationserver.sph.umich.edu/share/results/63761073c54672cda5b853ff39cb3f0f/chr_16.log https://imputationserver.sph.umich.edu/share/results/5f14b005acc9b7e009eccfe306ce6505/chr_17.log https://imputationserver.sph.umich.edu/share/results/66eaaa1c2c39a03489e469b25ddfbf20/chr_18.log https://imputationserver.sph.umich.edu/share/results/294dbd112c5440a80b8a35fd33ad4542/chr_19.log https://imputationserver.sph.umich.edu/share/results/2701f8b8da63bdf9e0b5525469d1a20b/chr_2.log https://imputationserver.sph.umich.edu/share/results/e7b0932d88b577f1303aaff9d97b6df5/chr_20.log https://imputationserver.sph.umich.edu/share/results/d8caf44de02401e5fbba3beaa4c5eeb1/chr_21.log https://imputationserver.sph.umich.edu/share/results/17300cdbf1b86e2ae6e9c36b50daead7/chr_22.log https://imputationserver.sph.umich.edu/share/results/248249066b6ffa7f943ea069df12691f/chr_3.log https://imputationserver.sph.umich.edu/share/results/63a9373cb64e596d59564ca470ff9679/chr_4.log https://imputationserver.sph.umich.edu/share/results/d2dde3fec0976ffa46eb34a5c1fa8258/chr_5.log https://imputationserver.sph.umich.edu/share/results/3f906b58d084a776144ab59e849bb6f5/chr_6.log https://imputationserver.sph.umich.edu/share/results/445b43f10f8237ae1d53947034fd9f51/chr_7.log https://imputationserver.sph.umich.edu/share/results/8404e9d7cc06023c508025b99bb2ef71/chr_8.log https://imputationserver.sph.umich.edu/share/results/54ae77fb7376418d39f310d969f07012/chr_9.log
for i in {1..22}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/Pakistani/Pakistani/Imputation/mturchin20/chr_${i}.zip -p'lklPN{3YAH2fTp'
done





#From http://csg.sph.umich.edu/abecasis/mach/tour/imputation.html
#r2 cutoff .3
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v African`; do
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

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v African`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
        ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
		
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3

	join <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr*_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.info.gz | awk '{ if ($7 > .3) { print $1 } }' | sed 's/:/_/g' | grep -v SNP | sort | uniq -u) <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr*_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.info.gz | awk '{ print $1 }' | sed 's/:/_/g' | grep -v SNP | sort | uniq -u) | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.info.r2gt3.noDups.ChrBPs.gz
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr*_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim | awk '{ print $1 "_" $4 }' | sort | uniq > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim.ChrBPs
	join <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.info.r2gt3.noDups.ChrBPs.gz | sort) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim.ChrBPs | sort) | sed 's/_/:/g' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim.ImptHRC.info.r2gt3.noDups.ChrBPs

done

#From https://www.biostars.org/p/46060/ & https://sourceforge.net/p/vcftools/mailman/message/29115811/
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);')`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
        ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
		
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3

	for i in {1..22}; do
		echo $i
		sbatch -t 72:00:00 --mem 2g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.error --comment "$ancestry1 $ancestry2 $i" <(echo -e '#!/bin/sh';
		echo -e "\nrm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.plinkTemp.tped /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.plinkTemp.tfam";)	
#		rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp* 
#		echo -e "\nvcftools --gzvcf /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.vcf.gz --plink-tped --snps /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim.ImptHRC.info.r2gt3.noDups.ChrBPs --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp"; \
#		echo -e "\nplink --tfile /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp --geno 0 --make-bed --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.plinkTemp";)
	done
	sleep 1
done 

#		echo -e "\nplink --tfile /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp --exclude /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.missnp --make-bed --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp"; \
#		echo -e "\nrm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.bed /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.bim /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.fam /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.tped /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.tfam"; \	

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v African`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`

        echo $ancestry1 $ancestry2 /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.MergeList.Vs2.txt

        for chr in {2..22}; do
                echo "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${chr}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.plinkTemp.bed /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${chr}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.plinkTemp.bim /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${chr}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.plinkTemp.fam" 
        done > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.MergeList.Vs2.txt

done
#	rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.MergeList.Vs2.txt

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v African`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`

        echo $ancestry1 $ancestry2 /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.MergeList.Vs2.txt
        sbatch -t 24:00:00 --mem 20g --qos=Normal -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.MergeList.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.MergeList.error --comment "$ancestry1 $ancestry2" <(echo -e '#!/bin/sh'; echo -e "\nplink --bfile /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr1_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.plinkTemp --merge-list /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.MergeList.Vs2.txt --recodeAD --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno"; \
	echo -e "cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw | perl -lane 'if (\$. == 1) { @vals1; for (my \$i = 6; \$i <= \$#F; \$i++) { if (\$F[\$i] =~ m/HET/) { \$PH = 1 } else { push(@vals1, \$i); } } } print join(\"\t\", @F[@vals1]);' | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz"; \
	echo -e "rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.edit.gz /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw";)

done
	
#rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.MergeList.*
#	rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.edit.gz /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw

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

        echo $ancestry1 $ancestry2 /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.pruned.MergeList.Vs2.txt

	rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.pruned.bed /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.pruned.bim /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.pruned.fam /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.pruned.log

done

#	gzip -f /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw
#	gzip -f /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.pruned.raw

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);')`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`

        echo $ancestry1 $ancestry2 /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.txt 

	join <(cat /users/mturchin/data/ukbiobank_jun17/mturchin/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.txt | awk '{ print $1 "_" $2 "\t" $0 }' | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.txt | awk '{ print $1 "_" $2 "\t" $0 }' | sort -k 1,1) | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23 }' | grep -v FID | R -q -e "Data1 <- read.table(file('stdin'), header=F); for (i in 3:6) { Data1 <- cbind(Data1, residuals(lm(Data1[,i] ~ Data1[,10] + Data1[,11] + Data1[,12] + Data1[,13] + Data1[,14] + Data1[,15] + Data1[,16] + Data1[,17] + Data1[,18] + Data1[,19], na.action=na.exclude))); }; write.table(Data1, file=\"\", quote=FALSE, row.name=FALSE, col.name=FALSE);" | grep -v \> | cat <(echo "FID IID Height BMI Waist Hip FID IID SEX ANCESTRY AGE PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 Height BMI Waist Hip") -  | perl -lane 'print $F[0], "\t", $F[1], "\t", join("\t", @F[$#F-3..$#F]);' > /users/mturchin/data/ukbiobank_jun17/mturchin/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.top10resids.txt 

done





##zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' | head -n 10
##join <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' | sort) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim |  
##cat file1 | R -q -e "Data1 <- read.table(file('stdin'), header=F); Data2 <- Data1[,1]; Data3 <- Data1[,2:ncol(Data1)]; colnames(Data2) <- c(\"SNPID\"); colnames(Data3) <- c(\"SNPID\", ....\"POS\"?....); Data4 <- merge(Data2, Data3, by=\"SNPID\"); write.table(Data4, quote=FALSE, row.name=FALSE, col.name=FALSE);" | grep -v \> | 
zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.gz | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' > /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs
cmp <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs | sed 's/_/ /g' | awk '{ print $1 }') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | awk '{ print $2 }') | wc
paste <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs | sed 's/_/ /g' | awk '{ print $1 }') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | awk '{ print $2 }') | awk '{ if ($1 == $2) { print $0 } } ' | wc

NonSyn
Exonic
ExonicPlus
ExonicPlus20kb

for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);')`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v African`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);')`; do
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
			ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
				
			echo $i $j $ancestry1 $ancestry2 $ancestry3

			if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses ]; then
				mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses
			fi
			if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath ]; then
				mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath
			fi
			if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i ]; then
				mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i
			fi
			if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k ]; then
				mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k
			fi
			if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm ]; then
				mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm
			fi
		done
	done
done

~for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);')`; do
~	ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
~	ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
~	ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
~		
~	echo $i $j $ancestry1 $ancestry2 $ancestry3
~
~	mv /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/FORCE /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath
~
~done
~for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);')`; do
~	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);')`; do
~		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);')`; do
~			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
~			ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
~			ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
~				
~			echo $i $j $ancestry1 $ancestry2 $ancestry3 $k
~
~			ls -lrt /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ | awk '{ print $9 }' | perl -slane 'if ($F[0] =~ m/(.*InterPath\.)($iBash.*)/) { my $val1 = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/" . $1 . $2; my $val2 = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/" . $1 . "vs1." . $2; system("mv $val1 $val2");}' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 -iBash=$i -kBash=$k
~			ls -lrt /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i | awk '{ print $9 }' | perl -slane 'if ($F[0] =~ m/(.*InterPath\.)($iBash.*)/) { my $val1 = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/$iBash/" . $1 . $2; my $val2 = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/$iBash/" . $1 . "vs1." . $2; system("mv $val1 $val2");}' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 -iBash=$i -kBash=$k
~			ls -lrt /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k | awk '{ print $9 }' | perl -slane 'if ($F[0] =~ m/(.*InterPath\.)($iBash.*)/) { my $val1 = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/$iBash/$kBash/" . $1 . $2; my $val2 = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/$iBash/$kBash/" . $1 . "vs1." . $2; system("mv $val1 $val2");}' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 -iBash=$i -kBash=$k
~
~		done
~	done
~done
~			
~#			ls -lrt /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ | awk '{ print $9 }' | perl -slane 'if ($F[0] =~ m/(.*)(FORCE)(.*)/) { my $val1 = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/" . $1 . $2 . $3; my $val2 = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/" . $1 . "InterPath" . $3; system("mv $val1 $val2");}' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 -iBash=$i -kBash=$k
~#			ls -lrt /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i | awk '{ print $9 }' | perl -slane 'if ($F[0] =~ m/(.*)(FORCE)(.*)/) { my $val1 = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/$iBash/" . $1 . $2 . $3; my $val2 = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/$iBash/" . $1 . "InterPath" . $3; system("mv $val1 $val2");}' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 -iBash=$i -kBash=$k
~#			ls -lrt /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k | awk '{ print $9 }' | perl -slane 'if ($F[0] =~ m/(.*)(FORCE)(.*)/) { my $val1 = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/$iBash/$kBash/" . $1 . $2 . $3; my $val2 = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/$iBash/$kBash/" . $1 . "InterPath" . $3; system("mv $val1 $val2");}' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 -iBash=$i -kBash=$k

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
        ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
		
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3

#	zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | head -n 1 | perl -lane 'print join("\n", @F);' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.SNPIDs
#	paste <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.SNPIDs | sed 's/_/ /g' | awk '{ print $1 }') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim | awk '{ print $2 }') | awk '{ if ($1 == $2) { print $0 } } ' | wc
#	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.SNPIDs | wc
#	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim | wc
#	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.SNPIDs | awk '{ print $1 "\t" NR }' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.SNPIDsRowPos
#	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw | perl -lane 'print join("\t", @F[0..6]);' | sed 's/_/ /g' | awk '{ print $1 "\t" $2 }' | R -q -e "Data1 <- read.table(file('stdin'), header=T); colnames(Data1) <- c(\"FID\", \"IID2\"); \
#	Data2 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/mturchin/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.txt\", header=T); \ 
#	Data3 <- merge(Data1, Data2, by=\"FID\"); \	
#	write.table(Data3[,c(1,3:ncol(Data3))], \"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt\", quote=FALSE, row.name=FALSE, col.name=TRUE);" 
#	paste <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw | perl -lane 'print join("\t", @F[0..6]);' | sed 's/_/ /g' | awk '{ print $1 "\t" $2 }') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt | awk '{ print $1 }') | awk '{ if ($1 != $3) { print $0 } }' | wc
#	rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw
	
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt | awk '{ print $1 "\t" $2 }' | R -q -e "Data1 <- read.table(file('stdin'), header=T); colnames(Data1) <- c(\"FID\", \"IID2\"); \ 
	Data2 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/mturchin/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.top10resids.txt\", header=T); \ 
	Data3 <- merge(Data1, Data2, by=\"FID\"); \
	write.table(Data3[,c(1,3:ncol(Data3))], \"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.flashpca.top10resids.txt\", quote=FALSE, row.name=FALSE, col.name=TRUE);"
	paste <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt | awk '{ print $1 "\t" $2 }') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.flashpca.top10resids.txt | awk '{ print $1 }') | awk '{ if ($1 != $3) { print $0 } }' | wc

done

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);')`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
		
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3

	rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.txt /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.txt /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim | awk '{ print $1 "\t" $4 "\t" $4 "\t" $5 "\t" $6 }' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat 
	perl /users/mturchin/Software/20180611_annovar/table_annovar.pl /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat /users/mturchin/data/mturchin/20180611_annovar/humandb/ -buildver hg19 -out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar -remove -protocol refGene -operation g -nastring NA
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | grep exonic | grep -v ncRNA | grep -v unknown | awk '{ print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $11 }' | grep NM_ | perl -lane 'my @info1 = split(/:/, $F[4]); if ($info1[3] =~ m/c\.(\D)\d+(\D)/) { if ($1 eq $2) { ($F[2], $F[3]) = ($F[3], $F[2]); } } print $F[0], "\t", $F[1], "\t", $F[1], "\t", $F[2], "\t", $F[3] ' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.exonic.AAFix.pre.txt
	cat <(join -v 1 <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | awk '{ print $1 "_" $2 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 }' | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.exonic.AAFix.pre.txt | awk '{ print $1 "_" $2 }' | sort) | perl -lane 'print join("\t", @F[1..$#F]);') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.exonic.AAFix.pre.txt) | grep -v Start > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.exonic.AAFix.txt
	perl /users/mturchin/Software/20180611_annovar/table_annovar.pl /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.exonic.AAFix.txt /users/mturchin/data/mturchin/20180611_annovar/humandb/ -buildver hg19 -out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix -remove -protocol refGene -operation g -nastring NA 
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.txt | perl -F\\t -lane 'if ($. == 1) { %hash1; } if ($hash1{$F[6]}) { push(@{$hash1{$F[6]}}, $F[0] . ":" . $F[1] . "\t" . $F[5] . "," . $F[7] . "," . $F[8]); } else { $hash1{$F[6]} = [($F[0] . ":" . $F[1] . "\t" . $F[5] . "," . $F[7] . "," . $F[8])]; } if (eof()) { foreach my $gene1 (keys %hash1) { foreach my $snp1 (@{$hash1{$gene1}}) { print $gene1, "\t", $snp1; } } };' | sed 's/ /_/g ' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.txt
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.txt | perl -lane 'my @genes1 = split(/;/, $F[0]); my @info1 = split(/,/, $F[2]); my @dists1 = split(/;/, $info1[1]); if (($F[0] =~ m/;/) && ($info1[0] eq "intergenic")) { print $genes1[0], "\t", $F[1], "\t", $info1[0] . "_upstream", ",", $dists1[0], ",", $info1[2]; print $genes1[1], "\t", $F[1], "\t", $info1[0] . "_downstream", ",", $dists1[1], ",", $info1[2]; } elsif ($F[0] =~ m/;/) { print $genes1[0], "\t", $F[1], "\t", $info1[0], ",", $info1[1], ",", $info1[2]; print $genes1[1], "\t", $F[1], "\t", $info1[0], ",", $info[1], ",", $info1[2]; } else { print join("\t", @F); }' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.txt
	join <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.SNPIDsRowPos | sed 's/_/ /g' | awk '{ print $1 "\t" $3 }' | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.txt | awk '{ print $2 "\t" $1 "\t" $3 }' | sort -k 1,1) | awk '{ print $3 "\t" $1 "\t" $4 "\t" $2 }' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt 

done
	
#	cat /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | grep exonic | grep -v ncRNA | grep -v unknown | awk '{ print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $11 }' | grep NM_ | perl -lane 'my @info1 = split(/:/, $F[4]); if ($info1[3] =~ m/c\.(\D)\d+(\D)/) { if ($1 eq $2) { ($F[2], $F[3]) = ($F[3], $F[2]); } } print $F[0], "\t", $F[1], "\t", $F[1], "\t", $F[2], "\t", $F[3] ' > /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.exonic.AAFix.pre.txt
#	cat <(join -v 1 <(cat /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | awk '{ print $1 "_" $2 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 }' | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.exonic.AAFix.pre.txt | awk '{ print $1 "_" $2 }' | sort) | perl -lane 'print join("\t", @F[1..$#F]);') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.exonic.AAFix.pre.txt) > /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.exonic.AAFix.txt

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v African | grep Cari`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
		
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3

	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | wc
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.txt | wc	
	join <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | awk '{ print $1 "_" $2 }' | sort) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.txt | awk '{ print $1 "_" $2 }' | sort) | wc

done

#nonsynonymous, exonic, exonic + intronic + UTR, exonic + intronic + UTR + upstream/downstream + 20kb away
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);')`; do
	ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`; ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`; echo $pheno1 $ancestry1 $ancestry2 $ancestry3;

	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep nonsynonymous | sort -k 1,1 | perl -lane 'if ($. == 1) { @gene1; push(@gene1, $F[0]); push(@gene1, $F[3]); } else { if ($F[0] ne $gene1[0]) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); @gene1 = (); push(@gene1, $F[0]); push(@gene1, $F[3]); } else { push(@gene1, $F[3]); } if (eof()) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); } };' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.NonSyn.txt
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep exonic | sort -k 1,1 | perl -lane 'if ($. == 1) { @gene1; push(@gene1, $F[0]); push(@gene1, $F[3]); } else { if ($F[0] ne $gene1[0]) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); @gene1 = (); push(@gene1, $F[0]); push(@gene1, $F[3]); } else { push(@gene1, $F[3]); } if (eof()) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); } };' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.Exonic.txt
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep -E 'exonic|intronic|UTR' | sort -k 1,1 | perl -lane 'if ($. == 1) { @gene1; push(@gene1, $F[0]); push(@gene1, $F[3]); } else { if ($F[0] ne $gene1[0]) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); @gene1 = (); push(@gene1, $F[0]); push(@gene1, $F[3]); } else { push(@gene1, $F[3]); } if (eof()) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); } };' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus.txt
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep -E 'exonic|intronic|UTR|upstream|downstream' | sort -k 1,1 | perl -lane 'my @info1 = split(/,/, $F[2]); if ($info1[0] =~ m/intergenic/) { my @dists1 = split(/=/, $info1[1]); if ($dists1[1] <= 20000) { print join("\t", @F); } } else { print join("\t", @F); }' | perl -lane 'if ($. == 1) { @gene1; push(@gene1, $F[0]); push(@gene1, $F[3]); } else { if ($F[0] ne $gene1[0]) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); @gene1 = (); push(@gene1, $F[0]); push(@gene1, $F[3]); } else { push(@gene1, $F[3]); } if (eof()) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); } };' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus20kb.txt
	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.NonSyn.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(",", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 2) { print join("\t", @F); }' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.NonSyn.txt
	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.Exonic.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(",", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 2) { print join("\t", @F); }' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.Exonic.txt
	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(",", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 2) { print join("\t", @F); }' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.ExonicPlus.txt
	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus20kb.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(",", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 2) { print join("\t", @F); }' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.ExonicPlus20kb.txt; done

#	for k in `cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | perl -F\\\t -lane 'print $F[6];' | sort | uniq`; do
#		cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | perl -F\\t -slane 'if ($F[6] eq $gene1) { print $gene1, "\t", $F[0], ":", $F[1]; }' -- -gene1=$k
#	done > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.txt

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -E 'African|British|Indian'`; do
	ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
	ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3
	
	for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);')`; do
		join <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/British/British/mturchin20/$i/ukb_chrAll_v2.British.QCed.reqDrop.QCed.dropRltvs.PCAdrop.noX.${i}.Trans.ADD.assoc.linear.clumped.gz | awk '{ print $1 "_" $4 }' | sort) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.SNPIDsRowPos | sed 's/_/ /g' | awk '{ print $1 "\t" $3 }' | sed 's/:/_/g' | sort -k 1,1) | awk '{ print $2 }' | perl -ane 'if (eof()) { print $F[0], "\n"; } else { print $F[0], "," }' | xargs echo $i $i >> /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.British.British.Phenos.Trans.ADD.assoc.linear.clumped.wRowPos.txt 
	done > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.British.British.Phenos.Trans.ADD.assoc.linear.clumped.wRowPos.txt 
done 

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -E 'African|British|Indian'`; do
	ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
	ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3

	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.British.British.Phenos.Trans.ADD.assoc.linear.clumped.wRowPos.txt | perl -lane 'my @vals1 = split(/,/, $F[$#F]); print $F[0], "\t", scalar(@vals1);'

done 

#From https://blog.rstudio.com/2016/03/29/feather/, https://blog.dominodatalab.com/the-r-data-i-o-shootout/, https://stackoverflow.com/questions/1727772/quickly-reading-very-large-tables-as-dataframes

module load R/3.3.1; for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep Indian`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
		
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3

	R -q -e "library("data.table"); library("feather"); \ 
	ptm <- proc.time(); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); print(proc.time() - ptm); \ 
	Data3.mean <- apply(Data3, 2, mean); Data3.sd <- apply(Data3, 2, sd); Data3.sd[which(Data3.sd==0)] <- 1; Data3 <- t((t(Data3)-Data3.mean)/Data3.sd); \
	ptm <- proc.time(); Data3.cov <- 1/nrow(Data3) * tcrossprod(as.matrix(Data3)); print(proc.time() - ptm); \
	write.table(Data3.cov, \"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.cov.txt\", quote=FALSE, col.name=FALSE, row.name=FALSE);"
done 

#Indn: 13253.066; 

#Data3.cov2 <- 1/nrow(Data3) * tcrossprod(scale(as.matrix(Data3), FALSE, FALSE))
#	ptm <- proc.time(); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); print(proc.time() - ptm); 
#	Data4 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane \'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);\'', header=T);
#	write_feather(Data3, \"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.edit.fthr\"); \
#	print(proc.time() - ptm); \
#	write_feather(as.data.frame(Data3.cov), \"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.cov.fthr\"); \
#	write.table(Data3.cov, "/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.cov.txt", quote=FALSE, col.name=FALSE, row.name=FALSE);"

module load R/3.3.1; for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep Height`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep 20kb`; do
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt | wc | awk '{ print $1 }'`	
#			NumPaths=4672
			echo $i $ancestry1 $ancestry2 $ancestry3 $k
			
			if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways ]; then
				mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways
			fi
			if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k ]; then
				mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k
			fi
			
			for (( PathNum=1; PathNum <= $NumPaths; PathNum=PathNum+80 )); do
				sbatch -t 72:00:00 --mem 50g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${k}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${k}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
				echo -e "\nR -q -e \"library(\\\"data.table\\\"); \ 
				Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); \
				Pathways <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\\\", header=F); \
				for (i in $PathNum:($PathNum+79)) { \
					Pathways.Regions <- as.numeric(as.character(unlist(strsplit(as.character(Pathways[i,3]), \\\",\\\")))); \
					Data3.temp <- as.matrix(Data3)[,c(Pathways.Regions)]; \
					write.table(Data3.temp, paste(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Regions.c2.${k}.Pathways\\\", i, \\\".txt\\\", sep=\\\"\\\"), quote=FALSE, col.name=TRUE, row.name=FALSE); \
					};\"") 
			done; sleep 2
		done; 
	done;
done;		

#			R -q -e "library("data.table"); library("feather"); \ 
#			ptm <- proc.time(); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); print(proc.time() - ptm); \ 
#			Pathways <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\", header=F); \
#			for (i in 1:nrow(Pathways)) { \
#				print(i); \
#				Pathways.Regions <- as.numeric(as.character(unlist(strsplit(as.character(Pathways[i,3]), \",\")))); \
#				write.table(Data3.temp, paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Regions.c2.${k}.Pathways\", i, \".txt\", sep=\"\"), quote=FALSE, col.name=TRUE, row.name=FALSE); \

#From https://unix.stackexchange.com/questions/55392/in-bash-is-it-possible-to-use-an-integer-variable-in-the-loop-control-of-a-for
#Afr, Brit4k, Indn -- Hght, BMI / NonSyn, Exon, Exon+, Exon+20
#Afr: 36g (now probably 8g?); Brit4k: 10g; Indn: 12g  

#Moving to the pre-processed X matrix setup
#cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.SingleRun.vs1.cpp /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.SingleRun.vs2.cpp
#Including 2 lines of code for explicit covariate inclusion (from https://github.com/lorinanthony/MAPIT/blob/master/OpenMP%20Version/MAPIT_OpenMP.cpp, lines 148 & 149, but get rid of 'b.col(q+1) = trans(X.row(i));' part (null space against the single SNP as was part of original MAPIT setup), & changing input for Z matrix as well)
#cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.SingleRun.vs2.cpp /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.SingleRun.vs2.wCovs.vs1.cpp
#Getting rid of Kc version #20180815 NOTE -- never ended up using this 'vs3' for things, got a more 'official' version edited by Lorin eventually
##cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.SingleRun.vs2.wCovs.vs1.cpp /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.SingleRun.vs3.wCovs.vs1.cpp
##mv /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.SingleRun.vs3.wCovs.vs1.cpp /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/OLD1.InterPath.mtEdits.SingleRun.vs3.wCovs.vs1.cpp

#Original Vs1 Run (where issues and wonky plots were observed)
module load R/3.3.1; sleep 25200; for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep -E 'Height|BMI'`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep -v Plus`; do
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			NumSNPs=`zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | head -n 1 | perl -ane 'print scalar(@F);'` 
			NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt | wc | awk '{ print $1 }'`	
#			NumPaths=2
			
			echo $i $ancestry1 $ancestry2 $ancestry3 $k; for (( PathNum=1; PathNum <= $NumPaths; PathNum=PathNum+10 )); do
				sbatch -t 72:00:00 --mem 10g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
				echo -e "\nR -q -e \"library(\\\"data.table\\\"); library(\\\"doParallel\\\"); library(\\\"Rcpp\\\"); library(\\\"RcppArmadillo\\\"); library(\\\"RcppParallel\\\"); sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.SingleRun.vs2.cpp\\\"); neg.is.na <- Negate(is.na); \
				Pathways.Regions <- list(); cores = detectCores(); InterPath.output <- list(); InterPath.output\\\$Est <- c(); InterPath.output\\\$Eigenvalues <- c(); InterPath.output\\\$PVE <- c(); \
				Pathways <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\\\", header=F); \
				for (i in $PathNum:($PathNum+9)) { \
					Pathways.Regions[[1]] <- as.numeric(as.character(unlist(strsplit(as.character(Pathways[i,3]), \\\",\\\")))); \
					Data3 <- fread(paste(\\\"cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Regions.c2.${k}.Pathways\\\", as.character(i), \\\".txt\\\", sep=\\\"\\\"), header=T); \ 
					Data3.mean <- apply(Data3, 2, mean); Data3.sd <- apply(Data3, 2, sd); Data3 <- t((t(Data3)-Data3.mean)/Data3.sd); \
					Data3.cov <- as.matrix(read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.cov.txt\\\", header=F)); \ 
					InterPath.output.temp <- list(); X <- Data3; X.cov <- Data3.cov; rm(Data3); rm(Data3.cov); \
					Y <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.flashpca.top10resids.txt\\\", header=T); Y.Pheno <- Y\\\$$i; \ 
					Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; X.Pheno.noNAs <- X[neg.is.na(Y.Pheno),]; X.cov.Pheno.noNAs <- X.cov[neg.is.na(Y.Pheno),neg.is.na(Y.Pheno)]; \
					if (length(Pathways.Regions[[1]]) > 1) { K <- 1/nrow(X.Pheno.noNAs) * tcrossprod(as.matrix(X.Pheno.noNAs)); \
				InterPath.output.temp <- InterPath(t(X.Pheno.noNAs),Y.Pheno.noNAs,as.matrix(X.cov.Pheno.noNAs),K,Pathways.Regions,nrow(X.Pheno.noNAs),as.numeric(as.character($NumSNPs)),cores=cores); InterPath.output\\\$Est <- c(InterPath.output\\\$Est, InterPath.output.temp\\\$Est); InterPath.output\\\$Eigenvalues <- cbind(InterPath.output\\\$Eigenvalues, InterPath.output.temp\\\$Eigenvalues); InterPath.output\\\$PVE <- c(InterPath.output\\\$PVE, InterPath.output.temp\\\$PVE); } else { InterPath.output\\\$Est <- c(InterPath.output\\\$Est, NA); InterPath.output\\\$Eigenvalues <- cbind(InterPath.output\\\$Eigenvalues, rep(NA, length(Y.Pheno.noNAs))); InterPath.output\\\$PVE <- c(InterPath.output\\\$PVE, NA);};};\
				write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.PCadj.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
				write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.PCadj.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
				write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.PCadj.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")
			done; sleep 2
		done; 
	done; 
done
	
#		Data3.cov <- 1/nrow(Data3) * tcrossprod(as.matrix(Data3)); \
#		X.sub <- X[,as.numeric(as.character(Pathways.Regions))]; \
#		X.sub.cov <- 1/nrow(X.sub) * tcrossprod(as.matrix(X.sub));
#				print(c(dim(X), dim(X.cov), length(Y.Height), \\\"nana\\\", length(neg.is.na(Y.Height)))); \
#				Y.Height.noNAs <- Y.Height[!is.na(Y.Height)]; X.Height.noNAs <- X[!is.na(Y.Height),]; X.cov.Height.noNAs <- X.cov[!is.na(Y.Height),!is.na(Y.Height)]; \
#		InterPath.output <- InterPath(t(X),Y.Height,as.matrix(X.cov),Pathways.Regions); \
#				eval(parse(text=paste(\\\"Y.Pheno <- Y\\\$\\\", \\\"$i\\\", sep=\\\"\\\"))); \ 
#				Y.Height <- Y\\\$Height; \
#				Pathways.Regions <- list(); Counter1 = 1; \
#				for (i in $PathNum:($PathNum+10)) { Pathways.Regions[[Counter1]] <- as.numeric(as.character(unlist(strsplit(as.character(Pathways[i,3]), \\\",\\\")))); Counter1 = Counter1 + 1; }; \	
#				cores = detectCores(); InterPath.output <- InterPath(t(X.Height.noNAs),Y.Height.noNAs,as.matrix(X.cov.Height.noNAs),Pathways.Regions,cores=cores); \
#		sourceCpp("/users/mturchin/LabMisc/RamachandranLab/InterPath/InterPath.cpp");
#		Pathways.Regions[[1]] <- c(1,2,3,4,5)
#		Pathways.Regions[[2]] <- c(4,7,8,1,3,4,6)
#		InterPath.output <- InterPath(t(X),Y.Height,as.matrix(X.cov),Pathways.Regions)	
#		InterPath.output <- InterPath(t(X),Y.Height[1:nrow(X)],Pathways.Regions)	
#		InterPath.output <- InterPath(t(X2),Y.Height,Pathways.Regions)	
#		ptm <- proc.time(); 
#		Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane \'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] =~ m/HET/) { $PH = 1 } else { push(@vals1, $i); } } } print join("t", @F[@vals1]);\'', header=T); 
#		print(proc.time() - ptm); 
#					AllPos1 <- c(); for (i in $PathNum:($PathNum+1)) { AllPos1 <- c(AllPos1, as.numeric(as.character(unlist(strsplit(as.character(Pathways[i,3]), \\\",\\\"))))); }; \
#					Data3.temp <- matrix(0, ncol=ncol(Data3), nrow=nrow(Data3)); Data3.temp[,AllPos1] <- Data3[,AllPos1]; Data3 <- Data3.temp; rm(Data3.temp); \
#					K <- 1/nrow(X.Pheno.noNAs) * tcrossprod(as.matrix(X.Pheno.noNAs[,Pathways.Regions[[1]]])); \
#					Y <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt\\\", header=T); Y.Pheno <- Y\\\$$i; \ 

#From https://stackoverflow.com/questions/8903239/how-to-calculate-time-difference-in-bash-script
for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep Hip`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);')`; do
			SECONDS=0;
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
			ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt | wc | awk '{ print $1 }'`	
#			NumPaths=20
			
			echo $i $ancestry1 $ancestry2 $ancestry3 $k
			for (( PathNum=1; PathNum <= $NumPaths; PathNum=PathNum+10 )); do
				R -q -e "library(\"CompQuadForm\"); neg.is.na <- Negate(is.na); \
				Pathways <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\", header=F); \ 
				InterPath.output.Est <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Paths${PathNum}.Est.txt\", header=F); \
				InterPath.output.Eigenvalues <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Paths${PathNum}.Eigenvalues.txt\", header=F); \
				InterPath.output.PVE <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Paths${PathNum}.PVE.txt\", header=F); \
				Results1 <- c(); Counter1 <- 1; for (i in $PathNum:($PathNum+9)) { \
					if (neg.is.na(InterPath.output.Est[Counter1,1])) { \ 
						Lambda <- sort(InterPath.output.Eigenvalues[,Counter1], decreasing=TRUE); \
						Davies.Output <- davies(InterPath.output.Est[Counter1,1], lambda=Lambda, acc=1e-8); \
						pVal <- 2*min(1-Davies.Output\$Qq, Davies.Output\$Qq); \
						Results1 <- rbind(Results1, c(as.character(Pathways[i,1]), InterPath.output.Est[Counter1,1], InterPath.output.PVE[Counter1,1], pVal)); \ 
						Counter1 = Counter1 + 1; \
					}; \
				}; write.table(Results1, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);"
			done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.AllPaths.Results.txt.pre.gz
			duration=$SECONDS; echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
		done;
	done;
done;

#	lambda = sort(vc.mod$Eigenvalues[,i],decreasing = T)
#	Davies_Method = davies(vc.mod$Est[i], lambda = lambda, acc=1e-8)
#	pvals[i] = 2*min(1-Davies_Method$Qq,Davies_Method$Qq)
#	names(pvals)[i] = names(vc.ts[i])
#	ptm <- proc.time(); ... write(proc.time() - ptm, stderr());"

#				R -q -e "library("CompQuadForm"); 
#				Pathways <- read.table("/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.ExonicPlus20kb.txt", header=F);  
#				InterPath.output.Est <- read.table("/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ExonicPlus20kb/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.Height.ExonicPlus20kb.Paths1.Est.txt", header=F); 
#				InterPath.output.Eigenvalues <- read.table("/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ExonicPlus20kb/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.Height.ExonicPlus20kb.Paths1.Eigenvalues.txt", header=F); 
#				InterPath.output.PVE <- read.table("/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ExonicPlus20kb/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.Height.ExonicPlus20kb.Paths1.PVE.txt", header=F); 
#				Results1 <- c(); 
#				Counter1 <- 1; for (i in 1:(1+10)) { 
#					Lambda <- sort(InterPath.output.Eigenvalues[,Counter1], decreasing=TRUE); 
#					Davies.Output <- davies(InterPath.output.Est[Counter1,1], lambda=Lambda, acc=1e-8); 
#					pVal <- 2*min(1-Davies.Output$Qq, Davies.Output$Qq); 
#					print(c(Pathways[i,1], InterPath.output.Est[Counter1,1], InterPath.output.PVE[Counter1,1], pVal));  
#					Results1 <- rbind(Results1, c(as.character(Pathways[i,1]), InterPath.output.Est[Counter1,1], InterPath.output.PVE[Counter1,1], pVal)); 
#					Counter1 = Counter1 + 1; 
#				};

#for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);')`; do
for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep Height`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
			ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			
			echo $i $ancestry1 $ancestry2 $ancestry3 $k

			R -q -e "library("RColorBrewer"); Data1 <- read.table("/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.NonSyn.AllPaths.Results.txt.pre.gz", header=F); 
			Data2 <- read.table("/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.Exonic.AllPaths.Results.txt.pre.gz", header=F); 
			Data3 <- read.table("/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.ExonicPlus.AllPaths.Results.txt.pre.gz", header=F); 
			Data4 <- read.table("/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.ExonicPlus20kb.AllPaths.Results.txt.pre.gz", header=F); 
			png("/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.AllStrats.AllPaths.Results.pre.QQplots.Recal.vs1.png", width=2000, height=2000, res=300); par(mar=c(5,5,4,2)); 
			xVals1 <- seq(1/nrow(Data1), 1, by=1/nrow(Data1)); xVals2 <- seq(1/nrow(Data2), 1, by=1/nrow(Data2)); xVals3 <- seq(1/nrow(Data3), 1, by=1/nrow(Data3)); xVals4 <- seq(1/nrow(Data4), 1, by=1/nrow(Data4)); 
			xlimMax <- max(c(-log10(xVals1), -log10(xVals2), -log10(xVals3), -log10(xVals4))); ylimMax <- max(c(-log10(Data1[,4]), -log10(Data2[,4]), -log10(Data3[,4]), -log10(Data4[,4]))); 
			plot(-log10(xVals1[order(xVals1, decreasing=TRUE)]), -log10(Data1[order(Data1[,4], decreasing=TRUE),4]), main="InterPath Prelim Results: $ancestry2 ($i)", xlab="-log10(Expected pValues)", ylab="-log10(Observed pValues)", xlim=c(0,max(-log10(xVals4[order(xVals4, decreasing=TRUE)]))), ylim=c(0,ylimMax), type="b", pch=16, col=brewer.pal(12, "Paired")[7], cex=1.5, cex.main=1.5, cex.axis=1.5, cex.lab=1.5); 
			points(-log10(xVals2[order(xVals2, decreasing=TRUE)]), -log10(Data2[order(Data2[,4], decreasing=TRUE),4]), type="b", pch=16, col=brewer.pal(12, "Paired")[3], cex=1.5);  
			points(-log10(xVals3[order(xVals3, decreasing=TRUE)]), -log10(Data3[order(Data3[,4], decreasing=TRUE),4]), type="b", pch=16, col=brewer.pal(12, "Paired")[9], cex=1.5);  
			points(-log10(xVals4[order(xVals4, decreasing=TRUE)]), -log10(Data4[order(Data4[,4], decreasing=TRUE),4]), type="b", pch=16, col=brewer.pal(12, "Paired")[1],  cex=1.5);  
			legend("topleft", c("NonSyn", "Exonic", "ExonicPlus", "ExonicPlus20kb"), pch=c(16,16,16,16), col=c(brewer.pal(12, "Paired")[7], brewer.pal(12, "Paired")[3], brewer.pal(12, "Paired")[9], brewer.pal(12, "Paired")[1]), bg="transparent"); 
			abline(0,1, col="BLACK"); 
			dev.off();"
	done;
done;

7, 3, 9, 1 
#PhenoColors <- c(brewer.pal(12, "Paired")[1], brewer.pal(12, "Paired")[3], brewer.pal(12, "Paired")[5], brewer.pal(12, "Paired")[7], brewer.pal(12, "Paired")[9], brewer.pal(12, "Paired")[12])
#			R -q -e "library(\"RColorBrewer\"); Data1 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.${i}.NonSyn.AllPaths.Results.txt.pre.gz\", header=F); \
#			pvals <- Data1[,4]; q = qchisq(pvals,df=1,lower.tail = FALSE); summary(q); gs = median(q); m = 0.4549; summary((q/gs)*m); pvals = pchisq((q/gs)*m,df=1,lower.tail = FALSE);  Data1[,4] <- pvals; 
#			pvals <- Data2[,4]; q = qchisq(pvals,df=1,lower.tail = FALSE); summary(q); gs = median(q); m = 0.4549; summary((q/gs)*m); pvals = pchisq((q/gs)*m,df=1,lower.tail = FALSE);  Data2[,4] <- pvals; 
#			pvals <- Data3[,4]; q = qchisq(pvals,df=1,lower.tail = FALSE); summary(q); gs = median(q); m = 0.4549; summary((q/gs)*m); pvals = pchisq((q/gs)*m,df=1,lower.tail = FALSE);  Data3[,4] <- pvals; 
#			pvals <- Data4[,4]; q = qchisq(pvals,df=1,lower.tail = FALSE); summary(q); gs = median(q); m = 0.4549; summary((q/gs)*m); pvals = pchisq((q/gs)*m,df=1,lower.tail = FALSE);  Data4[,4] <- pvals; 

#for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep Height`; do
#for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);')`; do
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -E 'British|Indian'`; do
		ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
		ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
		
		echo $i $ancestry1 $ancestry2 $ancestry3 $k

		R -q -e "library(\"RColorBrewer\"); Data1 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.AllPaths.Results.txt.pre.gz\", header=F); \
		Data2 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.Exonic.AllPaths.Results.txt.pre.gz\", header=F); \
		Data3 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/BMI/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.BMI.NonSyn.AllPaths.Results.txt.pre.gz\", header=F); \
		Data4 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/BMI/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.BMI.Exonic.AllPaths.Results.txt.pre.gz\", header=F); \
		png(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.HeightBMI.AllStrats.AllPaths.Results.pre.QQplots.vs1.png\", width=2000, height=2000, res=300);  par(mar=c(5,5,4,2)); \
		xVals1 <- seq(1/nrow(Data1), 1, by=1/nrow(Data1)); xVals2 <- seq(1/nrow(Data2), 1, by=1/nrow(Data2)); xVals3 <- seq(1/nrow(Data3), 1, by=1/nrow(Data3)); xVals4 <- seq(1/nrow(Data4), 1, by=1/nrow(Data4)); \
		xlimMax <- max(c(-log10(xVals1), -log10(xVals2), -log10(xVals3), -log10(xVals4))); ylimMax <- max(c(-log10(Data1[,4]), -log10(Data2[,4]), -log10(Data3[,4]), -log10(Data4[,4]))); \
		plot(-log10(xVals1[order(xVals1, decreasing=TRUE)]), -log10(Data1[order(Data1[,4], decreasing=TRUE),4]), main=\"InterPath Prelim Results: $ancestry2\", xlab=\"-log10(Expected pValues)\", ylab=\"-log10(Observed pValues)\", xlim=c(0,xlimMax), ylim=c(0,ylimMax), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[7], cex=1.5, cex.main=1.5, cex.axis=1.5, cex.lab=1.5); \
		points(-log10(xVals2[order(xVals2, decreasing=TRUE)]), -log10(Data2[order(Data2[,4], decreasing=TRUE),4]), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[3], cex=1.5); \ 
		points(-log10(xVals3[order(xVals3, decreasing=TRUE)]), -log10(Data3[order(Data3[,4], decreasing=TRUE),4]), type=\"b\", pch=17, col=brewer.pal(12, \"Paired\")[7], cex=1.5); \ 
		points(-log10(xVals4[order(xVals4, decreasing=TRUE)]), -log10(Data4[order(Data4[,4], decreasing=TRUE),4]), type=\"b\", pch=17, col=brewer.pal(12, \"Paired\")[3], cex=1.5); \ 
		legend(\"topleft\", c(\"Height: NonSyn\", \"Height: Exonic\", \"BMI: NonSyn\", \"BMI: Exonic\"), pch=c(16,16,17,17), col=c(brewer.pal(12, \"Paired\")[7], brewer.pal(12, \"Paired\")[3], brewer.pal(12, \"Paired\")[7], brewer.pal(12, \"Paired\")[3]), bg=\"transparent\"); \
		abline(0,1, col=\"BLACK\"); \
		dev.off();"
done;

#From MacBook Air
#mkdir /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath 
#mkdir /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Rnd1 
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.Height.*.AllPaths.Results.txt.pre.gz /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Rnd1/.
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/*/*/mturchin20/Analyses/InterPath/*AllStrats.AllPaths.Results.pre.QQplots.vs1.png /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Rnd1/.

for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);')`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);')`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep Exonic | grep ExonicPlus | grep 20kb`; do
			SECONDS=0;
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
			ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			
			echo $i $j $k $ancestry1 $ancestry2 $ancestry3 
			
			if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/PrevRuns ]; then
				mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/PrevRuns
			fi

			cd  /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k 
			tar -cvzf /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/PrevRuns/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.AllPaths.Est.tar.gz ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Paths*.Est.txt 
			tar -cvzf /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/PrevRuns/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.AllPaths.Eigenvalues.tar.gz ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Paths*.Eigenvalues.txt 
			tar -cvzf /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/PrevRuns/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.AllPaths.PVE.tar.gz ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Paths*.PVE.txt 

			duration=$SECONDS; echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
		done;
	done;
done;

for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);')`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);')`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);')`; do
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
			ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			
			echo $i $j $k $ancestry1 $ancestry2 $ancestry3 
		
			rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Paths*.Est.txt
			rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Paths*.Eigenvalues.txt
			rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Paths*.PVE.txt

		done;
	done;
done;

#GWAS Top Hits Pathway Setup
mkdir /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/slurm
mkdir /users/mturchin/data/ukbiobank_jun17/subsets/Indian/Indian/mturchin20/Analyses/InterPath/slurm
mkdir /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran4000/mturchin20/Analyses/InterPath/slurm

module load R/3.3.1; for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep Indian`; do
	ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
	NumSNPs=`zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | head -n 1 | perl -ane 'print scalar(@F);'` 
	echo $i $ancestry1 $ancestry2 $ancestry3 $k
	
	sbatch -t 72:00:00 --mem 120g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.British.British.Phenos.Trans.ADD.assoc.linear.clumped.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.British.British.Phenos.Trans.ADD.assoc.linear.clumped.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
	echo -e "\nR -q -e \"library(\\\"data.table\\\"); library(\\\"doParallel\\\"); library(\\\"Rcpp\\\"); library(\\\"RcppArmadillo\\\"); library(\\\"RcppParallel\\\"); sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.SingleRun.vs2.cpp\\\"); neg.is.na <- Negate(is.na); \
	Data1 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); \
	Pathways.Regions <- list(); cores = detectCores(); InterPath.output <- list(); InterPath.output\\\$Est <- c(); InterPath.output\\\$Eigenvalues <- c(); InterPath.output\\\$PVE <- c(); \
	Pathways <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.British.British.Phenos.Trans.ADD.assoc.linear.clumped.wRowPos.txt\\\", header=F); \
	for (i in 1:4) { \
		Pathways.Regions[[1]] <- as.numeric(as.character(unlist(strsplit(as.character(Pathways[i,3]), \\\",\\\")))); \
		Data3 <- as.matrix(Data1, nrow=nrow(Data1), ncol=ncol(Data1))[,Pathways.Regions[[1]]]; \ 
		Data3.mean <- apply(Data3, 2, mean); Data3.sd <- apply(Data3, 2, sd); Data3 <- t((t(Data3)-Data3.mean)/Data3.sd); \
		Data3.cov <- as.matrix(read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.cov.txt\\\", header=F)); \ 
		InterPath.output <- list(); X <- Data3; X.cov <- Data3.cov; rm(Data3); rm(Data3.cov); \
		Y <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt\\\", header=T); Pheno <- as.character(Pathways[i,1]); eval(parse(text=paste(\\\"Y.Pheno <- Y\\\$\\\", Pheno, sep=\\\"\\\"))); \ 
		Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; X.Pheno.noNAs <- X[neg.is.na(Y.Pheno),]; X.cov.Pheno.noNAs <- X.cov[neg.is.na(Y.Pheno),neg.is.na(Y.Pheno)]; \
		if (length(Pathways.Regions[[1]]) > 1) { K <- 1/nrow(X.Pheno.noNAs) * tcrossprod(as.matrix(X.Pheno.noNAs)); \
		InterPath.output <- InterPath(t(X.Pheno.noNAs),Y.Pheno.noNAs,as.matrix(X.cov.Pheno.noNAs),K,Pathways.Regions,nrow(X.Pheno.noNAs),as.numeric(as.character($NumSNPs)),cores=cores); \
		write.table(InterPath.output\\\$Est, paste(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.InterPath.vs1.British.British.\\\", as.character(Pheno), \\\".Trans.ADD.assoc.linear.clumped.Est.txt\\\", sep=\\\"\\\"), quote=FALSE, row.name=FALSE, col.name=FALSE); \
		write.table(InterPath.output\\\$Eigenvalues, paste(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.InterPath.vs1.British.British.\\\", as.character(Pheno), \\\".Trans.ADD.assoc.linear.clumped.Eigenvalues.txt\\\", sep=\\\"\\\"), quote=FALSE, row.name=FALSE, col.name=FALSE); \
		write.table(InterPath.output\\\$PVE, paste(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.InterPath.vs1.British.British.\\\", as.character(Pheno), \\\".Trans.ADD.assoc.linear.clumped.PVE.txt\\\", sep=\\\"\\\"), quote=FALSE, row.name=FALSE, col.name=FALSE); } else { InterPath.output\\\$Est <- c(InterPath.output\\\$Est, NA); InterPath.output\\\$Eigenvalues <- cbind(InterPath.output\\\$Eigenvalues, rep(NA, length(Y.Pheno.noNAs))); InterPath.output\\\$PVE <- c(InterPath.output\\\$PVE, NA);};};\"");
done

#		InterPath.output.temp <- list();...Force.output.temp <- InterPath(t(...
#		InterPath.output\\\$Est <- c(InterPath.output\\\$Est, InterPath.output.temp\\\$Est); InterPath.output\\\$Eigenvalues <- cbind(InterPath.output\\\$Eigenvalues, InterPath.output.temp\\\$Eigenvalues); InterPath.output\\\$PVE <- c(InterPath.output\\\$PVE, InterPath.output.temp\\\$PVE); \ 

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep Indian`; do
	ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
	ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
	NumSNPs=`zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | head -n 1 | perl -ane 'print scalar(@F);'` 
	echo $i $ancestry1 $ancestry2 $ancestry3 $k
	
	R -q -e "library(\"CompQuadForm\"); neg.is.na <- Negate(is.na); \
	Pathways <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.British.British.Phenos.Trans.ADD.assoc.linear.clumped.wRowPos.txt\", header=F); \
	Results1 <- c(); \ 
	for (i in 1:4) { \
		InterPath.output.Est <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.InterPath.vs1.British.British.\", as.character(Pathways[i,1]), \".Trans.ADD.assoc.linear.clumped.Est.txt\", sep=\"\"), header=F); \
		InterPath.output.Eigenvalues <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.InterPath.vs1.British.British.\", as.character(Pathways[i,1]), \".Trans.ADD.assoc.linear.clumped.Eigenvalues.txt\", sep=\"\"), header=F); \
		InterPath.output.PVE <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.InterPath.vs1.British.British.\", as.character(Pathways[i,1]), \".Trans.ADD.assoc.linear.clumped.PVE.txt\", sep=\"\"), header=F); \
		if (neg.is.na(InterPath.output.Est)) { \ 
			Lambda <- sort(InterPath.output.Eigenvalues[,1], decreasing=TRUE); \
			Davies.Output <- davies(InterPath.output.Est, lambda=Lambda, acc=1e-8); \
			pVal <- 2*min(1-Davies.Output\$Qq, Davies.Output\$Qq); \
			Results1 <- rbind(Results1, c(as.character(Pathways[i,1]), InterPath.output.Est, InterPath.output.PVE, pVal)); }; \ 
	}; write.table(Results1, file=\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.InterPath.vs1.British.British.AllPhenos.Trans.ADD.assoc.linear.clumped.Results.txt\", quote=FALSE, col.name=FALSE, row.name=FALSE);"
done;
	
#	done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.${i}.${k}.AllPaths.Results.txt.pre.gz
#		print(InterPath.output.Eigenvalues[,1]); \

#With Explicit Covariates Matrix Z setup
module load R/3.3.1; sleep 25200; for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep Height`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep NonSyn`; do
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`; NumSNPs=`zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | head -n 1 | perl -ane 'print scalar(@F);'` 
			NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt | wc | awk '{ print $1 }'`	
#			NumPaths=952
			echo $i $ancestry1 $ancestry2 $ancestry3 $k
			
			for (( PathNum=1; PathNum <= $NumPaths; PathNum=PathNum+10 )); do
				sbatch -t 72:00:00 --mem 12g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
				echo -e "\nR -q -e \"library(\\\"data.table\\\"); library(\\\"doParallel\\\"); library(\\\"Rcpp\\\"); library(\\\"RcppArmadillo\\\"); library(\\\"RcppParallel\\\"); sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.SingleRun.vs2.wCovs.vs1.cpp\\\"); neg.is.na <- Negate(is.na); \
				Y <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt\\\", header=T); Covars <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.txt\", header=T); Pathways <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\\\", header=F); \
				Pathways.Regions <- list(); cores = detectCores(); InterPath.output <- list(); InterPath.output\\\$Est <- c(); InterPath.output\\\$Eigenvalues <- c(); InterPath.output\\\$PVE <- c(); Y.Pheno <- Y\\\$$i; \
				for (i in $PathNum:($PathNum+9)) { \
					Data3 <- fread(paste(\\\"cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Regions.c2.${k}.Pathways\\\", as.character(i), \\\".txt\\\", sep=\\\"\\\"), header=T); Data3.cov <- as.matrix(read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.cov.txt\\\", header=F)); Pathways.Regions[[1]] <- as.numeric(as.character(unlist(strsplit(as.character(Pathways[i,3]), \\\",\\\")))); \ 
					Data3.mean <- apply(Data3, 2, mean); Data3.sd <- apply(Data3, 2, sd); Data3 <- t((t(Data3)-Data3.mean)/Data3.sd); \ 
					InterPath.output.temp <- list(); X <- Data3; X.cov <- Data3.cov; rm(Data3); rm(Data3.cov); Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; X.Pheno.noNAs <- X[neg.is.na(Y.Pheno),]; X.cov.Pheno.noNAs <- X.cov[neg.is.na(Y.Pheno),neg.is.na(Y.Pheno)]; Z <- Covars[neg.is.na(Y.Pheno),(ncol(Covars)-9):ncol(Covars)]; \
					if (length(Pathways.Regions[[1]]) > 1) { K <- 1/nrow(X.Pheno.noNAs) * tcrossprod(as.matrix(X.Pheno.noNAs)); \
						InterPath.output.temp <- InterPath(t(X.Pheno.noNAs),Y.Pheno.noNAs,as.matrix(X.cov.Pheno.noNAs),K,Z,Pathways.Regions,nrow(X.Pheno.noNAs),as.numeric(as.character($NumSNPs)),cores=cores); InterPath.output\\\$Est <- c(InterPath.output\\\$Est, InterPath.output.temp\\\$Est); InterPath.output\\\$Eigenvalues <- cbind(InterPath.output\\\$Eigenvalues, InterPath.output.temp\\\$Eigenvalues); InterPath.output\\\$PVE <- c(InterPath.output\\\$PVE, InterPath.output.temp\\\$PVE); } else { InterPath.output\\\$Est <- c(InterPath.output\\\$Est, NA); InterPath.output\\\$Eigenvalues <- cbind(InterPath.output\\\$Eigenvalues, rep(NA, length(Y.Pheno.noNAs))); InterPath.output\\\$PVE <- c(InterPath.output\\\$PVE, NA);};};\
				write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.wCov.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
				write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.wCov.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
				write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.wCov.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")
			done; sleep 2
		done; 
	done; 
done
	














~~~
#20180615
(InterPath) [  mturchin@node621  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.CompHRC.flipSNPs.rsIDs | wc
      0       0       0
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' | wc        
 379088  379088 4722679
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' | head -n 10
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
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | head -n 10
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
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' | head -n 10      1:729632_T                                                                                                                                                                                                                      1:752721_G                                                                                                                                                                                                                      1:754105_T                                                                                                                                                                                                                      1:756604_G                                                                                                                                                                                                                      1:759036_A                                                                                                                                                                                                                      1:761147_C
1:767096_G
1:768448_A
1:772927_T
1:779322_G
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' | tail -n 10
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
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | tail -n 10
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
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' > /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cmp /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | awk '{ print $2 }')
/users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs /dev/fd/63 differ: char 9, line 1
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$paste /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | awk '{ print $2 }') | head -n 10
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
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cmp <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs | sed 's/_/ /g' | awk '{ print $1 }') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | awk '{ print $2 }') | wc
      0       0       0
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$paste <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs | sed 's/_/ /g' | awk '{ print $1 }') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | awk '{ print $2 }') | awk '{ if ($1 == $2) { print $0 } } ' | wc
 373401  746802 8690114
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.SNPIDs | wc
 373401  373401 5091859
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim | wc
 373401 2240406 10930520
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | perl -F\\t -lane 'print $#F;' | sort | uniq -c             
 373402 9
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | head -n 10
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
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -ane 'print join("\n", @F[2..$#F]), "\n";' | wc
 444687  444687 2767870
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -ane 'print join("\n", @F[2..$#F]), "\n";' | sort | uniq | wc
  21095   21095  144468
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | perl -F\\t -lane 'print $F[6];' | sort | uniq | wc
  33331   33331  341748
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.IntergenicSplit.wRowPos.txt | sed 's/,/ /g' | awk '{ print $3 }' | sort | uniq -c
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
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.txt | wc                                                 
 556099 1663603 26951222
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.txt | sed 's/,/ /g' | awk '{ print $3 }' | sort | uniq -c
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
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep nonsynonymous | sort -k 1,1 | awk '{ print $1 }' | sort | uniq -c | awk '{ if ($1 > 1) { print $0 } } ' | wc
    503    1006    7288
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep nonsynonymous | sort -k 1,1 | awk '{ print $1 }' | sort | uniq -c | awk '{ if ($1 > 1) { print $0 } } ' | R -q -e "Data1 <- read.table(file('stdin'), header=F); table(Data1[,1]);"
> Data1 <- read.table(file('stdin'), header=F); table(Data1[,1]);

  2   3   4   5   6   7   9  10  11  17 
344  97  33  10   8   5   2   2   1   1 
> 
> 
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep nonsynonymous | sort -k 1,1 | head -n 10
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
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep -E 'exonic|intronic|UTR|upstream|downstream' | sort -k 1,1 | perl -lane 'my @info1 = split(/,/, $F[2]); if ($info1[0] =~ m/intergenic/) { my @dists1 = split(/=/, $info1[1]); if ($dists1[1] <= 20000) { print join("\t", @F); } } else { print join("\t", @F); }' | sed 's/,/ /g' | sed 's/=/ /'g | grep dist | awk '{ print $5 }' | sort | uniq -c | sort -g -k 2,2 | head -n 10
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
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep -E 'exonic|intronic|UTR|upstream|downstream' | sort -k 1,1 | perl -lane 'my @info1 = split(/,/, $F[2]); if ($info1[0] =~ m/intergenic/) { my @dists1 = split(/=/, $info1[1]); if ($dists1[1] <= 20000) { print join("\t", @F); } } else { print join("\t", @F); }' | sed 's/,/ /g' | sed 's/=/ /'g | grep dist | awk '{ print $5 }' | sort | uniq -c | sort -rg -k 2,2 | head -n 10
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
[  mturchin@node833  ~/LabMisc/RamachandranLab/InterPath]$for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
>         ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
>         ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
>         ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
>
>         echo $pheno1 $ancestry1 $ancestry2 $ancestry3
>
>         R -q -e "library(\"data.table\"); library(\"feather\"); \
>         ptm <- proc.time(); \
>         Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane \'if (\$. == 1) { @vals1; for (my \$i = 6; \$i <= \$#F; \$i++) { if (\$F[\$i] =~ m/HET/) { \$PH = 1 } else { push(@vals1, \$i); } } } print join(\"\t\", @F[@vals1]);\'', header=T); \
>         print(proc.time() - ptm); \
>         Data3.cov <- 1/nrow(Data3) * tcrossprod(as.matrix(Data3)); \
>         print(proc.time() - ptm); \
>         write_feather(Data3, \"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.edit.fthr\"); \
>         print(proc.time() - ptm); \
>         write_feather(as.data.frame(Data3.cov), \"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.edit.cov.fthr\"); \
>         print(proc.time() - ptm);"
> done
African African Afr
> library("data.table"); library("feather");         ptm <- proc.time();         Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane \'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] =~ m/HET/) { $PH = 1 } else { push(@vals1, $i); } } } print join("\t", @F[@vals1]);\'', header=T);         print(proc.time() - ptm);         Data3.cov <- 1/nrow(Data3) * tcrossprod(as.matrix(Data3));         print(proc.time() - ptm);         write_feather(Data3, "/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.edit.fthr");         print(proc.time() - ptm);         write_feather(as.data.frame(Data3.cov), "/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.edit.cov.fthr");         print(proc.time() - ptm);
   user  system elapsed
369.874   6.506 376.312
    user   system  elapsed
1094.043    9.185 1103.273
    user   system  elapsed
1097.186   14.257 1125.467
    user   system  elapsed
1097.442   14.340 1126.362
>
>
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath/Vs1]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.fam | head -n 10
1001592_1001592 1001592_1001592 0       0       0       -9
1002560_1002560 1002560_1002560 0       0       0       -9
1003036_1003036 1003036_1003036 0       0       0       -9
1004593_1004593 1004593_1004593 0       0       0       -9
1008167_1008167 1008167_1008167 0       0       0       -9
1009965_1009965 1009965_1009965 0       0       0       -9
1010953_1010953 1010953_1010953 0       0       0       -9
1012491_1012491 1012491_1012491 0       0       0       -9
1013297_1013297 1013297_1013297 0       0       0       -9
1013636_1013636 1013636_1013636 0       0       0       -9
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath/Vs1]$zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.gz | perl -lane 'print join("\t", @F[0..6]);' | sed 's/_/ /g' | awk '{ print $1 "\t" $2 }' | head -n 10
FID     IID
1001592 1001592
1002560 1002560
1003036 1003036
1004593 1004593
1008167 1008167
1009965 1009965
1010953 1010953
1012491 1012491
1013297 1013297
[  mturchin@node621  ~/LabMisc/RamachandranLab/InterPath]$zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.gz | perl -lane 'print join("\t", @F[0..6]);' | head -n 10
FID     IID     PAT     MAT     SEX     PHENOTYPE       1:729632_T
1001592_1001592 1001592_1001592 0       0       0       -9      0
1002560_1002560 1002560_1002560 0       0       0       -9      0
1003036_1003036 1003036_1003036 0       0       0       -9      0
1004593_1004593 1004593_1004593 0       0       0       -9      0
1008167_1008167 1008167_1008167 0       0       0       -9      0
1009965_1009965 1009965_1009965 0       0       0       -9      0
1010953_1010953 1010953_1010953 0       0       0       -9      0
1012491_1012491 1012491_1012491 0       0       0       -9      0
1013297_1013297 1013297_1013297 0       0       0       -9      1
[  mturchin@node621  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.Phenos.Transformed.txt | head -n 10
FID IID Height BMI Waist Hip
1001592 1001592 0.118664115101393 1.17841032424695 1.24353357925829 1.89450422410686
1002560 1002560 0.734183319747881 2.04290941388721 1.4416885105281 1.8805781156021
1003036 1003036 0.217969859995042 1.82585180750194 1.85728627141396 1.94552746926506
1004593 1004593 -0.0567068367051468 1.17620744874895 1.13275845670766 1.37836394899338
1008167 1008167 0.322764805585384 0.15175399426104 0.0479378611405927 0.0196024993412122
1009965 1009965 -1.17248489348928 1.28457721261652 1.28787234250449 1.42107036024099
1010953 1010953 -2.3492759185006 0.600994502398072 0.498626444501437 -0.166021390280289
1012491 1012491 1.0238031036417 -0.370391689731401 0.587105771979022 0.304346330032523
1013297 1013297 -0.157620896463327 0.509608883879002 0.526718205408012 0.130864454904108
[  mturchin@node621  ~/LabMisc/RamachandranLab/InterPath]$zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.gz | perl -lane 'print join("\t", @F[0..6]);' | tail -n 10
6008418_6008418 6008418_6008418 0       0       0       -9      1
6010753_6010753 6010753_6010753 0       0       0       -9      0
6011864_6011864 6011864_6011864 0       0       0       -9      0
6014349_6014349 6014349_6014349 0       0       0       -9      0
6018726_6018726 6018726_6018726 0       0       0       -9      0
6021711_6021711 6021711_6021711 0       0       0       -9      0
6023159_6023159 6023159_6023159 0       0       0       -9      0
6024552_6024552 6024552_6024552 0       0       0       -9      0
6025294_6025294 6025294_6025294 0       0       0       -9      0
6025455_6025455 6025455_6025455 0       0       0       -9      0
[  mturchin@node621  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.Phenos.Transformed.txt | tail -n 10
6008418 6008418 0.982149466245254 1.8348408914847 1.90215382198857 1.7746508525584
6010753 6010753 -2.14872618535448 -0.112070836661497 -0.101954153468391 -0.657232980242345
6011864 6011864 -1.05915623212098 1.73280000459313 1.4416885105281 1.55669008893538
6014349 6014349 -1.6506530741247 1.01323377099938 -0.382815612615648 0.436144212003557
6018726 6018726 1.82079001431327 1.21554668507235 1.38362419075144 1.62866290673721
6021711 6021711 0.591590115197903 -0.313927165407927 0.116665723858939 -0.393853203047272
6023159 6023159 0.180081946610816 0.67332717198737 0.222343460859681 1.30968233210111
6024552 6024552 0.761755798679409 1.00128624289528 1.29887607899065 1.42313443202277
6025294 6025294 -0.37373188626138 1.91986810870986 2.03468262436548 2.09352649294885
6025455 6025455 1.11598780644384 -0.60759044566529 -0.3396295383116 -0.404174127710831
[  mturchin@node621  ~/LabMisc/RamachandranLab/InterPath]$        paste <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.gz | perl -lane 'print join("\t", @F[0..6]);' | sed 's/_/ /g' | awk '{ print $1 "\t" $2 }') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.Phenos.Transformed.txt | awk '{ print $1 }') | head -n 10
FID     IID     FID
1001592 1001592 1001592
1002560 1002560 1002560
1003036 1003036 1003036
1004593 1004593 1004593
1008167 1008167 1008167
1009965 1009965 1009965
1010953 1010953 1010953
1012491 1012491 1012491
1013297 1013297 1013297
[  mturchin@node621  ~/LabMisc/RamachandranLab/InterPath]$        paste <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.gz | perl -lane 'print join("\t", @F[0..6]);' | sed 's/_/ /g' | awk '{ print $1 "\t" $2 }') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.Phenos.Transformed.txt | awk '{ print $1 }') | awk '{ if ($1 != $3) { print $0 } }' | wc
      0       0       0
> Data3[1:5,1:5]
     1:729632_T 1:752721_G 1:754105_T 1:756604_G 1:759036_A
[1,] -0.2454005  0.7235517 -0.2416615  0.2437967 -0.2393963
[2,] -0.2454005 -0.8359579 -0.2416615 -1.1756625 -0.2393963
[3,] -0.2454005  0.7235517 -0.2416615  1.6632559 -0.2393963
[4,] -0.2454005 -0.8359579 -0.2416615 -1.1756625 -0.2393963
[5,] -0.2454005 -0.8359579 -0.2416615 -1.1756625 -0.2393963
> Data3.cov[1:5,1:5]
           [,1]       [,2]        [,3]        [,4]       [,5]
[1,] 69.0008039  0.2584061 -0.11086390 -0.25001738  0.7030988
[2,]  0.2584061 75.3408524 -0.41715998 -0.21214926  0.5146039
[3,] -0.1108639 -0.4171600 67.43890405  0.01682419 -0.1411564
[4,] -0.2500174 -0.2121493  0.01682419 70.24424557 -0.1175144
[5,]  0.7030988  0.5146039 -0.14115639 -0.11751442 69.1273797
> table(Data3[,1])

-0.245400499390126   4.07364828987609
              2905                175
> table(Data3[,2])

-0.835957895660271  0.723551685583602   2.28306126682747
              1679               1151                250
> table(Data3[,3])

-0.241661491966026   4.13667612718316
              2910                170
> table(Data3[,4])

-1.17566247031684 0.243796725518467  1.66325592135378
             1074              1461               545
> table(Data3[,5])

-0.239396325134154   4.17581733602269
              2913                167
[  mturchin@node621  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | wc
 408348 4108101 29116909
[  mturchin@node621  ~/LabMisc/RamachandranLab/InterPath]$join <(cat /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | awk '{ print $1 "_" $2 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 }' | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.exonic.AAFix.pre.txt | awk '{ print $1 "_" $2 }' | sort) | wc                 24621  147726  897900
[  mturchin@node621  ~/LabMisc/RamachandranLab/InterPath]$join -v 1 <(cat /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | awk '{ print $1 "_" $2 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 }' | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.exonic.AAFix.pre.txt | awk '{ print $1 "_" $2 }' | sort) | wc
 383727 2302362 13998273
[  mturchin@node621  ~/LabMisc/RamachandranLab/InterPath]$R -q -e "408348 - 383727"
> 408348 - 383727
[1] 24621
> 
> 
[  mturchin@node621  ~/LabMisc/RamachandranLab/InterPath]$cat <(join -v 1 <(cat /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | awk '{ print $1 "_" $2 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 }' | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.exonic.AAFix.pre.txt | awk '{ print $1 "_" $2 }' | sort) | perl -lane 'print join("\t", @F[1..$#F]);') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.hg19_multianno.exonic.AAFix.pre.txt) | wc
 408348 2041740 10143725
[  mturchin@smp013  ~/LabMisc/RamachandranLab/InterPath/Vs1]$for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep Indian`; do
>         ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
>         ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
>
>         echo $pheno1 $ancestry1 $ancestry2 $ancestry3
>
>         R -q -e "library("data.table"); library("feather"); \
>         ptm <- proc.time(); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); print(proc.time() - ptm); \
>         Data3.mean <- apply(Data3, 2, mean); Data3.sd <- apply(Data3, 2, sd); Data3.sd[which(Data3.sd==0)] <- 1; Data3 <- t((t(Data3)-Data3.mean)/Data3.sd); \
>         ptm <- proc.time(); Data3.cov <- 1/nrow(Data3) * tcrossprod(as.matrix(Data3)); print(proc.time() - ptm); \
>         write.table(Data3.cov, \"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.cov.txt\", quote=FALSE, col.name=FALSE, row.name=FALSE);"
> done
Indian Indian
...
     user    system   elapsed
13206.894     0.022 13253.066
>
>
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -E 'African|British|Indian'`; do
>         ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
>         ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
>         echo $pheno1 $ancestry1 $ancestry2 $ancestry3
>
>         cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.British.British.Phenos.Trans.ADD.assoc.linear.clumped.wRowPos.txt | perl -lane 'my @vals1 = split(/,/, $F[$#F]); print $F[0], "\t", scalar(@vals1);'
>
> done
African African
Height  2184
BMI     472
Waist   317
Hip     473
British British.Ran4000
Height  1991
BMI     458
Waist   304
Hip     446
Indian Indian
Height  965
BMI     232
Waist   174
Hip     231
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -E 'African|British|Indian'`; do
>         ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`                                                                                                                                     >         ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
>         echo $pheno1 $ancestry1 $ancestry2 $ancestry3                                                                                                                                                                         > 
>         cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.British.British.Phenos.Trans.ADD.assoc.linear.clumped.wRowPos.txt | perl -lane 'my @vals1 = split(/,/, $F[$#F]); print $F[0], "\t", scalar(@vals1);'
> 
> done
African African
Height  965                                                                                                                                                                                                                     BMI     232
Waist   174
Hip     231
British British.Ran4000
Height  2184
BMI     472
Waist   317
Hip     473
Indian Indian
Height  1991
BMI     458
Waist   304
Hip     446
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$join <(cat /users/mturchin/data/ukbiobank_jun17/mturchin/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.txt | awk '{ print $1 "_" $2 "\t" $0 }' | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.txt | awk '{ print $1 "_" $2 "\t" $0 }' | sort -k 1,1) | head -n 10
1001592_1001592 1001592 1001592 0.118664115101393 1.17841032424695 1.24353357925829 1.89450422410686 1001592 1001592 1 African 45 0.04010937 -0.06472882 -0.01579029 -0.0028063 0.00672054 -0.01071403 -0.004161221 0.02206256 0.04384965 -0.04821377
1002560_1002560 1002560 1002560 0.734183319747881 2.04290941388721 1.4416885105281 1.8805781156021 1002560 1002560 0 African 49 -0.05246519 -0.09301313 0.03203676 0.01590448 -0.06281703 -0.01510241 0.05077683 0.004314775 -0.0129375 -0.0155569
1003036_1003036 1003036 1003036 0.217969859995042 1.82585180750194 1.85728627141396 1.94552746926506 1003036 1003036 0 African 41 0.05604095 0.01899383 0.004137209 -0.05324963 -0.005971676 -0.003017448 -0.004951464 -0.006294852 -0.01491296 0.004271876
1004593_1004593 1004593 1004593 -0.0567068367051468 1.17620744874895 1.13275845670766 1.37836394899338 1004593 1004593 0 African 50 0.04044754 0.04018205 0.004412988 0.051949 0.002875837 0.01332061 0.01493853 0.002691042 0.001538728 -0.008088887
1008167_1008167 1008167 1008167 0.322764805585384 0.15175399426104 0.0479378611405927 0.0196024993412122 1008167 1008167 0 African 59 0.03080848 -0.09631819 -0.02496982 0.01417764 0.02021924 -0.005031504 0.003485892 0.03258397 0.02134949 0.006985639
1009965_1009965 1009965 1009965 -1.17248489348928 1.28457721261652 1.28787234250449 1.42107036024099 1009965 1009965 0 African 58 0.05294918 0.05247505 0.0107574 0.03858143 0.004743263 0.01621968 0.02457712 -0.005881776 -0.01046 0.0151324
1010953_1010953 1010953 1010953 -2.3492759185006 0.600994502398072 0.498626444501437 -0.166021390280289 1010953 1010953 1 African 48 0.04929234 0.02053753 0.01234649 -0.04077462 -0.01188147 -0.0006558676 -0.009392189 -0.004736571 -0.003848354 0.01392407
1012491_1012491 1012491 1012491 1.0238031036417 -0.370391689731401 0.587105771979022 0.304346330032523 1012491 1012491 0 African 41 -0.415706 0.009725911 0.1694342 -0.02152109 0.124938 -0.06456354 0.04600945 0.01419852 0.009809398 -0.03911038
1013297_1013297 1013297 1013297 -0.157620896463327 0.509608883879002 0.526718205408012 0.130864454904108 1013297 1013297 1 African 46 0.03162186 0.04584261 -0.0113241 0.02832364 0.01017347 0.01959432 0.0002467427 0.006316162 0.01381216 -0.0051383
1013636_1013636 1013636 1013636 1.18596748140771 0.182411110585496 0.353722355972194 0.654152844441752 1013636 1013636 0 African 48 0.05394196 0.02237311 0.006092331 -0.04029114 -0.006966595 0.001310247 -0.0007259383 0.0001437942 -0.01494688 0.01508398
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$join <(cat /users/mturchin/data/ukbiobank_jun17/mturchin/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.txt | awk '{ print $1 "_" $2 "\t" $0 }' | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.txt | awk '{ print $1 "_" $2 "\t" $0 }' | sort -k 1,1) | wc
   3083   67826  769092
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/mturchin/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.txt | wc
   501731 3010386 44341746
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.txt | wc
   3083   46245  451340
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$join <(cat /users/mturchin/data/ukbiobank_jun17/mturchin/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.txt | awk '{ print $1 "_" $2 "\t" $0 }' | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.txt | awk '{ print $1 "_" $2 "\t" $0 }' | sort -k 1,1) | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23 }' | grep FID
FID     IID     Height  BMI     Waist   Hip     SEX     ANCESTRY        AGE     PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9     PC10
FID     IID     Height  BMI     Waist   Hip     SEX     ANCESTRY        AGE     PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9     PC10



~~~














