#!/bin/sh

###20180605 -- InterPath

#20180605
#TOC

 - 20180605: Vs1

#20180813
#See https://stackoverflow.com/questions/4089430/how-can-i-determine-the-url-that-a-local-git-repository-was-originally-cloned-fr, https://stackoverflow.com/questions/42830557/git-remote-add-origin-vs-remote-set-url-origin/42830632, https://conda.io/docs/user-guide/tasks/manage-environments.html#cloning-an-environment
#NOTE -- Made switch from `InterPath` to `InterPath`, via ':.,$ s/InterPath/InterPath/g', and manually follow-up with directories and such; made changes to git repo as well

#20181106
#Re-built website structure/platform using `workflowr` (see: https://jdblischak.github.io/workflowr/index.html)
R -q -e "library(\"workflowr\"); wflow_start(\"/users/mturchin/LabMisc/RamachandranLab/InterPath\", existing=TRUE, git=FALSE); wflow_build();"

#20180605
#Conda environment setup/details/etc
#See `/users/mturchin/RamachandranLab.CCV_General.Workbook.vs1.sh` for initial setup 

conda create -n InterPath
source activate InterPath
#20181003 NOTE -- with the system change from CentOS to RedHat, thinking I need to redo everything here; was getting some failures in the InterPath.cpp code and it may just be some fundamental issues. So just redoing everything here with the new name of InterPath2
##conda create -n InterPath2
##source activate InterPath2
#2018181003 NOTE -- actually realized I should probably just reinstall miniconda from scratch so, so made copy '/users/mturchin/data/mturchin/miniconda2RH', changed $PATH (manually on one screen window first just to experiment and make sure things work), and then created a new 'InterPath' environment under 'minicoda2RH'
conda create -n InterPath
source activate InterPath
#conda update conda
conda install armadillo eigen boost gcc
conda install R perl java-jdk
conda install git
conda install plink bedtools vcftools vcftools bwa samtools picard gatk imagemagick gnuplot eigensoft tabix
conda install r-base r-devtools r-knitr r-testthat r-cairo r-ashr r-rcolorbrewer r-essentials r-extrafont fonts-anaconda
#conda install eigen boost gcc
##20180311 (From MultiEthnicGWAS)
##conda install flashpca -- failed
#conda install eigen
##conda install spectra -- installed a Python package called spectra, not what was intended
#conda install boost
##conda install libgomp -- failed
#conda install gcc
conda install r-doParallel r-Rcpp r-RcppArmadillo r-RcppParallel r-CompQuadForm r-Matrix r-MASS r-truncnorm
#20180611 #20180814 NOTE -- not sure if I actually installed fortran this way, since I tried to reinstall using these same commands and conda said fortran was not available? Briefly looking online suggests that installing 'gcc' also installs 'gfortran', which is what you want?
#conda install armadillo fortran
conda install r-data.table r-bigmemory
conda install julia
##conda install r-coop -- failed, install via 'install.packages("coop")'
##conda install r-rbenchmark -- failed, install via 'install.packages("rbenchmark")'
##conda install r-cpgen -- failed, install via 'install.packages("cpgen")'
##conda install docker -- failed
#20180815 NOTE -- re-install for movement into 'InterPath' from 'FORCE'
#From https://stackoverflow.com/questions/42231764/how-can-i-rename-a-conda-environment (but 'cloning' didn't seem to really work in the end?....see comments from 'MultiEthnicGWAS' too)
###conda install R perl java-jdk
###conda install git plink bedtools vcftools vcftools bwa samtools picard gatk imagemagick gnuplot eigensoft tabix
###conda install r-base r-devtools 
###conda install r-knitr r-testthat r-cairo r-ashr 
###conda install r-rcolorbrewer r-essentials r-extrafont fonts-anaconda 
###conda install eigen boost gcc 
###conda install r-doParallel r-Rcpp r-RcppArmadillo r-RcppParallel r-CompQuadForm r-Matrix r-MASS r-truncnorm 
###conda install armadillo r-data.table r-bigmemory julia
conda install -c bioconda p7zip #20181108


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
#g++ Test1.cpp -o Test1.out  -O2 /gpfs/data/sramacha/mturchin/miniconda2/envs/FORCE/lib/libgfortran.so.4 /gpfs/data/sramacha/mturchin/miniconda2/envs/FORCE/lib/libarmadillo.so.8 -larmadillo -lopenblas -larpack -lgfortran -L/gpfs/data/sramacha/mturchin/miniconda2/envs/FORCE/lib
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
##UKBioBankPops=`echo "African;African;Afr British;British.Ran4000;Brit4k Caribbean;Caribbean;Carib Chinese;Chinese;Chi Indian;Indian;Indn Pakistani;Pakistani;Pkstn"`;
UKBioBankPops=`echo "African;African;Afr British;British.Ran4000;Brit4k British;British.Ran10000;Brit10k Caribbean;Caribbean;Carib Chinese;Chinese;Chi Indian;Indian;Indn Irish;Irish;Irish Pakistani;Pakistani;Pkstn"`;
##UKBioBankPops=`echo "African;African;Afr British;British.Ran4000;Brit4k British;British.Ran10000;Brit10k British;British.Ran100000;Brit100k Caribbean;Caribbean;Carib Chinese;Chinese;Chi Indian;Indian;Indn Irish;Irish;Irish Pakistani;Pakistani;Pkstn"`;
UKBioBankPops=`echo "African;African;Afr;472840 British;British.Ran4000;Brit4k;138503 British;British.Ran10000;Brit10k;9827442 Caribbean;Caribbean;Carib;328593 Chinese;Chinese;Chi;842743 Indian;Indian;Indn;549281 Irish;Irish;Irish;902143 Pakistani;Pakistani;Pkstn;232849"`;

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
sourceCpp("/users/mturchin/LabMisc/RamachandranLab/InterPath/InterPath.cpp"); Data3 <- X; set.seed(123); for (i in c(1e4, 1e5, 2e5, 3e5, ncol(Data3))) {
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

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep British | grep -v Ran4000 | grep -v Ran100000`; do
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
#mkdir /Volumes/NO\ NAME/British.Ran10000 /Volumes/NO\ NAME/British.Ran100000
#mkdir /Users/mturchin20/Documents/Work/LabMisc/Data/UKBioBank/subsets
#mkdir /Users/mturchin20/Documents/Work/LabMisc/Data/UKBioBank/subsets/British
#mkdir /Users/mturchin20/Documents/Work/LabMisc/Data/UKBioBank/subsets/British/British.Ran100000
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chr*_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/African 
scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/British/British/mturchin20/ukb_chr*_v2.British.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/British 
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran4000/mturchin20/ukb_chr*_v2.British.Ran4000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/British.Ran4000
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/mturchin20/ukb_chr*_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/Caribbean 
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/Chinese/Chinese/mturchin20/ukb_chr*_v2.Chinese.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/Chinese 
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/Indian/Indian/mturchin20/ukb_chr*_v2.Indian.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/Indian 
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/Irish/Irish/mturchin20/ukb_chr*_v2.Irish.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/Irish 
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/Pakistani/Pakistani/mturchin20/ukb_chr*_v2.Pakistani.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/Pakistani 
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/mturchin20/ukb_chr*_v2.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/British.Ran10000
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran100000/mturchin20/ukb_chr*_v2.British.Ran100000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Volumes/NO\ NAME/British.Ran100000
scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/mturchin20/ukb_chr*_v2.British.Ran10000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Users/mturchin20/Documents/Work/LabMisc/Data/UKBioBank/subsets/British/British.Ran10000/.
scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran100000/mturchin20/ukb_chr*_v2.British.Ran100000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.vcf.gz /Users/mturchin20/Documents/Work/LabMisc/Data/UKBioBank/subsets/British/British.Ran100000/.
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
cd /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20
wget -c https://imputationserver.sph.umich.edu/share/results/d82d556ca41e50a3758e51b6a7f866fd/qcreport.html https://imputationserver.sph.umich.edu/share/results/7ef9040846b1650a5597bb8efd84c024/chr_1.zip https://imputationserver.sph.umich.edu/share/results/bc75d166f1d5a482e295ed724c5bd26/chr_2.zip https://imputationserver.sph.umich.edu/share/results/350ab6b9f5f37a224074cc2a47e52dd5/chr_3.zip https://imputationserver.sph.umich.edu/share/results/39c306f89a177995422ecf7e77706e0c/chr_4.zip https://imputationserver.sph.umich.edu/share/results/d211ae5f352c88b95657191860ac9793/chr_5.zip https://imputationserver.sph.umich.edu/share/results/dd46ec80babb33c8d595a14a3bff524/chr_6.zip https://imputationserver.sph.umich.edu/share/results/fb9400ac2ee3688f54cd05fd88252b63/chr_7.zip https://imputationserver.sph.umich.edu/share/results/cd71894604446cb825bec1b0a9d1b625/chr_8.zip https://imputationserver.sph.umich.edu/share/results/a5da9db134f8db927fd26940e4adc903/statistics.txt https://imputationserver.sph.umich.edu/share/results/b5c899698ab2815f34b528ce89e3ba5c/chr_1.log https://imputationserver.sph.umich.edu/share/results/9d7b00e90a42b61502744b74affa43a4/chr_2.log https://imputationserver.sph.umich.edu/share/results/6b43388ddb97fd63f577a31b0560d47d/chr_3.log https://imputationserver.sph.umich.edu/share/results/bea1deac3e07229a118f2a2de5240cad/chr_4.log https://imputationserver.sph.umich.edu/share/results/d8ee2cd98ba1adec8764c2e1b42f8f9a/chr_5.log https://imputationserver.sph.umich.edu/share/results/4340a6b78d36edf6718602b8e59b4802/chr_6.log https://imputationserver.sph.umich.edu/share/results/d34e6eea778c274c8c899a87ee49df07/chr_7.log https://imputationserver.sph.umich.edu/share/results/6c31b717fa61388f3e8385575abe2269/chr_8.log
mv /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/qcreport.html /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/qcreport.1.html
mv /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/statistics.txt /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/statistics.1.txt
wget -c https://imputationserver.sph.umich.edu/share/results/13a898b3d02840022fbcc8875d64d2a5/qcreport.html https://imputationserver.sph.umich.edu/share/results/d14a7a1fef66ecd3602653ba90f9b329/chr_10.zip https://imputationserver.sph.umich.edu/share/results/99b852e5f5f1a5345d85f571c1ea81bc/chr_11.zip https://imputationserver.sph.umich.edu/share/results/a870e6d14085f675c81292bd5596dac0/chr_12.zip https://imputationserver.sph.umich.edu/share/results/85db580b996d62ec23c82930be426426/chr_13.zip https://imputationserver.sph.umich.edu/share/results/b2f58592769463fb757f9302aca126f2/chr_14.zip https://imputationserver.sph.umich.edu/share/results/e6d00bd40d02e79bfd63edba8269c4e8/chr_15.zip https://imputationserver.sph.umich.edu/share/results/caa037f3aa537f2dfb8f3c82d2f9b275/chr_16.zip https://imputationserver.sph.umich.edu/share/results/bbd1f5505f7384ca0b11561268cb75b3/chr_17.zip https://imputationserver.sph.umich.edu/share/results/c7b7a1f2541c72668616d8cc47f00dd6/chr_18.zip https://imputationserver.sph.umich.edu/share/results/5ed64b9b215fe261920d938c36a69c0e/chr_19.zip https://imputationserver.sph.umich.edu/share/results/1311f3619b6cd04fb1069734ef0197d2/chr_20.zip https://imputationserver.sph.umich.edu/share/results/9b7ac01118bdca168b98fe07bae9cf32/chr_21.zip https://imputationserver.sph.umich.edu/share/results/5827b86ca0c3d686642128a9ad63101d/chr_22.zip https://imputationserver.sph.umich.edu/share/results/a26979d43eb48e77497b7c2b67c4dc0d/chr_9.zip https://imputationserver.sph.umich.edu/share/results/65ed8223e2987d1f355596b800fcad58/statistics.txt https://imputationserver.sph.umich.edu/share/results/a879c466d9d24e4dc51d1fd6b3612a54/chr_10.log https://imputationserver.sph.umich.edu/share/results/e58b98643d14ccae90e6a9c8312d3b50/chr_11.log https://imputationserver.sph.umich.edu/share/results/b29fe8d4ca3036ec6bbd0fd5e3b4e378/chr_12.log https://imputationserver.sph.umich.edu/share/results/17808d27eff84db2a44d3d49b5307094/chr_13.log https://imputationserver.sph.umich.edu/share/results/5f9d85aee8b7469efc5dcd9b33c1ae6c/chr_14.log https://imputationserver.sph.umich.edu/share/results/2c6b4b0a1da69a5c851957aa0061e543/chr_15.log https://imputationserver.sph.umich.edu/share/results/c095af3835f1393b0f9257f7565f4736/chr_16.log https://imputationserver.sph.umich.edu/share/results/a1ae308042c45c0c16a6a511dc3912dd/chr_17.log https://imputationserver.sph.umich.edu/share/results/6a347d7c696b2d50b010031e2d7a5042/chr_18.log https://imputationserver.sph.umich.edu/share/results/b775216a1c24fbe3a0ebff48ec803e94/chr_19.log https://imputationserver.sph.umich.edu/share/results/6dda06399ada10ee7b5c25d93df37a10/chr_20.log https://imputationserver.sph.umich.edu/share/results/5e46ef2476dad9c2c9c64006079c0752/chr_21.log https://imputationserver.sph.umich.edu/share/results/3cd3f7a4586bf4c90c107044e41172c3/chr_22.log https://imputationserver.sph.umich.edu/share/results/ff9e7dfbf843bb3609201bfaf7603344/chr_9.log
mv /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/qcreport.html /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/qcreport.2.html
mv /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/statistics.txt /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/statistics.2.txt
for i in {1..8}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/chr_${i}.zip -p'T9JZ9bUhYaglF>'
done
for i in {9..22}; do
	7za x /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran10000/Imputation/mturchin20/chr_${i}.zip -p'wKq9ViltKGPt0'
done
#20181104 NOTE -- Michigan imputation server can only do up to 20k samples, so cannot for the moment do 100k Brit imputation through that route
#cd /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran100000/Imputation/mturchin20
#wget 
#for i in {1..22}; do
#	7za x /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran100000/Imputation/mturchin20/chr_${i}.zip -p'
#done
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
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v African | grep Ran10000`; do
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

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v African | grep Ran10000`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
        ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
		
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3

	join <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr*_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.info.gz | awk '{ if ($7 > .3) { print $1 } }' | sed 's/:/_/g' | grep -v SNP | sort | uniq -u) <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr*_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.info.gz | awk '{ print $1 }' | sed 's/:/_/g' | grep -v SNP | sort | uniq -u) | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.info.r2gt3.noDups.ChrBPs.gz
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chr*_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim | awk '{ print $1 "_" $4 }' | sort | uniq > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim.ChrBPs
	join <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.info.r2gt3.noDups.ChrBPs.gz | sort) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim.ChrBPs | sort) | sed 's/_/:/g' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim.ImptHRC.info.r2gt3.noDups.ChrBPs

done

#From https://www.biostars.org/p/46060/ & https://sourceforge.net/p/vcftools/mailman/message/29115811/; 20181108 From https://unix.stackexchange.com/questions/223778/how-to-run-an-infinite-loop-in-the-background
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep Ran10000`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
        ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`
		
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3

	for i in {7..10}; do
		echo $i
		sbatch -t 72:00:00 --mem 2g --account=ccmb-condo -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.slurm.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.error --comment "$ancestry1 $ancestry2 $i" <(echo -e '#!/bin/sh';
		echo -e "\nvcftools --gzvcf /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.vcf.gz --plink-tped --snps /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.bim.ImptHRC.info.r2gt3.noDups.ChrBPs --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp"; \
		echo -e "\nplink --tfile /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp --geno 0 --make-bed --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.plinkTemp";)
#		echo -e "\nrm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.plinkTemp.tped /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.plinkTemp.tfam";)	
#		rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp* 

	done 
	sleep 1
done 
#		echo -e "\nplink --tfile /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp --exclude /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.missnp --make-bed --out /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp"; \
#		echo -e "\nrm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.bed /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.bim /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.fam /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.tped /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chr${i}_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.plinkTemp.tfam"; \	

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v African | grep Ran10000`; do
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

cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chr1_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.plinkTemp.fam | awk '{ print $1 }' | sed 's/_/ /g' > /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.plinkTemp.FIDIIDs

#From https://stats.stackexchange.com/questions/11000/how-does-r-handle-missing-values-in-lm (for lm(..., na.action=na.exclude))
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);')`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`

        echo $ancestry1 $ancestry2 /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.txt 

	join <(cat /users/mturchin/data/ukbiobank_jun17/mturchin/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.txt | awk '{ print $1 "_" $2 "\t" $0 }' | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.txt | awk '{ print $1 "_" $2 "\t" $0 }' | sort -k 1,1) | awk '{ print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23 }' | grep -v FID | R -q -e "Data1 <- read.table(file('stdin'), header=F); for (i in 3:6) { Data1 <- cbind(Data1, residuals(lm(Data1[,i] ~ Data1[,10] + Data1[,11] + Data1[,12] + Data1[,13] + Data1[,14] + Data1[,15] + Data1[,16] + Data1[,17] + Data1[,18] + Data1[,19], na.action=na.exclude))); }; write.table(Data1, file=\"\", quote=FALSE, row.name=FALSE, col.name=FALSE);" | grep -v \> | grep -v Height | cat <(echo "FID IID Height BMI Waist Hip FID IID SEX ANCESTRY AGE PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 Height BMI Waist Hip") -  | perl -lane 'print $F[0], "\t", $F[1], "\t", join("\t", @F[$#F-3..$#F]);' > /users/mturchin/data/ukbiobank_jun17/mturchin/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.top10resids.txt 

done

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v African`; do
        ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
        ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`

        echo $ancestry1 $ancestry2 /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.txt 

	join <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.txt | awk '{ print $1 "_" $2 "\t" $0 }' | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt | awk '{ print $1 "_" $2 }' | sort) | grep -v ANCESTRY | cat <(echo "FID_IID FID IID SEX ANCESTRY AGE PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10") - | perl -lane 'print join("\t", @F[1..$#F]);' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.sort.ImptHRC.dose.100geno.raw.txt 

done

#From https://stackoverflow.com/questions/5774813/short-formula-call-for-many-variables-when-building-a-model
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v 'African' | grep Irish`; do
	for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep Plus`; do
		ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
		NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt | wc | awk '{ print $1 }'`
#		NumPaths=80
		echo $ancestry1 $ancestry2 $ancestry3 $k
                         	
		if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/phenos ]; then
			mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/phenos
		fi

		for (( PathNum=1; PathNum <= $NumPaths; PathNum=PathNum+10 )); do
			sbatch -t 72:00:00 --mem 31g --account=ccmb-condo -o temp1.output -e temp1.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
			echo -e "\nR -q -e \"library(\\\"data.table\\\"); \
			Data1 <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt\\\", header=T); Data3 <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\\\", header=F); neg.is.na <- Negate(is.na); \
			for (i in $PathNum:($PathNum+9)) { \
				if (i <= $NumPaths) { \ 
					Data2 <- fread(paste(\\\"zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Regions.c2.${k}.Pathways\\\", as.character(i), \\\".txt.gz\\\", sep=\\\"\\\"), header=T); \ 
					PhenoNew <- Data1[,c(1:2)]; \
					Data2.mean <- apply(Data2, 2, mean); Data2.sd <- apply(Data2, 2, sd); Data2.sd[which(Data2.sd==0)] <- 1; Data2 <- t((t(Data2)-Data2.mean)/Data2.sd); \
					for (j in 3:6) { \
						PhenoNew <- cbind(PhenoNew, residuals(lm(Data1[,j] ~ Data2, na.action=na.exclude))); \  
					}; \
					colnames(PhenoNew) <- colnames(Data1); \
					write.table(PhenoNew, file=paste(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.Regions.c2.${k}.Pathways\\\", i, \\\".txt\\\", sep=\\\"\\\"), quote=FALSE, row.name=FALSE, col.name=TRUE); \
				}; \
			};\""); 
		done;		
	done;
done;		
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v 'African' | grep Irish`; do
	for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep Plus`; do
		ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
		echo $ancestry1 $ancestry2 $ancestry3 $k

		gzip -f /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.*Regions.c2.${k}.Pathway*.txt         	
	done;
done;		

#Orig
#		R -q -e "library(\"data.table\"); Data1 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt\", header=T); Data2 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); Data3 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\", header=F); neg.is.na <- Negate(is.na); \
#		write.table(PhenoNew, file=paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.Regions.c2.${k}.Pathways\", i, \".txt\", sep=\"\"), quote=FALSE, row.name=FALSE, col.name=TRUE); \
#PCadj
#		R -q -e "library(\"data.table\"); Data1 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.flashpca.top10resids.txt\", header=T); Data2 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); Data3 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\", header=F); neg.is.na <- Negate(is.na); \			
#			write.table(PhenoNew, file=paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.top10resids.Regions.c2.${k}.Pathways\", i, \".txt\", sep=\"\"), quote=FALSE, row.name=FALSE, col.name=TRUE); \

#				print(i); \
#				Pathways.Regions <- as.numeric(as.character(unlist(strsplit(as.character(Data3[i,3]), \",\")))); \
#				Data2.sub <- data.frame(as.matrix(Data2)[,c(Pathways.Regions)]); \
#				Data2.sub.mean <- apply(Data2.sub, 2, mean); Data2.sub.sd <- apply(Data2.sub, 2, sd); Data2.sub.sd[which(Data2.sub.sd==0)] <- 1; Data2.sub <- t((t(Data2.sub)-Data2.sub.mean)/Data2.sub.sd); \
#					PhenoNew <- cbind(PhenoNew, residuals(lm(Data1[,j] ~ Data2.sub, na.action=na.exclude))); \  

#xnam <- paste("x", 1:25, sep="")
#(fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"))))
#lm(fmla, data = myData)
#		R -q -e "library(\"data.table\"); Data1 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/mturchin/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.African.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.top10resids.txt\", header=T); Data2 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); Data3 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.NonSyn.txt\", header=F); \
#			print(length(Pathways.Regions));
#				PhenoNew <- cbind(PhenoNew, residuals(lm(Data1[,j] ~ Data2.sub, na.action=na.exclude))); \  
#				PhenoNew <- cbind(PhenoNew, residuals(lm(as.formula(paste(\"Data1[,\", j, \"] ~ \", paste(paste(\"Data2.sub[,\", 1:ncol(Data2.sub), \"]\", sep=\"\"), collapse=\" + \"), sep=\"\")), na.action=na.exclude))); \
##				PhenoNew <- cbind(PhenoNew, residuals(lm(Y.Pheno.noNAs ~ X.Pheno.noNAs, na.action=na.exclude))); \  
##				PhenoNew <- cbind(PhenoNew, residuals(lm(as.formula(paste(\"Y.Pheno.noNAs ~ \", paste(paste(\"X.Pheno.noNAs[,\", 1:ncol(X.Pheno.noNAs), \"]\", sep=\"\"), collapse=\" + \"), sep=\"\")), na.action=na.exclude))); \

















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
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
	ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`; ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`; echo $pheno1 $ancestry1 $ancestry2 $ancestry3;

#	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep nonsynonymous | sort -k 1,1 | perl -lane 'if ($. == 1) { @gene1; push(@gene1, $F[0]); push(@gene1, $F[3]); } else { if ($F[0] ne $gene1[0]) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); @gene1 = (); push(@gene1, $F[0]); push(@gene1, $F[3]); } else { push(@gene1, $F[3]); } if (eof()) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); } };' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.NonSyn.txt
#	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep exonic | sort -k 1,1 | perl -lane 'if ($. == 1) { @gene1; push(@gene1, $F[0]); push(@gene1, $F[3]); } else { if ($F[0] ne $gene1[0]) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); @gene1 = (); push(@gene1, $F[0]); push(@gene1, $F[3]); } else { push(@gene1, $F[3]); } if (eof()) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); } };' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.Exonic.txt
#	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep -E 'exonic|intronic|UTR' | sort -k 1,1 | perl -lane 'if ($. == 1) { @gene1; push(@gene1, $F[0]); push(@gene1, $F[3]); } else { if ($F[0] ne $gene1[0]) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); @gene1 = (); push(@gene1, $F[0]); push(@gene1, $F[3]); } else { push(@gene1, $F[3]); } if (eof()) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); } };' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus.txt
#	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.txt | grep -E 'exonic|intronic|UTR|upstream|downstream' | sort -k 1,1 | perl -lane 'my @info1 = split(/,/, $F[2]); if ($info1[0] =~ m/intergenic/) { my @dists1 = split(/=/, $info1[1]); if ($dists1[1] <= 20000) { print join("\t", @F); } } else { print join("\t", @F); }' | perl -lane 'if ($. == 1) { @gene1; push(@gene1, $F[0]); push(@gene1, $F[3]); } else { if ($F[0] ne $gene1[0]) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); @gene1 = (); push(@gene1, $F[0]); push(@gene1, $F[3]); } else { push(@gene1, $F[3]); } if (eof()) { print $gene1[0], "\t", join(",", @gene1[1..$#gene1]); } };' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus20kb.txt
	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.NonSyn.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(",", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 2) { print join("\t", @F); }' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.NonSyn.txt
	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.Exonic.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(",", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 2) { print join("\t", @F); }' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.Exonic.txt
	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(",", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 2) { print join("\t", @F); }' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.ExonicPlus.txt
	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus20kb.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(",", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 2) { print join("\t", @F); }' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.ExonicPlus20kb.txt; done

#	for k in `cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | perl -F\\\t -lane 'print $F[6];' | sort | uniq`; do:w
uuuu
#		cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.txt | perl -F\\t -slane 'if ($F[6] eq $gene1) { print $gene1, "\t", $F[0], ":", $F[1]; }' -- -gene1=$k
#	done > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.bim.AnnovarFormat.TableAnnovar.hg19_multianno.GeneSNPs.txt

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);')`; do
	ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`; ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`; 
	echo $pheno1 $ancestry1 $ancestry2 $ancestry3;

	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.NonSyn.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; my@info2; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); push(@info2, $entry1); } } print $F[0], "\t", $F[1], "\t", join(",", @info2), "\t", join(";", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 3) { print join("\t", @F); }' | perl -lane 'my @vals1 = split(/;/, $F[$#F]); my $count1 = 1; my @print1; foreach my $entry1 (@vals1) { my @vals2 = split(/,/, $entry1); my @print2; for (my $entry2 = 0; $entry2 <= $#vals2; $entry2++) { push(@print2, $count1); $count1++ } push(@print1, join(",", @print2)); } $F[$#F] = join(";", @print1); print join("\t", @F);' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.geneVSpath.wGenes.NonSyn.txt 
	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.Exonic.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; my@info2; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); push(@info2, $entry1); } } print $F[0], "\t", $F[1], "\t", join(",", @info2), "\t", join(";", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 3) { print join("\t", @F); }' | perl -lane 'my @vals1 = split(/;/, $F[$#F]); my $count1 = 1; my @print1; foreach my $entry1 (@vals1) { my @vals2 = split(/,/, $entry1); my @print2; for (my $entry2 = 0; $entry2 <= $#vals2; $entry2++) { push(@print2, $count1); $count1++ } push(@print1, join(",", @print2)); } $F[$#F] = join(";", @print1); print join("\t", @F);' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.geneVSpath.wGenes.Exonic.txt 
	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; my@info2; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); push(@info2, $entry1); } } print $F[0], "\t", $F[1], "\t", join(",", @info2), "\t", join(";", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 3) { print join("\t", @F); }' | perl -lane 'my @vals1 = split(/;/, $F[$#F]); my $count1 = 1; my @print1; foreach my $entry1 (@vals1) { my @vals2 = split(/,/, $entry1); my @print2; for (my $entry2 = 0; $entry2 <= $#vals2; $entry2++) { push(@print2, $count1); $count1++ } push(@print1, join(",", @print2)); } $F[$#F] = join(";", @print1); print join("\t", @F);' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.geneVSpath.wGenes.ExonicPlus.txt 
	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus20kb.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; my@info2; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); push(@info2, $entry1); } } print $F[0], "\t", $F[1], "\t", join(",", @info2), "\t", join(";", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 3) { print join("\t", @F); }' | perl -lane 'my @vals1 = split(/;/, $F[$#F]); my $count1 = 1; my @print1; foreach my $entry1 (@vals1) { my @vals2 = split(/,/, $entry1); my @print2; for (my $entry2 = 0; $entry2 <= $#vals2; $entry2++) { push(@print2, $count1); $count1++ } push(@print1, join(",", @print2)); } $F[$#F] = join(";", @print1); print join("\t", @F);' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.geneVSpath.wGenes.ExonicPlus20kb.txt 
done

#	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.NonSyn.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(";", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 2) { print join("\t", @F); }' | perl -lane 'my @vals1 = split(/;/, $F[$#F]); my $count1 = 1; my @print1; foreach my $entry1 (@vals1) { my @vals2 = split(/,/, $entry1); my @print2; for (my $entry2 = 0; $entry2 <= $#vals2; $entry2++) { push(@print2, $count1); $count1++ } push(@print1, join(",", @print2)); } $F[$#F] = join(";", @print1); print join("\t", @F);' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.geneVSpath.NonSyn.txt 
#	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.Exonic.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(";", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 2) { print join("\t", @F); }' | perl -lane 'my @vals1 = split(/;/, $F[$#F]); my $count1 = 1; my @print1; foreach my $entry1 (@vals1) { my @vals2 = split(/,/, $entry1); my @print2; for (my $entry2 = 0; $entry2 <= $#vals2; $entry2++) { push(@print2, $count1); $count1++ } push(@print1, join(",", @print2)); } $F[$#F] = join(";", @print1); print join("\t", @F);' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.geneVSpath.Exonic.txt 
#	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(";", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 2) { print join("\t", @F); }' | perl -lane 'my @vals1 = split(/;/, $F[$#F]); my $count1 = 1; my @print1; foreach my $entry1 (@vals1) { my @vals2 = split(/,/, $entry1); my @print2; for (my $entry2 = 0; $entry2 <= $#vals2; $entry2++) { push(@print2, $count1); $count1++ } push(@print1, join(",", @print2)); } $F[$#F] = join(";", @print1); print join("\t", @F);' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.geneVSpath.ExonicPlus.txt 
#	cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.ExonicPlus20kb.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(";", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 2) { print join("\t", @F); }' | perl -lane 'my @vals1 = split(/;/, $F[$#F]); my $count1 = 1; my @print1; foreach my $entry1 (@vals1) { my @vals2 = split(/,/, $entry1); my @print2; for (my $entry2 = 0; $entry2 <= $#vals2; $entry2++) { push(@print2, $count1); $count1++ } push(@print1, join(",", @print2)); } $F[$#F] = join(";", @print1); print join("\t", @F);' > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.geneVSpath.ExonicPlus20kb.txt 

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
	write.table(Data3.cov, \"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.noFix.cov.txt\", quote=FALSE, col.name=FALSE, row.name=FALSE);"
done 

#20181028 NOTE -- the 'noFix' thing I believe was just to account for this single SNP that past previous filters in the Chinese subset that was causing the entire covariance matrix just to be NA; I had previously looked into this and only found 1 SNP that appeared fix in any population (and then subsequently realized this single SNP may have been messing up the original covariance matrix; running the above code produced a correct covariance matrix for Chinese subset). So just soft linking the other files to the new filename schema even though I don't believe anything would actually change between the two versions of the files
#20181028 NOTE -- there is an extremely small difference in the values if I rerun the covarianc calculations (like at the <10e-8 order of magnitude), but since I'm practically rerunning the start of the Vs2 code from the getgo now, I'm just going to rerun everything with the 'noFix' in the output, and go from there (eg rerun everything with the 'new', barely any different 'noFix' versions) 
##ln -s /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.cov.txt /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.noFix.cov.txt
##ln -s /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran4000/mturchin20/Analyses/InterPath/ukb_chrAll_v2.British.Ran4000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.cov.txt /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran4000/mturchin20/Analyses/InterPath/ukb_chrAll_v2.British.Ran4000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.noFix.cov.txt
##ln -s /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/mturchin20/Analyses/InterPath/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.cov.txt /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/mturchin20/Analyses/InterPath/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.noFix.cov.txt
##ln -s /users/mturchin/data/ukbiobank_jun17/subsets/Indian/Indian/mturchin20/Analyses/InterPath/ukb_chrAll_v2.Indian.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.cov.txt /users/mturchin/data/ukbiobank_jun17/subsets/Indian/Indian/mturchin20/Analyses/InterPath/ukb_chrAll_v2.Indian.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.noFix.cov.txt
##ln -s /users/mturchin/data/ukbiobank_jun17/subsets/Pakistani/Pakistani/mturchin20/Analyses/InterPath/ukb_chrAll_v2.Pakistani.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.cov.txt /users/mturchin/data/ukbiobank_jun17/subsets/Pakistani/Pakistani/mturchin20/Analyses/InterPath/ukb_chrAll_v2.Pakistani.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.noFix.cov.txt

#Indn: 13253.066; 

#Data3.cov2 <- 1/nrow(Data3) * tcrossprod(scale(as.matrix(Data3), FALSE, FALSE))
#	ptm <- proc.time(); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); print(proc.time() - ptm); 
#	Data4 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.raw.gz | perl -lane \'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);\'', header=T);
#	write_feather(Data3, \"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw.edit.fthr\"); \
#	print(proc.time() - ptm); \
#	write_feather(as.data.frame(Data3.cov), \"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.cov.fthr\"); \
#	write.table(Data3.cov, "/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.cov.txt", quote=FALSE, col.name=FALSE, row.name=FALSE);"
#	print(table(apply(Data3, 2, sd) == 0));"

module load R/3.3.1; for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep Height`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v -E 'African' | grep Irish`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);')`; do
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt | wc | awk '{ print $1 }'`	
#			NumPaths=1442
			echo $i $ancestry1 $ancestry2 $ancestry3 $k
			
			if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways ]; then
				mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways
			fi
			if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k ]; then
				mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k
			fi
			
			for (( PathNum=1; PathNum <= $NumPaths; PathNum=PathNum+80 )); do
				sbatch -t 72:00:00 --mem 65g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${k}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${k}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
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
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep Irish`; do
	for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep 20kb`; do
		ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
		echo $i $ancestry1 $ancestry2 $ancestry3 $k
	
		gzip -f /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Regions.c2.${k}.Pathways*.txt
	done;
done;

#			R -q -e "library("data.table"); library("feather"); \ 
#			ptm <- proc.time(); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); print(proc.time() - ptm); \ 
#			Pathways <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\", header=F); \
#			for (i in 1:nrow(Pathways)) { \
#				print(i); \
#				Pathways.Regions <- as.numeric(as.character(unlist(strsplit(as.character(Pathways[i,3]), \",\")))); \
#				write.table(Data3.temp, paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Regions.c2.${k}.Pathways\", i, \".txt\", sep=\"\"), quote=FALSE, col.name=TRUE, row.name=FALSE); \

#20181001 NOTE -- Was going to do the below setup for the 'Vs2 gene*pathway' runs, but each pathway covariance file post-gz was ~50Mb (+/- 25Mb) in African NonSyn, so overall just seemed like potentially creating 1-2Tb worth of information was not worth it; it was onlty taking 1-2 mins to create matrix per pathway so far too, so possibly fine to have this get created within code anyways
#20181001 NOTE -- Changing mind againm, going to keep this setup; however since it doesn't take particularly long to create these files, I will intend on deleting them after initial runs & results to preserve space; if I need to rerun and redo things, I can recreate them and rerun this code to do so
module load R/3.3.2; for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v -E 'African|Irish' | grep -E 'Paki'`; do
	for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);')`; do
        	ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
		NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt | wc | awk '{ print $1 }'`	
#		NumPaths=2
		echo $i $ancestry1 $ancestry2 $ancestry3 $k
			
		if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/cov ]; then
			mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/cov
		fi
			
		for (( PathNum=1; PathNum <= $NumPaths; PathNum=PathNum+79 )); do
			sbatch -t 72:00:00 --mem 20g --account=ccmb-condo -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/cov/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${k}.cov.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/cov/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${k}.cov.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
			echo -e "\nR -q -e \"library(\\\"data.table\\\"); \ 
			Pathways <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\\\", header=F); \
			Pathways.Check <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/c2.all.v6.1.wcp_comps.symbols.${ancestry2}.Regions.c2.${k}.txt\\\", header=F); \
			for (i in $PathNum:($PathNum+79)) { \
				if ((i <= $NumPaths) && (Pathways.Check[i,ncol(Pathways.Check)])) { Data3 <- fread(paste(\\\"zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Regions.c2.${k}.Pathways\\\", as.character(i), \\\".txt.gz\\\", sep=\\\"\\\"), header=T); \
				Data3.mean <- apply(Data3, 2, mean); Data3.sd <- apply(Data3, 2, sd); Data3.sd[which(Data3.sd==0)] <- 1; Data3 <- t((t(Data3)-Data3.mean)/Data3.sd); \
				Data3.cov <- 1/nrow(Data3) * tcrossprod(as.matrix(Data3)); \
				write.table(Data3.cov, paste(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/cov/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Regions.c2.${k}.Pathways\\\", as.character(i), \\\".cov.txt\\\", sep=\\\"\\\"), quote=FALSE, col.name=FALSE, row.name=FALSE);}; \
			};\"") 
		done
	done
done 
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v Irish | head -n 7 | tail -n 1`; do
	for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);')`; do
		ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
		echo $i $ancestry1 $ancestry2 $ancestry3 $k
	
		gzip -f /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/cov/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Regions.c2.${k}.Pathways*.cov.txt
	done;
done;
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | head -n 1`; do
	for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | head -n 1`; do
		ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
		echo $i $ancestry1 $ancestry2 $ancestry3 $k
	
		rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/cov/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Regions.c2.${k}.Pathways*.cov.txt.gz
	done;
done;

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

#Original Vs1 Run (where issues and wonky plots were observed); also used for runs of other code versions that can use the exact same code setup (aside from swapping in/out filenames, such as pheno filename or cpp code filename, and changing associated output filenames)
#NOTE -- copy and pasted `/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjDrop.cpp` from associated Slack channel and from Lorin's code posted on 2010903, then made file `/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjDrop.mtEdits.SingleRun.vs1.cpp` to begin edits
##cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjDrop.mtEdits.SingleRun.vs1.cpp /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjDrop.mtEdits.SingleRun.vs1.copy1.cpp ...copy2-5...
#20180912 -- cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjDrop.mtEdits.SingleRun.vs1.cpp /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjDrop.mtEdits.SingleRun.vs1.GjProj.vs1.cpp
#20180912 -- cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjDrop.mtEdits.SingleRun.vs1.cpp /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjDrop.mtEdits.SingleRun.vs1.GjOnly.vs1.cpp
module load R/3.3.1; sleep 28800; for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep -E 'Height|BMI'`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep -v Plus`; do
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			NumSNPs=`zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | head -n 1 | perl -ane 'print scalar(@F);'` 
			NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt | wc | awk '{ print $1 }'`	
#			NumPaths=2
			echo $i $ancestry1 $ancestry2 $ancestry3 $k; 
			
			LpCnt=1; LpCnt2=1; for (( PathNum=1; PathNum <= $NumPaths; PathNum=PathNum+10 )); do
			       sbatch -t 72:00:00 --mem 16g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.GjOnly.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.GjOnly.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
				echo -e "\nR -q -e \"library(\\\"data.table\\\"); library(\\\"doParallel\\\"); library(\\\"Rcpp\\\"); library(\\\"RcppArmadillo\\\"); library(\\\"RcppParallel\\\"); sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjDrop.mtEdits.SingleRun.vs1.GjOnly.vs1.cpp\\\"); neg.is.na <- Negate(is.na); \
				Pathways.Regions <- list(); cores = detectCores(); InterPath.output <- list(); InterPath.output\\\$Est <- c(); InterPath.output\\\$Eigenvalues <- c(); InterPath.output\\\$PVE <- c(); \
				Pathways <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\\\", header=F); \
				for (i in $PathNum:($PathNum+9)) { Y <- c(); Y.Pheno <- c(); Y.Pheno.noNAs <- c(); \
					if (i > nrow(Pathways)) { Pathways.Regions[[1]] <- 1; Y.Pheno.noNAs <- rep(NA, nrow(InterPath.output\\\$Eigenvalues)); } else { Pathways.Regions[[1]] <- as.numeric(as.character(unlist(strsplit(as.character(Pathways[i,3]), \\\",\\\")))); \
					Y <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt\\\", header=T); Y.Pheno <- Y\\\$$i; Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; }; \ 
					if (length(Pathways.Regions[[1]]) > 1) { Data3 <- fread(paste(\\\"zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Regions.c2.${k}.Pathways\\\", as.character(i), \\\".txt.gz\\\", sep=\\\"\\\"), header=T); \ 
					Data3.mean <- apply(Data3, 2, mean); Data3.sd <- apply(Data3, 2, sd); Data3 <- t((t(Data3)-Data3.mean)/Data3.sd); \
					Data3.cov <- as.matrix(read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.cov.txt\\\", header=F)); \ 
					InterPath.output.temp <- list(); X <- Data3; X.cov <- Data3.cov; rm(Data3); rm(Data3.cov); \
					X.Pheno.noNAs <- X[neg.is.na(Y.Pheno),]; X.cov.Pheno.noNAs <- X.cov[neg.is.na(Y.Pheno),neg.is.na(Y.Pheno)]; \
					K <- 1/nrow(X.Pheno.noNAs) * tcrossprod(as.matrix(X.Pheno.noNAs)); \
					InterPath.output.temp <- InterPath(t(X.Pheno.noNAs),Y.Pheno.noNAs,as.matrix(X.cov.Pheno.noNAs),K,Pathways.Regions,nrow(X.Pheno.noNAs),as.numeric(as.character($NumSNPs)),cores=cores); InterPath.output\\\$Est <- c(InterPath.output\\\$Est, InterPath.output.temp\\\$Est); InterPath.output\\\$Eigenvalues <- cbind(InterPath.output\\\$Eigenvalues, InterPath.output.temp\\\$Eigenvalues); InterPath.output\\\$PVE <- c(InterPath.output\\\$PVE, InterPath.output.temp\\\$PVE); } else { InterPath.output\\\$Est <- c(InterPath.output\\\$Est, NA); InterPath.output\\\$Eigenvalues <- cbind(InterPath.output\\\$Eigenvalues, rep(NA, length(Y.Pheno.noNAs))); InterPath.output\\\$PVE <- c(InterPath.output\\\$PVE, NA);}; PH <- 1; }; \
				write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
				write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
				write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")
			done; sleep 2
		done; 
	done; 
done
#				LpCnt2=$((LpCnt2+1)); if [ $LpCnt2 == 5 ] ; then LpCnt2=1; fi
#				LpCnt=LpCnt+1; if [ $LpCnt == 50 ] ; then LpCnt=1; sleep 90; fi
	
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

#Orig Setup
#			       sbatch -t 72:00:00 --mem 10g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
#				...sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.SingleRun.vs2.cpp\\\");...
#					Y <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt\\\", header=T); Y.Pheno <- Y\\\$$i; Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; }; \ 
#				write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")
#yTyAdj Setup 		(this is just moving from 'q(2) = as.scalar(yc.t()*yc)' to 'q(2) = as_scalar(yc.t()*(eye<mat>(n,n)-(b*btb_inv)*b.t())*yc);'; NOTE -- that all other versions inherently got this change to begin with as well (ie just incorporated it into everything as I was going through the other versions, so just made the change everywhere at once))
#			       sbatch -t 72:00:00 --mem 10g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.yTyAdj.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.yTyAdj.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
#				...sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.SingleRun.vs2.cpp\\\");...
#					Y <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt\\\", header=T); Y.Pheno <- Y\\\$$i; Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; }; \ 
#				write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.yTyAdj.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.yTyAdj.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.yTyAdj.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")
#PCadj setup
#			       sbatch -t 72:00:00 --mem 10g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.PCadj.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.PCadj.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
#				...sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.SingleRun.vs2.cpp\\\");...
#					Y <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.flashpca.top10resids.txt\\\", header=T); Y.Pheno <- Y\\\$$i; Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; }; \ 
#				write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.PCadj.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.PCadj.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.PCadj.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")
#GjDrop setup
#			       sbatch -t 72:00:00 --mem 10g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.GjDrop.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.GjDrop.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
#				...sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjDrop.mtEdits.SingleRun.vs1.cpp\\\");...
#					Y <- read.table(paste(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.Regions.c2.${k}.Pathways\\\", i, \\\".txt.gz\\\", sep=\\\"\\\"), header=T); Y.Pheno <- Y\\\$$i; \ 
#				write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")
~#GjProj setup
~#			       sbatch -t 72:00:00 --mem 10g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.GjProj.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.GjProj.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
~#				...sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjProj.mtEdits.SingleRun.vs1.cpp\\\");...
~#					Y <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.flashpca.top10resids.txt\\\", header=T); Y.Pheno <- Y\\\$$i; Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; }; \ 
~#				write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjProj.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
~#				write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjProj.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
~#				write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjProj.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")
#GjOnly setup
#			       sbatch -t 72:00:00 --mem 10g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.GjOnly.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.GjOnly.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
#				...sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjDrop.mtEdits.SingleRun.vs1.GjOnly.vs1.cpp\\\");...
#					Y <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt\\\", header=T); Y.Pheno <- Y\\\$$i; Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; }; \ 
#				write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")


#From https://stackoverflow.com/questions/8903239/how-to-calculate-time-difference-in-bash-script
for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep -E 'Height|BMI'`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep -v Plus`; do
			SECONDS=0;
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
			ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt | wc | awk '{ print $1 }'`	
#			NumPaths=20
			
			echo $i $ancestry1 $ancestry2 $ancestry3 $k
			for (( PathNum=1; PathNum <= $NumPaths; PathNum=PathNum+10 )); do
				R -q -e "library(\"CompQuadForm\"); neg.is.na <- Negate(is.na); \
				Pathways <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\", header=F); \ 
				InterPath.output.Est <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.yTyAdj.Paths${PathNum}.Est.txt\", header=F); \
				InterPath.output.Eigenvalues <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.yTyAdj.Paths${PathNum}.Eigenvalues.txt\", header=F); \
				InterPath.output.PVE <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.yTyAdj.Paths${PathNum}.PVE.txt\", header=F); \
				Results1 <- c(); Counter1 <- 1; for (i in $PathNum:($PathNum+9)) { \
					if (neg.is.na(InterPath.output.Est[Counter1,1])) { \ 
						Lambda <- sort(InterPath.output.Eigenvalues[,Counter1], decreasing=TRUE); \
						Davies.Output <- davies(InterPath.output.Est[Counter1,1], lambda=Lambda, acc=1e-8); \
						pVal <- 2*min(1-Davies.Output\$Qq, Davies.Output\$Qq); \
						Results1 <- rbind(Results1, c(as.character(Pathways[i,1]), InterPath.output.Est[Counter1,1], InterPath.output.PVE[Counter1,1], pVal)); \ 
						Counter1 = Counter1 + 1; \
					}; \
				}; write.table(Results1, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);"
			done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.yTyAdj.AllPaths.Results.txt.pre.gz
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

#From `/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Source.Simulations.vs2.R` onwards, copied/edited recal procedure from Isabella & Lorin
for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | head -n 1`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep -v Plus | head -n 1`; do
			SECONDS=0;
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
			ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			
			echo $i $ancestry1 $ancestry2 $ancestry3 $k
			R -q -e "library(\"CompQuadForm\"); neg.is.na <- Negate(is.na); \
			InterPath.AllPaths.output <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.PCadj.AllPaths.Results.txt.pre.gz\", header=F); \ 
			pVals <- InterPath.AllPaths.output[,ncol(InterPath.AllPaths.output)]; pVals[pVals >= 0]; \
			quants <- qchisq(pVals, df=1, lower.tail=FALSE); \
			quants.medn <- median(quants); m <- 0.4549; \
			pVals <- pchisq((quants/quants.medn) * m, df=1, lower.tail=FALSE); \
			InterPath.AllPaths.output[,ncol(InterPath.AllPaths.output)] <- pVals; \
			write.table(InterPath.AllPaths.output, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);" | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.PCadj.AllPaths.Results.Recal.txt.pre.gz


			write.table(Results1, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);"
			done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.PCadj.AllPaths.Results.txt.pre.gz
			duration=$SECONDS; echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
		done;
	done;
done;

#  pvals = davies.pvals[davies.pvals>=0]; summary(pvals)
#  q = qchisq(pvals,df=1,lower.tail = FALSE); summary(q)
#  gs = median(q); m = 0.4549; summary((q/gs)*m)
#  pvals = pchisq((q/gs)*m,df=1,lower.tail = FALSE)

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




#for l in `cat <(echo "yTyAdj GjDrop PCadj wCov GjDrop_wCov" | perl -lane 'print join("\n", @F);')`; do
#for l in `cat <(echo "GjOnly GjOnly_wCov" | perl -lane 'print join("\n", @F);')`; do
	for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep -E 'Height|BMI'`; do
		for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
			for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep -v Plus`; do
				SECONDS=0;
				ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
				ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
				
				echo $i $j $k $ancestry1 $ancestry2 $ancestry3 
				
				if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/PrevRuns ]; then
					mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/PrevRuns
				fi
				if [ ! -d /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurmn/PrevRuns ]; then
					mkdir /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/PrevRuns
				fi
	
				cd  /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k 
				tar -cvzf /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/PrevRuns/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.${l}.AllPaths.Est.tar.gz ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.${l}.Paths*.Est.txt 
				tar -cvzf /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/PrevRuns/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.${l}.AllPaths.Eigenvalues.tar.gz ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.${l}.Paths*.Eigenvalues.txt 
				tar -cvzf /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/PrevRuns/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.${l}.AllPaths.PVE.tar.gz ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.${l}.Paths*.PVE.txt 
	
				cd  /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm 
				tar -cvzf /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/PrevRuns/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.${l}.AllPaths.slurm.output.tar.gz ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.${l}.Pathways*.slurm.output 
				tar -cvzf /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/PrevRuns/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.${l}.AllPaths.slurm.error.tar.gz ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.${l}.Pathways*.slurm.error 
	
				duration=$SECONDS; echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
			done;
		done;
	done;
done;

#for l in `cat <(echo "yTyAdj GjDrop PCadj wCov GjDrop_wCov" | perl -lane 'print join("\n", @F);')`; do
for l in `cat <(echo "GjOnly GjOnly_wCov" | perl -lane 'print join("\n", @F);')`; do
	for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep -E 'Height|BMI'`; do
		for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
			for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep -v Plus`; do
				ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
				ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
				
				echo $i $j $k $ancestry1 $ancestry2 $ancestry3 
			
				rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.${l}.Paths*.Est.txt
				rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.${l}.Paths*.Eigenvalues.txt
				rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.${l}.Paths*.PVE.txt
	
				rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.${l}.Pathways*.slurm.output
				rm /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.${l}.Pathways*.slurm.error
	
			done;
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
#20180904 NOTE -- cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjDrop.mtEdits.SingleRun.vs1.cpp /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjDrop.mtEdits.SingleRun.vs1.wCovs.vs1.cpp
#20180912 -- cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjDrop.mtEdits.SingleRun.vs1.wCovs.vs1.cpp /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjOnly.mtEdits.SingleRun.vs1.wCovs.vs1.cpp 
module load R/3.3.1; sleep 25200; for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep -E 'Height|BMI' | grep BMI`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep -v Plus`; do
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`; NumSNPs=`zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | head -n 1 | perl -ane 'print scalar(@F);'` 
			NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt | wc | awk '{ print $1 }'`	
#			NumPaths=2
			echo $i $ancestry1 $ancestry2 $ancestry3 $k
			
			LpCnt=1; LpCnt2=1; for (( PathNum=1; PathNum <= $NumPaths; PathNum=PathNum+10 )); do
				sbatch -t 72:00:00 --mem 16g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.GjOnly_wCov.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.GjOnly_wCov.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
				echo -e "\nR -q -e \"library(\\\"data.table\\\"); library(\\\"doParallel\\\"); library(\\\"Rcpp\\\"); library(\\\"RcppArmadillo\\\"); library(\\\"RcppParallel\\\"); sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjOnly.mtEdits.SingleRun.vs1.wCovs.vs1.cpp\\\"); neg.is.na <- Negate(is.na); \
				Covars <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.sort.ImptHRC.dose.100geno.raw.txt\\\", header=T); Pathways <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\\\", header=F); \
				Pathways.Regions <- list(); cores = detectCores(); InterPath.output <- list(); InterPath.output\\\$Est <- c(); InterPath.output\\\$Eigenvalues <- c(); InterPath.output\\\$PVE <- c(); \ 
				for (i in $PathNum:($PathNum+9)) { Y <- c(); Y.Pheno <- c(); Y.Pheno.noNAs <- c(); \
					if (i > nrow(Pathways)) { Pathways.Regions[[1]] <- 1; Y.Pheno.noNAs <- rep(NA, nrow(InterPath.output\\\$Eigenvalues)); } else { Pathways.Regions[[1]] <- as.numeric(as.character(unlist(strsplit(as.character(Pathways[i,3]), \\\",\\\")))); \
					Y <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt\\\", header=T); Y.Pheno <- Y\\\$$i; Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; }; \ 
					if (length(Pathways.Regions[[1]]) > 1) { Data3 <- fread(paste(\\\"zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Regions.c2.${k}.Pathways\\\", as.character(i), \\\".txt.gz\\\", sep=\\\"\\\"), header=T); Data3.cov <- as.matrix(read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.cov.txt\\\", header=F)); \ 
					Data3.mean <- apply(Data3, 2, mean); Data3.sd <- apply(Data3, 2, sd); Data3 <- t((t(Data3)-Data3.mean)/Data3.sd); \ 
					InterPath.output.temp <- list(); X <- Data3; X.cov <- Data3.cov; rm(Data3); rm(Data3.cov); X.Pheno.noNAs <- X[neg.is.na(Y.Pheno),]; X.cov.Pheno.noNAs <- X.cov[neg.is.na(Y.Pheno),neg.is.na(Y.Pheno)]; Z <- Covars[neg.is.na(Y.Pheno),(ncol(Covars)-9):ncol(Covars)]; \
					G <- 1/nrow(X.Pheno.noNAs) * tcrossprod(as.matrix(X.Pheno.noNAs)); \
					InterPath.output.temp <- InterPath(t(X.Pheno.noNAs),Y.Pheno.noNAs,as.matrix(X.cov.Pheno.noNAs),G,t(as.matrix(Z)),Pathways.Regions,nrow(X.Pheno.noNAs),as.numeric(as.character($NumSNPs)),cores=cores); InterPath.output\\\$Est <- c(InterPath.output\\\$Est, InterPath.output.temp\\\$Est); InterPath.output\\\$Eigenvalues <- cbind(InterPath.output\\\$Eigenvalues, InterPath.output.temp\\\$Eigenvalues); InterPath.output\\\$PVE <- c(InterPath.output\\\$PVE, InterPath.output.temp\\\$PVE); } else { InterPath.output\\\$Est <- c(InterPath.output\\\$Est, NA); InterPath.output\\\$Eigenvalues <- cbind(InterPath.output\\\$Eigenvalues, rep(NA, length(Y.Pheno.noNAs))); InterPath.output\\\$PVE <- c(InterPath.output\\\$PVE, NA);};}; \
				write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly_wCov.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
				write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly_wCov.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
				write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly_wCov.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")
			done; sleep 2
		done; 
	done; 
done
#				LpCnt2=$((LpCnt2+1)); if [ $LpCnt2 == 5 ] ; then LpCnt2=1; fi
#				LpCnt=LpCnt+1; if [ $LpCnt == 100 ] ; then LpCnt=1; sleep 30; fi

#wCov setup
#				sbatch -t 72:00:00 --mem 10g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.wCov.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.wCov.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
#				...sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.mtEdits.SingleRun.vs2.wCovs.vs1.cpp\\\");...
#					Y <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt\\\", header=T); Y.Pheno <- Y\\\$$i; Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; }; \ 
#				write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.wCov.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.wCov.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.wCov.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")
#GjDrop + wCov setup
#				sbatch -t 72:00:00 --mem 10g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.GjDrop_wCov.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.GjDrop_wCov.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
#				...sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjDrop.mtEdits.SingleRun.vs1.wCovs.vs1.cpp\\\");...
#					Y <- read.table(paste(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.Regions.c2.${k}.Pathways\\\", i, \\\".txt.gz\\\", sep=\\\"\\\"), header=T); Y.Pheno <- Y\\\$$i; Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; }; \ 
#				write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop_wCov.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop_wCov.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop_wCov.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")
#GjOnly + wCov setup
#				sbatch -t 72:00:00 --mem 10g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.GjOnly_wCov.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.GjOnly_wCov.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
#				...sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.GjOnly.mtEdits.SingleRun.vs1.wCovs.vs1.cpp\\\");...
#					Y <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt\\\", header=T); Y.Pheno <- Y\\\$$i; Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; }; \ 
#				write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly_wCov.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly_wCov.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly_wCov.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")

#Vs1 Results Collection
for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep -E 'Height|BMI'`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep -v Plus`; do
			SECONDS=0;
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
			ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt | wc | awk '{ print $1 }'`	
#			NumPaths=20

			echo $i $ancestry1 $ancestry2 $ancestry3 $k
			for (( PathNum=1; PathNum <= $NumPaths; PathNum=PathNum+10 )); do
				R -q -e "library(\"CompQuadForm\"); neg.is.na <- Negate(is.na); \
				Pathways <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\", header=F); \ 
				InterPath.output.Est <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.Paths${PathNum}.Est.txt\", header=F); \
				InterPath.output.Eigenvalues <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.Paths${PathNum}.Eigenvalues.txt\", header=F); \
				InterPath.output.PVE <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.Paths${PathNum}.PVE.txt\", header=F); \
				Results1 <- c(); Counter1 <- 1; for (i in $PathNum:($PathNum+9)) { \
					if (neg.is.na(InterPath.output.Est[Counter1,1])) { \ 
						Lambda <- sort(InterPath.output.Eigenvalues[,Counter1], decreasing=TRUE); \
						Davies.Output <- davies(InterPath.output.Est[Counter1,1], lambda=Lambda, acc=1e-8); \
						pVal <- 2*min(1-Davies.Output\$Qq, Davies.Output\$Qq); \
						Results1 <- rbind(Results1, c(as.character(Pathways[i,1]), InterPath.output.Est[Counter1,1], InterPath.output.PVE[Counter1,1], pVal)); \ 
						Counter1 = Counter1 + 1; \
					}; \
				}; write.table(Results1, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);"
			done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.AllPaths.Results.txt.pre.gz
			duration=$SECONDS; echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
		done;
	done;
done;

#Orig
#				InterPath.output.Est <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Paths${PathNum}.Est.txt\", header=F); \
#				InterPath.output.Eigenvalues <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Paths${PathNum}.Eigenvalues.txt\", header=F); \
#				InterPath.output.PVE <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Paths${PathNum}.PVE.txt\", header=F); \
#			done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.AllPaths.Results.txt.pre.gz
#yTyAdj
#				InterPath.output.Est <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.yTyAdj.Paths${PathNum}.Est.txt\", header=F); \
#				InterPath.output.Eigenvalues <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.yTyAdj.Paths${PathNum}.Eigenvalues.txt\", header=F); \
#				InterPath.output.PVE <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.yTyAdj.Paths${PathNum}.PVE.txt\", header=F); \
#			done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.yTyAdj.AllPaths.Results.txt.pre.gz
#PCadj
#				InterPath.output.Est <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.PCadj.Paths${PathNum}.Est.txt\", header=F); \
#				InterPath.output.Eigenvalues <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.PCadj.Paths${PathNum}.Eigenvalues.txt\", header=F); \
#				InterPath.output.PVE <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.PCadj.Paths${PathNum}.PVE.txt\", header=F); \
#			done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.PCadj.AllPaths.Results.txt.pre.gz
#GjDrop
#				InterPath.output.Est <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop.Paths${PathNum}.Est.txt\", header=F); \
#				InterPath.output.Eigenvalues <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop.Paths${PathNum}.Eigenvalues.txt\", header=F); \
#				InterPath.output.PVE <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop.Paths${PathNum}.PVE.txt\", header=F); \
#			done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop.AllPaths.Results.txt.pre.gz
#wCov
#				InterPath.output.Est <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.wCov.Paths${PathNum}.Est.txt\", header=F); \
#				InterPath.output.Eigenvalues <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.wCov.Paths${PathNum}.Eigenvalues.txt\", header=F); \
#				InterPath.output.PVE <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.wCov.Paths${PathNum}.PVE.txt\", header=F); \
#			done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.wCov.AllPaths.Results.txt.pre.gz
#GjDrop + wCov
#				InterPath.output.Est <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop_wCov.Paths${PathNum}.Est.txt\", header=F); \
#				InterPath.output.Eigenvalues <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop_wCov.Paths${PathNum}.Eigenvalues.txt\", header=F); \
#				InterPath.output.PVE <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop_wCov.Paths${PathNum}.PVE.txt\", header=F); \
#			done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop_wCov.AllPaths.Results.txt.pre.gz
#GjOnly
#				InterPath.output.Est <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.Paths${PathNum}.Est.txt\", header=F); \
#				InterPath.output.Eigenvalues <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.Paths${PathNum}.Eigenvalues.txt\", header=F); \
#				InterPath.output.PVE <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.Paths${PathNum}.PVE.txt\", header=F); \
#			done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.AllPaths.Results.txt.pre.gz
#GjOnly + wCov
#				InterPath.output.Est <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly_wCov.Paths${PathNum}.Est.txt\", header=F); \
#				InterPath.output.Eigenvalues <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly_wCov.Paths${PathNum}.Eigenvalues.txt\", header=F); \
#				InterPath.output.PVE <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly_wCov.Paths${PathNum}.PVE.txt\", header=F); \
#			done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly_wCov.AllPaths.Results.txt.pre.gz

#2018090 NOTE -- mv /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop.AllPaths.Results.txt.pre.gz /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/OLD1.ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop.AllPaths.Results.txt.pre.gz
#2018090 NOTE -- mv /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop.AllPaths.Results.Recal.txt.pre.gz /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/OLD1.ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop.AllPaths.Results.Recal.txt.pre.gz

#From `/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Source.Simulations.vs2.R` onwards, copied/edited recal procedure from Isabella & Lorin
for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep -E 'Height|BMI'`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep -v Plus`; do
			SECONDS=0;
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
			ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			
			R -q -e "library(\"CompQuadForm\"); neg.is.na <- Negate(is.na); \
			InterPath.AllPaths.output <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.AllPaths.Results.txt.pre.gz\", header=F); \ 
			pVals <- InterPath.AllPaths.output[,ncol(InterPath.AllPaths.output)]; pVals.no0 <- pVals[pVals > 0]; \
			quants <- qchisq(pVals.no0, df=1, lower.tail=FALSE); \
			quants.medn <- median(quants); m <- 0.4549; \
			pVals.no0 <- pchisq((quants/quants.medn) * m, df=1, lower.tail=FALSE); \
			InterPath.AllPaths.output <- InterPath.AllPaths.output[pVals > 0,];  InterPath.AllPaths.output[,ncol(InterPath.AllPaths.output)] <- pVals.no0; \
			write.table(InterPath.AllPaths.output, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);" | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.AllPaths.Results.Recal.txt.pre.gz
			duration=$SECONDS; echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
		done;
	done;
done;

#  pvals = davies.pvals[davies.pvals>=0]; summary(pvals)
#  q = qchisq(pvals,df=1,lower.tail = FALSE); summary(q)
#  gs = median(q); m = 0.4549; summary((q/gs)*m)
#  pvals = pchisq((q/gs)*m,df=1,lower.tail = FALSE)

#Orig
#			InterPath.AllPaths.output <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.AllPaths.Results.txt.pre.gz\", header=F); \ 
#			write.table(InterPath.AllPaths.output, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);" | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.AllPaths.Results.Recal.txt.pre.gz
#yTyAdj
#			InterPath.AllPaths.output <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.yTyAdj.AllPaths.Results.txt.pre.gz\", header=F); \ 
#			write.table(InterPath.AllPaths.output, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);" | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.yTyAdj.AllPaths.Results.Recal.txt.pre.gz
#PCadj
#			InterPath.AllPaths.output <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.PCadj.AllPaths.Results.txt.pre.gz\", header=F); \ 
#			write.table(InterPath.AllPaths.output, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);" | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.PCadj.AllPaths.Results.Recal.txt.pre.gz
#GjDrop
#20180906 NOTE -- reran this under BMI with `...pVals[pVals > 0];...` because was getting so many 0s leading to complete breakdown of the procedure (just got 1s and NAs as a result); the lack of correction for population stratification is just pushing the results to be too extreme (eg 'p-values so low they become zero'), which we know not to be true, so just getting around this to be able to plot things properly later in R (not reflecting a belief at all in the results being true or reasonable)
#			InterPath.AllPaths.output <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop.AllPaths.Results.txt.pre.gz\", header=F); \ 
#			write.table(InterPath.AllPaths.output, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);" | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop.AllPaths.Results.Recal.txt.pre.gz
#wCov
#			InterPath.AllPaths.output <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.wCov.AllPaths.Results.txt.pre.gz\", header=F); \ 
#			write.table(InterPath.AllPaths.output, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);" | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.wCov.AllPaths.Results.Recal.txt.pre.gz
#GjDrop + wCov
#			InterPath.AllPaths.output <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop_wCov.AllPaths.Results.txt.pre.gz\", header=F); \ 
#			write.table(InterPath.AllPaths.output, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);" | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjDrop_wCov.AllPaths.Results.Recal.txt.pre.gz
#GjOnly
#			InterPath.AllPaths.output <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.AllPaths.Results.txt.pre.gz\", header=F); \ 
#			write.table(InterPath.AllPaths.output, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);" | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly.AllPaths.Results.Recal.txt.pre.gz
#GjOnly + wCov
#			InterPath.AllPaths.output <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly_wCov.AllPaths.Results.txt.pre.gz\", header=F); \ 
#			write.table(InterPath.AllPaths.output, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);" | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.GjOnly_wCov.AllPaths.Results.Recal.txt.pre.gz

mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses
mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd1AdditiveMdls

#R -q -e "library(\"RColorBrewer\"); DataTypes <- c(\"yTyAdj\", \"GjDrop\"); \
#R -q -e "library(\"RColorBrewer\"); DataTypes <- c(\"yTyAdj\", \"GjDrop\", \"PCadj\", \"wCov\", \"GjDrop_wCov\"); \ 
R -q -e "library(\"RColorBrewer\"); DataTypes <- c(\"yTyAdj\", \"GjOnly\", \"PCadj\", \"wCov\", \"GjOnly_wCov\"); \ 
R -q -e "library(\"RColorBrewer\"); DataTypes <- c(\"GjOnly\",  \"GjOnly_wCov\"); \
	png(\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd1AdditiveMdls/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.HeightBMI.NonSynExonic.Rnd1Vrsns.AllPaths.Results.vs2.png\", height=4000, width=4500, res=300); par(oma=c(1,1,4,14), mar=c(5,5,4,2), mfrow=c(2,2)); \ 
	for (i in DataTypes) { \
		Data1a <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.\", i, \".AllPaths.Results.txt.pre.gz\", sep=\"\"), header=F); \	
		Data2a <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.Exonic.\", i, \".AllPaths.Results.txt.pre.gz\", sep=\"\"), header=F); \	
		Data3a <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/BMI/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.BMI.NonSyn.\", i, \".AllPaths.Results.txt.pre.gz\", sep=\"\"), header=F); \
		Data4a <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/BMI/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.BMI.Exonic.\", i, \".AllPaths.Results.txt.pre.gz\", sep=\"\"), header=F); \
		Data1b <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.\", i, \".AllPaths.Results.Recal.txt.pre.gz\", sep=\"\"), header=F); \
		Data2b <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.Exonic.\", i, \".AllPaths.Results.Recal.txt.pre.gz\", sep=\"\"), header=F); \
		Data3b <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/BMI/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.BMI.NonSyn.\", i, \".AllPaths.Results.Recal.txt.pre.gz\", sep=\"\"), header=F); \
		Data4b <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/BMI/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.BMI.Exonic.\", i, \".AllPaths.Results.Recal.txt.pre.gz\", sep=\"\"), header=F); \	
		Data1a <- Data1a[Data1a[,4] > 0,]; Data2a <- Data2a[Data2a[,4] > 0,]; Data3a <- Data3a[Data3a[,4] > 0,]; Data4a <- Data4a[Data4a[,4] > 0,]; Data1b <- Data1b[Data1b[,4] > 0,]; Data2b <- Data2b[Data2b[,4] > 0,]; Data3b <- Data3b[Data3b[,4] > 0,]; Data4b <- Data4b[Data4b[,4] > 0,]; \
		xVals1a <- seq(1/nrow(Data1a), 1, by=1/nrow(Data1a)); xVals2a <- seq(1/nrow(Data2a), 1, by=1/nrow(Data2a)); xVals3a <- seq(1/nrow(Data3a), 1, by=1/nrow(Data3a)); xVals4a <- seq(1/nrow(Data4a), 1, by=1/nrow(Data4a)); \ 
		xVals1b <- seq(1/nrow(Data1b), 1, by=1/nrow(Data1b)); xVals2b <- seq(1/nrow(Data2b), 1, by=1/nrow(Data2b)); xVals3b <- seq(1/nrow(Data3b), 1, by=1/nrow(Data3b)); xVals4b <- seq(1/nrow(Data4b), 1, by=1/nrow(Data4b)); \
		xalimMax <- max(c(-log10(xVals1a), -log10(xVals2a), -log10(xVals3a), -log10(xVals4a))); yalimMax <- max(c(-log10(Data1a[,4]), -log10(Data2a[,4]), -log10(Data3a[,4]), -log10(Data4a[,4]))); \
		xblimMax <- max(c(-log10(xVals1b), -log10(xVals2b), -log10(xVals3b), -log10(xVals4b))); yblimMax <- max(c(-log10(Data1b[,4]), -log10(Data2b[,4]), -log10(Data3b[,4]), -log10(Data4b[,4]))); \
		\
		Title1 <- c(); if (i == \"yTyAdj\") { Title1 <- \"Orig\"; } else if (i == \"GjDrop\") { Title1 <- \"GjDrop\"; } else if (i == \"PCadj\") { Title1 <- \"PCRegress\"; } else if (i == \"wCov\") { Title1 <- \"PCProject\"; } else if (i == \"GjDrop_wCov\") { Title1 <- \"GjDrop+PCProject\"; } else if (i == \"GjOnly\") { Title1 <- \"GjOnly\"; } else if (i == \"GjOnly_wCov\") { Title1 <- \"GjOnly+PCProject\"; } else { Title1 <- \"Error1a\"; }; \	
		\
		plot(-log10(xVals1a[order(xVals1a, decreasing=TRUE)]), -log10(Data1a[order(Data1a[,4], decreasing=TRUE),4]), main=Title1, xlab=\"-log10(Expected pValues)\", ylab=\"-log10(Observed pValues)\", xlim=c(0,xalimMax), ylim=c(0,yalimMax), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[7], cex=1.5, cex.main=1.5, cex.axis=1.5, cex.lab=1.5); \
		points(-log10(xVals2a[order(xVals2a, decreasing=TRUE)]), -log10(Data2a[order(Data2a[,4], decreasing=TRUE),4]), type=\"b\", pch=17, col=brewer.pal(12, \"Paired\")[7], cex=1.5); \ 
		points(-log10(xVals3a[order(xVals3a, decreasing=TRUE)]), -log10(Data3a[order(Data3a[,4], decreasing=TRUE),4]), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[3], cex=1.5); \
		points(-log10(xVals4a[order(xVals4a, decreasing=TRUE)]), -log10(Data4a[order(Data4a[,4], decreasing=TRUE),4]), type=\"b\", pch=17, col=brewer.pal(12, \"Paired\")[3], cex=1.5); \
		abline(0,1, lwd=2, col=\"BLACK\"); \
		plot(-log10(xVals1b[order(xVals1b, decreasing=TRUE)]), -log10(Data1b[order(Data1b[,4], decreasing=TRUE),4]), main=paste(Title1, \".Recal\", sep=\"\"), xlab=\"-log10(Expected pValues)\", ylab=\"-log10(Observed pValues)\", xlim=c(0,xblimMax), ylim=c(0,yblimMax), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[7], cex=1.5, cex.main=1.5, cex.axis=1.5, cex.lab=1.5); \ 
		points(-log10(xVals2b[order(xVals2b, decreasing=TRUE)]), -log10(Data2b[order(Data2b[,4], decreasing=TRUE),4]), type=\"b\", pch=17, col=brewer.pal(12, \"Paired\")[7], cex=1.5); \ 
		points(-log10(xVals3b[order(xVals3b, decreasing=TRUE)]), -log10(Data3b[order(Data3b[,4], decreasing=TRUE),4]), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[3], cex=1.5); \
		points(-log10(xVals4b[order(xVals4b, decreasing=TRUE)]), -log10(Data4b[order(Data4b[,4], decreasing=TRUE),4]), type=\"b\", pch=17, col=brewer.pal(12, \"Paired\")[3], cex=1.5); \ 
		abline(0,1, lwd=2, col=\"BLACK\"); \
		\
	}; mtext(\"African\", line=-.75, outer=TRUE, cex=2.5); par(fig = c(0, 1, 0, 1), mfrow=c(1,1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot(0, 0, type = \"n\", bty = \"n\", xaxt = \"n\", yaxt = \"n\"); legend(\"topright\", c(\"Height\", \"BMI\", \"NonSyn\", \"Exonic\"), pch=c(NA, NA, 16, 17), text.col=c(brewer.pal(12, \"Paired\")[7], brewer.pal(12, \"Paired\")[3], \"BLACK\", \"BLACK\"), xpd=TRUE, inset=c(.027,.0385), bg=\"transparent\", cex=1.35, y.intersp=2); dev.off(); \
"

#		print(quantile(-log10(Data1a[,4]))); print(quantile(-log10(Data2a[,4]))); print(quantile(-log10(Data3a[,4]))); print(quantile(-log10(Data4a[,4]))); \
#		print(table(is.na(Data1a[,4]))); print(table(is.na(Data2a[,4]))); print(table(is.na(Data3a[,4]))); print(table(is.na(Data4a[,4]))); print(table(is.na(Data1b[,4]))); print(table(is.na(Data2b[,4]))); print(table(is.na(Data3b[,4]))); print(table(is.na(Data4b[,4]))); \
#R -q -e "vals1 <- rnorm(1000); vals2 <- rnorm(1000); png(\"nana.png\", height=2000, width=2000, res=300); plot(vals1, vals2, main=\"yoyo1\", xlab=\"nana1\", ylab=\"nana2\", cex=1, cex.main=20, cex.axis=3, cex.lab=3); dev.off();"
#	png(\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd1AdditiveMdls/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.HeightBMI.NonSynExonic.Rnd1Vrsns.AllPaths.Results.vs2.png\", height=10000, width=4250, res=300); par(oma=c(1,1,4,14), mar=c(5,5,4,2), mfrow=c(5,2)); \ 

#From MacBook Air
#mkdir /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Vs1
#mkdir /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Vs1/Analyses
#mkdir /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd1AdditiveMdls/
#scp -p  mturchin@ssh.ccv.brown.edu:/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd1AdditiveMdls/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.HeightBMI.NonSynExonic.Rnd1Vrsns.AllPaths.Results.vs1.png /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd1AdditiveMdls/.
#scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/pathways/nana.png /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd1AdditiveMdls/.
scp -p  mturchin@ssh.ccv.brown.edu:/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd1AdditiveMdls/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.HeightBMI.NonSynExonic.Rnd1Vrsns.AllPaths.Results.vs*.png /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd1AdditiveMdls/.


#join <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.PCadj.AllPaths.Results.txt.pre.gz | sort -g -k 4,4 | perl -lane 'if ($. == 1) { $count1 = 1; } print join("\t", @F), "\t", $count1; $count1++;' | sort -g -k 1,1 | awk '{ print $1 "\t" $5 }') <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.wCov.AllPaths.Results.txt.pre.gz | sort -g -k 4,4 | perl -lane 'if ($. == 1) { $count1 = 1; } print join("\t", @F), "\t", $count1; $count1++;' | sort -g -k 1,1 | awk '{ print $1 "\t" $5 }') 
join <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.PCadj.AllPaths.Results.txt.pre.gz | sort -g -k 4,4 | perl -lane 'if ($. == 1) { $count1 = 1; } print join("\t", @F), "\t", $count1; $count1++;' | sort -g -k 1,1 | awk '{ print $1 "\t" $5 }') <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.wCov.AllPaths.Results.txt.pre.gz | sort -g -k 4,4 | perl -lane 'if ($. == 1) { $count1 = 1; } print join("\t", @F), "\t", $count1; $count1++;' | sort -g -k 1,1 | awk '{ print $1 "\t" $5 }') | sort -g -k 2,2 | R -q -e "Data1 <- read.table(file('stdin'), header=F); Results1 <- c(); for (i in seq(1,nrow(Data1)-101,by=25)) { Results1 <- rbind(Results1, c(mean(Data1[i:(i+100),2] - Data1[i:(i+100),3]), mean(abs(Data1[i:(i+100),2] - Data1[i:(i+100),3])), sd(Data1[i:(i+100),2] - Data1[i:(i+100),3]), sd(abs(Data1[i:(i+100),2] - Data1[i:(i+100),3])))); }; print(Results1);" | vi -
join <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.PCadj.AllPaths.Results.txt.pre.gz | sort -g -k 4,4 | perl -lane 'if ($. == 1) { $count1 = 1; } print join("\t", @F), "\t", $count1; $count1++;' | sort -g -k 1,1 | awk '{ print $1 "\t" $5 }') <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.GjDrop.AllPaths.Results.txt.pre.gz | sort -g -k 4,4 | perl -lane 'if ($. == 1) { $count1 = 1; } print join("\t", @F), "\t", $count1; $count1++;' | sort -g -k 1,1 | awk '{ print $1 "\t" $5 }') | sort -g -k 2,2 | R -q -e "Data1 <- read.table(file('stdin'), header=F); Results1 <- c(); for (i in seq(1,nrow(Data1)-101,by=25)) { Results1 <- rbind(Results1, c(mean(Data1[i:(i+100),2] - Data1[i:(i+100),3]), mean(abs(Data1[i:(i+100),2] - Data1[i:(i+100),3])), sd(Data1[i:(i+100),2] - Data1[i:(i+100),3]), sd(abs(Data1[i:(i+100),2] - Data1[i:(i+100),3])))); }; print(Results1);" | vi -

join <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.PCadj.AllPaths.Results.txt.pre.gz | sort -g -k 4,4 | perl -lane 'if ($. == 1) { $count1 = 1; } print join("\t", @F), "\t", $count1; $count1++;' | sort -g -k 1,1 | awk '{ print $1 "\t" $5 }') <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.wCov.AllPaths.Results.txt.pre.gz | sort -g -k 4,4 | perl -lane 'if ($. == 1) { $count1 = 1; } print join("\t", @F), "\t", $count1; $count1++;' | sort -g -k 1,1 | awk '{ print $1 "\t" $5 }') | sort -g -k 2,2 | R -q -e "Data1 <- read.table(file('stdin'), header=F); png(\"/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.PCadj.AllPaths.Results.Comp.wCov.vs1.png\", width=2000, height=2000, res=300); par(mar=c(5,5,4,2)); plot(seq(1,nrow(Data1),by=1), Data1[,3], main=\"African Height NonSyn -- PCRegress vs. PCProject\", xlab=\"PCRegress Gene Order\", ylab=\"Ranking\", type=\"l\", lwd=1, col=\"RED\", cex=1, cex.main=1, cex.lab=1, cex.axis=1); lines(seq(1,nrow(Data1),by=1), Data1[,2], lwd=1, col=\"BLUE\", cex=1); legend(\"topleft\", c(\"PCRegress\", \"PCProject\"), lty=c(1,1), col=c(\"BLUE\", \"RED\"), cex=1); dev.off();"
join <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.PCadj.AllPaths.Results.txt.pre.gz | sort -g -k 4,4 | perl -lane 'if ($. == 1) { $count1 = 1; } print join("\t", @F), "\t", $count1; $count1++;' | sort -g -k 1,1 | awk '{ print $1 "\t" $5 }') <(zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.yTyAdj.AllPaths.Results.txt.pre.gz | sort -g -k 4,4 | perl -lane 'if ($. == 1) { $count1 = 1; } print join("\t", @F), "\t", $count1; $count1++;' | sort -g -k 1,1 | awk '{ print $1 "\t" $5 }') | sort -g -k 2,2 | R -q -e "Data1 <- read.table(file('stdin'), header=F); png(\"/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.PCadj.AllPaths.Results.Comp.yTyAdj.vs1.png\", width=2000, height=2000, res=300); par(mar=c(5,5,4,2)); plot(seq(1,nrow(Data1),by=1), Data1[,3], main=\"African Height NonSyn -- PCRegress vs. Original\", xlab=\"PCRegress Gene Order\", ylab=\"Ranking\", type=\"l\", lwd=1, col=\"RED\", cex=1, cex.main=1, cex.lab=1, cex.axis=1); lines(seq(1,nrow(Data1),by=1), Data1[,2], lwd=1, col=\"BLUE\", cex=1); legend(\"topleft\", c(\"PCRegress\", \"Original\"), lty=c(1,1), col=c(\"BLUE\", \"RED\"), cex=1); dev.off();"

scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.PCadj.AllPaths.Results.Comp.wCov.vs1.png /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd1AdditiveMdls/.
scp -p mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.PCadj.AllPaths.Results.Comp.*.vs1.png /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd1AdditiveMdls/.











#20180820
#Vs2 Runs (moving towards using epistatic interaction models) 

#pathway*remaining genome
#NOTE -- copy and pasted `/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Source.Vs2.cpp` & `/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Source.Simulations.Vs2.R` from associated Slack channel and from Lorin's code posted on 20180731
#cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Source.Vs2.cpp /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Vs2.GjDrop.mtEdits.SingleRun.vs1.wCovs.vs1.cpp
#cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Vs2.GjDrop.mtEdits.SingleRun.vs1.wCovs.vs1.cpp /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Vs2.GjDrop.mtEdits.SingleRun.vs1.wCovs.vs1.GG.cpp
module load R/3.3.1; sleep 25200; for i in `cat <(echo "Height;1254 BMI;58923 Waist;49281 Hip;37485" | perl -lane 'print join("\n", @F);') | grep -E 'Height|BMI' | head -n 1`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v Irish | head -n 1 | tail -n 1`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep -v Plus | head -n 1`; do
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`; NumSNPs=`zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | head -n 1 | perl -ane 'print scalar(@F);'`; Pheno1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; PhenoSeed1=`echo $i | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`; AncSeed1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[3];'`
			NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt | wc | awk '{ print $1 }'`	
#			NumPaths=2
			echo $Pheno1 $ancestry1 $ancestry2 $ancestry3 $k
			
			LpCnt=1; LpCnt2=1; for (( PathNum=1; PathNum <= $NumPaths; PathNum=PathNum+10 )); do
	       sbatch -t 36:00:00 --mem 16g --account ccmb-condo -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$Pheno1/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.Vs2.GjDrop_wCov_GK_perm1.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$Pheno1/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.Vs2.GjDrop_wCov_GK_perm1.Pathways${PathNum}.slurm.error --comment "$Pheno1 $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
				echo -e "\nR -q -e \"library(\\\"data.table\\\"); library(\\\"doParallel\\\"); library(\\\"Rcpp\\\"); library(\\\"RcppArmadillo\\\"); library(\\\"RcppParallel\\\"); sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Vs2.GjDrop.mtEdits.SingleRun.vs1.wCovs.vs1.cpp\\\"); neg.is.na <- Negate(is.na); neg.is.true <- Negate(isTRUE); \
				Covars <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.sort.ImptHRC.dose.100geno.raw.txt\\\", header=T); Pathways <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\\\", header=F); Pathways.Check <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/c2.all.v6.1.wcp_comps.symbols.${ancestry2}.Regions.c2.${k}.txt\\\", header=F); Y.Check <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.Regions.c2.${k}.Pathways1.txt.gz\\\", header=T); Y.Check.Pheno <- Y.Check\\\$$Pheno1; Y.Check.Pheno.noNAs <- Y.Check.Pheno[neg.is.na(Y.Check.Pheno)]; Y.Seed <- $AncSeed1 + $PhenoSeed1; \
				Pathways.Regions <- list(); cores = detectCores(); InterPath.output <- list(); InterPath.output\\\$Est <- c(); InterPath.output\\\$Eigenvalues <- c(); InterPath.output\\\$PVE <- c(); \ 
				for (i in $PathNum:($PathNum+9)) { set.seed(Y.Seed); Y <- c(); Y.Pheno <- c(); Y.Pheno.noNAs <- c(); \
					if (i > nrow(Pathways)) { Pathways.Regions[[1]] <- 1; Y.Pheno.noNAs <- rep(NA, nrow(InterPath.output\\\$Eigenvalues)); } else if (neg.is.true(Pathways.Check[i,ncol(Pathways.Check)])) { Pathways.Regions[[1]] <- 1; Y.Pheno.noNAs <- rep(NA, length(Y.Check.Pheno.noNAs)); } else { Pathways.Regions[[1]] <- as.numeric(as.character(unlist(strsplit(as.character(Pathways[i,3]), \\\",\\\")))); \
				      Y <- read.table(paste(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.Regions.c2.${k}.Pathways\\\", i, \\\".txt.gz\\\", sep=\\\"\\\"), header=T); Y.Pheno <- Y\\\$$Pheno1; Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; Y.Pheno.noNAs <- sample(Y.Pheno.noNAs);}; \ 
					if (length(Pathways.Regions[[1]]) > 1) { Data3 <- fread(cmd=paste(\\\"zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Regions.c2.${k}.Pathways\\\", as.character(i), \\\".txt.gz\\\", sep=\\\"\\\"), header=T); Data3.cov <- as.matrix(read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.noFix.cov.txt\\\", header=F)); \ 
					Data3.mean <- apply(Data3, 2, mean); Data3.sd <- apply(Data3, 2, sd); Data3 <- t((t(Data3)-Data3.mean)/Data3.sd); \ 
					InterPath.output.temp <- list(); X <- Data3; X.cov <- Data3.cov; rm(Data3); rm(Data3.cov); X.Pheno.noNAs <- X[neg.is.na(Y.Pheno),]; X.cov.Pheno.noNAs <- X.cov[neg.is.na(Y.Pheno),neg.is.na(Y.Pheno)]; Z <- Covars[neg.is.na(Y.Pheno),(ncol(Covars)-9):ncol(Covars)]; \
					G <- 1/nrow(X.Pheno.noNAs) * tcrossprod(as.matrix(X.Pheno.noNAs)); \
					InterPath.output.temp <- InterPath(t(X.Pheno.noNAs),Y.Pheno.noNAs,as.matrix(X.cov.Pheno.noNAs),G,t(as.matrix(Z)),Pathways.Regions,nrow(X.Pheno.noNAs),as.numeric(as.character($NumSNPs)),cores=cores); InterPath.output\\\$Est <- c(InterPath.output\\\$Est, InterPath.output.temp\\\$Est); InterPath.output\\\$Eigenvalues <- cbind(InterPath.output\\\$Eigenvalues, InterPath.output.temp\\\$Eigenvalues); InterPath.output\\\$PVE <- c(InterPath.output\\\$PVE, InterPath.output.temp\\\$PVE); } else { InterPath.output\\\$Est <- c(InterPath.output\\\$Est, NA); InterPath.output\\\$Eigenvalues <- cbind(InterPath.output\\\$Eigenvalues, rep(NA, length(Y.Pheno.noNAs))); InterPath.output\\\$PVE <- c(InterPath.output\\\$PVE, NA);};}; \
			       write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$Pheno1/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK_perm1.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$Pheno1/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK_perm1.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$Pheno1/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK_perm1.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")
			done; sleep 2
		done; 
	done; 
done
#				LpCnt2=$((LpCnt2+1)); if [  $LpCnt2 == 5 ] ; then LpCnt2=1; fi
#				LpCnt=LpCnt+1; if [ $LpCnt == 100 ] ; then LpCnt=1; sleep 30; fi

#cat /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran4000/mturchin20/Analyses/InterPath/*/*/slurm/ukb_chrAll_v2.British.Ran4000.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.*.Vs2.GjDrop_wCov_GK.Pathways4721.slurm.error | vi -

#GjDrop + wCov + GK setup
#				sbatch -t 72:00:00 --mem 10g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.Vs2.GjDrop_wCov_GK.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.Vs2.GjDrop_wCov_GK.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
#				...sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Vs2.GjDrop.mtEdits.SingleRun.vs1.wCovs.vs1.cpp\\\");...
#					Y <- read.table(paste(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.Regions.c2.${k}.Pathways\\\", i, \\\".txt.gz\\\", sep=\\\"\\\"), header=T); Y.Pheno <- Y\\\$$i; Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; }; \ 
#				write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")
##GjDrop + wCov + GG setup
##				sbatch -t 72:00:00 --mem 10g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.Vs2.GjDrop_wCov_GG.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.Vs2.GjDrop_wCov_GG.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
##				...sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Vs2.GjDrop.mtEdits.SingleRun.vs1.wCovs.vs1.GG.cpp\\\");...
##					Y <- read.table(paste(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.Regions.c2.${k}.Pathways\\\", i, \\\".txt.gz\\\", sep=\\\"\\\"), header=T); Y.Pheno <- Y\\\$$i; Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; }; \ 
##				write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GG.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
##				write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GG.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
##				write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GG.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")
#GjDrop + wCov + GK + perm1 setup
#				sbatch -t 72:00:00 --mem 10g -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.Vs2.GjDrop_wCov_GK_perm1.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.Vs2.GjDrop_wCov_GK_perm1.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
#				...sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Vs2.GjDrop.mtEdits.SingleRun.vs1.wCovs.vs1.cpp\\\");...
#				      Y <- read.table(paste(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.Regions.c2.${k}.Pathways\\\", i, \\\".txt.gz\\\", sep=\\\"\\\"), header=T); Y.Pheno <- Y\\\$$Pheno1; Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; Y.Pheno.noNAs <- sample(Y.Pheno.noNAs);}; \ 
#				write.table(InterPath.output\\\$Est, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK_perm1.Paths${PathNum}.Est.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$Eigenvalues, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK_perm1.Paths${PathNum}.Eigenvalues.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE); \
#				write.table(InterPath.output\\\$PVE, \\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK_perm1.Paths${PathNum}.PVE.txt\\\", quote=FALSE, row.name=FALSE, col.name=FALSE);\"")


#Vs2 Error File Hack Check
mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/DataRunChecks

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v Irish`; do
	echo $j
	for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);')`; do
		echo $i
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);')`; do
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
			ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			
#			ls -lrt /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/*error | grep _GK | awk '{ print $5 }' | sort | uniq -c 
			ls -lrt /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ | grep -v -E "Oct 30|Oct 31|Nov  1" | grep error$ | grep _GK | awk '{ print $5 }' | sort | uniq -c 

		done 
	done
done > /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/DataRunChecks/Vs2ErrorFileHackChecks.ErrorFileSizes.vs1.txt
> /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/DataRunChecks/Vs2ErrorFileHackChecks.ErrorFileSizes.vs1.txt
		
##			ls -lrt /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/*error | grep _GK | awk '{ print $5 }' | sort | uniq -c | xargs echo $k

#Vs2 Results Collection
for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep -v -E 'Height|BMI'`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v Irish`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep Plus`; do
			SECONDS=0;
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
			ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt | wc | awk '{ print $1 }'`	
#			NumPaths=20

			echo $i $ancestry1 $ancestry2 $ancestry3 $k
			for (( PathNum=1; PathNum <= $NumPaths; PathNum=PathNum+10 )); do
				R -q -e "library(\"CompQuadForm\"); neg.is.na <- Negate(is.na); \
				Pathways <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\", header=F); \ 
				InterPath.output.Est <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.Paths${PathNum}.Est.txt\", header=F); \
				InterPath.output.Eigenvalues <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.Paths${PathNum}.Eigenvalues.txt\", header=F); \
				InterPath.output.PVE <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.Paths${PathNum}.PVE.txt\", header=F); \
				Results1 <- c(); Counter1 <- 1; for (i in $PathNum:($PathNum+9)) { \
					if (i <= $NumPaths) { if (neg.is.na(InterPath.output.Est[Counter1,1])) { \ 
						Lambda <- sort(InterPath.output.Eigenvalues[,Counter1], decreasing=TRUE); \
						Davies.Output <- davies(InterPath.output.Est[Counter1,1], lambda=Lambda, acc=1e-8); \
						pVal <- 2*min(1-Davies.Output\$Qq, Davies.Output\$Qq); \
						Results1 <- rbind(Results1, c(as.character(Pathways[i,1]), InterPath.output.Est[Counter1,1], InterPath.output.PVE[Counter1,1], pVal)); \ 
						Counter1 = Counter1 + 1; \
					} else { Results1 <- rbind(Results1, c(as.character(Pathways[i,1]), rep(NA, 3))); }; \
				};}; write.table(Results1, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);"
			done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.AllPaths.Results.txt.pre.gz
			duration=$SECONDS; echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
		done;
	done;
done;

#GjDrop + wCov + GK setup
#				InterPath.output.Est <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.Paths${PathNum}.Est.txt\", header=F); \
#				InterPath.output.Eigenvalues <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.Paths${PathNum}.Eigenvalues.txt\", header=F); \
#				InterPath.output.PVE <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.Paths${PathNum}.PVE.txt\", header=F); \
#			done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.AllPaths.Results.txt.pre.gz
#GjDrop + wCov + GG setup
#				InterPath.output.Est <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GG.Paths${PathNum}.Est.txt\", header=F); \
#				InterPath.output.Eigenvalues <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GG.Paths${PathNum}.Eigenvalues.txt\", header=F); \
#				InterPath.output.PVE <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GG.Paths${PathNum}.PVE.txt\", header=F); \
#			done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GG.AllPaths.Results.txt.pre.gz

~for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep -E 'Height|BMI'`; do
~	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
~		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep -v Plus`; do
~			SECONDS=0;
~			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
~			ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
~			
~			R -q -e "library(\"CompQuadForm\"); neg.is.na <- Negate(is.na); \
~			InterPath.AllPaths.output <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GG.AllPaths.Results.txt.pre.gz\", header=F); \ 
~			pVals <- InterPath.AllPaths.output[,ncol(InterPath.AllPaths.output)]; pVals.no0 <- pVals[pVals > 0]; \
~			quants <- qchisq(pVals.no0, df=1, lower.tail=FALSE); \
~			quants.medn <- median(quants); m <- 0.4549; \
~			pVals.no0 <- pchisq((quants/quants.medn) * m, df=1, lower.tail=FALSE); \
~			InterPath.AllPaths.output <- InterPath.AllPaths.output[pVals > 0,];  InterPath.AllPaths.output[,ncol(InterPath.AllPaths.output)] <- pVals.no0; \
~			write.table(InterPath.AllPaths.output, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);" | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GG.AllPaths.Results.Recal.txt.pre.gz
~			duration=$SECONDS; echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
~		done;
~	done;
~done;
~
~#GjDrop + wCov + GK setup
~#			InterPath.AllPaths.output <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.AllPaths.Results.txt.pre.gz\", header=F); \ 
~#			write.table(InterPath.AllPaths.output, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);" | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.AllPaths.Results.Recal.txt.pre.gz
~#GjDrop + wCov + GG setup
~#			InterPath.AllPaths.output <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GG.AllPaths.Results.txt.pre.gz\", header=F); \ 
~#			write.table(InterPath.AllPaths.output, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);" | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GG.AllPaths.Results.Recal.txt.pre.gz

mkdir /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls

R -q -e "library(\"RColorBrewer\"); DataTypes <- c(\"GjDrop_wCov_GK\", \"GjDrop_wCov_GG\"); \ 
	png(\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.HeightBMI.NonSynExonic.Rnd2Vrsns.AllPaths.Results.vs1.png\", height=4000, width=4500, res=300); par(oma=c(1,1,4,14), mar=c(5,5,4,2), mfrow=c(2,2)); \ 
	for (i in DataTypes) { \
		Data1a <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.Vs2.\", i, \".AllPaths.Results.txt.pre.gz\", sep=\"\"), header=F); \	
		Data2a <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.Exonic.Vs2.\", i, \".AllPaths.Results.txt.pre.gz\", sep=\"\"), header=F); \	
		Data3a <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/BMI/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.BMI.NonSyn.Vs2.\", i, \".AllPaths.Results.txt.pre.gz\", sep=\"\"), header=F); \
		Data4a <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/BMI/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.BMI.Exonic.Vs2.\", i, \".AllPaths.Results.txt.pre.gz\", sep=\"\"), header=F); \
		Data1b <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.Vs2.\", i, \".AllPaths.Results.Recal.txt.pre.gz\", sep=\"\"), header=F); \
		Data2b <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.Exonic.Vs2.\", i, \".AllPaths.Results.Recal.txt.pre.gz\", sep=\"\"), header=F); \
		Data3b <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/BMI/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.BMI.NonSyn.Vs2.\", i, \".AllPaths.Results.Recal.txt.pre.gz\", sep=\"\"), header=F); \
		Data4b <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/BMI/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.BMI.Exonic.Vs2.\", i, \".AllPaths.Results.Recal.txt.pre.gz\", sep=\"\"), header=F); \	
		Data1a <- Data1a[Data1a[,4] > 0,]; Data2a <- Data2a[Data2a[,4] > 0,]; Data3a <- Data3a[Data3a[,4] > 0,]; Data4a <- Data4a[Data4a[,4] > 0,]; Data1b <- Data1b[Data1b[,4] > 0,]; Data2b <- Data2b[Data2b[,4] > 0,]; Data3b <- Data3b[Data3b[,4] > 0,]; Data4b <- Data4b[Data4b[,4] > 0,]; \
		xVals1a <- seq(1/nrow(Data1a), 1, by=1/nrow(Data1a)); xVals2a <- seq(1/nrow(Data2a), 1, by=1/nrow(Data2a)); xVals3a <- seq(1/nrow(Data3a), 1, by=1/nrow(Data3a)); xVals4a <- seq(1/nrow(Data4a), 1, by=1/nrow(Data4a)); \ 
		xVals1b <- seq(1/nrow(Data1b), 1, by=1/nrow(Data1b)); xVals2b <- seq(1/nrow(Data2b), 1, by=1/nrow(Data2b)); xVals3b <- seq(1/nrow(Data3b), 1, by=1/nrow(Data3b)); xVals4b <- seq(1/nrow(Data4b), 1, by=1/nrow(Data4b)); \
		xalimMax <- max(c(-log10(xVals1a), -log10(xVals2a), -log10(xVals3a), -log10(xVals4a))); yalimMax <- max(c(-log10(Data1a[,4]), -log10(Data2a[,4]), -log10(Data3a[,4]), -log10(Data4a[,4]))); \
		xblimMax <- max(c(-log10(xVals1b), -log10(xVals2b), -log10(xVals3b), -log10(xVals4b))); yblimMax <- max(c(-log10(Data1b[,4]), -log10(Data2b[,4]), -log10(Data3b[,4]), -log10(Data4b[,4]))); \
		\
		Title1 <- c(); if (i == \"yTyAdj\") { Title1 <- \"Orig\"; } else if (i == \"GjDrop_wCov_GK\") { Title1 <- \"G*K + PCProject\"; } else if (i == \"GjDrop_wCov_GG\") { Title1 <- \"G*G + PCProject\"; } else {  Title1 <- \"Error1a\"; }; \ 
		\
		plot(-log10(xVals1a[order(xVals1a, decreasing=TRUE)]), -log10(Data1a[order(Data1a[,4], decreasing=TRUE),4]), main=Title1, xlab=\"-log10(Expected pValues)\", ylab=\"-log10(Observed pValues)\", xlim=c(0,xalimMax), ylim=c(0,yalimMax), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[7], cex=1.5, cex.main=1.5, cex.axis=1.5, cex.lab=1.5); \
		points(-log10(xVals2a[order(xVals2a, decreasing=TRUE)]), -log10(Data2a[order(Data2a[,4], decreasing=TRUE),4]), type=\"b\", pch=17, col=brewer.pal(12, \"Paired\")[7], cex=1.5); \ 
		points(-log10(xVals3a[order(xVals3a, decreasing=TRUE)]), -log10(Data3a[order(Data3a[,4], decreasing=TRUE),4]), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[3], cex=1.5); \
		points(-log10(xVals4a[order(xVals4a, decreasing=TRUE)]), -log10(Data4a[order(Data4a[,4], decreasing=TRUE),4]), type=\"b\", pch=17, col=brewer.pal(12, \"Paired\")[3], cex=1.5); \
		abline(0,1, lwd=2, col=\"BLACK\"); \
		plot(-log10(xVals1b[order(xVals1b, decreasing=TRUE)]), -log10(Data1b[order(Data1b[,4], decreasing=TRUE),4]), main=paste(Title1, \".Recal\", sep=\"\"), xlab=\"-log10(Expected pValues)\", ylab=\"-log10(Observed pValues)\", xlim=c(0,xblimMax), ylim=c(0,yblimMax), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[7], cex=1.5, cex.main=1.5, cex.axis=1.5, cex.lab=1.5); \ 
		points(-log10(xVals2b[order(xVals2b, decreasing=TRUE)]), -log10(Data2b[order(Data2b[,4], decreasing=TRUE),4]), type=\"b\", pch=17, col=brewer.pal(12, \"Paired\")[7], cex=1.5); \ 
		points(-log10(xVals3b[order(xVals3b, decreasing=TRUE)]), -log10(Data3b[order(Data3b[,4], decreasing=TRUE),4]), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[3], cex=1.5); \
		points(-log10(xVals4b[order(xVals4b, decreasing=TRUE)]), -log10(Data4b[order(Data4b[,4], decreasing=TRUE),4]), type=\"b\", pch=17, col=brewer.pal(12, \"Paired\")[3], cex=1.5); \ 
		abline(0,1, lwd=2, col=\"BLACK\"); \
		\
	}; mtext(\"African\", line=-.75, outer=TRUE, cex=2.5); par(fig = c(0, 1, 0, 1), mfrow=c(1,1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot(0, 0, type = \"n\", bty = \"n\", xaxt = \"n\", yaxt = \"n\"); legend(\"topright\", c(\"Height\", \"BMI\", \"NonSyn\", \"Exonic\"), pch=c(NA, NA, 16, 17), text.col=c(brewer.pal(12, \"Paired\")[7], brewer.pal(12, \"Paired\")[3], \"BLACK\", \"BLACK\"), xpd=TRUE, inset=c(.0425,.1), bg=\"transparent\", cex=1.35, y.intersp=2); dev.off(); \
"

#inset values for # of rows:
#2 rows = c(
#5 rows = c(.027,.0385)

#From MacBook Air
#mkdir /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/
scp -p  mturchin@ssh.ccv.brown.edu:/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.HeightBMI.NonSynExonic.Rnd2Vrsns.AllPaths.Results.vs1.png /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/.

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v Irish`; do
	echo $j
	for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);')`; do
		echo $i
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);')`; do
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
			ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`

#			zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.AllPaths.Results.txt.pre.gz | awk '{ if (($4 != 0) && ($4 <= 5e-5)) { print $0 } }' | wc | awk '{ print $1 }'	
			NoZero=`zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.AllPaths.Results.txt.pre.gz | awk '{ if ($4 != 0) { print $0 } }' | wc | awk '{ print $1 }'`
			NoNA=`zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.AllPaths.Results.txt.pre.gz | awk '{ if ($4 != "NA") { print $0 } }' | wc | awk '{ print $1 }'`
			NoZeroAndNA=`zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_GK.AllPaths.Results.txt.pre.gz | awk '{ if (($4 != 0) && ($4 != "NA")) { print $0 } }' | wc | awk '{ print $1 }'`
			echo $NoZero $NoNA $NoZeroAndNA

		done 
	done
done > /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/Rnd2AdditiveMdls.Vs2.GK.totalPaths.counts.vs1.txt
/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/Rnd2AdditiveMdls.Vs2.GK.sigPaths.counts.vs1.txt

R -q -e "library(\"RColorBrewer\"); UKBioBankPops <- c(\"African;African\",\"British;British.Ran4000\",\"Caribbean;Caribbean\",\"Chinese;Chinese\",\"Indian;Indian\",\"Pakistani;Pakistani\"); DataTypes <- c(\"GjDrop_wCov_GK\"); \
	png(\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/ukb_chrAll_v2.AllPops.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.AllPhenos.AllTypes.Rnd2Vrsns.AllPaths.Results.vs2.png\", height=12000, width=8500, res=300); par(oma=c(1,1,4,14), mar=c(5,5,4,2), mfrow=c(6,4)); \ 
	for (j in UKBioBankPops) { ancestry1 = strsplit(j, \";\")[[1]][1]; ancestry2 = strsplit(j, \";\")[[1]][2]; \	
		for (i in DataTypes) { \
			for (k in c(\"Height\", \"BMI\", \"Waist\", \"Hip\")) { \
				Data1 <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/\", ancestry1, \"/\", ancestry2, \"/mturchin20/Analyses/InterPath/\", k, \"/ukb_chrAll_v2.\", ancestry2, \".QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.\", k, \".NonSyn.Vs2.\", i, \".AllPaths.Results.txt.pre.gz\", sep=\"\"), header=F); \	
				Data2 <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/\", ancestry1, \"/\", ancestry2, \"/mturchin20/Analyses/InterPath/\", k, \"/ukb_chrAll_v2.\", ancestry2, \".QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.\", k, \".Exonic.Vs2.\", i, \".AllPaths.Results.txt.pre.gz\", sep=\"\"), header=F); \	
				Data3 <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/\", ancestry1, \"/\", ancestry2, \"/mturchin20/Analyses/InterPath/\", k, \"/ukb_chrAll_v2.\", ancestry2, \".QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.\", k, \".ExonicPlus.Vs2.\", i, \".AllPaths.Results.txt.pre.gz\", sep=\"\"), header=F); \	
				Data4 <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/\", ancestry1, \"/\", ancestry2, \"/mturchin20/Analyses/InterPath/\", k, \"/ukb_chrAll_v2.\", ancestry2, \".QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.\", k, \".ExonicPlus20kb.Vs2.\", i, \".AllPaths.Results.txt.pre.gz\", sep=\"\"), header=F); \	
				Data1 <- Data1[!is.na(Data1[,4]) & Data1[,4] != 0,]; Data2 <- Data2[!is.na(Data2[,4]) & Data2[,4] != 0,]; Data3 <- Data3[!is.na(Data3[,4]) & Data3[,4] != 0,]; Data4 <- Data4[!is.na(Data4[,4]) & Data4[,4] != 0,]; \  
				xVals1 <- seq(1/nrow(Data1), 1, by=1/nrow(Data1)); xVals2 <- seq(1/nrow(Data2), 1, by=1/nrow(Data2)); xVals3 <- seq(1/nrow(Data3), 1, by=1/nrow(Data3)); xVals4 <- seq(1/nrow(Data4), 1, by=1/nrow(Data4)); \ 
				xlimMax <- max(c(-log10(xVals1), -log10(xVals2), -log10(xVals3), -log10(xVals4))); ylimMax <- max(c(-log10(Data1[,4]), -log10(Data2[,4]), -log10(Data3[,4]), -log10(Data4[,4]))); \
				plot(-log10(xVals1[order(xVals1, decreasing=TRUE)]), -log10(Data1[order(Data1[,4], decreasing=TRUE),4]), main=paste(ancestry2, \": \", k, sep=\"\"), xlab=\"-log10(Expected pValues)\", ylab=\"-log10(Observed pValues)\", xlim=c(0,max(-log10(xVals4[order(xVals4, decreasing=TRUE)]))), ylim=c(0,ylimMax), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[7], cex=1.5, cex.main=1.5, cex.axis=1.5, cex.lab=1.5); \ 
				points(-log10(xVals2[order(xVals2, decreasing=TRUE)]), -log10(Data2[order(Data2[,4], decreasing=TRUE),4]), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[3], cex=1.5); \
				points(-log10(xVals3[order(xVals3, decreasing=TRUE)]), -log10(Data3[order(Data3[,4], decreasing=TRUE),4]), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[9], cex=1.5); \
				points(-log10(xVals4[order(xVals4, decreasing=TRUE)]), -log10(Data4[order(Data4[,4], decreasing=TRUE),4]), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[1],  cex=1.5); \ 
				abline(0,1, lwd=2, col=\"BLACK\"); \
				legend(\"topleft\", c(\"NonSyn\", \"Exonic\", \"ExonicPlus\", \"ExonicPlus20kb\"), pch=c(16,16,16,16), col=c(brewer.pal(12, \"Paired\")[7], brewer.pal(12, \"Paired\")[3], brewer.pal(12, \"Paired\")[9], brewer.pal(12, \"Paired\")[1]), bg=\"transparent\", cex=1.5); \ 
			}; \
		}; \
	}; \
"
	
#				print(quantile(-log10(Data1[,4]))); \
#				print(quantile(-log10(Data2[,4]))); \
#				print(quantile(-log10(Data3[,4]))); \
#				print(quantile(-log10(Data4[,4]))); \
#				print(ylimMax); \
#	legend(\"topright\", c(\"NonSyn\", \"Exonic\", \"ExonicPlus\", \"ExonicPlus20kb\"), pch=c(16, 16, 16, 16), text.col=c(brewer.pal(12, \"Paired\")[7], brewer.pal(12, \"Paired\")[3], brewer.pal(12, \"Paired\")[9], brewer.pal(12, \"Paired\")[1]), xpd=TRUE, inset=c(.0425,.1), bg=\"transparent\", cex=1.35, y.intersp=2); dev.off(); \

#From MacBook Air
scp -p  mturchin@ssh.ccv.brown.edu:/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/ukb_chrAll_v2.AllPops.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.AllPhenos.AllTypes.Rnd2Vrsns.AllPaths.Results.vs2.png /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/.
#/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.HeightBMI.NonSynExonic.Rnd2Vrsns.AllPaths.Results.vs2.png 

ln -s /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/ukb_chrAll_v2.AllPops.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.AllPhenos.AllTypes.Rnd2Vrsns.AllPaths.Results.vs2.png /users/mturchin/LabMisc/RamachandranLab/InterPath/images/20181102.InterPath.Vs2.GK.Analysis.MainResults.vs1.png

#For Poster
for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);')`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
			ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			
			echo $i $ancestry1 $ancestry2 $ancestry3 $k

			R -q -e "library(\"RColorBrewer\"); Data1 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.NonSyn.Vs2.GjDrop_wCov_GK.AllPaths.Results.txt.pre.gz\", header=F); \ 
			Data2 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.Exonic.Vs2.GjDrop_wCov_GK.AllPaths.Results.txt.pre.gz\", header=F); \
			Data3 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.ExonicPlus.Vs2.GjDrop_wCov_GK.AllPaths.Results.txt.pre.gz\", header=F); \
			Data4 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.ExonicPlus20kb.Vs2.GjDrop_wCov_GK.AllPaths.Results.txt.pre.gz\", header=F); \ 
			Data1 <- Data1[!is.na(Data1[,4]) & Data1[,4] != 0,]; Data2 <- Data2[!is.na(Data2[,4]) & Data2[,4] != 0,]; Data3 <- Data3[!is.na(Data3[,4]) & Data3[,4] != 0,]; Data4 <- Data4[!is.na(Data4[,4]) & Data4[,4] != 0,]; \ 
			png(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.AllStrats.Vs2.GjDrop_wCov_GK.AllPaths.Results.pre.QQplots.vs1.png\", width=2000, height=2000, res=300); par(mar=c(5,5,4,2)); \ 
			xVals1 <- seq(1/nrow(Data1), 1, by=1/nrow(Data1)); xVals2 <- seq(1/nrow(Data2), 1, by=1/nrow(Data2)); xVals3 <- seq(1/nrow(Data3), 1, by=1/nrow(Data3)); xVals4 <- seq(1/nrow(Data4), 1, by=1/nrow(Data4)); \ 
			xlimMax <- max(c(-log10(xVals1), -log10(xVals2), -log10(xVals3), -log10(xVals4))); ylimMax <- max(c(-log10(Data1[,4]), -log10(Data2[,4]), -log10(Data3[,4]), -log10(Data4[,4]))); \
			plot(-log10(xVals1[order(xVals1, decreasing=TRUE)]), -log10(Data1[order(Data1[,4], decreasing=TRUE),4]), main=\"InterPath Prelim Results: $ancestry2 ($i)\", xlab=\"-log10(Expected pValues)\", ylab=\"-log10(Observed pValues)\", xlim=c(0,max(-log10(xVals4[order(xVals4, decreasing=TRUE)]))), ylim=c(0,ylimMax), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[7], cex=1.5, cex.main=1.5, cex.axis=1.5, cex.lab=1.5); \ 
			points(-log10(xVals2[order(xVals2, decreasing=TRUE)]), -log10(Data2[order(Data2[,4], decreasing=TRUE),4]), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[3], cex=1.5); \
			points(-log10(xVals3[order(xVals3, decreasing=TRUE)]), -log10(Data3[order(Data3[,4], decreasing=TRUE),4]), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[9], cex=1.5); \
			points(-log10(xVals4[order(xVals4, decreasing=TRUE)]), -log10(Data4[order(Data4[,4], decreasing=TRUE),4]), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[1],  cex=1.5); \
			legend(\"topleft\", c(\"NonSyn\", \"Exonic\", \"ExonicPlus\", \"ExonicPlus20kb\"), pch=c(16,16,16,16), col=c(brewer.pal(12, \"Paired\")[7], brewer.pal(12, \"Paired\")[3], brewer.pal(12, \"Paired\")[9], brewer.pal(12, \"Paired\")[1]), bg=\"transparent\"); \ 
			abline(0,1, col=\"BLACK\"); \ 
			dev.off();"
	done;
done;
for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -E 'British|Indian'`; do
		ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
		ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
		
		echo $i $ancestry1 $ancestry2 $ancestry3 $k

		R -q -e "library(\"RColorBrewer\"); Data1 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.Vs2.GjDrop_wCov_GK.AllPaths.Results.txt.pre.gz\", header=F); \
		Data2 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/Height/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.ExonicPlus.Vs2.GjDrop_wCov_GK.AllPaths.Results.txt.pre.gz\", header=F); \
		Data3 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/BMI/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.BMI.NonSyn.Vs2.GjDrop_wCov_GK.AllPaths.Results.txt.pre.gz\", header=F); \
		Data4 <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/BMI/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.BMI.ExonicPlus.Vs2.GjDrop_wCov_GK.AllPaths.Results.txt.pre.gz\", header=F); \
		Data1 <- Data1[!is.na(Data1[,4]) & Data1[,4] != 0,]; Data2 <- Data2[!is.na(Data2[,4]) & Data2[,4] != 0,]; Data3 <- Data3[!is.na(Data3[,4]) & Data3[,4] != 0,]; Data4 <- Data4[!is.na(Data4[,4]) & Data4[,4] != 0,]; \ 
		png(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.HeightBMI.AllStrats.Vs2.GjDrop_wCov_GK.AllPaths.Results.pre.QQplots.vs1.png\", width=2000, height=2000, res=300);  par(mar=c(5,5,4,2)); \
		xVals1 <- seq(1/nrow(Data1), 1, by=1/nrow(Data1)); xVals2 <- seq(1/nrow(Data2), 1, by=1/nrow(Data2)); xVals3 <- seq(1/nrow(Data3), 1, by=1/nrow(Data3)); xVals4 <- seq(1/nrow(Data4), 1, by=1/nrow(Data4)); \
		xlimMax <- max(c(-log10(xVals1), -log10(xVals2), -log10(xVals3), -log10(xVals4))); ylimMax <- max(c(-log10(Data1[,4]), -log10(Data2[,4]), -log10(Data3[,4]), -log10(Data4[,4]))); \
		plot(-log10(xVals1[order(xVals1, decreasing=TRUE)]), -log10(Data1[order(Data1[,4], decreasing=TRUE),4]), main=\"InterPath Prelim Results: $ancestry1\", xlab=\"-log10(Expected pValues)\", ylab=\"-log10(Observed pValues)\", xlim=c(0,xlimMax), ylim=c(0,ylimMax), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[7], cex=1.5, cex.main=1.5, cex.axis=1.5, cex.lab=1.5); \
		points(-log10(xVals2[order(xVals2, decreasing=TRUE)]), -log10(Data2[order(Data2[,4], decreasing=TRUE),4]), type=\"b\", pch=16, col=brewer.pal(12, \"Paired\")[9], cex=1.5); \ 
		points(-log10(xVals3[order(xVals3, decreasing=TRUE)]), -log10(Data3[order(Data3[,4], decreasing=TRUE),4]), type=\"b\", pch=17, col=brewer.pal(12, \"Paired\")[7], cex=1.5); \ 
		points(-log10(xVals4[order(xVals4, decreasing=TRUE)]), -log10(Data4[order(Data4[,4], decreasing=TRUE),4]), type=\"b\", pch=17, col=brewer.pal(12, \"Paired\")[9], cex=1.5); \ 
		legend(\"topleft\", c(\"Height: NonSyn\", \"Height: ExonicPlus\", \"BMI: NonSyn\", \"BMI: ExonicPlus\"), pch=c(16,16,17,17), col=c(brewer.pal(12, \"Paired\")[7], brewer.pal(12, \"Paired\")[9], brewer.pal(12, \"Paired\")[7], brewer.pal(12, \"Paired\")[9]), bg=\"transparent\"); \
		abline(0,1, col=\"BLACK\"); \
		dev.off();"
done;

#On MacBook Pro
#scp -p  mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/*/*/mturchin20/Analyses/InterPath/*GjDrop_wCov_GK.AllPaths.Results.pre.QQplots.vs1.png /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/.

R -q -e "library(\"RColorBrewer\"); UKBioBankPops <- c(\"African;African\",\"British;British.Ran4000\",\"Caribbean;Caribbean\",\"Chinese;Chinese\",\"Indian;Indian\",\"Pakistani;Pakistani\"); DataTypes <- c(\"GjDrop_wCov_GK\"); \
	for (j in UKBioBankPops) { ancestry1 = strsplit(j, \";\")[[1]][1]; ancestry2 = strsplit(j, \";\")[[1]][2]; \	
		png(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/\", ancestry1, \"/\", ancestry2, \"/mturchin20/Analyses/InterPath/ukb_chrAll_v2.\", ancestry2, \".QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.AllPhenos.GjDrop_wCov_GK.Rnd2Vrsns.AllPaths.pValHists.vs1.png\", sep=\"\"), height=8000, width=8000, res=300); par(oma=c(1,12,8,1), mar=c(5,5,4,2), mfcol=c(4,4)); \
		for (i in DataTypes) { \
			for (k in c(\"Height\", \"BMI\", \"Waist\", \"Hip\")) { \
				Data1 <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/\", ancestry1, \"/\", ancestry2, \"/mturchin20/Analyses/InterPath/\", k, \"/ukb_chrAll_v2.\", ancestry2, \".QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.\", k, \".NonSyn.Vs2.\", i, \".AllPaths.Results.txt.pre.gz\", sep=\"\"), header=F); \	
				Data2 <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/\", ancestry1, \"/\", ancestry2, \"/mturchin20/Analyses/InterPath/\", k, \"/ukb_chrAll_v2.\", ancestry2, \".QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.\", k, \".Exonic.Vs2.\", i, \".AllPaths.Results.txt.pre.gz\", sep=\"\"), header=F); \	
				Data3 <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/\", ancestry1, \"/\", ancestry2, \"/mturchin20/Analyses/InterPath/\", k, \"/ukb_chrAll_v2.\", ancestry2, \".QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.\", k, \".ExonicPlus.Vs2.\", i, \".AllPaths.Results.txt.pre.gz\", sep=\"\"), header=F); \	
				Data4 <- read.table(paste(\"/users/mturchin/data/ukbiobank_jun17/subsets/\", ancestry1, \"/\", ancestry2, \"/mturchin20/Analyses/InterPath/\", k, \"/ukb_chrAll_v2.\", ancestry2, \".QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.\", k, \".ExonicPlus20kb.Vs2.\", i, \".AllPaths.Results.txt.pre.gz\", sep=\"\"), header=F); \	
				Data1 <- Data1[!is.na(Data1[,4]) & Data1[,4] != 0,]; Data2 <- Data2[!is.na(Data2[,4]) & Data2[,4] != 0,]; Data3 <- Data3[!is.na(Data3[,4]) & Data3[,4] != 0,]; Data4 <- Data4[!is.na(Data4[,4]) & Data4[,4] != 0,]; \ 
				Data1 <- Data1[!is.na(Data1[,4]),]; Data1[Data1[,4] == 0,4] <- 3e-8; Data2 <- Data2[!is.na(Data2[,4]),]; Data2[Data2[,4] == 0,4] <- 3e-8; Data3 <- Data3[!is.na(Data3[,4]),]; Data3[Data3[,4] == 0,4] <- 3e-8; Data4 <- Data4[!is.na(Data4[,4]),]; Data4[Data4[,4] == 0,4] <- 3e-8; \
				hist(Data1[,4], main=\"\", xlab=\"pValues\", breaks=25, cex=1.5, cex.main=1.5, cex.axis=1.5, cex.lab=1.5); \
				hist(Data2[,4], main=\"\", xlab=\"pValues\", breaks=25, cex=1.5, cex.main=1.5, cex.axis=1.5, cex.lab=1.5); \
				hist(Data3[,4], main=\"\", xlab=\"pValues\", breaks=25, cex=1.5, cex.main=1.5, cex.axis=1.5, cex.lab=1.5); \
				hist(Data4[,4], main=\"\", xlab=\"pValues\", breaks=25, cex=1.5, cex.main=1.5, cex.axis=1.5, cex.lab=1.5); \
			}; \
		}; mtext(\"Height\", side=3, outer=TRUE, at=.125, cex=3); mtext(\"BMI\", side=3, outer=TRUE, at=.385, cex=3); mtext(\"Waist\", side=3, outer=TRUE, at=.63, cex=3); mtext(\"Hip\", side=3, outer=TRUE, at=.885, cex=3); \
		mtext(\"NonSyn\", side=2, line=4, outer=TRUE, at=.13, cex=3); mtext(\"Exonic\", side=2, line=4, outer=TRUE, at=.38, cex=3); mtext(\"ExonicPlus\", side=2, line=4, outer=TRUE, at=.63, cex=3); mtext(\"ExonicPlus20kb\", side=2, line=4, outer=TRUE, at=.88, cex=3); dev.off(); \
	}; \
"

#				print(c(min(Data1[,4]), min(Data2[,4]), min(Data3[,4]), min(Data4[,4]))); \
	
#On MacBook Pro
#scp -p  mturchin@ssh.ccv.brown.edu:/users/mturchin/data/ukbiobank_jun17/subsets/*/*/mturchin20/Analyses/InterPath/*AllPhenos.*.Rnd2Vrsns.AllPaths.pValHists.vs1.png /Users/mturchin20/Documents/Work/LabMisc/RamachandranLab/InterPath/Vs1/Analyses/Rnd2AdditiveMdls/.

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v Irish`; do
	ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
#	ln -s /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.AllPhenos.GjDrop_wCov_GK.Rnd2Vrsns.AllPaths.pValHists.vs1.png /users/mturchin/LabMisc/RamachandranLab/InterPath/images/pValHists/20181108.InterPath.Vs2.GK.Analysis.pValHists.$ancestry2.vs1.png 
	cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.fam | wc
done








#NOTE -- once again logged into http://software.broadinstitute.org/gsea/msigdb/collections.jsp, manually downloaded files to MacBook Pro, and then scp'ed onto Oscar server here
#cd /users/mturchin/data/mturchin/Broad/MSigDB
#From MackBook Pro
#mkdir /Users/mturchin20/Documents/Work/LabMisc/Data/MSigDB/
#scp -p /Users/mturchin20/Documents/Work/LabMisc/Data/MSigDB/c2.c*gmt mturchin@ssh.ccv.brown.edu:/users/mturchin/data/mturchin/Broad/MSigDB/.

join <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.v6.2.symbols.gmt | awk '{ print $1 }' | sort) <(cat <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.biocarta.v6.2.symbols.gmt | awk '{ print $1 "\tBIOCARTA" }') <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.kegg.v6.2.symbols.gmt | awk '{ print $1 "\tKEGG" }') <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.reactome.v6.2.symbols.gmt | awk '{ print $1 "\tREACTOME" }') | sort -k 1,1) | awk '{ print $2 }' | sort | uniq -c

cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | awk '{ print $1 "\t" $2 }' > /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.First2Cols.gmt

join -a 1 -a 2 -e "FALSE" -o 0 1.1 2.1 <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | awk '{ print $1 }' | sort) <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.v6.2.symbols.gmt | awk '{ print $1 }' | sort) | awk '{ print $1 "\t" $3 }' | perl -lane 'if ($F[1] ne "FALSE") { $F[1] = "TRUE"; } print join("\t", @F);' | R -q -e "Data1 <- read.table(file('stdin'), header=F); Data2 <- read.table(\"/users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.First2Cols.gmt\", header=F); colnames(Data1) <- c(\"PATH\", \"CP\"); Data2 <- cbind(Data2, seq(1,nrow(Data2))); colnames(Data2) <- c(\"PATH\", \"CP\", \"ORDER\"); Data3 <- merge(Data2, Data1, by=\"PATH\"); Data3 <- Data3[order(Data3[3], decreasing=FALSE),]; write.table(Data3, \"/users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.wcp_comps.symbols.gmt\", quote=FALSE, row.name=FALSE, col.name=FALSE);"

cmp <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | awk '{ print $1 }') <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.wcp_comps.symbols.gmt | awk '{ print $1 }') | wc 

for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);')`; do
	ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`; ancestry3=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[2];'`; echo $pheno1 $ancestry1 $ancestry2 $ancestry3;

	join <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.wcp_comps.symbols.gmt | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.NonSyn.txt | awk '{ print $1 }' | sort) | sort -g -k 3,3 > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/c2.all.v6.1.wcp_comps.symbols.${ancestry2}.Regions.c2.NonSyn.txt
	join <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.wcp_comps.symbols.gmt | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.Exonic.txt | awk '{ print $1 }' | sort) | sort -g -k 3,3 > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/c2.all.v6.1.wcp_comps.symbols.${ancestry2}.Regions.c2.Exonic.txt
	join <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.wcp_comps.symbols.gmt | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.ExonicPlus.txt | awk '{ print $1 }' | sort) | sort -g -k 3,3 > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/c2.all.v6.1.wcp_comps.symbols.${ancestry2}.Regions.c2.ExonicPlus.txt
	join <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.wcp_comps.symbols.gmt | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.ExonicPlus20kb.txt | awk '{ print $1 }' | sort) | sort -g -k 3,3 > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/c2.all.v6.1.wcp_comps.symbols.${ancestry2}.Regions.c2.ExonicPlus20kb.txt

done







737649               63      batch  mturchin ccmb-condo 2018-11-02T12:57:17 2018-11-03T02:54:11   04:45:28          8     FAILED      1:0                 BMI Caribbean Caribbean Exonic 171
737649.batch      batch                      ccmb-condo 2018-11-03T02:18:30 2018-11-03T02:54:11   04:45:28          8     FAILED      1:0
741120               63      batch  mturchin ccmb-condo 2018-11-03T01:12:45 2018-11-03T06:34:50   00:31:12          8     FAILED      1:0         Height Caribbean Caribbean ExonicPlus 1291
741120.batch      batch                      ccmb-condo 2018-11-03T06:30:56 2018-11-03T06:34:50   00:31:12          8     FAILED      1:0
741994               63      batch  mturchin ccmb-condo 2018-11-03T01:13:03 2018-11-03T11:12:14   03:45:52          8     FAILED      1:0             BMI Caribbean Caribbean ExonicPlus 551
741994.batch      batch                      ccmb-condo 2018-11-03T10:44:00 2018-11-03T11:12:14   03:45:52          8     FAILED      1:0
742448               63      batch  mturchin ccmb-condo 2018-11-03T01:13:12 2018-11-03T13:27:28   01:19:12          8     FAILED      1:0         BMI Caribbean Caribbean ExonicPlus20kb 341
742448.batch      batch                      ccmb-condo 2018-11-03T13:17:34 2018-11-03T13:27:28   01:19:12          8     FAILED      1:0
742493               63      batch  mturchin ccmb-condo 2018-11-03T01:13:13 2018-11-03T14:16:53   01:32:48          8     FAILED      1:0         BMI Caribbean Caribbean ExonicPlus20kb 791 
742493.batch      batch                      ccmb-condo 2018-11-03T14:05:17 2018-11-03T14:16:53   01:32:48          8     FAILED      1:0                                                    
748446               63      batch  mturchin ccmb-condo 2018-11-03T15:45:45 2018-11-03T16:43:33   00:13:44          8     FAILED      1:0              Height Chinese Chinese ExonicPlus 691                                                  748446.batch      batch                      ccmb-condo 2018-11-03T16:41:50 2018-11-03T16:43:33   00:13:44          8     FAILED      1:0
748922               63      batch  mturchin ccmb-condo 2018-11-03T15:45:54 2018-11-03T16:59:13   00:13:44          8     FAILED      1:0          Height Chinese Chinese ExonicPlus20kb 711
748922.batch      batch                      ccmb-condo 2018-11-03T16:57:30 2018-11-03T16:59:13   00:13:44          8     FAILED      1:0                                                                                                     749342               63      batch  mturchin ccmb-condo 2018-11-03T15:46:02 2018-11-03T17:11:27   00:25:04          8     FAILED      1:0                 BMI Chinese Chinese ExonicPlus 171
749342.batch      batch                      ccmb-condo 2018-11-03T17:08:19 2018-11-03T17:11:27   00:25:04          8     FAILED      1:0                                                                                                     749367               63      batch  mturchin ccmb-condo 2018-11-03T15:46:02 2018-11-03T17:12:43   00:11:12          8     FAILED      1:0                 BMI Chinese Chinese ExonicPlus 421
749367.batch      batch                      ccmb-condo 2018-11-03T17:11:19 2018-11-03T17:12:43   00:11:12          8     FAILED      1:0
749434               63      batch  mturchin ccmb-condo 2018-11-03T15:46:03 2018-11-03T17:18:27   00:12:40          8     FAILED      1:0                BMI Chinese Chinese ExonicPlus 1091
749434.batch      batch                      ccmb-condo 2018-11-03T17:16:52 2018-11-03T17:18:27   00:12:40          8     FAILED      1:0                                                                                                     749445               63      batch  mturchin ccmb-condo 2018-11-03T15:46:03 2018-11-03T17:20:07   00:18:40          8     FAILED      1:0                BMI Chinese Chinese ExonicPlus 1201                                                  749445.batch      batch                      ccmb-condo 2018-11-03T17:17:47 2018-11-03T17:20:07   00:18:40          8     FAILED      1:0                                                                                                     749801               63      batch  mturchin ccmb-condo 2018-11-03T15:46:11 2018-11-03T17:25:57   00:17:28          8     FAILED      1:0              BMI Chinese Chinese ExonicPlus20kb 11                                                  749801.batch      batch                      ccmb-condo 2018-11-03T17:23:46 2018-11-03T17:25:57   00:17:28          8     FAILED      1:0                                                                                                     
778985               63      batch  mturchin ccmb-condo 2018-11-07T11:36:35 2018-11-10T09:19:29 1-02:14:56          8     FAILED      1:0               Height Indian Indian ExonicPlus 1071 
778985.batch      batch                      ccmb-condo 2018-11-10T06:02:37 2018-11-10T09:19:29 1-02:14:56          8     FAILED      1:0                                                    


#gene*pathway (gene*G)
#cp -p /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Vs2.GjDrop.mtEdits.SingleRun.vs1.wCovs.vs1.cpp /users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Vs2.geneG.mtEdits.SingleRun.vs1.wCovs.vs1.cpp
module load R/3.4.3_mkl; sleep 25200; for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep -E 'Height|BMI'`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v Irish | head -n 6 | tail -n 1`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep -v Plus`; do
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`; NumSNPs=`zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | head -n 1 | perl -ane 'print scalar(@F);'` 
			NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt | wc | awk '{ print $1 }'`	
#			NumPaths=2
			echo $i $ancestry1 $ancestry2 $ancestry3 $k
		
			LpCnt=1; LpCnt2=1; for (( PathNum=1; PathNum <= $NumPaths; PathNum=PathNum+10 )); do
			       sbatch -t 36:00:00 -n 8 -N 1-1 --mem 81g --account ccmb-condo -o /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.Vs2.GjDrop_wCov_geneG.Pathways${PathNum}.slurm.output -e /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/slurm/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.c2.Exonic.${i}.${k}.Vs2.GjDrop_wCov_geneG.Pathways${PathNum}.slurm.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
				echo -e "R -q -e \"library(\\\"data.table\\\"); library(\\\"doParallel\\\"); library(\\\"Rcpp\\\"); library(\\\"RcppArmadillo\\\"); library(\\\"RcppParallel\\\"); sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Vs1/InterPath.Vs2.geneG.mtEdits.SingleRun.vs1.wCovs.vs1.cpp\\\"); neg.is.na <- Negate(is.na); neg.is.true <- Negate(isTRUE); \
				Covars <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.sort.ImptHRC.dose.100geno.raw.txt\\\", header=T); Pathways <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.geneVSpath.${k}.txt\\\", header=F); Pathways.Check <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/c2.all.v6.1.wcp_comps.symbols.${ancestry2}.Regions.c2.${k}.txt\\\", header=F); Pathways.Regions <- list(); cores = detectCores(); \ 
				for (i in $PathNum:($PathNum+9)) { Y <- c(); Y.Pheno <- c(); Y.Pheno.noNAs <- c(); InterPath.output <- list(); SkipFlag1 <- 0; \
					if (i <= nrow(Pathways)) { if (Pathways.Check[i,ncol(Pathways.Check)]) { Pathways.Regions <- lapply(strsplit(unlist(strsplit(as.character(Pathways[i,3]), \\\";\\\")),\\\",\\\"), function(x) { return(as.numeric(as.character(x))); }); if (length(Pathways.Regions) > 1) { \ 
						Y <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt\\\", header=T); Y.Pheno <- Y\\\$$i; Y.Pheno.noNAs <- Y.Pheno[neg.is.na(Y.Pheno)]; \ 
						Data3 <- fread(cmd=paste(\\\"zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Regions.c2.${k}.Pathways\\\", as.character(i), \\\".txt.gz\\\", sep=\\\"\\\"), header=T); \
						Data3.cov <- as.matrix(fread(cmd=paste(\\\"zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/pathways/$k/cov/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.Regions.c2.${k}.Pathways\\\", as.character(i), \\\".cov.txt.gz\\\", sep=\\\"\\\"), header=F)); \
						Data3.mean <- apply(Data3, 2, mean); Data3.sd <- apply(Data3, 2, sd); Data3 <- t((t(Data3)-Data3.mean)/Data3.sd); \ 
						X <- Data3; X.cov <- Data3.cov; rm(Data3); rm(Data3.cov); X.Pheno.noNAs <- X[neg.is.na(Y.Pheno),]; X.cov.Pheno.noNAs <- X.cov[neg.is.na(Y.Pheno),neg.is.na(Y.Pheno)]; Z <- Covars[neg.is.na(Y.Pheno),(ncol(Covars)-9):ncol(Covars)]; \
						Y.Pheno.GjResids <- c(); for (j in 1:length(Pathways.Regions)) { Y.Pheno.GjResids <- cbind(Y.Pheno.GjResids, residuals(lm(Y.Pheno ~ X[,Pathways.Regions[[j]]], na.action=na.exclude))); }; Y.Pheno.GjResids.noNAs <- Y.Pheno.GjResids[neg.is.na(Y.Pheno),]; \
						InterPath.output <- InterPath(t(X.Pheno.noNAs),Y.Pheno.GjResids.noNAs,as.matrix(X.cov.Pheno.noNAs),t(as.matrix(Z)),Pathways.Regions,cores=cores); \ 
					} else { SkipFlag1 <- 1; }; } else { SkipFlag1 <- 1; }; if (SkipFlag1 == 1) { InterPath.output\\\$Est <- NA; InterPath.output\\\$Eigenvalues <- NA; InterPath.output\\\$PVE <- NA; }; \
					write.table(InterPath.output\\\$Est, paste(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_geneG.Paths\\\", as.character(i), \\\".Est.txt\\\", sep=\\\"\\\"), quote=FALSE, row.name=FALSE, col.name=FALSE); write.table(InterPath.output\\\$Eigenvalues, paste(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_geneG.Paths\\\", as.character(i), \\\".Eigenvalues.txt\\\", sep=\\\"\\\"), quote=FALSE, row.name=FALSE, col.name=FALSE); write.table(InterPath.output\\\$PVE, paste(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_geneG.Paths\\\", as.character(i), \\\".PVE.txt\\\", sep=\\\"\\\"), quote=FALSE, row.name=FALSE, col.name=FALSE); \
				};};\"")
			done; sleep 2
		done; 
	done; 
done
				
#				echo -e "\nmodule load gcc; export C_INCLUDE_PATH=$C_INCLUDE_PATH:/users/mturchin/data/mturchin/miniconda2RH/envs/InterPath/gcc/include/c++/parallel/:/users/mturchin/data/mturchin/miniconda2RH/envs/InterPath/gcc/include/c++/tr1; export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/users/mturchin/data/mturchin/miniconda2RH/envs/InterPath/gcc/include/c++/parallel/:/users/mturchin/data/mturchin/miniconda2RH/envs/InterPath/gcc/include/c++/tr1; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/users/mturchin/data/mturchin/miniconda2RH/envs/InterPath/gcc/include/c++/parallel/:/users/mturchin/data/mturchin/miniconda2RH/envs/InterPath/gcc/include/c++/tr1"; 
##						...print(c(dim(X), dim(X.cov), length(neg.is.na(Y.Pheno)), Pathways.Regions, summary(Pathways.Regions), dim(Pathways.Regions)));...
#						...print(c(dim(X), dim(X.cov), length(neg.is.na(Y.Pheno)), Pathways.Regions));...

#R -q -e "library(\"data.table\"); library(\"doParallel\"); library(\"Rcpp\"); library(\"RcppArmadillo\"); library(\"RcppParallel\"); sourceCpp(\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Test2.cpp\"); \
#	val1 <- \"1,2,3,4;2,3,4,5;3,4,5,6,7\"; \
#	val2 <- unlist(strsplit(val1, \";\")); \
#	print(val2); \
#	val3 <- lapply(strsplit(val2, \",\"), function(x) { return(as.numeric(as.character(x))); }); \
#	val4 <- Test2(val3); \
#	print(val4); \
#"
#sbatch -t 72:00:00 --mem 20g --account ccmb-condo -o Test2.out -e Test2.error --comment "$i $ancestry1 $ancestry2 $k $PathNum" <(echo -e '#!/bin/sh';
#	echo -e "\nR -q -e \"library(\\\"data.table\\\"); library(\\\"doParallel\\\"); library(\\\"Rcpp\\\"); library(\\\"RcppArmadillo\\\"); library(\\\"RcppParallel\\\"); sourceCpp(\\\"/users/mturchin/LabMisc/RamachandranLab/InterPath/Test2.cpp\\\"); \
#        val1 <- \\\"1,2,3,4;2,3,4,5;3,4,5,6,7\\\"; \
#        val2 <- unlist(strsplit(val1, \\\";\\\")); \
#        print(val2); \
#        val3 <- lapply(strsplit(val2, \\\",\\\"), function(x) { return(as.numeric(as.character(x))); }); \
#        val4 <- Test2(val3); \
#        print(val4);\"")


 Pathways.Check <- read.table(\\\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/c2.all.v6.1.wcp_comps.symbols.${ancestry2}.Regions.c2.${k}.txt\\\", header=F); 

for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep -E "Height|BMI"`; do
	for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v Irish | grep African`; do
		for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep Plus`; do
			SECONDS=0;
			ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
			ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
			NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt | wc | awk '{ print $1 }'`	
#			NumPaths=20

			echo $i $ancestry1 $ancestry2 $ancestry3 $k
			for (( PathNum=1; PathNum <= $NumPaths; PathNum=PathNum+1 )); do
				CheckFlag1=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/c2.all.v6.1.wcp_comps.symbols.${ancestry2}.Regions.c2.${k}.txt | head -n 10 | perl -sane 'if ($. == $PathNum2) { print $F[$#F]; }' -- -PathNum2=$PathNum`
				if [ $CheckFlag1 == "TRUE" ] ; then
					R -q -e "library(\"CompQuadForm\"); neg.is.na <- Negate(is.na); \
					Pathways.wGenes <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.geneVSpath.wGenes.${k}.txt\", header=F); \ 
					InterPath.output.Est <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_geneG.Paths${PathNum}.Est.txt\", header=F); \
					InterPath.output.Eigenvalues <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_geneG.Paths${PathNum}.Eigenvalues.txt\", header=F); \
					InterPath.output.PVE <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/$k/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_geneG.Paths${PathNum}.PVE.txt\", header=F); \
					Results1 <- c(); pVals1 <- c(); Genes <- strsplit(Pathways.wGenes[$PathNum,3], ",")[[1]]; for (i in 1:ncol(InterPath.output.Eigenvalues)) { \
						if (neg.is.na(InterPath.output.Est[i,1])) { \ 
							Lambda <- sort(InterPath.output.Eigenvalues[,i], decreasing=TRUE); \
							Davies.Output <- davies(InterPath.output.Est[i,1], lambda=Lambda, acc=1e-8); \
							pVal <- 2*min(1-Davies.Output\$Qq, Davies.Output\$Qq); \
							Results1 <- c(Results1, paste(c(as.character(Genes[i]), InterPath.output.Est[i,1], InterPath.output.PVE[i,1], pVal), collapse=\",\"); 
							pVals1 <- c(pVals, pVal); \ 
						} else { Results1 <- c(Results1, paste(c(as.character(Genes[i]), rep(NA,3)), collapse=\",\")); pVals1 <- c(pVals1, NA); }; \
					}; \ 
					Chisq.Stat <- -2 * sum(log(Results1), na.rm=TRUE); Chisq.pVal <- "PH"; if (Chisq.Stat == 0) { Chisq.pVal <- -9; } else { Chisq.pVal <- pchisq(Chisq.Stat, df=2*length(Results1), lower.tail=FALSE); Results1 <- c(as.character(Pathways.Genes[$PathNum,1]), Chisq.pval, paste(Results1, collapse=\";\")); }; \ 
					write.table(Results1, file=\"\", quote=FALSE, col.name=FALSE, row.name=FALSE);"
				fi
			done | grep -v \> | gzip > /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/$i/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.${i}.${k}.Vs2.GjDrop_wCov_geneG.AllPaths.Results.txt.pre.gz
			duration=$SECONDS; echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
		done;
	done;
done;

#					Pathways <- read.table(\"/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt\", header=F); \ 























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
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.raw | perl -lane 'if ($. == 1) { @vals1; for (my $i = 6; $i <= $#F; $i++) { if ($F[$i] !~ m/HET/) { push(@vals1, $i); } } } print join("\t", @F[@vals1]);' | head -n 1 | perl -lane 'print join("\n", @F);' | head -n 10      
1:729632_T                                                                                                                                                                                                                      
1:752721_G                                                                                                                                                                                                                      
1:754105_T                                                                                                                                                                                                                      
1:756604_G                                                                                                                                                                                                                      
1:759036_A                                                                                                                                                                                                                      
1:761147_C
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
#20180820
(InterPath) [  mturchin@login002  ~]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chr1_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.plinkTemp.fam | head -n 10
1001592_1001592 1001592_1001592 0 0 0 -9
1002560_1002560 1002560_1002560 0 0 0 -9
1003036_1003036 1003036_1003036 0 0 0 -9
1004593_1004593 1004593_1004593 0 0 0 -9
1008167_1008167 1008167_1008167 0 0 0 -9
1009965_1009965 1009965_1009965 0 0 0 -9
1010953_1010953 1010953_1010953 0 0 0 -9
1012491_1012491 1012491_1012491 0 0 0 -9
1013297_1013297 1013297_1013297 0 0 0 -9
1013636_1013636 1013636_1013636 0 0 0 -9
(InterPath) [  mturchin@login002  ~]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.plinkTemp.FIDIIDs | head -n 10
1001592 1001592
1002560 1002560
1003036 1003036
1004593 1004593
1008167 1008167
1009965 1009965
1010953 1010953
1012491 1012491
1013297 1013297
1013636 1013636
(InterPath) [  mturchin@login002  ~]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chr1_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.plinkTemp.fam | tail -n 10
6008418_6008418 6008418_6008418 0 0 0 -9
6010753_6010753 6010753_6010753 0 0 0 -9
6011864_6011864 6011864_6011864 0 0 0 -9
6014349_6014349 6014349_6014349 0 0 0 -9
6018726_6018726 6018726_6018726 0 0 0 -9
6021711_6021711 6021711_6021711 0 0 0 -9
6023159_6023159 6023159_6023159 0 0 0 -9
6024552_6024552 6024552_6024552 0 0 0 -9
6025294_6025294 6025294_6025294 0 0 0 -9
6025455_6025455 6025455_6025455 0 0 0 -9
(InterPath) [  mturchin@login002  ~]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.plinkTemp.FIDIIDs | tail -n 10
6008418 6008418
6010753 6010753
6011864 6011864
6014349 6014349
6018726 6018726
6021711 6021711
6023159 6023159
6024552 6024552
6025294 6025294
6025455 6025455
(InterPath) [  mturchin@login002  ~]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.txt | head -n 10
FID     IID     SEX     ANCESTRY        AGE     PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9     PC10
1001592 1001592 1 African 45 0.04010937 -0.06472882 -0.01579029 -0.0028063 0.00672054 -0.01071403 -0.004161221 0.02206256 0.04384965 -0.04821377
1002560 1002560 0 African 49 -0.05246519 -0.09301313 0.03203676 0.01590448 -0.06281703 -0.01510241 0.05077683 0.004314775 -0.0129375 -0.0155569
1003036 1003036 0 African 41 0.05604095 0.01899383 0.004137209 -0.05324963 -0.005971676 -0.003017448 -0.004951464 -0.006294852 -0.01491296 0.004271876
1004593 1004593 0 African 50 0.04044754 0.04018205 0.004412988 0.051949 0.002875837 0.01332061 0.01493853 0.002691042 0.001538728 -0.008088887
1008167 1008167 0 African 59 0.03080848 -0.09631819 -0.02496982 0.01417764 0.02021924 -0.005031504 0.003485892 0.03258397 0.02134949 0.006985639
1009965 1009965 0 African 58 0.05294918 0.05247505 0.0107574 0.03858143 0.004743263 0.01621968 0.02457712 -0.005881776 -0.01046 0.0151324
1010953 1010953 1 African 48 0.04929234 0.02053753 0.01234649 -0.04077462 -0.01188147 -0.0006558676 -0.009392189 -0.004736571 -0.003848354 0.01392407
1012491 1012491 0 African 41 -0.415706 0.009725911 0.1694342 -0.02152109 0.124938 -0.06456354 0.04600945 0.01419852 0.009809398 -0.03911038
1013297 1013297 1 African 46 0.03162186 0.04584261 -0.0113241 0.02832364 0.01017347 0.01959432 0.0002467427 0.006316162 0.01381216 -0.0051383
(InterPath) [  mturchin@login002  ~]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt | head -n 10
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
(InterPath) [  mturchin@login002  ~]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.sort.ImptHRC.dose.100geno.raw.txt | head -n 10
FID     IID     SEX     ANCESTRY        AGE     PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9     PC10
1001592 1001592 1       African 45      0.04010937      -0.06472882     -0.01579029     -0.0028063      0.00672054      -0.01071403     -0.004161221    0.02206256      0.04384965      -0.04821377
1002560 1002560 0       African 49      -0.05246519     -0.09301313     0.03203676      0.01590448      -0.06281703     -0.01510241     0.05077683      0.004314775     -0.0129375      -0.0155569
1003036 1003036 0       African 41      0.05604095      0.01899383      0.004137209     -0.05324963     -0.005971676    -0.003017448    -0.004951464    -0.006294852    -0.01491296     0.004271876
1004593 1004593 0       African 50      0.04044754      0.04018205      0.004412988     0.051949        0.002875837     0.01332061      0.01493853      0.002691042     0.001538728     -0.008088887
1008167 1008167 0       African 59      0.03080848      -0.09631819     -0.02496982     0.01417764      0.02021924      -0.005031504    0.003485892     0.03258397      0.02134949      0.006985639
1009965 1009965 0       African 58      0.05294918      0.05247505      0.0107574       0.03858143      0.004743263     0.01621968      0.02457712      -0.005881776    -0.01046        0.0151324
1010953 1010953 1       African 48      0.04929234      0.02053753      0.01234649      -0.04077462     -0.01188147     -0.0006558676   -0.009392189    -0.004736571    -0.003848354    0.01392407
1012491 1012491 0       African 41      -0.415706       0.009725911     0.1694342       -0.02152109     0.124938        -0.06456354     0.04600945      0.01419852      0.009809398     -0.03911038
1013297 1013297 1       African 46      0.03162186      0.04584261      -0.0113241      0.02832364      0.01017347      0.01959432      0.0002467427    0.006316162     0.01381216      -0.0051383
(InterPath) [  mturchin@login002  ~]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt | head -n 10
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
(InterPath) [  mturchin@login002  ~]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.Phenos.Transformed.txt | tail -n 10
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
(InterPath) [  mturchin@login002  ~]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.sort.ImptHRC.dose.100geno.raw.txt | head -n 10
FID     IID     SEX     ANCESTRY        AGE     PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9     PC10
1001592 1001592 1       African 45      0.04010937      -0.06472882     -0.01579029     -0.0028063      0.00672054      -0.01071403     -0.004161221    0.02206256      0.04384965      -0.04821377
1002560 1002560 0       African 49      -0.05246519     -0.09301313     0.03203676      0.01590448      -0.06281703     -0.01510241     0.05077683      0.004314775     -0.0129375      -0.0155569
1003036 1003036 0       African 41      0.05604095      0.01899383      0.004137209     -0.05324963     -0.005971676    -0.003017448    -0.004951464    -0.006294852    -0.01491296     0.004271876
1004593 1004593 0       African 50      0.04044754      0.04018205      0.004412988     0.051949        0.002875837     0.01332061      0.01493853      0.002691042     0.001538728     -0.008088887
1008167 1008167 0       African 59      0.03080848      -0.09631819     -0.02496982     0.01417764      0.02021924      -0.005031504    0.003485892     0.03258397      0.02134949      0.006985639
1009965 1009965 0       African 58      0.05294918      0.05247505      0.0107574       0.03858143      0.004743263     0.01621968      0.02457712      -0.005881776    -0.01046        0.0151324
1010953 1010953 1       African 48      0.04929234      0.02053753      0.01234649      -0.04077462     -0.01188147     -0.0006558676   -0.009392189    -0.004736571    -0.003848354    0.01392407
1012491 1012491 0       African 41      -0.415706       0.009725911     0.1694342       -0.02152109     0.124938        -0.06456354     0.04600945      0.01419852      0.009809398     -0.03911038
1013297 1013297 1       African 46      0.03162186      0.04584261      -0.0113241      0.02832364      0.01017347      0.01959432      0.0002467427    0.006316162     0.01381216      -0.0051383
(InterPath) [  mturchin@login002  ~]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.sort.ImptHRC.dose.100geno.raw.txt | tail -n 10
6008418 6008418 0       African 53      -0.05630276     -0.08667724     0.07282399      0.03140086      -0.1160429      -0.05296765     0.04631429      -0.09186453     0.001412289     -0.07666399
6010753 6010753 1       African 44      0.01078697      -0.1036524      -0.004505102    0.02005242      -0.02876383     -0.01757235     0.04441196      0.0574596       -0.0127933      0.02429221
6011864 6011864 0       African 49      0.03284634      -0.1009702      -0.02701615     0.01283153      0.04012466      -0.0006517096   -0.01007781     -0.01583565     0.01192305      0.008451749
6014349 6014349 1       African 49      0.01818009      0.03198973      0.0005825137    0.05307932      -0.01055752     -0.05562142     -0.07428205     0.01018994      -0.01069136     -0.01446091
6018726 6018726 1       African 64      0.05402393      0.02063491      0.01087056      -0.03448014     -0.007878176    -0.00722104     -0.01175685     -0.001874121    0.001628825     0.01807416
6021711 6021711 0       African 43      -0.2581659      -0.02012493     -0.09674245     0.009629706     0.01726518      -0.004650545    0.06122912      -0.001842116    -0.05859597     0.004980408
6023159 6023159 1       African 55      0.0560072       0.04890589      0.01140521      0.04719828      0.003705662     0.02530067      0.02440856      0.009025501     0.008926544     -0.005566456
6024552 6024552 0       African 57      -0.04452432     -0.09559        0.06703958      0.01933989      -0.1068556      -0.0386704      0.06353684      -0.04127274     -0.01821965     -0.03980357
6025294 6025294 1       African 50      0.06248499      0.0175956       0.009511427     -0.05038411     -0.005810936    -0.003879298    -0.005019973    -0.006737438    -0.006394164    -0.001246517
6025455 6025455 0       African 45      0.02788425      -0.1059086      -0.02180886     0.0008563429    0.03900867      0.006091423     -0.002064998    0.006942791     0.01616838      0.01242587
#20180903
> val1 <- rnorm(100); val2 <- matrix(runif(400), ncol=4);
> val3 <- residuals(lm(val1 ~ val2[,1] + val2[,2] + val2[,3] + val2[,4]))
> val4 <- residuals(lm(val1 ~ val2))
> val5 <- residuals(lm(val1 ~ val2[,1] + val2[,2] + val2[,3]))
> head(val3)
          1           2           3           4           5           6 
-0.05020331 -0.73260196 -0.13801556 -0.25875213 -0.99948778  0.47538947 
> head(val4)
          1           2           3           4           5           6 
-0.05020331 -0.73260196 -0.13801556 -0.25875213 -0.99948778  0.47538947 
> head(val5)
          1           2           3           4           5           6 
 0.18811731 -0.41527515 -0.09983512 -0.48137387 -1.07903715  0.65783540 
> tail(val3)
        95         96         97         98         99        100 
-0.5934427  0.8438430 -0.7437926  1.3315325 -0.3405165 -0.3413527 
> tail(val4)
        95         96         97         98         99        100 
-0.5934427  0.8438430 -0.7437926  1.3315325 -0.3405165 -0.3413527 
> tail(val5)
        95         96         97         98         99        100 
-0.6934394  0.8826347 -0.8245094  1.2992763 -0.2049817 -0.1608099 
> val6 <- residuals(lm(val1 ~ val2[,1:3]))
> head(val6)
          1           2           3           4           5           6 
 0.18811731 -0.41527515 -0.09983512 -0.48137387 -1.07903715  0.65783540 
> tail(val6)
        95         96         97         98         99        100 
-0.6934394  0.8826347 -0.8245094  1.2992763 -0.2049817 -0.1608099 
> val7 <- residuals(lm(val1 ~ val2[,1] + val2[,2] + val2[,4]))
> val8 <- residuals(lm(val1 ~ val2[,c(1,2,4)]))
> head(val7)
         1          2          3          4          5          6 
-0.1823634 -0.9146229 -0.1843958 -0.1075598 -1.1088176  0.5349682 
> head(val8)
         1          2          3          4          5          6 
-0.1823634 -0.9146229 -0.1843958 -0.1075598 -1.1088176  0.5349682 
> tail(val7)
        95         96         97         98         99        100 
-0.5425506  0.6321150 -0.6125844  1.4050724 -0.4181078 -0.4992774 
> tail(val8)
        95         96         97         98         99        100 
-0.5425506  0.6321150 -0.6125844  1.4050724 -0.4181078 -0.4992774 
> val9 <- residuals(lm(val1 ~ val2[,c(1,2,4)] + val2[,3]))
> head(val9)
          1           2           3           4           5           6 
-0.05020331 -0.73260196 -0.13801556 -0.25875213 -0.99948778  0.47538947 
> tail(val9)
        95         96         97         98         99        100 
-0.5934427  0.8438430 -0.7437926  1.3315325 -0.3405165 -0.3413527 
[  mturchin@node628  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/pathways/NonSyn/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.African.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.top10resids.Regions.c2.NonSyn.Pathways1.txt | head -n 10                                                                                                                                 FID IID Height BMI Waist Hip
1001592 1001592 0.660046734924629 0.667870305200669 0.890437815929841 1.47104585550913
1002560 1002560 0.97363198221407 1.60211450476216 1.21698301519892 1.46147770852491
1003036 1003036 0.508991800824005 1.27776666161257 1.51693104735486 1.60736094764931
1004593 1004593 0.329534625835474 0.783817569779739 0.760495818081324 1.12900240388704
1008167 1008167 0.609206833203642 -0.479244107773219 -0.413384326092569 -0.558963570021779
1009965 1009965 -0.747667906428213 0.564834821853729 0.733468458454059 0.975724390362749
1010953 1010953 -2.06750348760394 0.178256145584831 0.229447265843397 -0.446644353383066
1012491 1012491 1.27790750003048 -0.238885128056914 0.297054272398364 0.0774529312095399
1013297 1013297 0.307365483483242 -0.0119302971908837 0.114153564877426 -0.197082377656367
[  mturchin@node628  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/pathways/NonSyn/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.African.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.top10resids.Regions.c2.NonSyn.Pathways1.vs2.txt | head -n 10                                                                                                                             FID IID Height BMI Waist Hip
1001592 1001592 0.660046734924629 0.667870305200669 0.890437815929841 1.47104585550913
1002560 1002560 0.97363198221407 1.60211450476216 1.21698301519892 1.46147770852491
1003036 1003036 0.508991800824005 1.27776666161257 1.51693104735486 1.60736094764931
1004593 1004593 0.329534625835474 0.783817569779739 0.760495818081324 1.12900240388704
1008167 1008167 0.609206833203642 -0.479244107773219 -0.413384326092569 -0.558963570021779
1009965 1009965 -0.747667906428213 0.564834821853729 0.733468458454059 0.975724390362749
1010953 1010953 -2.06750348760394 0.178256145584831 0.229447265843397 -0.446644353383066
1012491 1012491 1.27790750003048 -0.238885128056914 0.297054272398364 0.0774529312095399
1013297 1013297 0.307365483483242 -0.0119302971908837 0.114153564877426 -0.197082377656367
[  mturchin@node628  ~/LabMisc/RamachandranLab/InterPath]$cmp /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/pathways/NonSyn/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.African.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.top10resids.Regions.c2.NonSyn.Pathways1.txt /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/pathways/NonSyn/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.African.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.top10resids.Regions.c2.NonSyn.Pathways1.vs2.txt
[  mturchin@node628  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/pathways/NonSyn/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.African.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.top10resids.Regions.c2.NonSyn.Pathways2.txt | head -n 10                                                                                                                                 FID IID Height BMI Waist Hip
1001592 1001592 0.452990769988482 0.624487204431815 0.892737344746463 1.47026976079333
1002560 1002560 1.12619382895407 1.69124245117424 1.24380673306564 1.53225389899892
1003036 1003036 0.444246408577252 1.31749267082513 1.50573794504002 1.5788883412859
1004593 1004593 0.486991613008768 0.731220012039863 0.846776769357214 1.17836524913718
1008167 1008167 0.676996618769845 -0.493141251819213 -0.394245470803315 -0.559822537672775
1009965 1009965 -0.581451462902888 0.708042543837926 0.830359741377277 1.15300861586114
1010953 1010953 -2.07700639110897 0.134408339448983 0.171235567613378 -0.535308609792916
1012491 1012491 1.25814007336742 -0.183503705345015 0.299560144669338 0.127134940613715
1013297 1013297 0.396874312766688 0.0195285334172241 0.203584244244566 -0.116202448370685
[  mturchin@node628  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/pathways/NonSyn/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.African.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.top10resids.Regions.c2.NonSyn.Pathways2.vs2.txt | head -n 10
FID IID Height BMI Waist Hip
1001592 1001592 0.452990769988482 0.624487204431815 0.892737344746463 1.47026976079333
1002560 1002560 1.12619382895407 1.69124245117424 1.24380673306564 1.53225389899892
1003036 1003036 0.444246408577252 1.31749267082513 1.50573794504002 1.5788883412859
1004593 1004593 0.486991613008768 0.731220012039863 0.846776769357214 1.17836524913718
1008167 1008167 0.676996618769845 -0.493141251819213 -0.394245470803315 -0.559822537672775
1009965 1009965 -0.581451462902888 0.708042543837926 0.830359741377277 1.15300861586114
1010953 1010953 -2.07700639110897 0.134408339448983 0.171235567613378 -0.535308609792916
1012491 1012491 1.25814007336742 -0.183503705345015 0.299560144669338 0.127134940613715
1013297 1013297 0.396874312766688 0.0195285334172241 0.203584244244566 -0.116202448370685
[  mturchin@node628  ~/LabMisc/RamachandranLab/InterPath]$cmp /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/pathways/NonSyn/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.African.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.top10resids.Regions.c2.NonSyn.Pathways2.txt /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/pathways/NonSyn/phenos/ukb9200.2017_8_WinterRetreat.Phenos.Transformed.Edit.African.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.top10resids.Regions.c2.NonSyn.Pathways2.vs2.txt
[  mturchin@node628  ~/LabMisc/RamachandranLab/InterPath]$
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/NonSyn/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.wCov.vs2.Paths1.Eigenvalues.txt | head -n 10
-0.00224732301252183 -0.00227160099446316 -0.0024347731369674 -0.00162181320507432 -0.00256814857989494 -0.00196876443482982 -0.00228593547401008 -0.00256158964828816 -0.00214861920508991 -0.00161672091579026
-0.00223519321917511 -0.00225632345751242 -0.00242293311417706 -0.00161209465757868 -0.00255293430744794 -0.00195957601310138 -0.00226983470769512 -0.0025449265594011 -0.00213527190548288 -0.00160707387155447
-0.00222754937210879 -0.0022490132441315 -0.00241277641608698 -0.00160542669140689 -0.002543237649315 -0.00195290651068728 -0.00226195691915617 -0.00253636407178965 -0.00212823010852434 -0.00160081884032346
-0.00222298870267528 -0.00224594010180176 -0.00241135433492922 -0.00160475493194552 -0.00254107522493046 -0.00195053447630417 -0.00226071814320222 -0.00253135183838703 -0.00212601863941403 -0.00159932168357168
-0.00222022883597646 -0.00224216115948836 -0.00240656247188272 -0.0016004317218279 -0.00253635447857526 -0.00194747053580779 -0.00225635592060953 -0.00252888491464825 -0.00212340349527207 -0.00159648959828272
-0.00221777436210844 -0.00224160413509536 -0.00240338056223393 -0.00159933914474707 -0.00253268178877434 -0.00194480191968743 -0.00225362019041798 -0.00252667850334964 -0.00212031323548973 -0.00159572629709841
-0.00221195189800478 -0.00223480067671544 -0.00239944226448667 -0.00159577406950688 -0.00252823335156869 -0.0019414549839269 -0.00224928123680466 -0.00252178329415863 -0.00211600790537148 -0.00159178171069788
-0.00220841370025964 -0.00223041042868392 -0.0023936618520097 -0.00159279940777228 -0.00252186492123962 -0.0019373203090957 -0.00224368965256377 -0.00251614137093512 -0.00210990703304562 -0.00158759759902796
-0.00220458363636256 -0.00222735401712495 -0.0023910184085639 -0.0015917898588491 -0.00252055655843352 -0.0019343851391708 -0.00224255208215947 -0.0025136251445472 -0.00210716972245337 -0.00158684720324457
-0.00220408189363105 -0.00222436301145285 -0.00238771443048796 -0.00158836455987046 -0.00251482631265684 -0.00193259465799658 -0.00223817803450726 -0.00250875781813836 -0.00210485138110701 -0.00158365445917131
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/NonSyn/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.wCov.Paths1.Eigenvalues.txt | head -n 10
-0.00224732301252183 -0.00227160099446316 -0.0024347731369674 -0.00162181320507432 -0.00256814857989494 -0.00196876443482982 -0.00228593547401008 -0.00256158964828816 -0.00214861920508991 -0.00161672091579026
-0.00223519321917511 -0.00225632345751242 -0.00242293311417706 -0.00161209465757868 -0.00255293430744794 -0.00195957601310138 -0.00226983470769512 -0.0025449265594011 -0.00213527190548288 -0.00160707387155447
-0.00222754937210879 -0.0022490132441315 -0.00241277641608698 -0.00160542669140689 -0.002543237649315 -0.00195290651068728 -0.00226195691915617 -0.00253636407178965 -0.00212823010852434 -0.00160081884032346
-0.00222298870267528 -0.00224594010180176 -0.00241135433492922 -0.00160475493194552 -0.00254107522493046 -0.00195053447630417 -0.00226071814320222 -0.00253135183838703 -0.00212601863941403 -0.00159932168357168
-0.00222022883597646 -0.00224216115948836 -0.00240656247188272 -0.0016004317218279 -0.00253635447857526 -0.00194747053580779 -0.00225635592060953 -0.00252888491464825 -0.00212340349527207 -0.00159648959828272
-0.00221777436210844 -0.00224160413509536 -0.00240338056223393 -0.00159933914474707 -0.00253268178877434 -0.00194480191968743 -0.00225362019041798 -0.00252667850334964 -0.00212031323548973 -0.00159572629709841
-0.00221195189800478 -0.00223480067671544 -0.00239944226448667 -0.00159577406950688 -0.00252823335156869 -0.0019414549839269 -0.00224928123680466 -0.00252178329415863 -0.00211600790537148 -0.00159178171069788
-0.00220841370025964 -0.00223041042868392 -0.0023936618520097 -0.00159279940777228 -0.00252186492123962 -0.0019373203090957 -0.00224368965256377 -0.00251614137093512 -0.00210990703304562 -0.00158759759902796
-0.00220458363636256 -0.00222735401712495 -0.0023910184085639 -0.0015917898588491 -0.00252055655843352 -0.0019343851391708 -0.00224255208215947 -0.0025136251445472 -0.00210716972245337 -0.00158684720324457
-0.00220408189363105 -0.00222436301145285 -0.00238771443048796 -0.00158836455987046 -0.00251482631265684 -0.00193259465799658 -0.00223817803450726 -0.00250875781813836 -0.00210485138110701 -0.00158365445917131
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cmp /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/NonSyn/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.wCov.Paths1.Eigenvalues.txt /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/NonSyn/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.wCov.vs2.Paths1.Eigenvalues.txt
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/NonSyn/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.wCov.vs2.Paths11.Eigenvalues.txt | head -n 10
-0.00228608331978058 -0.00258198952258782 -0.00279802874670921 -0.00262466524406116 -0.00248663716192726 -0.00240548609945551 -0.00242827685241754 -0.00240065357457377 -0.00176678425726935 -0.00243112720642938
-0.00227315343221865 -0.00256595870300828 -0.00278722797232656 -0.0026116139891195 -0.0024718407167853 -0.002390253260126 -0.00241060738560444 -0.00238511073716972 -0.00175563730556512 -0.00241345535909982
-0.00226418260534757 -0.00255747664513899 -0.00277453436994319 -0.00260161829929734 -0.00246191185200117 -0.00238158407703431 -0.00240347419748708 -0.00237610187260281 -0.00174888564427693 -0.00240636267844654
-0.00226304149305798 -0.00255327764399049 -0.00277134170218369 -0.0025997581796038 -0.00246089440250557 -0.0023808241537985 -0.00240151926246741 -0.00237539098231764 -0.00174788338296937 -0.00240395925317263
-0.00225941908479925 -0.00254976156185724 -0.0027664630663484 -0.00259566991570418 -0.00245787808736097 -0.0023771382329633 -0.0023976047864874 -0.00237008009760211 -0.00174430914515081 -0.00239816196171044
-0.0022536132984761 -0.00254530625991795 -0.00275730213336582 -0.00259014405561033 -0.00245570443512917 -0.00237423190526185 -0.00239452927339241 -0.00236639859141501 -0.00174341538477266 -0.00239454538647804
-0.00225129039218352 -0.00254097666969169 -0.00275517270853164 -0.00258458901613077 -0.00244765136226556 -0.00236882016399134 -0.00238849925573827 -0.00236224550152693 -0.00173979230486835 -0.0023922908458085
-0.00224544816692075 -0.0025365729842384 -0.00275198983601283 -0.0025792576081589 -0.00244377774265599 -0.00236347153282933 -0.00238491281308742 -0.00235840643312811 -0.00173455981925537 -0.00238524215493502
-0.00224309359594719 -0.00253419767756771 -0.00274661165714535 -0.00257754695442635 -0.00244129123450553 -0.00236283207574842 -0.00238207344804981 -0.00235453365426134 -0.0017342933743257 -0.00238356422598178
-0.00224041067587359 -0.00253061267473722 -0.0027438798370123 -0.00257337460062414 -0.00243729246645276 -0.00235699613281001 -0.00237917919954816 -0.00235182867693919 -0.00173010260298487 -0.00237880665259901
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/NonSyn/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.wCov.Paths11.Eigenvalues.txt | head -n 10
-0.00228608331978058 -0.00258198952258782 -0.00279802874670921 -0.00262466524406116 -0.00248663716192726 -0.00240548609945551 -0.00242827685241754 -0.00240065357457377 -0.00176678425726935 -0.00243112720642938
-0.00227315343221865 -0.00256595870300828 -0.00278722797232656 -0.0026116139891195 -0.0024718407167853 -0.002390253260126 -0.00241060738560444 -0.00238511073716972 -0.00175563730556512 -0.00241345535909982
-0.00226418260534757 -0.00255747664513899 -0.00277453436994319 -0.00260161829929734 -0.00246191185200117 -0.00238158407703431 -0.00240347419748708 -0.00237610187260281 -0.00174888564427693 -0.00240636267844654
-0.00226304149305798 -0.00255327764399049 -0.00277134170218369 -0.0025997581796038 -0.00246089440250557 -0.0023808241537985 -0.00240151926246741 -0.00237539098231764 -0.00174788338296937 -0.00240395925317263
-0.00225941908479925 -0.00254976156185724 -0.0027664630663484 -0.00259566991570418 -0.00245787808736097 -0.0023771382329633 -0.0023976047864874 -0.00237008009760211 -0.00174430914515081 -0.00239816196171044
-0.0022536132984761 -0.00254530625991795 -0.00275730213336582 -0.00259014405561033 -0.00245570443512917 -0.00237423190526185 -0.00239452927339241 -0.00236639859141501 -0.00174341538477266 -0.00239454538647804
-0.00225129039218352 -0.00254097666969169 -0.00275517270853164 -0.00258458901613077 -0.00244765136226556 -0.00236882016399134 -0.00238849925573827 -0.00236224550152693 -0.00173979230486835 -0.0023922908458085
-0.00224544816692075 -0.0025365729842384 -0.00275198983601283 -0.0025792576081589 -0.00244377774265599 -0.00236347153282933 -0.00238491281308742 -0.00235840643312811 -0.00173455981925537 -0.00238524215493502
-0.00224309359594719 -0.00253419767756771 -0.00274661165714535 -0.00257754695442635 -0.00244129123450553 -0.00236283207574842 -0.00238207344804981 -0.00235453365426134 -0.0017342933743257 -0.00238356422598178
-0.00224041067587359 -0.00253061267473722 -0.0027438798370123 -0.00257337460062414 -0.00243729246645276 -0.00235699613281001 -0.00237917919954816 -0.00235182867693919 -0.00173010260298487 -0.00237880665259901
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cmp /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/NonSyn/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.wCov.Paths11.Eigenvalues.txt /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/NonSyn/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.NonSyn.wCov.vs2.Paths11.Eigenvalues.txt
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/Exonic/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.Exonic.wCov.vs2.Paths1.Eigenvalues.txt | head -n 10
-0.00208529573453502 -0.00237479138083705 -0.00222615547924653 -0.00154980220060829 -0.00268282606722285 -0.00210913897692833 -0.00235266208092273 -0.00250796247042682 -0.00224517673197677 -0.00163515053163698
-0.00207461295297694 -0.00235906070100884 -0.00221653685859553 -0.00154034681787756 -0.00266806324872499 -0.00209910991627643 -0.00233556089475806 -0.00249164426239387 -0.00223094264446737 -0.00162527029468844
-0.00206602041550254 -0.0023514909668878 -0.00220646725191893 -0.00153424704940354 -0.00265527524315673 -0.00209115592546522 -0.0023285428022549 -0.00248236468863702 -0.00222413085416501 -0.00161861484504148
-0.00206353534153473 -0.00234775632126711 -0.00220391873044517 -0.00153366680867276 -0.00265323972947756 -0.00209007346345023 -0.00232700401859095 -0.00247784460993444 -0.0022217701212163 -0.00161747290861791
-0.0020612491393431 -0.00234493730751815 -0.00220106622770982 -0.00152929957039372 -0.0026502001927936 -0.00208726975029471 -0.00232264239473506 -0.00247594809699733 -0.00221984348892846 -0.00161424631218674
-0.00205722225051175 -0.00234220375391898 -0.00219776494258054 -0.00152796039870802 -0.00264497288461055 -0.00208320223340403 -0.00232007748186762 -0.0024732256551667 -0.00221578172632891 -0.00161277208586346
-0.00205203475983023 -0.00233698755189239 -0.00219355241017397 -0.0015248141104731 -0.00264053281297807 -0.00207916034871658 -0.00231560283894312 -0.00246863098666736 -0.00220923566742202 -0.0016095779508007
-0.00205057714763665 -0.00233244361799973 -0.00218931079460097 -0.00152185513086955 -0.0026351218336122 -0.00207608855232234 -0.00231002869531123 -0.00246269434438078 -0.00220546433434159 -0.00160572687426547
-0.00204504484404881 -0.00232838651325776 -0.00218767112684689 -0.00152139004026205 -0.00263224799646068 -0.00207312939228383 -0.00230849271344094 -0.00246067435125828 -0.0022013170678739 -0.00160485423810082
-0.00204419905699197 -0.00232590289434806 -0.00218343486793937 -0.0015179308475812 -0.00262751711858506 -0.00206993302448462 -0.00230369620406778 -0.00245588746708482 -0.00219840971476974 -0.00160166028238099
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/Exonic/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.Exonic.wCov.Paths1.Eigenvalues.txt | head -n 10
-0.00208529573453502 -0.00237479138083705 -0.00222615547924653 -0.00154980220060829 -0.00268282606722285 -0.00210913897692833 -0.00235266208092273 -0.00250796247042682 -0.00224517673197677 -0.00163515053163698
-0.00207461295297694 -0.00235906070100884 -0.00221653685859553 -0.00154034681787756 -0.00266806324872499 -0.00209910991627643 -0.00233556089475806 -0.00249164426239387 -0.00223094264446737 -0.00162527029468844
-0.00206602041550254 -0.0023514909668878 -0.00220646725191893 -0.00153424704940354 -0.00265527524315673 -0.00209115592546522 -0.0023285428022549 -0.00248236468863702 -0.00222413085416501 -0.00161861484504148
-0.00206353534153473 -0.00234775632126711 -0.00220391873044517 -0.00153366680867276 -0.00265323972947756 -0.00209007346345023 -0.00232700401859095 -0.00247784460993444 -0.0022217701212163 -0.00161747290861791
-0.0020612491393431 -0.00234493730751815 -0.00220106622770982 -0.00152929957039372 -0.0026502001927936 -0.00208726975029471 -0.00232264239473506 -0.00247594809699733 -0.00221984348892846 -0.00161424631218674
-0.00205722225051175 -0.00234220375391898 -0.00219776494258054 -0.00152796039870802 -0.00264497288461055 -0.00208320223340403 -0.00232007748186762 -0.0024732256551667 -0.00221578172632891 -0.00161277208586346
-0.00205203475983023 -0.00233698755189239 -0.00219355241017397 -0.0015248141104731 -0.00264053281297807 -0.00207916034871658 -0.00231560283894312 -0.00246863098666736 -0.00220923566742202 -0.0016095779508007
-0.00205057714763665 -0.00233244361799973 -0.00218931079460097 -0.00152185513086955 -0.0026351218336122 -0.00207608855232234 -0.00231002869531123 -0.00246269434438078 -0.00220546433434159 -0.00160572687426547
-0.00204504484404881 -0.00232838651325776 -0.00218767112684689 -0.00152139004026205 -0.00263224799646068 -0.00207312939228383 -0.00230849271344094 -0.00246067435125828 -0.0022013170678739 -0.00160485423810082
-0.00204419905699197 -0.00232590289434806 -0.00218343486793937 -0.0015179308475812 -0.00262751711858506 -0.00206993302448462 -0.00230369620406778 -0.00245588746708482 -0.00219840971476974 -0.0016016602(InterPath) [  mturchin@login002  (InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cmp /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/Exonic/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.Exonic.wCov.Paths1.Eigenvalues.txt /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/Height/Exonic/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.Height.Exonic.wCov.vs2.Paths1.Eigenvalues.txt
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$
(InterPath) [  mturchin@login002  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/pathways]$for i in `cat <(echo "Height BMI Waist Hip" | perl -lane 'print join("\n", @F);') | grep -E 'Height|BMI'`; do
>         for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep African`; do
>                 for k in `cat <(echo "NonSyn Exonic ExonicPlus ExonicPlus20kb" | perl -lane 'print join("\n", @F);') | grep -v Plus`; do
>                         SECONDS=0;
>                         ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
>                         ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
>                         NumPaths=`cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.${k}.txt | wc | awk '{ print $1 }'`
...
Height African African NonSyn
8 minutes and 37 seconds elapsed.
Height African African Exonic
8 minutes and 54 seconds elapsed.
BMI African African NonSyn
Warning message:
In davies(InterPath.output.Est[Counter1, 1], lambda = Lambda, acc = 1e-08) :
  Consider playing with 'lim' or 'acc'.
Warning messages:
1: In davies(InterPath.output.Est[Counter1, 1], lambda = Lambda, acc = 1e-08) :
  Consider playing with 'lim' or 'acc'.
2: In davies(InterPath.output.Est[Counter1, 1], lambda = Lambda, acc = 1e-08) :
  Consider playing with 'lim' or 'acc'.
Warning message:
In davies(InterPath.output.Est[Counter1, 1], lambda = Lambda, acc = 1e-08) :
  Consider playing with 'lim' or 'acc'.
18 minutes and 2 seconds elapsed.
BMI African African Exonic
Warning message:
In davies(InterPath.output.Est[Counter1, 1], lambda = Lambda, acc = 1e-08) :
  Consider playing with 'lim' or 'acc'.
Warning message:
In davies(InterPath.output.Est[Counter1, 1], lambda = Lambda, acc = 1e-08) :
  Consider playing with 'lim' or 'acc'.
Warning message:
In davies(InterPath.output.Est[Counter1, 1], lambda = Lambda, acc = 1e-08) :
  Consider playing with 'lim' or 'acc'.
Warning message:
In davies(InterPath.output.Est[Counter1, 1], lambda = Lambda, acc = 1e-08) :
  Consider playing with 'lim' or 'acc'.
17 minutes and 55 seconds elapsed.
[  mturchin@node646  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/pathways]$for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);')`; do
>         ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`
>         ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
>
>         echo $pheno1 $ancestry1 $ancestry2 $ancestry3
>
>         R -q -e "library("data.table"); library("feather"); \
>         ptm <- proc.time(); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); print(proc.time() - ptm); \
>         print(table(apply(Data3, 2, sd) == 0));"
> done
African African
> library(data.table); library(feather);         ptm <- proc.time(); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/African/African/Imputation/mturchin20/ukb_chrAll_v2.African.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); print(proc.time() - ptm);         print(table(apply(Data3, 2, sd) == 0));
   user  system elapsed
 47.167   2.482  49.625

 FALSE
226304
>
>
British British.Ran4000
> library(data.table); library(feather);         ptm <- proc.time(); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran4000/Imputation/mturchin20/ukb_chrAll_v2.British.Ran4000.QCed.reqDrop.QCed.dropRltvs
.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); print(proc.time() - ptm);         print(table(apply(Data3, 2, sd) == 0));
   user  system elapsed
143.384   7.470 150.762

 FALSE
598811
>
>
Caribbean Caribbean
> library(data.table); library(feather);         ptm <- proc.time(); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/Imputation/mturchin20/ukb_chrAll_v2.Caribbean.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); print(proc.time() - ptm);         print(table(apply(Data3, 2, sd) == 0));
   user  system elapsed
 96.876   5.250 102.064

 FALSE
408347
>
>
Chinese Chinese
> library(data.table); library(feather);         ptm <- proc.time(); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/Chinese/Chinese/Imputation/mturchin20/ukb_chrAll_v2.Chinese.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); print(proc.time() - ptm);         print(table(apply(Data3, 2, sd) == 0));
   user  system elapsed
 39.085   2.168  41.227

 FALSE   TRUE
343957      1
>
>
Indian Indian
> library(data.table); library(feather);         ptm <- proc.time(); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/Indian/Indian/Imputation/mturchin20/ukb_chrAll_v2.Indian.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); print(proc.time() - ptm);         print(table(apply(Data3, 2, sd) == 0));
   user  system elapsed
152.501   8.778 161.179

 FALSE   TRUE
504569      1
>
>
Irish Irish
> library(data.table); library(feather);         ptm <- proc.time(); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/Irish/Irish/Imputation/mturchin20/ukb_chrAll_v2.Irish.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); print(proc.time() - ptm);         print(table(apply(Data3, 2, sd) == 0));
   user  system elapsed
324.863  21.356 346.013
Killed
Pakistani Pakistani
> library(data.table); library(feather);         ptm <- proc.time(); Data3 <- fread('zcat /users/mturchin/data/ukbiobank_jun17/subsets/Pakistani/Pakistani/Imputation/mturchin20/ukb_chrAll_v2.Pakistani.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz', header=T); print(proc.time() - ptm);         print(table(apply(Data3, 2, sd) == 0));
   user  system elapsed
 61.909   3.438  65.309

 FALSE
515520
>
>
[  mturchin@node646  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/pathways]$zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | R -q -e "Data1 <- read.table(file('stdin'), header=T); Data1[1:20,apply(Data1, 2, sd) == 0];"
> Data1 <- read.table(file('stdin'), header=T); Data1[1:20,apply(Data1, 2, sd) == 0];
 [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
>
>
[  mturchin@node646  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/pathways]$zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | R -q -e "Data1 <- read.table(file('stdin'), header=T); which(apply(Data1, 2, sd) == 0);"
> Data1 <- read.table(file('stdin'), header=T); which(apply(Data1, 2, sd) == 0);
^[[A^[[BX20.256763_0
      323883
>
>
[  mturchin@node646  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/pathways]$zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | wc
   1424 489796192 983595428
[  mturchin@node646  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/pathways]$zcat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.raw.edit.gz | perl -lane 'print $#F;' | sort | uniq -c
   1424 343957
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$ls -lrt /users/mturchin/data/ukbiobank_jun17/subsets/*/*/mturchin20/Analyses/InterPath/ukb_chrAll_v2.*QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.sort.ImptHRC.dose.100geno.raw.txt
-rw-r--r-- 1 mturchin sramacha  451131 Aug 20 14:27 /users/mturchin/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/ukb_chrAll_v2.African.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.sort.ImptHRC.dose.100geno.raw.txt
-rw-r--r-- 1 mturchin sramacha  622212 Sep 24 00:23 /users/mturchin/data/ukbiobank_jun17/subsets/British/British.Ran4000/mturchin20/Analyses/InterPath/ukb_chrAll_v2.British.Ran4000.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.sort.ImptHRC.dose.100geno.raw.txt
-rw-r--r-- 1 mturchin sramacha  543149 Sep 24 00:23 /users/mturchin/data/ukbiobank_jun17/subsets/Caribbean/Caribbean/mturchin20/Analyses/InterPath/ukb_chrAll_v2.Caribbean.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.sort.ImptHRC.dose.100geno.raw.txt
-rw-r--r-- 1 mturchin sramacha  208766 Sep 24 00:23 /users/mturchin/data/ukbiobank_jun17/subsets/Chinese/Chinese/mturchin20/Analyses/InterPath/ukb_chrAll_v2.Chinese.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.sort.ImptHRC.dose.100geno.raw.txt
-rw-r--r-- 1 mturchin sramacha  760914 Sep 24 00:23 /users/mturchin/data/ukbiobank_jun17/subsets/Indian/Indian/mturchin20/Analyses/InterPath/ukb_chrAll_v2.Indian.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.sort.ImptHRC.dose.100geno.raw.txt
-rw-r--r-- 1 mturchin sramacha 1703525 Sep 24 00:23 /users/mturchin/data/ukbiobank_jun17/subsets/Irish/Irish/mturchin20/Analyses/InterPath/ukb_chrAll_v2.Irish.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.sort.ImptHRC.dose.100geno.raw.txt
-rw-r--r-- 1 mturchin sramacha  239832 Sep 24 00:23 /users/mturchin/data/ukbiobank_jun17/subsets/Pakistani/Pakistani/mturchin20/Analyses/InterPath/ukb_chrAll_v2.Pakistani.QCed.pruned.QCed.dropRltvs.noX.PCAdrop.flashpca.pcs.wFullCovars.sort.ImptHRC.dose.100geno.raw.txt
(InterPath) [  mturchin@login002  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/BMI/Exonic/slurm]$join <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.v6.2.symbols.gmt | awk '{ print $1 }' | sort) <( \
cat <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.biocarta.v6.2.symbols.gmt | awk '{ print $1 "\tBIOCARTA" }') \
<(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.kegg.v6.2.symbols.gmt | awk '{ print $1 "\tKEGG" }') \
<(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.reactome.v6.2.symbols.gmt | awk '{ print $1 "\tREACTOME" }') | sort -k 1,1) | awk '{ print $2 }' | sort | uniq -c
    217 BIOCARTA
    186 KEGG
    674 REACTOME
(InterPath) [  mturchin@login002  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/BMI/Exonic/slurm]$join -v 1 <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.v6.2.symbols.gmt | awk '{ print $1 }' | sort) <( \
cat <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.biocarta.v6.2.symbols.gmt | awk '{ print $1 "\tBIOCARTA" }') \
<(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.kegg.v6.2.symbols.gmt | awk '{ print $1 "\tKEGG" }') \
<(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.reactome.v6.2.symbols.gmt | awk '{ print $1 "\tREACTOME" }') | sort -k 1,1) | wc
    252     252    5562
(InterPath) [  mturchin@login002  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/BMI/Exonic/slurm]$cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.biocarta.v6.2.symbols.gmt | wc
    217    4964   47380
(InterPath) [  mturchin@login002  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/BMI/Exonic/slurm]$cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.kegg.v6.2.symbols.gmt | wc
    186   13247  100017
(InterPath) [  mturchin@login002  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/BMI/Exonic/slurm]$cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.reactome.v6.2.symbols.gmt | wc
    674   38949  330737
(InterPath) [  mturchin@login002  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/BMI/Exonic/slurm]$join -v 1 <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.v6.2.symbols.gmt | awk '{ print $1 }' | sort) <( \
cat <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.biocarta.v6.2.symbols.gmt | awk '{ print $1 "\tBIOCARTA" }') \
<(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.kegg.v6.2.symbols.gmt | awk '{ print $1 "\tKEGG" }') \
<(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.reactome.v6.2.symbols.gmt | awk '{ print $1 "\tREACTOME" }') | sort -k 1,1) | head -n 10
NABA_BASEMENT_MEMBRANES
NABA_COLLAGENS
NABA_CORE_MATRISOME
NABA_ECM_AFFILIATED
NABA_ECM_GLYCOPROTEINS
NABA_ECM_REGULATORS
NABA_MATRISOME
NABA_MATRISOME_ASSOCIATED
NABA_PROTEOGLYCANS
NABA_SECRETED_FACTORS
(InterPath) [  mturchin@login002  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/BMI/Exonic/slurm]$cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.v6.2.symbols.gmt | wc
   1329   70280  575692
(InterPath) [  mturchin@login002  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/BMI/Exonic/slurm]$R -q -e "252 + 217 + 186 + 674"
> 252 + 217 + 186 + 674
[1] 1329
> 
> 
(InterPath) [  mturchin@login002  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/BMI/Exonic/slurm]$join -a 1 -a 2 -e "FALSE" -o 0 1.1 2.1 <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | awk '{ print $1 }' | sort) <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.cp.v6.2.symbols.gmt | awk '{ print $1 }' | sort) | awk '{ print $1 "\t" $3 }' | perl -lane 'if ($F[1] ne "FALSE") { $F[1] = "TRUE"; } print join("\t", @F);' | awk '{ print $2 }' | sort | uniq -c
   3409 FALSE
   1329 TRUE
(InterPath) [  mturchin@login002  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/BMI/Exonic/slurm]$cmp <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | awk '{ print $1 }') <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.wcp_comps.symbols.gmt | awk '{ print $1 }') | wc
      0       0       0
(InterPath) [  mturchin@login002  ~/data/ukbiobank_jun17/subsets/African/African/mturchin20/Analyses/InterPath/BMI/Exonic/slurm]$cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.wcp_comps.symbols.gmt | perl -lane 'print $F[$#F];' | sort | uniq -c
   3409 FALSE
   1329 TRUE
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.wcp_comps.symbols.gmt | wc
   4738   18952  592322
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.NonSyn.txt | wc
   4584   13752 1794113
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$join <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.wcp_comps.symbols.gmt | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.NonSyn.txt | awk '{ print $1 }' | sort) | sort -g -k 3,3 | wc
   4584   18336  573245
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$join <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.wcp_comps.symbols.gmt | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.NonSyn.txt | awk '{ print $1 }' | sort) | sort -g -k 3,3 | head -n 10
KEGG_GLYCOLYSIS_GLUCONEOGENESIS http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_GLYCOLYSIS_GLUCONEOGENESIS 1 TRUE
KEGG_CITRATE_CYCLE_TCA_CYCLE http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_CITRATE_CYCLE_TCA_CYCLE 2 TRUE
KEGG_PENTOSE_PHOSPHATE_PATHWAY http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_PENTOSE_PHOSPHATE_PATHWAY 3 TRUE
KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS 4 TRUE
KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM 5 TRUE
KEGG_GALACTOSE_METABOLISM http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_GALACTOSE_METABOLISM 6 TRUE
KEGG_ASCORBATE_AND_ALDARATE_METABOLISM http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_ASCORBATE_AND_ALDARATE_METABOLISM 7 TRUE
KEGG_FATTY_ACID_METABOLISM http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_FATTY_ACID_METABOLISM 8 TRUE
KEGG_STEROID_BIOSYNTHESIS http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_STEROID_BIOSYNTHESIS 9 TRUE
KEGG_PRIMARY_BILE_ACID_BIOSYNTHESIS http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_PRIMARY_BILE_ACID_BIOSYNTHESIS 10 TRUE
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cmp <(join <(cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.wcp_comps.symbols.gmt | sort -k 1,1) <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.NonSyn.txt | awk '{ print $1 }' | sort) | sort -g -k 3,3 | awk '{ print $1 }') <(cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.c2.NonSyn.txt | awk '{ print $1 }') | wc
      0       0       0
(InterPath) [  mturchin@login002  ~/LabMisc/RamachandranLab/InterPath]$cat /users/mturchin/data/mturchin/Broad/MSigDB/c2.all.v6.1.symbols.gmt | perl -slane 'if ($. == 1) { $input_file = "/users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1b/$ancestry2b/mturchin20/Analyses/InterPath/ukb_chrAll_v2.$ancestry2b.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.bim.AnnovarFormat.TableAnnovar.AAFix.hg19_multianno.GeneSNPs.SemiColonSplit.wRowPos.Regions.NonSyn.txt"; %hash1; open( my $input_fh, "<", $input_file ) || die "Cannot open $input_file: $!"; while(my @row = split(/\s+/, <$input_fh>)) { chomp @row; $hash1{$row[0]} = \$row[1]; } close($input_fh); } my @info1; foreach my $entry1 (@F[2..$#F]) { if ($hash1{$entry1}) { push(@info1, ${$hash1{$entry1}}); } } print $F[0], "\t", $F[1], "\t", join(";", @info1);' -- -ancestry1b=$ancestry1 -ancestry2b=$ancestry2 | perl -lane 'if ($#F == 2) { print join("\t", @F); }' | perl -lane 'my @vals1 = split(/;/, $F[$#F]); my $count1 = 1; my @print1; foreach my $entry1 (@vals1) { my @vals2 = split(/,/, $entry1); my @print2; for (my $entry2 = 0; $entry2 <= $#vals2; $entry2++) { push(@print2, $count1); $count1++ } push(@print1, join(",", @print2)); } $F[$#F] = join(";", @print1); print join("\t", @F);' | head -n 10
KEGG_GLYCOLYSIS_GLUCONEOGENESIS http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_GLYCOLYSIS_GLUCONEOGENESIS 1,2;3;4;5;6;7;8,9;10,11,12;13,14;15;16;17,18;19;20;21,22,23;24,25,26;27,28;29;30;31;32,33;34;35;36,37,38;39,40,41;42;43,44,45;46,47;48,49,50;51
KEGG_CITRATE_CYCLE_TCA_CYCLE    http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_CITRATE_CYCLE_TCA_CYCLE    1;2,3,4;5,6;7;8;9;10;11;12,13,14;15;16;17;18,19,20
KEGG_PENTOSE_PHOSPHATE_PATHWAY  http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_PENTOSE_PHOSPHATE_PATHWAY  1;2;3;4;5;6;7;8,9,10;11,12;13,14,15;16;17
KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS   http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS   1;2;3,4,5;6,7,8,9;10,11,12,13,14;15;16;17,18,19;20;21;22;23;24,25;26;27;28,29;30
KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM    http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM    1;2;3;4;5,6,7;8;9;10;11,12,13;14;15,16,17;18,19;20;21
KEGG_GALACTOSE_METABOLISM       http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_GALACTOSE_METABOLISM       1,2;3;4;5;6;7;8;9,10;11,12,13;14;15,16;17,18,19,20,21,22,23,24;25,26;27,28,29;30;31,32,33,34,35,36,37;38
KEGG_ASCORBATE_AND_ALDARATE_METABOLISM  http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_ASCORBATE_AND_ALDARATE_METABOLISM  1;2;3,4,5;6,7,8;9,10,11,12;13,14,15,16,17;18,19,20;21;22;23;24;25,26;27;28;29,30;31
KEGG_FATTY_ACID_METABOLISM      http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_FATTY_ACID_METABOLISM      1;2,3,4;5,6,7;8,9;10,11;12;13;14;15;16,17,18;19;20;21,22;23,24;25,26,27;28;29,30;31;32;33,34;35,36,37;38,39;40,41,42,43;44,45;46,47;48
KEGG_STEROID_BIOSYNTHESIS       http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_STEROID_BIOSYNTHESIS       1,2;3,4,5,6,7;8,9;10;11;12,13;14;15;16,17
KEGG_PRIMARY_BILE_ACID_BIOSYNTHESIS     http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_PRIMARY_BILE_ACID_BIOSYNTHESIS     1,2;3,4;5;6,7;8,9,10;11;12;13,14,15,16,17;18,19,20
#20181110
(InterPath) [  mturchin@login003  ~/LabMisc/RamachandranLab/InterPath]$for j in `cat <(echo $UKBioBankPops | perl -lane 'print join("\n", @F);') | grep -v Irish`; do
>         ancestry1=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[0];'`; ancestry2=`echo $j | perl -ane 'my @vals1 = split(/;/, $F[0]); print $vals1[1];'`
> #       ln -s /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/mturchin20/Analyses/InterPath/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.Regions.Exonic.c2.InterPath.vs1.AllPhenos.GjDrop_wCov_GK.Rnd2Vrsns.AllPaths.pValHists.vs1.png /users/mturchin/LabMisc/RamachandranLab/InterPath/images/pValHists/20181108.InterPath.Vs2.GK.Analysis.pValHists.$ancestry2.vs1.png
>         cat /users/mturchin/data/ukbiobank_jun17/subsets/$ancestry1/$ancestry2/Imputation/mturchin20/ukb_chrAll_v2.${ancestry2}.QCed.reqDrop.QCed.dropRltvs.PCAdrop.sort.ImptHRC.dose.100geno.fam | wc
> done
   3080   18480  126280
   3824   22944  156784
   3641   21846  149281
   1423    8538   58343
   5158   30948  211478
   1597    9582   65477



~~~














