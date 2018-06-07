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









