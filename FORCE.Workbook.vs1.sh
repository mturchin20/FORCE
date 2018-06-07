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

conda install r-base
conda install r-devtools
conda install r-knitr
conda install r-testthat
conda install r-cairo
conda install r-ashr
conda install r-rcolorbrewer
##20180311
##conda install flashpca -- failed
#conda install eigen
##conda install spectra -- installed a Python package called spectra, not what was intended
#conda install boost
##conda install libgomp -- failed
#conda install gcc
#20180315
conda install r-essentials
conda install -c anaconda fonts-anaconda
conda install -c bioconda r-extrafont
#20180525
conda install tabix






