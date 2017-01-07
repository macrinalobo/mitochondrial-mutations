#pmi information


rm(list=ls())

AD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/ADindRatPMI.txt',header=TRUE)
NAD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/NADindRatPMI.txt',header=TRUE)

numAD = nrow(AD_info)
numNAD = nrow(NAD_info)

t.test(as.numeric(AD_info[,4]),as.numeric(NAD_info[,4]),alternative="two.sided",mu=0,paired=FALSE,var.equal=TRUE,conf.level=0.95)