rm(list=ls())

AD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/ADind_withrare_thresh3.txt',header=TRUE)
NAD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/NADind_withrare_thresh3.txt',header=TRUE)

numAD = nrow(AD_info)
numNAD = nrow(NAD_info)


#rare muta and AD
t.test(as.numeric(AD_info[,12]),as.numeric(NAD_info[,12]),alternative="two.sided",mu=0,paired=FALSE,var.equal=TRUE,conf.level=0.95)