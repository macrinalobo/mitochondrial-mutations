
rm(list=ls())

AD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/ADind_withrare_thresh3.txt',header=TRUE)
NAD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/NADind_withrare_thresh3.txt',header=TRUE)

numAD = nrow(AD_info)
numNAD = nrow(NAD_info)

totAge_AD = sum(as.numeric(AD_info[,12]))
totAge_NAD = sum(as.numeric(NAD_info[,12]))
avgAge_AD = mean(as.numeric(AD_info[,12]))
avgAge_NAD = mean(as.numeric(NAD_info[,12]))
