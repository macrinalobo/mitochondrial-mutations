#stats num_mitochondria

rm(list=ls())

AD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/AD_mitCount.txt',header=TRUE)
NAD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/NAD_mitCount.txt',header=TRUE)

numAD = nrow(AD_info)
numNAD = nrow(NAD_info)

tot_AD = sum(as.numeric(AD_info[,106]))
mean_AD = mean(as.numeric(AD_info[,106]))
tot_NAD = sum(as.numeric(NAD_info[,106]))
mean_NAD = mean(as.numeric(NAD_info[,106]))