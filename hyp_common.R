#hypothessi testing for each common

rm(list=ls())

AD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/ADind_common3.txt')

NAD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/NADind_common3.txt')

colnames(AD_info) = AD_info[1,]
colnames(NAD_info) = NAD_info[1,]

numAD = nrow(AD_info)
numNAD = nrow(NAD_info)

AD_info = AD_info[2:numAD,]
NAD_info = NAD_info[2:numNAD,]

rownames(AD_info) = c(1:(numAD-1))
rownames(NAD_info) = c(1:(numNAD-1))

numAD = numAD - 1
numNAD = numNAD - 1

numCommon = dim(AD_info)[2]


for(count in 1:numCommon){
  print(t.test(AD_info[,count],NAD_info[,count],alternative="two.sided",mu=0,paired=FALSE,var.equal=TRUE,conf.level=0.95))
  print(count)
  readline()
}


