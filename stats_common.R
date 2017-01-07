rm(list=ls())

AD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/ADind_common3.txt',header=TRUE)
NAD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/NADind_common3.txt',header=TRUE)

numAD = nrow(AD_info)
numNAD = nrow(NAD_info)
#age
numCommon = dim(AD_info)[2]
avg_AD = rep(0,numCommon)
avg_NAD = rep(0,numCommon)
tot_AD = rep(0,numCommon)
tot_NAD = rep(0,numCommon)

for(count in 1:numCommon){
  tot_AD[count] = sum(as.numeric(AD_info[,count]))
  tot_NAD[count] = sum(as.numeric(NAD_info[,count]))
  avg_AD[count] = mean(as.numeric(AD_info[,count]))
  #stdAD = sd(as.numeric(AD_info[,2]))
  avg_NAD[count] = mean(as.numeric(NAD_info[,count]))
  #stdAge_NAD = sd(as.numeric(NAD_info[,2]))
}