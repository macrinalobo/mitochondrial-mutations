#stats

rm(list=ls())

AD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/ADindRat.txt',header=TRUE)
NAD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/NADindRat.txt',header=TRUE)

numAD = nrow(AD_info)
numNAD = nrow(NAD_info)
#age
totAge_AD = sum(as.numeric(AD_info[,2]))
totAge_NAD = sum(as.numeric(NAD_info[,2]))
avgAge_AD = mean(as.numeric(AD_info[,2]))
stdAge_AD = sd(as.numeric(AD_info[,2]))
avgAge_NAD = mean(as.numeric(NAD_info[,2]))
stdAge_NAD = sd(as.numeric(NAD_info[,2]))

#number male, female
#avgMale_AD = mean(as.numeric(AD_info[,3]))
#avgFemale_AD = (numAD - sum(as.numeric(AD_info[,3]))) / numAD

#avgMale_NAD = mean(as.numeric(NAD_info[,3]))
#avgFemale_NAD = (numNAD - sum(as.numeric(NAD_info[,3]))) / numNAD

#gender
numMale_AD = sum(as.numeric(AD_info[,3]))
numFemale_AD = numAD - numMale_AD
numMale_NAD = sum(as.numeric(NAD_info[,3]))
numFemale_NAD = numNAD - numMale_NAD
meanMaleAD = numMale_AD/(numMale_AD+numMale_NAD)
meanFemaleAD = numFemale_AD/(numFemale_AD+numFemale_NAD)
meanMaleNAD = numMale_NAD/(numMale_AD+numMale_NAD)
meanFemaleNAD = numFemale_NAD/(numFemale_AD+numFemale_NAD)
#stdMaleAD = 

#0/1, 1/1 genotypes
num01 = sum(as.numeric(AD_info[,6]))
num11 = sum(as.numeric(AD_info[,7]))
mean01 = mean(as.numeric(AD_info[,6]))
mean11 = mean(as.numeric(AD_info[,7]))

num01_NAD = sum(as.numeric(NAD_info[,6]))
num11_NAD = sum(as.numeric(NAD_info[,7]))
mean01_NAD = mean(as.numeric(NAD_info[,6]))
mean11_NAD = mean(as.numeric(NAD_info[,7]))

#allele ratios
ratTot_AD = sum(as.numeric(AD_info[,9]))
ratTot_NAD = sum(as.numeric(NAD_info[,9]))
ratMean_AD = mean(as.numeric(AD_info[,9]))
ratMean_NAD = mean(as.numeric(NAD_info[,9]))


