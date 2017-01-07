#correlation between number of mitochondria and each kind of mutation - not considering the effect of AD and NAD

#install.packages('Hmisc') 
#installed from terminal via
#sudo apt-get install r-cran-hmisc
library(Hmisc)
rm(list=ls())

AD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/AD_mitCount.txt',header=TRUE)
NAD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/NAD_mitCount.txt',header=TRUE)
info_tot = rbind(AD_info,NAD_info)
numAD = nrow(AD_info)
numNAD = nrow(NAD_info)


#removing completely empty columns
info_tot = info_tot[ , !apply(info_tot==0,2,all)]

correln = rep(0,105)

for(counter in 1:105){
  rcorr(info_tot[,106],info_tot[,counter])
  #print(rcorr(info_tot[,106],info_tot[,counter])$p.value)
  readline()
}





print(correln)