#hyp testing new version
rm(list=ls())

AD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/ADindRat.txt',header=TRUE)
NAD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/NADindRat.txt',header=TRUE)

numAD = nrow(AD_info)
numNAD = nrow(NAD_info)

#age and AD
t.test(as.numeric(AD_info[,2]),as.numeric(NAD_info[,2]),alternative="two.sided",mu=0,paired=FALSE,var.equal=TRUE,conf.level=0.95)
boxplot(as.numeric(AD_info[,2]), as.numeric(NAD_info[,2]), ylab="age",names=c("AD","healthy"),main="Age span and disease status")

#total number of mutations and AD
totmut_AD = as.numeric(AD_info[,6]) + as.numeric(AD_info[,7])
totmut_NAD = as.numeric(NAD_info[,6]) + as.numeric(NAD_info[,7])

t.test(totmut_AD,totmut_NAD,alternative="two.sided",mu=0,paired=FALSE,var.equal=TRUE,conf.level=0.95)



#old DP4 and AD
t.test(as.numeric(AD_info[,8]),as.numeric(NAD_info[,8]),alternative="two.sided",mu=0,paired=FALSE,var.equal=TRUE,conf.level=0.95)



#new DP4 and AD
t.test(as.numeric(AD_info[,9]),as.numeric(NAD_info[,9]),alternative="two.sided",mu=0,paired=FALSE,var.equal=TRUE,conf.level=0.95)


#Dp4 without log
t.test(as.numeric(AD_info[,10]),as.numeric(NAD_info[,10]),alternative="two.sided",mu=0,paired=FALSE,var.equal=TRUE,conf.level=0.95)

#DP4 without eps in den - BAD!!! (eps in denom is  required)
t.test(as.numeric(AD_info[,11]),as.numeric(NAD_info[,11]),alternative="two.sided",mu=0,paired=FALSE,var.equal=TRUE,conf.level=0.95)




#total number of mutations as comptuted from log DP4 values
totMutVec = append(as.numeric(AD_info[,9]),as.numeric(NAD_info[,9]))
totAgeVec = append(as.numeric(AD_info[,2]),as.numeric(NAD_info[,2]))
#chisq.test(totMutVec,totAgeVec)
cor(totMutVec,totAgeVec)

#normalization
mutVecNorm = (totMutVec - mean(totMutVec))/sd(totMutVec)
ageVecNorm = (totAgeVec - mean(totAgeVec))/sd(totAgeVec)

plot(mutVecNorm,ageVecNorm)
cor(mutVecNorm,ageVecNorm)

#trying t-test after normalization of new DP4 and AD which gave best t-test results - no change

dp4Norm_AD = (as.numeric(AD_info[,9]) - mean(totMutVec))/sd(totMutVec)
dp4Norm_NAD = (as.numeric(NAD_info[,9]) - mean(totMutVec))/sd(totMutVec)
t.test(dp4Norm_AD,dp4Norm_NAD,alternative="two.sided",mu=0,paired=FALSE,var.equal=TRUE,conf.level=0.95)

#gender and AD - Fisher
challenge.df = matrix(c(15,31,24,40), nrow = 2)
fisher.test(challenge.df)





