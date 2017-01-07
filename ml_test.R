#clustering, logistic regression and PCA

#building the data frame
#chosen dimensions - age, gender, allele ratio (log sum)
#110 rows, 3 columns
rm(list=ls())

AD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/ADindRat.txt',header=TRUE)
NAD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/NADindRat.txt',header=TRUE)

numAD = nrow(AD_info)
numNAD = nrow(NAD_info)

#combine AD_info and NAD_info
dataIP = rbind(AD_info,NAD_info)
dataIPTrunc = subset(dataIP,select=c(2,3,9))
dimnames(dataIPTrunc) = list(c(1:110),c('age','gender','allelic_ratio'))
dataIPTrunc = apply(dataIPTrunc,c(1,2),as.numeric)
mydata1 = scale(dataIPTrunc[,1])
mydata2 = scale(dataIPTrunc[,3])
mydata = cbind(mydata1,dataIPTrunc[,2],mydata2)
mydata <-data.frame(mydata)
#mydata$X2 = as.factor(mydata$X2)
indices = c(rep(1,numAD),rep(0,numNAD))
#plotting data
#install.packages("scatterplot3d")
library(scatterplot3d)
colnames(mydata) = c('age','gender','allele_ratio')

# data
DF <- data.frame(mydata$age,mydata$gender,mydata$allele_ratio,group = indices)

# create the plot
s3d <- with(DF, scatterplot3d(mydata$age, mydata$gender, mydata$allele_ratio, color = c(rep('green',numAD),rep('blue',numNAD)), pch = 19))

#legend
#legend(s3d$(mydata$age, mydata$gender, mydata$allele_ratiomydata$age, mydata$gender, mydata$allele_ratio.convert(0.5, 0.7, 0.5), pch = 19, yjust=0,legend = levels(DF$group), col = seq_along(levels(DF$group)))

#mydata = cbind(mydata,indices)
#logistic regression
fit.logit <- glm(indices~mydata$age+mydata$gender+mydata$allele_ratio,family="binomial")
fit.logit$coefficients

#probit regression
fit.probit <- glm(indices~mydata$age+mydata$gender+mydata$allele_ratio,family=binomial(link="probit"))
fit.probit$coefficients


#with dividing data into train and test sets (70-30)
numTrain = ceiling(0.7*110)
numTest = 110 - numTrain

#randomly generate nuTrain indices for training
indTrain = sample(1:110,numTrain,replace=FALSE)
mydataTrain = mydata[indTrain,]
mydataTrain1 = mydataTrain
indTest = setdiff(1:110, indTrain)

#rownames(mydataTest) = c(1:numTest)

#logisti cregression
fit.logit2 <- glm(indices[indTrain]~mydataTrain$age+mydataTrain$gender+mydataTrain$allele_ratio,family="binomial")
fit.logit2$coefficients

#probit regression
fit.probit2 <- glm(indices[indTrain]~mydataTrain$age+mydataTrain$gender+mydataTrain$allele_ratio,family=binomial(link="probit"))
fit.probit2$coefficients

summary(fit.logit2)
summary(fit.probit2)


#testing

#testing data is named mydataTrain for convienece to predict function
mydataTrain = mydata[indTest,]
predicted <- predict(fit.logit2,mydataTrain,type='response') 
#rescaling to [0.1]
#predicted = (predicted - min(predicted))/(max(predicted) - min(predicted))

indxNAD = names(which(predicted<=0.5))
indxAD = names(which(predicted>0.5))

tp = 0
fp = 0
tn = 0
fn = 0

#actual OP for the predcition of AD
actual_AD = indices[as.numeric(rownames(mydataTrain[indxAD,]))]
actual_NAD = indices[as.numeric(rownames(mydataTrain[indxNAD,]))]

tp = sum(actual_AD)
tn = length(actual_NAD) - sum(actual_NAD)
fp = length(actual_AD) - sum(actual_AD)
fn = sum(actual_NAD)
f1 = (2*tp)/(2*tp + fp + fn)
accuracy = (tp + tn)/numTest
error_rate = (1-accuracy)
tpr = tp/(tp+fn)
fpr = fp/(fp+tn)
precision = tp/(tp+fp)

print(accuracy)



#plotting logistic
#too many dimensions; dont' plot






#don't think there's need for clustering; simple plotting will do for this application
# #clustering
# # Determine number of clusters
# wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
# for (i in 2:15) {
#   wss[i] <- sum(kmeans(mydata,centers=i)$withinss)
# }
# plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
# 
# #cluster; with 2 clusters
# fit <- kmeans(mydata, 2) # 2 cluster solution
# # get cluster means
# aggregate(mydata,by=list(fit$cluster),FUN=mean)
# # append cluster assignment
# mydata <- data.frame(mydata, fit$cluster) 
# 

