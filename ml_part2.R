#ml_part2

#logit and probit with LOOCV

#building the data frame
#chosen dimensions - age, gender, allele ratio (log sum)
#110 rows, 3 columns
rm(list=ls())

AD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/ADind_rare_thresh3.txt',header=TRUE)
NAD_info = read.table('/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/NADind_rare_thresh3.txt',header=TRUE)

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
colnames(mydata) = c('age','gender','allele_ratio')
#mydata$X2 = as.factor(mydata$X2)
indices = c(rep(1,numAD),rep(0,numNAD))


#logit + LOOCV
dataSz = numAD + numNAD

#write separately for first iteration

thresh = seq(0,1,0.1)
numThresh = length(thresh)
tp = rep(0,numThresh)
fp = rep(0,numThresh)
tn = rep(0,numThresh)
fn = rep(0,numThresh)

for(count1 in 1:numThresh){
  threshVal = thresh[count1]
  predVec = rep(0,dataSz)
  
  #first element separately
  indTrain = indices[2:dataSz]
  myDataTrain = mydata[2:dataSz,]
  
  #fit.logit <- glm(indTrain~myDataTrain$age+myDataTrain$gender+myDataTrain$allele_ratio,family="binomial")
  fit.probit <- glm(indTrain~myDataTrain$age+myDataTrain$gender+myDataTrain$allele_ratio,family=binomial(link="probit"))
  #fit.logit$coefficients
  myDataTrain1 = myDataTrain
  #name test set as "myDataTrain"
  #iterating over different values of log reg threshold
  myDataTrain = mydata[1,]
  expOP = indices[1]
  predicted <- predict(fit.probit,myDataTrain,type='response')
  if(predicted>threshVal){
    predVec[1] = 1
  }else{
    predVec[1] = 0
  }
    
  
  for(count in 2:dataSz){
    #train on all except count
    indTrain = c(indices[1:(count-1)],indices[(count+1):dataSz])
    myDataTrain = rbind(mydata[1:(count-1),],mydata[(count+1):dataSz,])
    #fit.logit <- glm(indTrain~myDataTrain$age+myDataTrain$gender+myDataTrain$allele_ratio,family="binomial")
    fit.probit <- glm(indTrain~myDataTrain$age+myDataTrain$gender+myDataTrain$allele_ratio,family=binomial(link="probit"))
    #fit.logit$coefficients
    myDataTrain1 = myDataTrain
    #name test set as "myDataTrain"
    #iterating over different values of log reg threshold
    myDataTrain = mydata[count,]
    expOP = indices[count]
    predicted <- predict(fit.probit,myDataTrain,type='response')
    if(predicted>threshVal)
      predVec[count] = 1
    else
      predVec[count] = 0
    #rescaling to [0.1]
    #predicted = (predicted - min(predicted))/(max(predicted) - min(predicted))
  }
  
  #computing stats for this threshold value
  fp[count1] = length(which((indices - predVec) < 0))
  fn[count1] = length(which((indices-predVec)>0))
  #tn[count1] = length(which((indices-predVec)==0)) - 
  tp[count1] = numAD - fn[count1]
  tn[count1] = numNAD - fp[count1]
  
  
}

#ROC plot
tpr = tp/(tp+fn)
fpr = fp/(fp+tn)

plot(fpr,tpr)
