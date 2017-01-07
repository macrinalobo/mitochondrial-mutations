#for project

rm(list=ls())

setwd("~/Documents/BioinfCoure/project/ROSMAP_Mito")
pathDir = '/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/MT_'


sampleInfo = '/home/macrina/Documents/BioinfCoure/project/ROSMAP_SamplesInfoForMito.txt'
sampleData = read.table(sampleInfo,header=TRUE)
data_ext = '_Direct_d500k.vcf'
numData = nrow(sampleData)

numPos = 16569
numFiles = 146
#ADFile = 'ADInfo.csv'
AD_info = matrix(0,numFiles,11)
#AD_df = data.frame(sample_name)
#NADFile = 'NADInfo.csv'
NAD_info = matrix(0,numFiles,11)
#the size will be reduced later
numAD = 0
numNAD = 0
nameMatCol = c('sample_name','age','gender','pmi','0','0/1','1/1','allele_ratio','al_rat_ind','alNo_log','al_no_eps')
eps = 0.001#to prevent denom of dp4 ratio going to infinity
sampleDataMat = as.matrix(sampleData)
#nameInfo = as.character(sampleData$Samples)
#disIndo = as.character(sampleData$AD)
#genInfo = as.
checker = 0
for(counter in 1:numData){#change to numData
  vcf_name = sampleDataMat[counter,1]
  #print(vcf_name)
  #readline()
  
  vcf_name = paste(pathDir,vcf_name,data_ext,sep="")
  if (file.exists(vcf_name)){
    #print(vcf_name)
    #readline()
    info = file.info(vcf_name)
    if(info$size < 3000){#empty vcf file; size < 3kB since only headers are present; exact numbr might vary if not mandatory to list all positions
      print("no data")
      checker = checker + 1
      next()
    }
    vcf_df<-read.table(vcf_name,header=FALSE)
    
    if(dim(vcf_df)[1]<numPos){
      checker = checker + 1
      print("not all position info. hence ignore")
      next()
    }
    AD_stat = sampleDataMat[counter,5]
    if(is.na(AD_stat)){#no AD information given
      print("no AD info")
      checker = checker + 1
      next()
    }
    
    
    genotype_individual = as.character(vcf_df[,10])
    
    
    ##for better performance switch to lapply
    
    #position_individual = as.integer(vcf_df[,2])
    #vcf_dfMat = as.matrix(vcf_df)
    
    #temp = strsplit(as.character(vcf_dfMat[1:numPos,"V8"]),"DP4=")
    #lapply(vcf_dMat, function){dp4 = unlist(strsplit(unlist(strsplit(unlist(temp[[2]]))[[2]],";MQ")),","))}
    #dp4 = unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(vcf_dfMat[1:numPos,"V8"]),"DP4="))[[2]],";MQ")),","))
    dp4Num = 0
    dp4Den = 0
    dp4Ind = 0
    dp4No_log = 0
    dp4_no_eps = 0
    for(count2 in 1:numPos){
      temp = strsplit(as.character(vcf_df[[count2,"V8"]]),"DP4=")
      if(length(temp[[1]])==1){#no dp4 field; means no mutations; ignore position
        #print("bad2")
        checker = checker + 1
        next()
      }
      #dp4 = unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(vcf_df[[count2,"V8"]]),"DP4="))[[2]],";MQ")),","))
      dp4 = unlist(strsplit(unlist(strsplit(unlist(temp)[[2]],";MQ")),","))
      dp4 = as.integer(dp4)#number of mutated positions in genotype
      tempNum = dp4[1]+dp4[2]
      tempDen = dp4[3] + dp4[4]
      # print(tempNum)
      # print(tempDen)
      # readline()
      dp4Num = dp4Num + tempNum
      dp4Den = dp4Den + tempDen
      dp4Ind = dp4Ind + log2((tempNum+eps)/(tempDen+eps))
      dp4No_log = dp4No_log + (tempNum/(tempDen+eps))
      if(tempDen!=0){
        dp4_no_eps = dp4_no_eps + log2((tempNum+0.01)/tempDen) #only consider mutated positiions
      }
        
    }
    #all rows don't have dp4 field!
    dp4_ratio = dp4Num/(dp4Den + eps)
    #dp4_ratio = sum((dp4[1:numPos,1]+dp4[1:numPos,2]))/sum((dp4[,3]+dp4[,4]+eps))
    #dp4Sum = sum(log2(dp4_ratio))
    # print(dp4Ind)
    # readline()
    if(AD_stat==' 1'){
      numAD = numAD + 1
      AD_info[numAD,1] = sampleDataMat[counter,1]
      AD_info[numAD,2] = sampleDataMat[counter,2]
      AD_info[numAD,3] = sampleDataMat[counter,3]
      AD_info[numAD,4] = sampleDataMat[counter,4]
      AD_info[numAD,6]=length(grep('0/1',genotype_individual))
      AD_info[numAD,7]=length(grep('1/1',genotype_individual)) 
      AD_info[numAD,5]=length(which('0'==genotype_individual))#0 occurs alone; can't grep since string may contain other 0s
      AD_info[numAD,8] = log2(dp4_ratio)
      AD_info[numAD,9] = dp4Ind
      AD_info[numAD,10] = dp4No_log
      AD_info[numAD,11] = dp4_no_eps
    }
    if(AD_stat==' 0'){
      numNAD = numNAD + 1
      NAD_info[numNAD,1] = sampleDataMat[counter,1]
      NAD_info[numNAD,2] = sampleDataMat[counter,2]
      NAD_info[numNAD,3] = sampleDataMat[counter,3]
      NAD_info[numNAD,4] = sampleDataMat[counter,4]
      NAD_info[numNAD,6]=length(grep('0/1',genotype_individual))
      NAD_info[numNAD,7]=length(grep('1/1',genotype_individual)) 
      NAD_info[numNAD,5]=length(which('0'==genotype_individual))#0 occurs alone; can't grep since string may contain other 0s
      NAD_info[numNAD,8] = log2(dp4_ratio)
      NAD_info[numNAD,9] = dp4Ind
      NAD_info[numNAD,10] = dp4No_log
      NAD_info[numNAD,11] = dp4_no_eps
      
    }
    print("good")
  }
}

if(numAD!=0){
  AD_info = AD_info[1:numAD,]
  #assign names and save
  dimnames(AD_info)=list(1:numAD,nameMatCol)
  #save into csvs
  #AD_info=write.table('AD.csv",header=TRUE)
  
  write.table(AD_info,'ADindRatFin.txt',row.names =FALSE,col.names = TRUE)
}

if(numNAD!=0){
  NAD_info=NAD_info[1:numNAD,]
  
  
  dimnames(NAD_info)=list(1:numNAD,nameMatCol)
  
  write.table(NAD_info,'NADindRatFin.txt',row.names =FALSE,col.names = TRUE)
  #write.table(NAD_info,file="NAD.txt")
}  

