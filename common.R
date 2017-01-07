
#common mutations

#rare_common_part2

#go through each .vcf file and count number of rare mutations in it. also check if AD or not; creates AD_info and NAD_info file with raremutations

#copied and vairied test1.R

#for project

##add shor tmodification for saving number of arre mutations per cell in the AD_info table


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
commonAD = matrix(0,numFiles,23)
#AD_df = data.frame(sample_name)
#NADFile = 'NADInfo.csv'
commonNAD = matrix(0,numFiles,23)
#the size will be reduced later
numAD = 0
numNAD = 0
nameMatCol = c(as.character(1:23))
eps = 0.001#to prevent denom/num of dp4 ratio going to infinity
sampleDataMat = as.matrix(sampleData)
#nameInfo = as.character(sampleData$Samples)
#disIndo = as.character(sampleData$AD)
#genInfo = as.
checker = 0

#load rare mutation info
common_mutInfo = as.matrix(read.table('commonmutation_info3.txt'))
common_expanded = which(common_mutInfo>=3,arr.ind=TRUE)
common_expanded = common_expanded[order(as.numeric(rownames(common_expanded))),]
mutThresh = 3
numCommon = dim(common_expanded)[1]
#do an if-else for absence of file

#numRare = c(0,0) #AD, not AD
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
    # readline()
    
    # 
    indx_01 = grep('0/1',genotype_individual)
    indx_11 = grep('1/1',genotype_individual)
    # #print(counter)
    # 
    # #check if rare
    # 
    # 
    # #print(indx_11)
    # #print(indx_01)
    # #readline()
    
    
    idxTemp_01 = as.character(intersect(as.numeric(rownames(common_mutInfo)),indx_01))
    
    idxTemp_11 = as.character(intersect(as.numeric(rownames(common_mutInfo)),indx_11))
    
    if(AD_stat==' 1'){
      numAD = numAD + 1
      if(length(idxTemp_01)>0){
        for(count5 in 1:length(idxTemp_01)){
          if(common_expanded[idxTemp_01[count5],2]==1){
            commonAD[numAD,common_expanded[idxTemp_01[count5],1]] = 1
          }
        }
      }
      
      if(length(idxTemp_11)>0){
        for(count5 in 1:length(idxTemp_11)){
          if(common_expanded[idxTemp_11[count5],2]==2){
            commonAD[numAD,common_expanded[idxTemp_11[count5],1]] = 1
          }
        }
      }
      
    }
    
    if(AD_stat==' 0'){
      numNAD = numNAD + 1
      if(length(idxTemp_01)>0){
        for(count5 in 1:length(idxTemp_01)){
          if(common_expanded[idxTemp_01[count5],2]==1){
            commonNAD[numNAD,common_expanded[idxTemp_01[count5],1]] = 1
          }
        }
      }
      if(length(idxTemp_11)>0){
        for(count5 in 1:length(idxTemp_11)){
          if(common_expanded[idxTemp_11[count5],2]==2){
            commonNAD[numNAD,common_expanded[idxTemp_11[count5],1]] = 1
          }
        }
        
      }
    }
    
    
    #print("good")
    print(counter)
    #readline()
  }
}

##can do but chose not to
#remove columns with all zeros
#commonAD = commonAD[ , !apply(commonAD==0,2,all)]
#commonNAD = commonNAD[ , !apply(commonNAD==0,2,all)]
###

if(numAD!=0){
  commonAD = commonAD[1:numAD,]
  #assign names and save
  dimnames(commonAD)=list(1:numAD,rownames(common_expanded))
  #AD_info=write.table('AD.csv",header=TRUE)
  
  write.table(commonAD,'ADind_common3.txt',row.names =FALSE,col.names = TRUE)
}

if(numNAD!=0){
  commonNAD=commonNAD[1:numNAD,]
  
  
  dimnames(commonNAD)=list(1:numNAD,rownames(common_expanded))
  
  write.table(commonNAD,'NADind_common3.txt',row.names =FALSE,col.names = TRUE)
  #write.table(NAD_info,file="NAD.txt")
}



