#purpose: does number of mitochondria and genotype have any relation?
#for each .vcf file a table with columns comprising all mutations + one extra column indicating avereagenumber of reads (gives estimate of number of mitochondria in the cell);
#average number of reads = (total number of reads=sum of all DP4 values)/number of positons

#also test average number of reads inAD/NAD



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

#the size will be reduced later
numAD = 0
numNAD = 0
sampleDataMat = as.matrix(sampleData)
#nameInfo = as.character(sampleData$Samples)
#disIndo = as.character(sampleData$AD)
#genInfo = as.
checker = 0

#load rare mutation info
mutInfo = as.matrix(read.table('mutation_info.txt'))
mutInfo_expanded = which(mutInfo>0,arr.ind = TRUE)
numMut = dim(mutInfo_expanded)[1]
AD_info = matrix(0,numFiles,(numMut+1))
#AD_df = data.frame(sample_name)
#NADFile = 'NADInfo.csv'
NAD_info = matrix(0,numFiles,(numMut+1))


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

   
    
    
    idxTemp_01 = as.character(intersect(as.numeric(rownames(mutInfo)),indx_01))
    
    idxTemp_11 = as.character(intersect(as.numeric(rownames(mutInfo)),indx_11))
    depth = 0
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
      depth = depth + dp4[1]+dp4[2]+dp4[3]+dp4[4]
    }
    
    depth = depth/numPos
    
    if(AD_stat==' 1'){
      numAD = numAD + 1
      if(length(idxTemp_01)>0){
        for(count5 in 1:length(idxTemp_01)){
          if(mutInfo_expanded[idxTemp_01[count5],2]==1){
            AD_info[numAD,mutInfo_expanded[idxTemp_01[count5],1]] = 1
          }
        }
      
      }

      if(length(idxTemp_11)>0){
        for(count5 in 1:length(idxTemp_11)){
          if(mutInfo_expanded[idxTemp_11[count5],2]==2){
            AD_info[numAD,mutInfo_expanded[idxTemp_11[count5],1]] = 1
          }
        }
      }
      AD_info[numAD,(numMut+1)] = depth
    }
    
    if(AD_stat==' 0'){
      numNAD = numNAD + 1
      if(length(idxTemp_01)>0){de
        for(count5 in 1:length(idxTemp_01)){
          if(mutInfo_expanded[idxTemp_01[count5],2]==1){
            NAD_info[numNAD,mutInfo_expanded[idxTemp_01[count5],1]] = 1
          }
        }
      }
      if(length(idxTemp_11)>0){
        for(count5 in 1:length(idxTemp_11)){
          if(mutInfo_expanded[idxTemp_11[count5],2]==2){
            NAD_info[numNAD,mutInfo_expanded[idxTemp_11[count5],1]] = 1
          }
        }
          
      }
      NAD_info[numNAD,(numMut+1)] = depth
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
  AD_info = AD_info[1:numAD,]
  #assign names and save
  dimnames(AD_info)=list(1:numAD,c(rownames(mutInfo_expanded),'num_mit'))
  #AD_info=write.table('AD.csv",header=TRUE)
  
  write.table(AD_info,'AD_mitCount.txt',row.names =FALSE,col.names = TRUE)
}

if(numNAD!=0){
  NAD_info=NAD_info[1:numNAD,]
  
  
  dimnames(NAD_info)=list(1:numNAD,c(rownames(mutInfo_expanded),'num_mit'))
  
  write.table(NAD_info,'NAD_mitCount.txt',row.names =FALSE,col.names = TRUE)
  #write.table(NAD_info,file="NAD.txt")
}



