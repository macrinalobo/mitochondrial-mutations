#rare and common mutations.
#observed mutation counts: 0/0, 0/1, 1/1 are used and not dp4 counts

rm(list=ls())

setwd("~/Documents/BioinfCoure/project/ROSMAP_Mito")
#pathDir = '/home/macrina/Documents/BioinfCoure/project/ROSMAP_Mito/MT_'

fileNames = list.files(pattern='*.vcf',full.names = TRUE)

# sampleInfo = '/home/macrina/Documents/BioinfCoure/project/ROSMAP_SamplesInfoForMito.txt'
# sampleData = read.table(sampleInfo,header=TRUE)
#data_ext = '_Direct_d500k.vcf'
#numData = nrow(sampleData)

numPos = 16569
numFiles = length(fileNames)

#matrix showing position names along rows and 0/0, 0/1, 1/1, allelic_ratios along columns
genotype_info = matrix(0,nrow=numPos,ncol=3)
dimnames(genotype_info) = list(c(1:numPos),c('0/1','1/1','allelic_ratio'))

checker = 0

numSelect = 0 #number of files slected for this analysis

fileNames = list.files(pattern='*.vcf',full.names = TRUE)

for(count in 1:numFiles){
  #read each .vcf file
  vcf_name = fileNames[count]
  #print(vcf_name)
  #readline()
  
  #even vcf files which don't have their AD information in .txt file are considered here 
  #since our goal is only to identify rare and common mutations
  
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
  
  genotype_individual = as.character(vcf_df[,10])
  
  genotype_info[grep('0/1',genotype_individual),'0/1'] = genotype_info[grep('0/1',genotype_individual),'0/1'] + 1
  genotype_info[grep('1/1',genotype_individual),'1/1'] = genotype_info[grep('1/1',genotype_individual),'1/1'] + 1
  
    
  numSelect = numSelect + 1
}

#rare mutation defined at 1%
#numThresh = ceiling(numSelect/100)
numThresh = 3
#positions of rare mutation

#removing rows with no mutations
row_sub = apply(genotype_info, 1, function(row) any(row !=0 ))
gen_inf_reduced = genotype_info[row_sub,1:2]

write.table(gen_inf_reduced,'mutation_info.txt',row.names =TRUE,col.names = TRUE)

#only rows with rare mutations
row_sub = unique(rownames(which((gen_inf_reduced<=numThresh & gen_inf_reduced > 0),arr.ind = TRUE)))
gen_rareMut = gen_inf_reduced[row_sub,]
#row_sub = apply(gen_inf_reduced, 1, function(row) any(row !=0 && row <= 2))
#gen_rareMut = gen_inf_reduced[row_sub,]
write.table(gen_rareMut,'raremutation_info_thresh3.txt',row.names =TRUE,col.names = TRUE)

#only preserving rows with raere mutations i.e. preserve rows with at least 1 (any) element <= 2
#row_sub = apply(gen_inf_reduced, 1, function(row) any(row !=0 ))

#only rows with common mutations
row_sub = unique(rownames(which((gen_inf_reduced>numThresh),arr.ind = TRUE)))
gen_rareMut = gen_inf_reduced[row_sub,]
#row_sub = apply(gen_inf_reduced, 1, function(row) any(row !=0 && row <= 2))
#gen_rareMut = gen_inf_reduced[row_sub,]
write.table(gen_rareMut,'commonmutation_info3.txt',row.names =TRUE,col.names = TRUE)


#############################################



