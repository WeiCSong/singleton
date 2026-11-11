library(pgenlibr)
library(parallel)
library(data.table)

trio=fread("file/trio.txt",data.table=F)
samplelist=fread("/dssg/home/acct-bioxsyy/bioxsyy-user2/ukb/file/test.sample",data.table=F)
samplelist=samplelist[-1,]

for (chr in 1:22){
pseudo=fread(paste(c("singleton/c",chr,".pseudosingleton.list"),collapse=""),data.table=F,head=F)
pseudo=pseudo[,1]
singleton=fread(paste(c("singleton/c",chr,".masked.singleton"),collapse=""),data.table=F)
sample=fread(paste(c("trio/trio.",chr,".psam"),collapse=""),data.table=F)
var=fread(paste(c("trio/trio.",chr,".pvar"),collapse=""),data.table=F)
var=var[,3,drop=F]


pgen=NewPgen(paste(c("trio/trio.",chr,".pgen"),collapse=""))
chunksize=500000
n=GetVariantCt(pgen)
chunk=ceiling(n/chunksize)
stlist=lapply(trio[,1],function(id){
  offidx <- which(samplelist[,1] == id)
  offst <- singleton[which(singleton[,5] == offidx),]
  return(offst)
})
rm(singleton)
gc()
for(i in 1:chunk){
  offset <- chunksize * (i - 1)
  start <- offset + 1
  end <- min(offset + chunksize, n)
  mat <- ReadIntList(pgen, start:end)
  #maskcol=which(apply(mat,2,max)==2 | colSums(mat,na.rm=F)>2)
  

  
  for(j in 1:nrow(trio)) {
    offspring <- which(sample[,1] == trio[j,1])
    father <- which(sample[,1] == trio[j,2])
    mother <- which(sample[,1] == trio[j,3])
    offst <- stlist[[j]]
    
    l <- paste("DRAGEN", offst[,1], sep=":")
    l <- paste(l, offst[,2], sep=":")
    l <- paste(l, offst[,3], sep=":")
    offst <- paste(l, offst[,4], sep=":")
    
    int <- mat[c(offspring, father, mother), ]
    prob <- which(int[1,] == 1)
    denovo <- prob[which(int[2,prob] + int[3,prob] == 0)]
    denovo <- var[denovo + offset, 1]
    denovo <- setdiff(denovo, "masked")
    dnm_in_st <- length(intersect(denovo, offst))
    dnm_out_st <- length(setdiff(denovo, offst))
    
    inhere <- prob[which(int[2,prob] + int[3,prob] == 1)]
    inhere <- var[inhere + offset, 1]
    inhere <- setdiff(inhere,"masked")
    inh_in_st <- length(intersect(inhere, pseudo))
    
    int <- data.frame(
      chr=chr,
      chunk=i,
      offset=offset,
      offspring = trio[j,1],
      dnm_in_st = dnm_in_st,
      dnm_out_st = dnm_out_st,
      inh_in_st = inh_in_st
    )
  fwrite(int,"trio.stat",sep="\t",append=T)
  }
  gc()
  print(i)
}
}

####################################################
stat=fread("trio.stat",data.table=F,fill=T)
library(dplyr)
library(ggplot2)
stat$group=paste(stat[,1],stat[,2],sep="_")
test <- stat %>%
  group_by(group) %>%
  summarise(
   meandos = mean(dnm_out_st, na.rm = TRUE)
  )
mask=test[which(test[,2]>1),1]
mask=data.frame(mask)[,1]

stat=stat[which(!stat$group %in% mask),]

result <- stat %>%
  group_by(offspring) %>%
  summarise(
    sum_value1 = sum(dnm_in_st, na.rm = TRUE),
    sum_value2 = sum(inh_in_st, na.rm = TRUE)
  ) %>%
  mutate(
    total = sum_value1 + sum_value2,
    ratio = sum_value1 / total  # 计算比例
  )

result_long <- result %>%
  select(offspring, sum_value1, sum_value2) %>%
  pivot_longer(
    cols = c(sum_value1, sum_value2),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = factor(variable, 
                      levels = c("sum_value2", "sum_value1"),
                      labels = c("inheritST", "denovoST"))  
  )
result_long$offspring=factor(result_long$offspring,levels=offlevel)
