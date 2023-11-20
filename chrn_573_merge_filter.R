#!/usr/bin/Rscript

#######################
# filtering
#######################

#Enter chr # 1-20 (or run as an array) 
CHROMNUM=2
name = read.table('pedigree_ordered_sample_list.txt')  
NUMIDS=dim(name)[1]

#merging then filter on merge (below)
#######################
# extract info from ASM files
#######################
name=as.matrix(name)

SNP=numeric()
CpG=numeric()
SNPlist=numeric()
CpGlist=numeric()
count=1

# split up the ASM outputs by chromosome

for(j in 1:NUMIDS){
  str=paste(name[j],'.chr.',CHROMNUM,'.R1_CGmap.PASS2.DP5.ASM',sep='') 
  data=read.table(str,header=T)
  N=nrow(data)
  
  for(k in 2:N)
  {
    SNP[count]=as.numeric(data[k,2])
    CpG[count]=as.numeric(data[k,6])
    count=count+1
  }
}

sSNP=unique(sort(SNP))
N=length(sSNP)
for(k in 1:N){
  idx=which(SNP==sSNP[k])
  CpGtmp=unique(sort(CpG[idx]))
  n1=length(CpGtmp)
  n2=length(SNPlist)
  SNPlist[(n2+1):(n2+n1)]=sSNP[k]
  CpGlist[(n2+1):(n2+n1)]=CpGtmp
}

r1=matrix(0,ncol=NUMIDS,nrow=length(CpGlist))
r2=r1
y1=r1
y2=r1
genotype=r1

for(j in 1:NUMIDS){
  print(j)
  str=paste(name[j],'.chr.',CHROMNUM,'.R1_CGmap.PASS2.DP5.ASM',sep='')
  if(file.exists(str))
  {
    data=read.table(str)
    data=as.matrix(data)
    N=nrow(data)
  }else{
    next
  }
  for(k in 2:N)
  {
    snp=as.numeric(data[k,2])
    cpg=as.numeric(data[k,6])
    idx=which(SNPlist==snp)
    idx1=which(CpGlist[idx]==cpg)
    if(data[k,3]==data[k,4]&&data[k,3]==data[k,5])
    {
      tmp=strsplit(data[k,7],'-')
      y1[idx[idx1],j]=as.numeric(tmp[[1]][1])
      y2[idx[idx1],j]=0
      r1[idx[idx1],j]=as.numeric(tmp[[1]][1])+as.numeric(tmp[[1]][2])
      r2[idx[idx1],j]=0
      genotype[idx[idx1],j]=0
    }else if(data[k,3]!=data[k,4]&&data[k,3]!=data[k,5]){
      tmp=strsplit(data[k,7],'-')
      y1[idx[idx1],j]=as.numeric(tmp[[1]][1])
      y2[idx[idx1],j]=0
      r1[idx[idx1],j]=as.numeric(tmp[[1]][1])+as.numeric(tmp[[1]][2])
      r2[idx[idx1],j]=0
      genotype[idx[idx1],j]=2
    }else{
      tmp=strsplit(data[k,7],'-')
      tmp2=strsplit(data[k,8],'-')
      y1[idx[idx1],j]=as.numeric(tmp[[1]][1])
      y2[idx[idx1],j]=as.numeric(tmp2[[1]][1])
      r1[idx[idx1],j]=as.numeric(tmp[[1]][1])+as.numeric(tmp[[1]][2])
      r2[idx[idx1],j]=as.numeric(tmp2[[1]][1])+as.numeric(tmp2[[1]][2])
      genotype[idx[idx1],j]=1
    }
  } }

str2=paste('data_chr',CHROMNUM,'.RData',sep='')
save(SNPlist,CpGlist,y1,y2,r1,r2,genotype,file=str2)

load(str2)
r1m=r1
r2m=r2
y1m=y1
y2m=y2
genotypem=genotype

ym=y1m+y2m
rm=r1m+r2m

N=nrow(genotypem)
NUMIDS=ncol(genotypem)
#######################
# filtering
#######################

# save SNP and CpGlist after each filter step (to see what filtering does)
# Following [49], we excluded CpG sites (i) that were measured in less than 20 individuals, (ii) where methylation levels fell below 10% or above 90% in at least 90% of measured individuals, (iii) that had a mean read depth less than 5, or (iv) that were paired with a SNP with MAF < 0.05 across individuals for whom DNA methylation estimates were available. To avoid potential mapping bias, we also excluded CpG sites with apparent differences in methylation levels between reference and alternate alleles that were larger than 0.6. 

# FILTER 1 - remove sites not covered in at least 10% of IDs(0.1)

num=numeric()
for(i in 1:N){
  num[i]=length(which(rm[i,]>0))
}
idx=which(num < NUMIDS*0.1)
y1m=y1m[-idx,]
y2m=y2m[-idx,]
r1m=r1m[-idx,]
r2m=r2m[-idx,]
genotypem=genotypem[-idx,]
ym=ym[-idx,]
rm=rm[-idx,]
CpGlist=CpGlist[-idx]
SNPlist=SNPlist[-idx]
N=nrow(genotypem)

#dimensions - what are you losing how many SNPs are lost
# quality control

Filt1_CpG_count <- length(CpGlist)
Filt1_SNP_count <- length(SNPlist)

print(Filt1_CpG_count)
print(Filt1_SNP_count)

write.table(CpGlist,paste('May10_FILT1_CpGs_chr',CHROMNUM,'.txt',sep=''),row.names=F,sep='\t')
write.table(SNPlist,paste('May10_FILT1_SNPs_chr',CHROMNUM,'.txt',sep=''),row.names=F,sep='\t')

# FILTER 2 - remove sites where methylation levels fall below 10% or above 90% in at least 90% of measured individuals
num=numeric()
num1=numeric()
num2=numeric()
for(i in 1:N){
  idx=which(rm[i,]>0)
  ratio=ym[i,idx]/rm[i,idx]
  num[i]=length(idx)
  num1[i]=length(which(ratio<0.1))
  num2[i]=length(which(ratio>0.9))
}
ratio1=num1/num
ratio2=num2/num
idx1=which(ratio1>0.9)
idx2=which(ratio2>0.9)
idx=union(idx1,idx2)
y1m=y1m[-idx,]
y2m=y2m[-idx,]
r1m=r1m[-idx,]
r2m=r2m[-idx,]
genotypem=genotypem[-idx,]
ym=ym[-idx,]
rm=rm[-idx,]
CpGlist=CpGlist[-idx]
SNPlist=SNPlist[-idx]
N=nrow(genotypem)


Filt2_CpG_count <- length(CpGlist)
Filt2_SNP_count <- length(SNPlist)

print(Filt2_CpG_count)
print(Filt2_SNP_count)

write.table(CpGlist,paste('May10_FILT2__test_CpGs_chr',CHROMNUM,'.txt',sep=''),row.names=F,sep='\t')
write.table(SNPlist,paste('May10_FILT2_SNPs_test_chr',CHROMNUM,'.txt',sep=''),row.names=F,sep='\t')

# FILTER 3 - remove sites with mean read depth less than 5
a_rdp=numeric()
for(i in 1:N){
  a_rdp[i]=sum(rm[i,])/length(which(rm[i,]>0))
}
idx=which(a_rdp<5)
y1m=y1m[-idx,]
y2m=y2m[-idx,]
r1m=r1m[-idx,]
r2m=r2m[-idx,]
genotypem=genotypem[-idx,]
ym=ym[-idx,]
rm=rm[-idx,]
CpGlist=CpGlist[-idx]
SNPlist=SNPlist[-idx]
N=nrow(genotypem)

Filt3_CpG_count <- length(CpGlist)
Filt3_SNP_count <- length(SNPlist)

print(Filt3_CpG_count)
print(Filt3_SNP_count)

write.table(CpGlist,paste('May10_test_FILT3_CpGs_chr',CHROMNUM,'.txt',sep=''),row.names=F,sep='\t')
write.table(SNPlist,paste('May10_test_FILT3_SNPs_chr',CHROMNUM,'.txt',sep=''),row.names=F,sep='\t')

# FILTER 4 - that were paired with a SNP with MAF < 0.05 across individuals for whom DNA methylation estimates were available
# 
maf=numeric()
for(i in 1:N){
  idx=which(rm[i,]>0)
  maf[i]=(length(which(genotypem[i,idx]==1))+length(which(genotypem[i,idx]==2))*2)/length(idx)/2
}
idx=union(which(maf<0.05),which(maf>0.95))
y1m=y1m[-idx,]
y2m=y2m[-idx,]
r1m=r1m[-idx,]
r2m=r2m[-idx,]
genotypem=genotypem[-idx,]
ym=ym[-idx,]
rm=rm[-idx,]
CpGlist=CpGlist[-idx]
SNPlist=SNPlist[-idx]
N=nrow(genotypem)

Filt4_CpG_count <- length(CpGlist)
Filt4_SNP_count <- length(SNPlist)

print(Filt4_CpG_count)
print(Filt4_SNP_count)

write.table(CpGlist,paste('May10_FILT4_test_CpGs_chr',CHROMNUM,'.txt',sep=''),row.names=F,sep='\t')
write.table(SNPlist,paste('May10_FILT4_test_SNPs_chr',CHROMNUM,'.txt',sep=''),row.names=F,sep='\t')

geno<-list()
geno[[1]]<-matrix(0,ncol=N,nrow=NUMIDS)
geno[[2]]<-matrix(0,ncol=N,nrow=NUMIDS)
for(i in 1:N){
  for(j in 1:NUMIDS)
  {
    if(is.na(genotypem[i,j]))
    {
      geno[[2]][j,i]=0/0
      geno[[1]][j,i]=0/0
    }else if(genotypem[i,j]==1){
      geno[[2]][j,i]=1
    }else if(genotypem[i,j]==2){
      geno[[2]][j,i]=1
      geno[[1]][j,i]=1
    }
  }
}
names(geno) <- c('hap1', 'hap2')
data<-list()
data[[1]]<-matrix(0,ncol=N,nrow=NUMIDS)
data[[2]]<-matrix(0,ncol=N,nrow=NUMIDS)
data[[3]]<-matrix(0,ncol=N,nrow=NUMIDS)
data[[4]]<-matrix(0,ncol=N,nrow=NUMIDS)
data[[5]]<-matrix(0,ncol=N,nrow=NUMIDS)
data[[6]]<-matrix(0,ncol=N,nrow=NUMIDS)
data[[1]]<-t(rm)
data[[2]]<-t(ym)
data[[3]]<-t(r1m)
data[[4]]<-t(r2m)
data[[5]]<-t(y1m)
data[[6]]<-t(y2m)
names(data) <- c('r', 'y', 'r1', 'r2', 'y1', 'y2')

save.image(file = "merge_filter_573_chr2.RData")
