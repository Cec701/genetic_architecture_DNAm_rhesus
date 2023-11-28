#!/usr/bin/Rscript

#get MAF for each SNP-CpG pair in results
#load RData for each separately  

library(data.table)
#complete dataset
fdr_corrected_res_all_chr <- fread("fdr_corrected_chr_res_all.txt")

#chr1
fdr_corrected_res_all_chr1 <- subset(fdr_corrected_res_all_chr, chr =="1")
load("merge_filter_573_chr1_aa.RData")

maf_test_aa <- maf[-idx]
maf_test_aa <- as.data.frame(maf_test_aa)
maf_test_aa$CpG <- as.vector(CpGlist)
maf_test_aa$SNP <- as.vector(SNPlist)

#load("merge_filter_573_chr1_ab.RData")

maf_test_ab <- maf[-idx]
maf_test_ab <- as.data.frame(maf_test_ab)
maf_test_ab$CpG <- as.vector(CpGlist)
maf_test_ab$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr1_ac.RData")

maf_test_ac <- maf[-idx]
maf_test_ac <- as.data.frame(maf_test_ac)
maf_test_ac$CpG <- as.vector(CpGlist)
maf_test_ac$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr1_ad.RData")

maf_test_ad <- maf[-idx]
maf_test_ad <- as.data.frame(maf_test_ad)
maf_test_ad$CpG <- as.vector(CpGlist)
maf_test_ad$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr1_ae.RData")

maf_test_ae <- maf[-idx]
maf_test_ae <- as.data.frame(maf_test_ae)
maf_test_ae$CpG <- as.vector(CpGlist)
maf_test_ae$SNP <- as.vector(SNPlist)

colnames(maf_test_aa) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_ab) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_ac) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_ad) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_ae) <- c("maf_test", "CpG", "SNP")

maf_test_chr1 <- rbind(maf_test_aa,maf_test_ab, maf_test_ac, maf_test_ad, maf_test_ae)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr1, maf_test_chr1, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr1.txt")

#chr2 
load("merge_filter_573_chr2.RData")

fdr_corrected_res_all_chr2 <- subset(fdr_corrected_res_all_chr, chr =="2")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr2, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr2.txt")

#chr3
fdr_corrected_res_all_chr3 <- subset(fdr_corrected_res_all_chr, chr =="3")
load("merge_filter_573_chr3_aa.RData")

maf_test_aa <- maf[-idx]
maf_test_aa <- as.data.frame(maf_test_aa)
maf_test_aa$CpG <- as.vector(CpGlist)
maf_test_aa$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr3_ab.RData")

maf_test_ab <- maf[-idx]
maf_test_ab <- as.data.frame(maf_test_ab)
maf_test_ab$CpG <- as.vector(CpGlist)
maf_test_ab$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr3_ac.RData")

maf_test_ac <- maf[-idx]
maf_test_ac <- as.data.frame(maf_test_ac)
maf_test_ac$CpG <- as.vector(CpGlist)
maf_test_ac$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr3_ad.RData")

maf_test_ad <- maf[-idx]
maf_test_ad <- as.data.frame(maf_test_ad)
maf_test_ad$CpG <- as.vector(CpGlist)
maf_test_ad$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr3_ae.RData")

maf_test_ae <- maf[-idx]
maf_test_ae <- as.data.frame(maf_test_ae)
maf_test_ae$CpG <- as.vector(CpGlist)
maf_test_ae$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr3_af.RData")

maf_test_af <- maf[-idx]
maf_test_af <- as.data.frame(maf_test_af)
maf_test_af$CpG <- as.vector(CpGlist)
maf_test_af$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr3_ag.RData")

maf_test_ag <- maf[-idx]
maf_test_ag <- as.data.frame(maf_test_ag)
maf_test_ag$CpG <- as.vector(CpGlist)
maf_test_ag$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr3_ah.RData")

maf_test_ah <- maf[-idx]
maf_test_ah <- as.data.frame(maf_test_ah)
maf_test_ah$CpG <- as.vector(CpGlist)
maf_test_ah$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr3_ai.RData")

maf_test_ai <- maf[-idx]
maf_test_ai <- as.data.frame(maf_test_ai)
maf_test_ai$CpG <- as.vector(CpGlist)
maf_test_ai$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr3_aj.RData")

maf_test_aj <- maf[-idx]
maf_test_aj <- as.data.frame(maf_test_aj)
maf_test_aj$CpG <- as.vector(CpGlist)
maf_test_aj$SNP <- as.vector(SNPlist)

colnames(maf_test_aa) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_ab) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_ac) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_ad) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_ae) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_af) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_ag) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_ah) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_ai) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_aj) <- c("maf_test", "CpG", "SNP")

maf_test_chr3 <- rbind(maf_test_aa,maf_test_ab, maf_test_ac, maf_test_ad, maf_test_ae, maf_test_af, maf_test_ag, maf_test_ah, maf_test_ai, maf_test_aj)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr3, maf_test_chr3, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr3.txt")

#4
fdr_corrected_res_all_chr4 <- subset(fdr_corrected_res_all_chr, chr =="4")
load("merge_filter_573_chr4.RData")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr4, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr4.txt")

#5
fdr_corrected_res_all_chr5 <- subset(fdr_corrected_res_all_chr, chr =="5")
load("merge_filter_573_chr5.RData")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr5, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr5.txt")

#6
fdr_corrected_res_all_chr6 <- subset(fdr_corrected_res_all_chr, chr =="6")
load("merge_filter_573_chr6.RData")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr6, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr6.txt")

#7
fdr_corrected_res_all_chr7 <- subset(fdr_corrected_res_all_chr, chr =="7")
load("merge_filter_573_chr7.RData")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr7, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr7.txt")

#8
fdr_corrected_res_all_chr8 <- subset(fdr_corrected_res_all_chr, chr =="8")
load("merge_filter_573_chr8.RData")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr8, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr8.txt")

#9
fdr_corrected_res_all_chr9 <- subset(fdr_corrected_res_all_chr, chr =="9")
load("merge_filter_573_chr9.RData")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr9, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr9.txt")

#10
fdr_corrected_res_all_chr10 <- subset(fdr_corrected_res_all_chr, chr =="10")
load("merge_filter_573_chr10.RData")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr10, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr10.txt")

#11
fdr_corrected_res_all_chr11 <- subset(fdr_corrected_res_all_chr, chr =="11")
load("merge_filter_573_chr11.RData")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr11, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr11.txt")

#12
fdr_corrected_res_all_chr12 <- subset(fdr_corrected_res_all_chr, chr =="12")
load("merge_filter_573_chr12.RData")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr12, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr12.txt")

#13
fdr_corrected_res_all_chr13 <- subset(fdr_corrected_res_all_chr, chr =="13")
load("merge_filter_573_chr13.RData")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr13, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr13.txt")

#14
fdr_corrected_res_all_chr14 <- subset(fdr_corrected_res_all_chr, chr =="14")
load("merge_filter_573_chr14.RData")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr14, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr14.txt")

#15
fdr_corrected_res_all_chr15 <- subset(fdr_corrected_res_all_chr, chr =="15")
load("merge_filter_573_chr15.RData")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr15, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr15.txt")

#16
fdr_corrected_res_all_chr16 <- subset(fdr_corrected_res_all_chr, chr =="16")
load("merge_filter_573_chr16.RData")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr16, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr16.txt")

#17
fdr_corrected_res_all_chr17 <- subset(fdr_corrected_res_all_chr, chr =="17")
load("merge_filter_573_chr17.RData")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr17, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr17.txt")

#18
fdr_corrected_res_all_chr18 <- subset(fdr_corrected_res_all_chr, chr =="18")
load("merge_filter_chr18_573.RData")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr18, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr18.txt")

#19
fdr_corrected_res_all_chr19 <- subset(fdr_corrected_res_all_chr, chr =="19")
load("merge_filter_573_chr19_aa.RData")

maf_test_aa <- maf[-idx]
maf_test_aa <- as.data.frame(maf_test_aa)
maf_test_aa$CpG <- as.vector(CpGlist)
maf_test_aa$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr19_ab.RData")

maf_test_ab <- maf[-idx]
maf_test_ab <- as.data.frame(maf_test_ab)
maf_test_ab$CpG <- as.vector(CpGlist)
maf_test_ab$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr19_ac.RData")

maf_test_ac <- maf[-idx]
maf_test_ac <- as.data.frame(maf_test_ac)
maf_test_ac$CpG <- as.vector(CpGlist)
maf_test_ac$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr19_ad.RData")

maf_test_ad <- maf[-idx]
maf_test_ad <- as.data.frame(maf_test_ad)
maf_test_ad$CpG <- as.vector(CpGlist)
maf_test_ad$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr19_ae.RData")

maf_test_ae <- maf[-idx]
maf_test_ae <- as.data.frame(maf_test_ae)
maf_test_ae$CpG <- as.vector(CpGlist)
maf_test_ae$SNP <- as.vector(SNPlist)

colnames(maf_test_aa) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_ab) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_ac) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_ad) <- c("maf_test", "CpG", "SNP")
colnames(maf_test_ae) <- c("maf_test", "CpG", "SNP")

maf_test_chr19 <- rbind(maf_test_aa,maf_test_ab, maf_test_ac, maf_test_ad, maf_test_ae)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr19, maf_test_chr19, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr19.txt")

#20
fdr_corrected_res_all_chr20 <- subset(fdr_corrected_res_all_chr, chr =="20")
load("merge_filter_573_chr20.RData")

maf_test <- maf[-idx]
maf_test <- as.data.frame(maf_test)
maf_test$CpG <- as.vector(CpGlist)
maf_test$SNP <- as.vector(SNPlist)

library(dplyr)
maf_merge <- merge(fdr_corrected_res_all_chr20, maf_test, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(maf_merge, file="fdr_corrected_res_maf_chr20.txt")

#sum the MAF in bins
fdr_corrected_all_maf_chr1 <- read.table("fdr_corrected_res_maf_chr1.txt", header=TRUE)
fdr_corrected_all_maf_chr2 <- read.table("fdr_corrected_res_maf_chr2.txt", header=TRUE)
fdr_corrected_all_maf_chr3 <- read.table("fdr_corrected_res_maf_chr3.txt", header=TRUE)
fdr_corrected_all_maf_chr4 <- read.table("fdr_corrected_res_maf_chr4.txt", header=TRUE)
fdr_corrected_all_maf_chr5 <- read.table("fdr_corrected_res_maf_chr5.txt", header=TRUE)
fdr_corrected_all_maf_chr6 <- read.table("fdr_corrected_res_maf_chr6.txt", header=TRUE)
fdr_corrected_all_maf_chr7 <- read.table("fdr_corrected_res_maf_chr7.txt", header=TRUE)
fdr_corrected_all_maf_chr8 <- read.table("fdr_corrected_res_maf_chr8.txt", header=TRUE)
fdr_corrected_all_maf_chr9 <- read.table("fdr_corrected_res_maf_chr9.txt", header=TRUE)
fdr_corrected_all_maf_chr10 <- read.table("fdr_corrected_res_maf_chr10.txt", header=TRUE)
fdr_corrected_all_maf_chr11 <- read.table("fdr_corrected_res_maf_chr11.txt", header=TRUE)
fdr_corrected_all_maf_chr12 <- read.table("fdr_corrected_res_maf_chr12.txt", header=TRUE)
fdr_corrected_all_maf_chr13 <- read.table("fdr_corrected_res_maf_chr13.txt", header=TRUE)
fdr_corrected_all_maf_chr14 <- read.table("fdr_corrected_res_maf_chr14.txt", header=TRUE)
fdr_corrected_all_maf_chr15 <- read.table("fdr_corrected_res_maf_chr15.txt", header=TRUE)
fdr_corrected_all_maf_chr16 <- read.table("fdr_corrected_res_maf_chr16.txt", header=TRUE)
fdr_corrected_all_maf_chr17 <- read.table("fdr_corrected_res_maf_chr17.txt", header=TRUE)
fdr_corrected_all_maf_chr18 <- read.table("fdr_corrected_res_maf_chr18.txt", header=TRUE)
fdr_corrected_all_maf_chr19 <- read.table("fdr_corrected_res_maf_chr19.txt", header=TRUE)
fdr_corrected_all_maf_chr20 <- read.table("fdr_corrected_res_maf_chr20.txt", header=TRUE)

fdr_corrected_res_all_w_maf <- rbind(fdr_corrected_all_maf_chr1, fdr_corrected_all_maf_chr2, fdr_corrected_all_maf_chr3, fdr_corrected_all_maf_chr4, fdr_corrected_all_maf_chr5,
                                     fdr_corrected_all_maf_chr6, fdr_corrected_all_maf_chr7, fdr_corrected_all_maf_chr8, fdr_corrected_all_maf_chr9, fdr_corrected_all_maf_chr10, fdr_corrected_all_maf_chr11, fdr_corrected_all_maf_chr12,fdr_corrected_all_maf_chr13,
                                     fdr_corrected_all_maf_chr14, fdr_corrected_all_maf_chr15,fdr_corrected_all_maf_chr16, fdr_corrected_all_maf_chr17, fdr_corrected_all_maf_chr18, 
                                     fdr_corrected_all_maf_chr19, fdr_corrected_all_maf_chr20)

write.table(fdr_corrected_res_all_w_maf, file="fdr_corrected_res_all_w_maf.txt")

# Number of loci
sum(fdr_corrected_res_all_w_maf$maf_test <= 0.1)
sum(fdr_corrected_res_all_w_maf$maf_test > 0.1 & fdr_corrected_res_all_w_maf$maf_test <= 0.2)
sum(fdr_corrected_res_all_w_maf$maf_test > 0.2 & fdr_corrected_res_all_w_maf$maf_test <= 0.3)
sum(fdr_corrected_res_all_w_maf$maf_test > 0.3 & fdr_corrected_res_all_w_maf$maf_test <= 0.4)
sum(fdr_corrected_res_all_w_maf$maf_test > 0.4 & fdr_corrected_res_all_w_maf$maf_test <= 0.5)
sum(fdr_corrected_res_all_w_maf$maf_test > 0.5 & fdr_corrected_res_all_w_maf$maf_test <= 0.6)
sum(fdr_corrected_res_all_w_maf$maf_test > 0.6 & fdr_corrected_res_all_w_maf$maf_test <= 0.7)
sum(fdr_corrected_res_all_w_maf$maf_test > 0.7 & fdr_corrected_res_all_w_maf$maf_test <= 0.8)
sum(fdr_corrected_res_all_w_maf$maf_test > 0.8 & fdr_corrected_res_all_w_maf$maf_test <= 0.9)

#Number of significant meQTL
maf_merge_res_0.05 <- subset(fdr_corrected_res_all_w_maf, fdr_FDR_IMAGE_res_all < 0.05)

sum(maf_merge_res_0.05$maf_test <= 0.1)
sum(maf_merge_res_0.05$maf_test > 0.1 & maf_merge_res_0.05$maf_test <= 0.2)
sum(maf_merge_res_0.05$maf_test > 0.2 & maf_merge_res_0.05$maf_test <= 0.3)
sum(maf_merge_res_0.05$maf_test > 0.3 & maf_merge_res_0.05$maf_test <= 0.4)
sum(maf_merge_res_0.05$maf_test > 0.4 & maf_merge_res_0.05$maf_test <= 0.5)
sum(maf_merge_res_0.05$maf_test > 0.5 & maf_merge_res_0.05$maf_test <= 0.6)
sum(maf_merge_res_0.05$maf_test > 0.6 & maf_merge_res_0.05$maf_test <= 0.7)
sum(maf_merge_res_0.05$maf_test > 0.7 & maf_merge_res_0.05$maf_test <= 0.8)
sum(maf_merge_res_0.05$maf_test > 0.8 & maf_merge_res_0.05$maf_test <= 0.9)
