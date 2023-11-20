#!/usr/bin/Rscript

#mean read count for each SNP-CpG pair in results
#and mean methylation level for each CpG site 
# load RData for each separately 

# methylation levels for each cpg site
# ym[i,]/rm[i,]

library(data.table)
#complete dataset with maf
fdr_corrected_res_all_chr <- fread("fdr_corrected_res_all_w_maf.txt")

#chr1
fdr_corrected_res_all_chr1 <- subset(fdr_corrected_res_all_chr, chr =="1")
load("merge_filter_573_chr1_aa.RData")

#ok I mean this works
test_chr1_aa_meth_levels <- ym/rm

sample_names <- read.table("pedigree_ordered_sample_list.txt", header=FALSE)
colnames(test_chr1_aa_meth_levels) <- sample_names$V1

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_aa <- as.data.frame(mean_methylation_list)
mean_meth_read_aa$mean_read <- mean_read_count_list
mean_meth_read_aa$CpG <- as.vector(CpGlist)
mean_meth_read_aa$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr1_ab.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_ab <- as.data.frame(mean_methylation_list)
mean_meth_read_ab$mean_read <- mean_read_count_list
mean_meth_read_ab$CpG <- as.vector(CpGlist)
mean_meth_read_ab$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr1_ac.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_ac <- as.data.frame(mean_methylation_list)
mean_meth_read_ac$mean_read <- mean_read_count_list
mean_meth_read_ac$CpG <- as.vector(CpGlist)
mean_meth_read_ac$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr1_ad.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_ad <- as.data.frame(mean_methylation_list)
mean_meth_read_ad$mean_read <- mean_read_count_list
mean_meth_read_ad$CpG <- as.vector(CpGlist)
mean_meth_read_ad$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr1_ae.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_ae <- as.data.frame(mean_methylation_list)
mean_meth_read_ae$mean_read <- mean_read_count_list
mean_meth_read_ae$CpG <- as.vector(CpGlist)
mean_meth_read_ae$SNP <- as.vector(SNPlist)

colnames(mean_meth_read_aa) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_ab) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_ac) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_ad) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_ae) <- c("mean_methylation", "mean_read", "CpG", "SNP")

mean_meth_read_chr1 <- rbind(mean_meth_read_aa,mean_meth_read_ab, mean_meth_read_ac, mean_meth_read_ad, mean_meth_read_ae)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr1, mean_meth_read_chr1, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr1.txt")

#chr2 
load("merge_filter_573_chr2.RData")

fdr_corrected_res_all_chr2 <- subset(fdr_corrected_res_all_chr, chr =="2")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr2, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr2.txt")

#chr3
fdr_corrected_res_all_chr3 <- subset(fdr_corrected_res_all_chr, chr =="3")
load("merge_filter_573_chr3_aa.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_aa <- as.data.frame(mean_methylation_list)
mean_meth_read_aa$mean_read <- mean_read_count_list
mean_meth_read_aa$CpG <- as.vector(CpGlist)
mean_meth_read_aa$SNP <- as.vector(SNPlist)


load("merge_filter_573_chr3_ab.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_ab <- as.data.frame(mean_methylation_list)
mean_meth_read_ab$mean_read <- mean_read_count_list
mean_meth_read_ab$CpG <- as.vector(CpGlist)
mean_meth_read_ab$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr3_ac.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_ac <- as.data.frame(mean_methylation_list)
mean_meth_read_ac$mean_read <- mean_read_count_list
mean_meth_read_ac$CpG <- as.vector(CpGlist)
mean_meth_read_ac$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr3_ad.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_ad <- as.data.frame(mean_methylation_list)
mean_meth_read_ad$mean_read <- mean_read_count_list
mean_meth_read_ad$CpG <- as.vector(CpGlist)
mean_meth_read_ad$SNP <- as.vector(SNPlist)


load("merge_filter_573_chr3_ae.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_ae <- as.data.frame(mean_methylation_list)
mean_meth_read_ae$mean_read <- mean_read_count_list
mean_meth_read_ae$CpG <- as.vector(CpGlist)
mean_meth_read_ae$SNP <- as.vector(SNPlist)


load("merge_filter_573_chr3_af.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_af <- as.data.frame(mean_methylation_list)
mean_meth_read_af$mean_read <- mean_read_count_list
mean_meth_read_af$CpG <- as.vector(CpGlist)
mean_meth_read_af$SNP <- as.vector(SNPlist)


load("merge_filter_573_chr3_ag.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_ag <- as.data.frame(mean_methylation_list)
mean_meth_read_ag$mean_read <- mean_read_count_list
mean_meth_read_ag$CpG <- as.vector(CpGlist)
mean_meth_read_ag$SNP <- as.vector(SNPlist)


load("merge_filter_573_chr3_ah.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_ah <- as.data.frame(mean_methylation_list)
mean_meth_read_ah$mean_read <- mean_read_count_list
mean_meth_read_ah$CpG <- as.vector(CpGlist)
mean_meth_read_ah$SNP <- as.vector(SNPlist)


load("merge_filter_573_chr3_ai.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_ai <- as.data.frame(mean_methylation_list)
mean_meth_read_ai$mean_read <- mean_read_count_list
mean_meth_read_ai$CpG <- as.vector(CpGlist)
mean_meth_read_ai$SNP <- as.vector(SNPlist)


load("merge_filter_573_chr3_aj.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_aj <- as.data.frame(mean_methylation_list)
mean_meth_read_aj$mean_read <- mean_read_count_list
mean_meth_read_aj$CpG <- as.vector(CpGlist)
mean_meth_read_aj$SNP <- as.vector(SNPlist)

colnames(mean_meth_read_aa) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_ab) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_ac) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_ad) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_ae) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_af) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_ag) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_ah) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_ai) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_aj) <- c("mean_methylation", "mean_read", "CpG", "SNP")

mean_meth_read_chr3 <- rbind(mean_meth_read_aa,mean_meth_read_ab, mean_meth_read_ac, mean_meth_read_ad, mean_meth_read_ae,
mean_meth_read_af, mean_meth_read_ag, mean_meth_read_ah, mean_meth_read_ai, mean_meth_read_aj)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr3, mean_meth_read_chr3, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr3.txt")

#4
fdr_corrected_res_all_chr4 <- subset(fdr_corrected_res_all_chr, chr =="4")
load("merge_filter_573_chr4.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr4, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr4.txt")

#5
fdr_corrected_res_all_chr5 <- subset(fdr_corrected_res_all_chr, chr =="5")
load("merge_filter_573_chr5.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr5, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr5.txt")

#6
fdr_corrected_res_all_chr6 <- subset(fdr_corrected_res_all_chr, chr =="6")
load("merge_filter_573_chr6.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr6, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr6.txt")

#7
fdr_corrected_res_all_chr7 <- subset(fdr_corrected_res_all_chr, chr =="7")
load("merge_filter_573_chr7.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr7, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr7.txt")

#8
fdr_corrected_res_all_chr8 <- subset(fdr_corrected_res_all_chr, chr =="8")
load("merge_filter_573_chr8.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr8, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr8.txt")

#9
fdr_corrected_res_all_chr9 <- subset(fdr_corrected_res_all_chr, chr =="9")
load("merge_filter_573_chr9.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr9, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr9.txt")

#10
fdr_corrected_res_all_chr10 <- subset(fdr_corrected_res_all_chr, chr =="10")
load("merge_filter_573_chr10.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr10, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr10.txt")

#11
fdr_corrected_res_all_chr11 <- subset(fdr_corrected_res_all_chr, chr =="11")
load("merge_filter_573_chr11.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr11, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr11.txt")

#12
fdr_corrected_res_all_chr12 <- subset(fdr_corrected_res_all_chr, chr =="12")
load("merge_filter_573_chr12.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr12, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr12.txt")


#13
fdr_corrected_res_all_chr13 <- subset(fdr_corrected_res_all_chr, chr =="13")
load("merge_filter_573_chr13.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr13, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr13.txt")

#14
fdr_corrected_res_all_chr14 <- subset(fdr_corrected_res_all_chr, chr =="14")
load("merge_filter_573_chr14.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)


library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr14, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr14.txt")

#15
fdr_corrected_res_all_chr15 <- subset(fdr_corrected_res_all_chr, chr =="15")
load("merge_filter_573_chr15.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr15, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr15.txt")

#16
fdr_corrected_res_all_chr16 <- subset(fdr_corrected_res_all_chr, chr =="16")
load("merge_filter_573_chr16.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr16, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr16.txt")

#17
fdr_corrected_res_all_chr17 <- subset(fdr_corrected_res_all_chr, chr =="17")
load("merge_filter_573_chr17.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr17, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr17.txt")

#18
fdr_corrected_res_all_chr18 <- subset(fdr_corrected_res_all_chr, chr =="18")
load("merge_filter_chr18_573.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr18, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr18.txt")

#19
fdr_corrected_res_all_chr19 <- subset(fdr_corrected_res_all_chr, chr =="19")
load("merge_filter_573_chr19_aa.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_aa <- as.data.frame(mean_methylation_list)
mean_meth_read_aa$mean_read <- mean_read_count_list
mean_meth_read_aa$CpG <- as.vector(CpGlist)
mean_meth_read_aa$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr19_ab.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_ab <- as.data.frame(mean_methylation_list)
mean_meth_read_ab$mean_read <- mean_read_count_list
mean_meth_read_ab$CpG <- as.vector(CpGlist)
mean_meth_read_ab$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr19_ac.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_ac <- as.data.frame(mean_methylation_list)
mean_meth_read_ac$mean_read <- mean_read_count_list
mean_meth_read_ac$CpG <- as.vector(CpGlist)
mean_meth_read_ac$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr19_ad.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_ad <- as.data.frame(mean_methylation_list)
mean_meth_read_ad$mean_read <- mean_read_count_list
mean_meth_read_ad$CpG <- as.vector(CpGlist)
mean_meth_read_ad$SNP <- as.vector(SNPlist)

load("merge_filter_573_chr19_ae.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

mean_meth_read_ae <- as.data.frame(mean_methylation_list)
mean_meth_read_ae$mean_read <- mean_read_count_list
mean_meth_read_ae$CpG <- as.vector(CpGlist)
mean_meth_read_ae$SNP <- as.vector(SNPlist)

colnames(mean_meth_read_aa) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_ab) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_ac) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_ad) <- c("mean_methylation", "mean_read", "CpG", "SNP")
colnames(mean_meth_read_ae) <- c("mean_methylation", "mean_read", "CpG", "SNP")

mean_meth_read_chr19 <- rbind(mean_meth_read_aa,mean_meth_read_ab, mean_meth_read_ac, mean_meth_read_ad, mean_meth_read_ae)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr19, mean_meth_read_chr19, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr19.txt")

#20
fdr_corrected_res_all_chr20 <- subset(fdr_corrected_res_all_chr, chr =="20")
load("merge_filter_573_chr20.RData")

missing <- c()
n <- c()
mean_methylation <- c()
missing_reads <- c()
x <- c()
mean_read_count <- c()
mean_methylation_list <- c()
mean_read_count_list <- c()
  
for (i in 1:nrow(ym)){
  missing <- sum(is.nan(ym[i,]/rm[i,]) == TRUE)
  n <- 573 - missing 
  mean_methylation <- sum(ym[i,]/rm[i,], na.rm = TRUE)/n
  missing_reads <- sum(rm[i,] == "0")
  x <- 573 - missing_reads
  mean_read_count <- sum(rm[i,])/x
  mean_methylation_list <- c(mean_methylation_list, mean_methylation)
  mean_read_count_list <- c(mean_read_count_list, mean_read_count)
}

rownames_testing <- paste('1', CpGlist, sep='_')
rownames(test_chr1_aa_meth_levels) <- rownames_testing



mean_meth_read <- as.data.frame(mean_methylation_list)
mean_meth_read$mean_read <- mean_read_count_list
mean_meth_read$CpG <- as.vector(CpGlist)
mean_meth_read$SNP <- as.vector(SNPlist)

library(dplyr)
mean_meth_read_merge <- merge(fdr_corrected_res_all_chr20, mean_meth_read, by.x = c("CpG", "SNP"), by.y = c("CpG", "SNP"))
write.table(mean_meth_read_merge, file="fdr_corrected_res_mean_read_chr20.txt")