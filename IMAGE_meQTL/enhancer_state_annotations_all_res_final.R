#!/usr/bin/Rscript

# Human PBMC Chromatin State Annotations

library(data.table)
fdr_corrected_res_annotate_all <- fread(file ="fdr_corrected_image_res_annotated_allchr.txt")
colnames(fdr_corrected_res_annotate_all[,1]) <- "V1.a"
colnames(fdr_corrected_res_annotate_all) <- make.unique(colnames(fdr_corrected_res_annotate_all))

# load chromatin state annotations (converted from human to macaque genome coordinates)
hglift_enhancers <- read.table("hglft_genome_274c5_4a3980.bed", header = FALSE, sep = "\t")
colnames(hglift_enhancers) <- c("chr", "enhancer_start", "enhancer_end")
#subset by chromosome

enhancers_macaque_chr1 <- subset(hglift_enhancers,chr =='chr1')
colnames(enhancers_macaque_chr1) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr2 <- subset(hglift_enhancers,chr =='chr2')
colnames(enhancers_macaque_chr2) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr3 <- subset(hglift_enhancers,chr =='chr3')
colnames(enhancers_macaque_chr3) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr4 <- subset(hglift_enhancers,chr =='chr4')
colnames(enhancers_macaque_chr4) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr5 <- subset(hglift_enhancers,chr =='chr5')
colnames(enhancers_macaque_chr5) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr6 <- subset(hglift_enhancers,chr =='chr6')
colnames(enhancers_macaque_chr6) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr7 <- subset(hglift_enhancers,chr =='chr7')
colnames(enhancers_macaque_chr7) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr8 <- subset(hglift_enhancers,chr =='chr8')
colnames(enhancers_macaque_chr8) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr9 <- subset(hglift_enhancers,chr =='chr9')
colnames(enhancers_macaque_chr9) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr10 <- subset(hglift_enhancers,chr =='chr10')
colnames(enhancers_macaque_chr10) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr11 <- subset(hglift_enhancers,chr =='chr11')
colnames(enhancers_macaque_chr11) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr12 <- subset(hglift_enhancers,chr =='chr12')
colnames(enhancers_macaque_chr12) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr13 <- subset(hglift_enhancers,chr =='chr13')
colnames(enhancers_macaque_chr13) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr14 <- subset(hglift_enhancers,chr =='chr14')
colnames(enhancers_macaque_chr14) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr15 <- subset(hglift_enhancers,chr =='chr15')
colnames(enhancers_macaque_chr15) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr16 <- subset(hglift_enhancers,chr =='chr16')
colnames(enhancers_macaque_chr16) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr17 <- subset(hglift_enhancers,chr =='chr17')
colnames(enhancers_macaque_chr17) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr18 <- subset(hglift_enhancers,chr =='chr18')
colnames(enhancers_macaque_chr18) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr19 <- subset(hglift_enhancers,chr =='chr19')
colnames(enhancers_macaque_chr19) <- c("chr", "enhancer_start", "enhancer_end", "state")
enhancers_macaque_chr20 <- subset(hglift_enhancers,chr =='chr20')
colnames(enhancers_macaque_chr20) <- c("chr", "enhancer_start", "enhancer_end", "state")

fdr_corrected_res_annotate_all$enhancer <- "NO"
fdr_corrected_res_annotate_all$enhancer_state <- "NO"

enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr1)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "1" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr1[i,]$enhancer_start : enhancers_macaque_chr1[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr1[i,]$state
  }
}

#chr2
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr2)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "2" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr2[i,]$enhancer_start : enhancers_macaque_chr2[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr2[i,]$state
  }
}

#chr3
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr3)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "3" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr3[i,]$enhancer_start : enhancers_macaque_chr3[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr3[i,]$state
  }
}

#chr4
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr4)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "4" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr4[i,]$enhancer_start : enhancers_macaque_chr4[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr4[i,]$state
  }
}

#chr5
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr5)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "5" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr5[i,]$enhancer_start : enhancers_macaque_chr5[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr5[i,]$state
  }
}

#chr6
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr6)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "6" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr6[i,]$enhancer_start : enhancers_macaque_chr6[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr6[i,]$state
  }
}

#chr7
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr7)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "7" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr7[i,]$enhancer_start : enhancers_macaque_chr7[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr7[i,]$state
  }
}

#chr8
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr8)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "8" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr8[i,]$enhancer_start : enhancers_macaque_chr8[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr8[i,]$state
  }
}

#chr9
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr9)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "9" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr9[i,]$enhancer_start : enhancers_macaque_chr9[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr9[i,]$state
  }
}

#chr10
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr10)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "10" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr10[i,]$enhancer_start : enhancers_macaque_chr10[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr10[i,]$state
  }
}

#chr11
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr11)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "11" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr11[i,]$enhancer_start : enhancers_macaque_chr11[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr11[i,]$state
  }
}

#chr12
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr12)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "12" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr12[i,]$enhancer_start : enhancers_macaque_chr12[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr12[i,]$state
  }
}

#chr13
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr13)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "13" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr13[i,]$enhancer_start : enhancers_macaque_chr13[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr13[i,]$state
  }
}

#chr14
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr14)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "14" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr14[i,]$enhancer_start : enhancers_macaque_chr14[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr14[i,]$state
  }
}

#chr15
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr15)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "15" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr15[i,]$enhancer_start : enhancers_macaque_chr15[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr15[i,]$state
  }
}

#chr16
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr16)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "16" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr16[i,]$enhancer_start : enhancers_macaque_chr16[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr16[i,]$state
  }
}

#chr17
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr17)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "17" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr17[i,]$enhancer_start : enhancers_macaque_chr17[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr17[i,]$state
  }
}

#chr18
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr18)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "18" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr18[i,]$enhancer_start : enhancers_macaque_chr18[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr18[i,]$state
  }
}

#chr19
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr19)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "19" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr19[i,]$enhancer_start : enhancers_macaque_chr19[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr19[i,]$state
  }
}

#chr20
enhancer_index <- c()
for (i in 1:nrow(enhancers_macaque_chr20)) {
  enhancer_index <- which(fdr_corrected_res_annotate_all$chr == "20" & fdr_corrected_res_annotate_all$CpG %in% (enhancers_macaque_chr20[i,]$enhancer_start : enhancers_macaque_chr20[i,]$enhancer_end) == TRUE)
  if  (sum(enhancer_index) > 0){
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer <- "YES"
    fdr_corrected_res_annotate_all[enhancer_index,]$enhancer_state <- enhancers_macaque_chr20[i,]$state
  }
}

write.table(fdr_corrected_res_annotate_all, file="fdr_corrected_image_res_annotated_allchr_wchromhmm.txt")


#odds ratio chromhmm
fdr_corrected_res_annotate_wchrom <- read.table("fdr_corrected_image_res_annotated_allchr_wchromhmm.txt", header= TRUE)
#subset data
fdr_corrected_res_annotate_wchrom_0.05 <- subset(fdr_corrected_res_annotate_wchrom, fdr_FDR_IMAGE_res_all < 0.05)
#and the sites that didn't pass threshold dataframe
fdr_corrected_res_annotate_wchrom_greater_0.05 <- subset(fdr_corrected_res_annotate_wchrom, fdr_FDR_IMAGE_res_all > 0.05)


#test for enrichment in specific chrom states
#1
enhancer_1_TssA_contingency_all <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state == "1_TssA"), sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state == "1_TssA")))
colnames(enhancer_1_TssA_contingency_all) <- "yes"
enhancer_1_TssA_contingency_all$no <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state != "1_TssA"),sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state != "1_TssA")))

fisher_test_enhancer_1_TssA_2side_all <- fisher.test(enhancer_1_TssA_contingency_all, alternative = "two.sided")
fisher_test_enhancer_1_TssA_2side_all$p.value
#pvalue 0
fisher.test(enhancer_1_TssA_contingency_all, alternative = "greater")
#2.2e-16
fisher.test(enhancer_1_TssA_contingency_all, alternative = "less")
#1
log2(fisher_test_enhancer_1_TssA_2side_all$estimate)
# 1.333831 
log2(fisher_test_enhancer_1_TssA_2side_all$conf.int)
# 1.304812 1.362986

#2
enhancer_2_TssAFlnk_contingency_all <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state == "2_TssAFlnk"), sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state == "2_TssAFlnk")))
colnames(enhancer_2_TssAFlnk_contingency_all) <- "yes"
enhancer_2_TssAFlnk_contingency_all$no <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state != "2_TssAFlnk"),sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state != "2_TssAFlnk")))

fisher_test_enhancer_2_TssAFlnk_2side_all <- fisher.test(enhancer_2_TssAFlnk_contingency_all, alternative = "two.sided")
fisher_test_enhancer_2_TssAFlnk_2side_all$p.value
#pvalue 0
fisher.test(enhancer_2_TssAFlnk_contingency_all, alternative = "greater")
#2.2e-16
fisher.test(enhancer_2_TssAFlnk_contingency_all, alternative = "less")
#1
log2(fisher_test_enhancer_2_TssAFlnk_2side_all$estimate)
#1.365902 
log2(fisher_test_enhancer_2_TssAFlnk_2side_all$conf.int)
#1.343348 1.388560

#3
enhancer_3_TxFlnk_contingency_all <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state == "3_TxFlnk"), sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state == "3_TxFlnk")))
colnames(enhancer_3_TxFlnk_contingency_all) <- "yes"
enhancer_3_TxFlnk_contingency_all$no <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state != "3_TxFlnk"),sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state != "3_TxFlnk")))

fisher_test_enhancer_3_TxFlnk_2side_all <- fisher.test(enhancer_3_TxFlnk_contingency_all, alternative = "two.sided")
fisher_test_enhancer_3_TxFlnk_2side_all$p.value
#pvalue 1.308716e-08
fisher.test(enhancer_3_TxFlnk_contingency_all, alternative = "greater")
#6.735e-09
fisher.test(enhancer_3_TxFlnk_contingency_all, alternative = "less")
#1
log2(fisher_test_enhancer_3_TxFlnk_2side_all$estimate)
# 0.5240784 
log2(fisher_test_enhancer_3_TxFlnk_2side_all$conf.int)
#0.3377847 0.7133249

#4
enhancer_4_Tx_contingency_all <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state == "4_Tx"), sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state == "4_Tx")))
colnames(enhancer_4_Tx_contingency_all) <- "yes"
enhancer_4_Tx_contingency_all$no <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state != "4_Tx"),sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state != "4_Tx")))

fisher_test_enhancer_4_Tx_2side_all <- fisher.test(enhancer_4_Tx_contingency_all, alternative = "two.sided")
fisher_test_enhancer_4_Tx_2side_all$p.value
#pvalue 0
fisher.test(enhancer_4_Tx_contingency_all, alternative = "greater")
#1
fisher.test(enhancer_4_Tx_contingency_all, alternative = "less")
#2.2e-16
log2(fisher_test_enhancer_4_Tx_2side_all$estimate)
#-0.9568238 
log2(fisher_test_enhancer_4_Tx_2side_all$conf.int)
#-0.9847425 -0.9288051

#5
enhancer_5_TxWk_contingency_all <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state == "5_TxWk"), sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state == "5_TxWk")))
colnames(enhancer_5_TxWk_contingency_all) <- "yes"
enhancer_5_TxWk_contingency_all$no <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state != "5_TxWk"),sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state != "5_TxWk")))

fisher_test_enhancer_5_TxWk_2side_all <- fisher.test(enhancer_5_TxWk_contingency_all, alternative = "two.sided")
fisher_test_enhancer_5_TxWk_2side_all$p.value
#pvalue 0
fisher.test(enhancer_5_TxWk_contingency_all, alternative = "greater")
#1
fisher.test(enhancer_5_TxWk_contingency_all, alternative = "less")
#2.2e-16
log2(fisher_test_enhancer_5_TxWk_2side_all$estimate)
#-0.7416501 
log2(fisher_test_enhancer_5_TxWk_2side_all$conf.int)
#-0.7666937 -0.7166285

#6
enhancer_6_EnhG_contingency_all <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state == "6_EnhG"), sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state == "6_EnhG")))
colnames(enhancer_6_EnhG_contingency_all) <- "yes"
enhancer_6_EnhG_contingency_all$no <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state != "6_EnhG"),sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state != "6_EnhG")))

fisher_test_enhancer_6_EnhG_2side_all <- fisher.test(enhancer_6_EnhG_contingency_all, alternative = "two.sided")
fisher_test_enhancer_6_EnhG_2side_all$p.value
#pvalue 7.423634e-16
fisher.test(enhancer_6_EnhG_contingency_all, alternative = "greater")
#1
fisher.test(enhancer_6_EnhG_contingency_all, alternative = "less")
#4.215e-16
log2(fisher_test_enhancer_6_EnhG_2side_all$estimate)
#-0.4080501 
log2(fisher_test_enhancer_6_EnhG_2side_all$conf.int)
#-0.5069316 -0.3089172

#7
enhancer_7_Enh_contingency_all <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state == "7_Enh"), sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state == "7_Enh")))
colnames(enhancer_7_Enh_contingency_all) <- "yes"
enhancer_7_Enh_contingency_all$no <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state != "7_Enh"),sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state != "7_Enh")))

fisher_test_enhancer_7_Enh_2side_all <- fisher.test(enhancer_7_Enh_contingency_all, alternative = "two.sided")
fisher_test_enhancer_7_Enh_2side_all$p.value
#pvalue 0.00305413
fisher.test(enhancer_7_Enh_contingency_all, alternative = "greater")
#0.001559
fisher.test(enhancer_7_Enh_contingency_all, alternative = "less")
# 0.9985
log2(fisher_test_enhancer_7_Enh_2side_all$estimate)
#0.08047034
log2(fisher_test_enhancer_7_Enh_2side_all$conf.int)
#0.02693053 0.13410156

#8
enhancer_8_ZNF_Rpts_contingency_all <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state == "8_ZNF/Rpts"), sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state == "8_ZNF/Rpts")))
colnames(enhancer_8_ZNF_Rpts_contingency_all) <- "yes"
enhancer_8_ZNF_Rpts_contingency_all$no <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state != "8_ZNF/Rpts"),sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state != "8_ZNF/Rpts")))

fisher_test_enhancer_8_ZNF_Rpts_2side_all <- fisher.test(enhancer_8_ZNF_Rpts_contingency_all, alternative = "two.sided")
fisher_test_enhancer_8_ZNF_Rpts_2side_all$p.value
#pvalue 0.0002046815
fisher.test(enhancer_8_ZNF_Rpts_contingency_all, alternative = "greater")
#0.0001077
fisher.test(enhancer_8_ZNF_Rpts_contingency_all, alternative = "less")
#0.9999
log2(fisher_test_enhancer_8_ZNF_Rpts_2side_all$estimate)
#  0.318204 
log2(fisher_test_enhancer_8_ZNF_Rpts_2side_all$conf.int)
#0.1471691 0.4915439

#9
enhancer_9_Het_contingency_all <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state == "9_Het"), sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state == "9_Het")))
colnames(enhancer_9_Het_contingency_all) <- "yes"
enhancer_9_Het_contingency_all$no <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state != "9_Het"),sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state != "9_Het")))

fisher_test_enhancer_9_Het_2side_all <- fisher.test(enhancer_9_Het_contingency_all, alternative = "two.sided")
fisher_test_enhancer_9_Het_2side_all$p.value
#pvalue 0.004459047
fisher.test(enhancer_9_Het_contingency_all, alternative = "greater")
#0.9979
fisher.test(enhancer_9_Het_contingency_all, alternative = "less")
#0.002248
log2(fisher_test_enhancer_9_Het_2side_all$estimate)
#-0.1005624 
log2(fisher_test_enhancer_9_Het_2side_all$conf.int)
#-0.16973344 -0.03113691

#10
enhancer_10_TssBiv_contingency_all <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state == "10_TssBiv"), sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state == "10_TssBiv")))
colnames(enhancer_10_TssBiv_contingency_all) <- "yes"
enhancer_10_TssBiv_contingency_all$no <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state != "10_TssBiv"),sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state != "10_TssBiv")))

fisher_test_enhancer_10_TssBiv_2side_all <- fisher.test(enhancer_10_TssBiv_contingency_all, alternative = "two.sided")
fisher_test_enhancer_10_TssBiv_2side_all$p.value
#pvalue 0
fisher.test(enhancer_10_TssBiv_contingency_all, alternative = "greater")
#2.2e-16
fisher.test(enhancer_10_TssBiv_contingency_all, alternative = "less")
#1
log2(fisher_test_enhancer_10_TssBiv_2side_all$estimate)
# 1.309642 
log2(fisher_test_enhancer_10_TssBiv_2side_all$conf.int)
#1.269139 1.350392

#11
enhancer_11_BivFlnk_contingency_all <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state == "11_BivFlnk"), sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state == "11_BivFlnk")))
colnames(enhancer_11_BivFlnk_contingency_all) <- "yes"
enhancer_11_BivFlnk_contingency_all$no <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state != "11_BivFlnk"),sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state != "11_BivFlnk")))

fisher_test_enhancer_11_BivFlnk_2side_all <- fisher.test(enhancer_11_BivFlnk_contingency_all, alternative = "two.sided")
fisher_test_enhancer_11_BivFlnk_2side_all$p.value
#pvalue 0
fisher.test(enhancer_11_BivFlnk_contingency_all, alternative = "greater")
#2.2e-16
fisher.test(enhancer_11_BivFlnk_contingency_all, alternative = "less")
#1
log2(fisher_test_enhancer_11_BivFlnk_2side_all$estimate)
#1.03818 
log2(fisher_test_enhancer_11_BivFlnk_2side_all$conf.int)
#0.9972854 1.0792835

#12
enhancer_12_EnhBiv_contingency_all <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state == "12_EnhBiv"), sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state == "12_EnhBiv")))
colnames(enhancer_12_EnhBiv_contingency_all) <- "yes"
enhancer_12_EnhBiv_contingency_all$no <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state != "12_EnhBiv"),sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state != "12_EnhBivk")))

fisher_test_enhancer_12_EnhBiv_2side_all <- fisher.test(enhancer_12_EnhBiv_contingency_all, alternative = "two.sided")
fisher_test_enhancer_12_EnhBiv_2side_all$p.value
#pvalue 1.603517e-198
fisher.test(enhancer_12_EnhBiv_contingency_all, alternative = "greater")
#2.2e-16
fisher.test(enhancer_12_EnhBiv_contingency_all, alternative = "less")
#1
log2(fisher_test_enhancer_12_EnhBiv_2side_all$estimate)
# 0.9852801 
log2(fisher_test_enhancer_12_EnhBiv_2side_all$conf.int)
#0.9168556 1.0541680

#13
enhancer_13_ReprPC_contingency_all <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state == "13_ReprPC"), sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state == "13_ReprPC")))
colnames(enhancer_13_ReprPC_contingency_all) <- "yes"
enhancer_13_ReprPC_contingency_all$no <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state != "13_ReprPC"),sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state != "13_ReprPC")))

fisher_test_enhancer_13_ReprPC_2side_all <- fisher.test(enhancer_13_ReprPC_contingency_all, alternative = "two.sided")
fisher_test_enhancer_13_ReprPC_2side_all$p.value
#pvalue 0
fisher.test(enhancer_13_ReprPC_contingency_all, alternative = "greater")
#2.2e-16
fisher.test(enhancer_13_ReprPC_contingency_all, alternative = "less")
#1
log2(fisher_test_enhancer_13_ReprPC_2side_all$estimate)
#0.7861157 
log2(fisher_test_enhancer_13_ReprPC_2side_all$conf.int)
#0.7561968 0.8160342

#14
enhancer_14_ReprPCWk_contingency_all <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state == "14_ReprPCWk"), sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state == "14_ReprPCWk")))
colnames(enhancer_14_ReprPCWk_contingency_all) <- "yes"
enhancer_14_ReprPCWk_contingency_all$no <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state != "14_ReprPCWk"),sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state != "14_ReprPCWk")))

fisher_test_enhancer_14_ReprPCWk_2side_all <- fisher.test(enhancer_14_ReprPCWk_contingency_all, alternative = "two.sided")
fisher_test_enhancer_14_ReprPCWk_2side_all$p.value
#pvalue 0
fisher.test(enhancer_14_ReprPCWk_contingency_all, alternative = "greater")
#1
fisher.test(enhancer_14_ReprPCWk_contingency_all, alternative = "less")
#2.2e-16
log2(fisher_test_enhancer_14_ReprPCWk_2side_all$estimate)
#-0.5457106
log2(fisher_test_enhancer_14_ReprPCWk_2side_all$conf.int)
#-0.5672917 -0.5240542

#15
enhancer_15_Quies_contingency_all <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state == "15_Quies"), sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state == "15_Quies")))
colnames(enhancer_15_Quies_contingency_all) <- "yes"
enhancer_15_Quies_contingency_all$no <- as.data.frame(c(sum(fdr_corrected_res_annotate_wchrom_0.05$enhancer_state != "15_Quies"),sum(fdr_corrected_res_annotate_wchrom_greater_0.05$enhancer_state != "15_Quies")))

fisher_test_enhancer_15_Quies_2side_all <- fisher.test(enhancer_15_Quies_contingency_all, alternative = "two.sided")
fisher_test_enhancer_15_Quies_2side_all$p.value
#pvalue 0
fisher.test(enhancer_15_Quies_contingency_all, alternative = "greater")
#1
fisher.test(enhancer_15_Quies_contingency_all, alternative = "less")
#2.2e-16
log2(fisher_test_enhancer_15_Quies_2side_all$estimate)
#-1.073735 
log2(fisher_test_enhancer_15_Quies_2side_all$conf.int)
#-1.088692 -1.058922

all_res_chromhmm_ann_cpg_log2_odds <- read.table("all_chr_res_chromhmm_logodds.txt", header=T, sep='\t')

ggplot(all_res_chromhmm_ann_cpg_log2_odds, aes(x = log2odds, y = type)) +
  geom_point(size = 4, aes(color=type)) + geom_errorbar(aes(xmax = log2upper, xmin = log2lower, color=type))
