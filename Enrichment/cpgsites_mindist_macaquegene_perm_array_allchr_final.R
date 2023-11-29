#!/usr/bin/Rscript
###---
#all chr genes cpg site within 50kb (nearest gene)
#permutations as array (get permuted proportions)
#loaded from chr minimum distance table
library(dplyr)
library(data.table)

#load
filterbymindist_chr1 <- read.table("chr_1_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr2 <- read.table("chr_2_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr3 <- read.table("chr_3_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr4 <- read.table("chr_4_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr5 <- read.table("chr_5_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr6 <- read.table("chr_6_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr7 <- read.table("chr_7_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr8 <- read.table("chr_8_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr9 <- read.table("chr_9_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr10 <- read.table("chr_10_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr11 <- read.table("chr_11_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr12 <- read.table("chr_12_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr13 <- read.table("chr_13_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr14 <- read.table("chr_14_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr15 <- read.table("chr_15_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr16 <- read.table("chr_16_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr17 <- read.table("chr_17_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr18 <- read.table("chr_18_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr19 <- read.table("chr_19_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)
filterbymindist_chr20 <- read.table("chr_20_cpg_TSS_50kb_nearestgene_.txt", header = TRUE)


filterbymindist_allchr <- rbind(filterbymindist_chr1, filterbymindist_chr2, filterbymindist_chr3, filterbymindist_chr4, filterbymindist_chr5,
filterbymindist_chr6, filterbymindist_chr7, filterbymindist_chr8, filterbymindist_chr9, filterbymindist_chr10,
filterbymindist_chr11, filterbymindist_chr12, filterbymindist_chr13, filterbymindist_chr14, filterbymindist_chr15, 
filterbymindist_chr16, filterbymindist_chr17, filterbymindist_chr18, filterbymindist_chr19, filterbymindist_chr20)

#Permute
filterbymindist_allchr$qvalue <- sample(filterbymindist_allchr$qvalue, replace = FALSE)

#Input permuted data to get the proportions for each gene
# group by gene 
by_gene <- filterbymindist_allchr %>% group_by(gene)

#each gene gets a number of sig and number of nonsig
chr_cpg_TSS_50kb_sig_counts_perm <- by_gene %>% summarise(sig_count=sum((qvalue < 0.05) == TRUE ), non_sig_count =sum((qvalue > 0.05) == TRUE ))

# Filter out any genes that have less than 10 tested sites because it will sway the proportino too much
less_than_10_index <- which(chr_cpg_TSS_50kb_sig_counts_perm$sig_count + chr_cpg_TSS_50kb_sig_counts_perm$non_sig_count < 10)

#remove those 
chr_cpg_TSS_50kb_sig_counts_perm_greater10 <- chr_cpg_TSS_50kb_sig_counts_perm[-less_than_10_index,]
#save total sites tested near gene
chr_cpg_TSS_50kb_sig_counts_perm_greater10$total_sites <- chr_cpg_TSS_50kb_sig_counts_perm_greater10$sig_count + chr_cpg_TSS_50kb_sig_counts_perm_greater10$non_sig_count
#get the proportion significant
chr_cpg_TSS_50kb_sig_counts_perm_greater10$proportion_sig <- chr_cpg_TSS_50kb_sig_counts_perm_greater10$sig_count/chr_cpg_TSS_50kb_sig_counts_perm_greater10$total_sites

#grab the array id (100 permutations, array=1-100)
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
# coerce the value to an integer
n <- as.numeric(slurm_arrayid)

#then save the permuted proportions -
write.table(chr_cpg_TSS_50kb_sig_counts_perm_greater10, file = paste("allchr_cpg_TSS_50kb_sig_counts_prop_permutation", n, ".txt", sep = "_"))

