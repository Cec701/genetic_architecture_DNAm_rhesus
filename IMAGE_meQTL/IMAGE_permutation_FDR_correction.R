#!/usr/bin/Rscript


install.packages("cobs", repos="http://cran.us.r-project.org")
install.packages("glmnet", repos="http://cran.us.r-project.org")
library(devtools)
install_github("janderson94/Revolution")
library(Revolution)

#Load IMAGE results
#chr1

image_res_chr1_aa <- read.table("24June_IMAGE_results_chr1_aa.txt", header = TRUE)

chr1_aa_CpG_SNP_index <- image_res_chr1_aa$Loc

chr1_aa_CpG <- read.table("24June_IMAGE_CpGs_chr1_aa.txt", header = TRUE)
chr1_aa_SNP <- read.table("24June_IMAGE_SNPs_chr1_aa.txt", header = TRUE)

image_res_chr1_aa$CpG <- chr1_aa_CpG[chr1_aa_CpG_SNP_index,]
image_res_chr1_aa$SNP <- chr1_aa_SNP[chr1_aa_CpG_SNP_index,]

image_res_chr1_ab <- read.table("24June_IMAGE_results_chr1_ab.txt", header = TRUE)

chr1_ab_CpG_SNP_index <- image_res_chr1_ab$Loc

chr1_ab_CpG <- read.table("24June_IMAGE_CpGs_chr1_ab.txt", header = TRUE)
chr1_ab_SNP <- read.table("24June_IMAGE_SNPs_chr1_ab.txt", header = TRUE)

image_res_chr1_ab$CpG <- chr1_ab_CpG[chr1_ab_CpG_SNP_index,]
image_res_chr1_ab$SNP <- chr1_ab_SNP[chr1_ab_CpG_SNP_index,]

image_res_chr1_ac <- read.table("24June_IMAGE_results_chr1_ac.txt", header = TRUE)

chr1_ac_CpG_SNP_index <- image_res_chr1_ac$Loc

chr1_ac_CpG <- read.table("24June_IMAGE_CpGs_chr1_ac.txt", header = TRUE)
chr1_ac_SNP <- read.table("24June_IMAGE_SNPs_chr1_ac.txt", header = TRUE)

image_res_chr1_ac$CpG <- chr1_ac_CpG[chr1_ac_CpG_SNP_index,]
image_res_chr1_ac$SNP <- chr1_ac_SNP[chr1_ac_CpG_SNP_index,]

image_res_chr1_ad <- read.table("24June_IMAGE_results_chr1_ad.txt", header = TRUE)

chr1_ad_CpG_SNP_index <- image_res_chr1_ad$Loc

chr1_ad_CpG <- read.table("24June_IMAGE_CpGs_chr1_ad.txt", header = TRUE)
chr1_ad_SNP <- read.table("24June_IMAGE_SNPs_chr1_ad.txt", header = TRUE)

image_res_chr1_ad$CpG <- chr1_ad_CpG[chr1_ad_CpG_SNP_index,]
image_res_chr1_ad$SNP <- chr1_ad_SNP[chr1_ad_CpG_SNP_index,]

image_res_chr1_ae <- read.table("24June_IMAGE_results_chr1_ae.txt", header = TRUE)

chr1_ae_CpG_SNP_index <- image_res_chr1_ae$Loc

chr1_ae_CpG <- read.table("24June_IMAGE_CpGs_chr1_ae.txt", header = TRUE)
chr1_ae_SNP <- read.table("24June_IMAGE_SNPs_chr1_ae.txt", header = TRUE)

image_res_chr1_ae$CpG <- chr1_ae_CpG[chr1_ae_CpG_SNP_index,]
image_res_chr1_ae$SNP <- chr1_ae_SNP[chr1_ae_CpG_SNP_index,]

image_res_chr1 <- rbind(image_res_chr1_aa, image_res_chr1_ab, image_res_chr1_ac, image_res_chr1_ad, image_res_chr1_ae)
image_res_chr1$chr <- "1"

#chr2

image_res_chr2 <- read.table("24May_IMAGE_results_chr2.txt", header=TRUE)
image_res_chr2$chr <- "2"

chr2_CpG_SNP_index <- image_res_chr2$Loc

#add cpg and snp list to image_res_chr2
chr2_CpG_image <- read.table("24May_IMAGE_CpGs_chr2.txt", header=T)
chr2_SNP_image <- read.table("24May_IMAGE_SNPs_chr2.txt", header=T)

image_res_chr2$CpG <- chr2_CpG_image[chr2_CpG_SNP_index,]
image_res_chr2$SNP <- chr2_SNP_image[chr2_CpG_SNP_index,]

#chr3
image_res_chr3_aa <- read.table("24June_IMAGE_results_chr3_aa.txt", header = TRUE)

chr3_aa_CpG_SNP_index <- image_res_chr3_aa$Loc

chr3_aa_CpG <- read.table("24June_IMAGE_CpGs_chr3_aa.txt", header = TRUE)
chr3_aa_SNP <- read.table("24June_IMAGE_SNPs_chr3_aa.txt", header = TRUE)

image_res_chr3_aa$CpG <- chr3_aa_CpG[chr3_aa_CpG_SNP_index,]
image_res_chr3_aa$SNP <- chr3_aa_SNP[chr3_aa_CpG_SNP_index,]

image_res_chr3_ab <- read.table("24June_IMAGE_results_chr3_ab.txt", header = TRUE)

chr3_ab_CpG_SNP_index <- image_res_chr3_ab$Loc

chr3_ab_CpG <- read.table("24June_IMAGE_CpGs_chr3_ab.txt", header = TRUE)
chr3_ab_SNP <- read.table("24June_IMAGE_SNPs_chr3_ab.txt", header = TRUE)

image_res_chr3_ab$CpG <- chr3_ab_CpG[chr3_ab_CpG_SNP_index,]
image_res_chr3_ab$SNP <- chr3_ab_SNP[chr3_ab_CpG_SNP_index,]

image_res_chr3_ac <- read.table("24June_IMAGE_results_chr3_ac.txt", header = TRUE)

chr3_ac_CpG_SNP_index <- image_res_chr3_ac$Loc

chr3_ac_CpG <- read.table("24June_IMAGE_CpGs_chr3_ac.txt", header = TRUE)
chr3_ac_SNP <- read.table("24June_IMAGE_SNPs_chr3_ac.txt", header = TRUE)

image_res_chr3_ac$CpG <- chr3_ac_CpG[chr3_ac_CpG_SNP_index,]
image_res_chr3_ac$SNP <- chr3_ac_SNP[chr3_ac_CpG_SNP_index,]

image_res_chr3_ad <- read.table("24June_IMAGE_results_chr3_ad.txt", header = TRUE)

chr3_ad_CpG_SNP_index <- image_res_chr3_ad$Loc

chr3_ad_CpG <- read.table("24June_IMAGE_CpGs_chr3_ad.txt", header = TRUE)
chr3_ad_SNP <- read.table("24June_IMAGE_SNPs_chr3_ad.txt", header = TRUE)

image_res_chr3_ad$CpG <- chr3_ad_CpG[chr3_ad_CpG_SNP_index,]
image_res_chr3_ad$SNP <- chr3_ad_SNP[chr3_ad_CpG_SNP_index,]

image_res_chr3_ae <- read.table("24June_IMAGE_results_chr3_ae.txt", header = TRUE)

chr3_ae_CpG_SNP_index <- image_res_chr3_ae$Loc

chr3_ae_CpG <- read.table("24June_IMAGE_CpGs_chr3_ae.txt", header = TRUE)
chr3_ae_SNP <- read.table("24June_IMAGE_SNPs_chr3_ae.txt", header = TRUE)

image_res_chr3_ae$CpG <- chr3_ae_CpG[chr3_ae_CpG_SNP_index,]
image_res_chr3_ae$SNP <- chr3_ae_SNP[chr3_ae_CpG_SNP_index,]

image_res_chr3_af <- read.table("24June_IMAGE_results_chr3_af.txt", header = TRUE)

chr3_af_CpG_SNP_index <- image_res_chr3_af$Loc

chr3_af_CpG <- read.table("24June_IMAGE_CpGs_chr3_af.txt", header = TRUE)
chr3_af_SNP <- read.table("24June_IMAGE_SNPs_chr3_af.txt", header = TRUE)

image_res_chr3_af$CpG <- chr3_af_CpG[chr3_af_CpG_SNP_index,]
image_res_chr3_af$SNP <- chr3_af_SNP[chr3_af_CpG_SNP_index,]

image_res_chr3_ag <- read.table("24June_IMAGE_results_chr3_ag.txt", header = TRUE)

chr3_ag_CpG_SNP_index <- image_res_chr3_ag$Loc

chr3_ag_CpG <- read.table("24June_IMAGE_CpGs_chr3_ag.txt", header = TRUE)
chr3_ag_SNP <- read.table("24June_IMAGE_SNPs_chr3_ag.txt", header = TRUE)

image_res_chr3_ag$CpG <- chr3_ag_CpG[chr3_ag_CpG_SNP_index,]
image_res_chr3_ag$SNP <- chr3_ag_SNP[chr3_ag_CpG_SNP_index,]

image_res_chr3_ah <- read.table("24June_IMAGE_results_chr3_ah.txt", header = TRUE)

chr3_ah_CpG_SNP_index <- image_res_chr3_ah$Loc

chr3_ah_CpG <- read.table("24June_IMAGE_CpGs_chr3_ah.txt", header = TRUE)
chr3_ah_SNP <- read.table("24June_IMAGE_SNPs_chr3_ah.txt", header = TRUE)

image_res_chr3_ah$CpG <- chr3_ah_CpG[chr3_ah_CpG_SNP_index,]
image_res_chr3_ah$SNP <- chr3_ah_SNP[chr3_ah_CpG_SNP_index,]

image_res_chr3_ai <- read.table("24June_IMAGE_results_chr3_ai.txt", header = TRUE)

chr3_ai_CpG_SNP_index <- image_res_chr3_ai$Loc

chr3_ai_CpG <- read.table("24June_IMAGE_CpGs_chr3_ai.txt", header = TRUE)
chr3_ai_SNP <- read.table("24June_IMAGE_SNPs_chr3_ai.txt", header = TRUE)

image_res_chr3_ai$CpG <- chr3_ai_CpG[chr3_ai_CpG_SNP_index,]
image_res_chr3_ai$SNP <- chr3_ai_SNP[chr3_ai_CpG_SNP_index,]

image_res_chr3_aj <- read.table("24June_IMAGE_results_chr3_aj.txt", header = TRUE)

chr3_aj_CpG_SNP_index <- image_res_chr3_aj$Loc

chr3_aj_CpG <- read.table("24June_IMAGE_CpGs_chr3_aj.txt", header = TRUE)
chr3_aj_SNP <- read.table("24June_IMAGE_SNPs_chr3_aj.txt", header = TRUE)

image_res_chr3_aj$CpG <- chr3_aj_CpG[chr3_aj_CpG_SNP_index,]
image_res_chr3_aj$SNP <- chr3_aj_SNP[chr3_aj_CpG_SNP_index,]

image_res_chr3 <- rbind(image_res_chr3_aa, image_res_chr3_ab, image_res_chr3_ac, image_res_chr3_ad, image_res_chr3_ae, 
image_res_chr3_af, image_res_chr3_ag, image_res_chr3_ah, image_res_chr3_ai, image_res_chr3_aj)
image_res_chr3$chr <- "3"

#chr4
image_res_chr4 <- read.table("24May_IMAGE_results_chr4.txt", header=TRUE)
image_res_chr4$chr <- "4"

chr4_CpG_SNP_index <- image_res_chr4$Loc

chr4_CpG_image <- read.table("24May_IMAGE_CpGs_chr4.txt", header=T)
chr4_SNP_image <- read.table("24May_IMAGE_SNPs_chr4.txt", header=T)

image_res_chr4$CpG <- chr4_CpG_image[chr4_CpG_SNP_index,]
image_res_chr4$SNP <- chr4_SNP_image[chr4_CpG_SNP_index,]

image_res_chr5 <- read.table("24May_IMAGE_results_chr5.txt", header=TRUE)
image_res_chr5$chr <- "5"

chr5_CpG_SNP_index <- image_res_chr5$Loc

chr5_CpG_image <- read.table("24May_IMAGE_CpGs_chr5.txt", header=T)
chr5_SNP_image <- read.table("24May_IMAGE_SNPs_chr5.txt", header=T)

image_res_chr5$CpG <- chr5_CpG_image[chr5_CpG_SNP_index,]
image_res_chr5$SNP <- chr5_SNP_image[chr5_CpG_SNP_index,]

image_res_chr6 <- read.table("24May_IMAGE_results_chr6.txt", header=TRUE)
image_res_chr6$chr <- "6"

chr6_CpG_SNP_index <- image_res_chr6$Loc

chr6_CpG_image <- read.table("24May_IMAGE_CpGs_chr6.txt", header=T)
chr6_SNP_image <- read.table("24May_IMAGE_SNPs_chr6.txt", header=T)

image_res_chr6$CpG <- chr6_CpG_image[chr6_CpG_SNP_index,]
image_res_chr6$SNP <- chr6_SNP_image[chr6_CpG_SNP_index,]

image_res_chr7 <- read.table("24May_IMAGE_results_chr7.txt", header=TRUE)
image_res_chr7$chr <- "7"

chr7_CpG_SNP_index <- image_res_chr7$Loc

chr7_CpG_image <- read.table("24May_IMAGE_CpGs_chr7.txt", header=T)
chr7_SNP_image <- read.table("24May_IMAGE_SNPs_chr7.txt", header=T)

image_res_chr7$CpG <- chr7_CpG_image[chr7_CpG_SNP_index,]
image_res_chr7$SNP <- chr7_SNP_image[chr7_CpG_SNP_index,]

image_res_chr8 <- read.table("20May_IMAGE_results_chr8.txt", header=TRUE)
image_res_chr8$chr <- "8"

chr8_CpG_SNP_index <- image_res_chr8$Loc

chr8_CpG_image <- read.table("20May_IMAGE_CpGs_chr8.txt", header=T)
chr8_SNP_image <- read.table("20May_IMAGE_SNPs_chr8.txt", header=T)

image_res_chr8$CpG <- chr8_CpG_image[chr8_CpG_SNP_index,]
image_res_chr8$SNP <- chr8_SNP_image[chr8_CpG_SNP_index,]

image_res_chr9 <- read.table("24May_IMAGE_results_chr9.txt", header=TRUE)
image_res_chr9$chr <- "9"

chr9_CpG_SNP_index <- image_res_chr9$Loc

chr9_CpG_image <- read.table("24May_IMAGE_CpGs_chr9.txt", header=T)
chr9_SNP_image <- read.table("24May_IMAGE_SNPs_chr9.txt", header=T)

image_res_chr9$CpG <- chr9_CpG_image[chr9_CpG_SNP_index,]
image_res_chr9$SNP <- chr9_SNP_image[chr9_CpG_SNP_index,]

image_res_chr10 <- read.table("24May_IMAGE_results_chr10.txt", header=TRUE)
image_res_chr10$chr <- "10"

chr10_CpG_SNP_index <- image_res_chr10$Loc

chr10_CpG_image <- read.table("24May_IMAGE_CpGs_chr10.txt", header=T)
chr10_SNP_image <- read.table("24May_IMAGE_SNPs_chr10.txt", header=T)

image_res_chr10$CpG <- chr10_CpG_image[chr10_CpG_SNP_index,]
image_res_chr10$SNP <- chr10_SNP_image[chr10_CpG_SNP_index,]

image_res_chr11 <- read.table("24May_IMAGE_results_chr11.txt", header=TRUE)
image_res_chr11$chr <- "11"

chr11_CpG_SNP_index <- image_res_chr11$Loc

chr11_CpG_image <- read.table("24May_IMAGE_CpGs_chr11.txt", header=T)
chr11_SNP_image <- read.table("24May_IMAGE_SNPs_chr11.txt", header=T)

image_res_chr11$CpG <- chr11_CpG_image[chr11_CpG_SNP_index,]
image_res_chr11$SNP <- chr11_SNP_image[chr11_CpG_SNP_index,]

image_res_chr12 <- read.table("24May_IMAGE_results_chr12.txt", header=TRUE)
image_res_chr12$chr <- "12"

chr12_CpG_SNP_index <- image_res_chr12$Loc

chr12_CpG_image <- read.table("24May_IMAGE_CpGs_chr12.txt", header=T)
chr12_SNP_image <- read.table("24May_IMAGE_SNPs_chr12.txt", header=T)

image_res_chr12$CpG <- chr12_CpG_image[chr12_CpG_SNP_index,]
image_res_chr12$SNP <- chr12_SNP_image[chr12_CpG_SNP_index,]

image_res_chr13 <- read.table("4June_IMAGE_results_chr13.txt", header=TRUE)
image_res_chr13$chr <- "13"

chr13_CpG_SNP_index <- image_res_chr13$Loc

chr13_CpG_image <- read.table("4June_IMAGE_CpGs_chr13.txt", header=T)
chr13_SNP_image <- read.table("4June_IMAGE_SNPs_chr13.txt", header=T)

image_res_chr13$CpG <- chr13_CpG_image[chr13_CpG_SNP_index,]
image_res_chr13$SNP <- chr13_SNP_image[chr13_CpG_SNP_index,]

image_res_chr14 <- read.table("4June_IMAGE_results_chr14.txt", header=TRUE)
image_res_chr14$chr <- "14"

chr14_CpG_SNP_index <- image_res_chr14$Loc

chr14_CpG_image <- read.table("4June_IMAGE_CpGs_chr14.txt", header=T)
chr14_SNP_image <- read.table("4June_IMAGE_SNPs_chr14.txt", header=T)

image_res_chr14$CpG <- chr14_CpG_image[chr14_CpG_SNP_index,]
image_res_chr14$SNP <- chr14_SNP_image[chr14_CpG_SNP_index,]

image_res_chr15 <- read.table("24May_IMAGE_results_chr15.txt", header=TRUE)
image_res_chr15$chr <- "15"

chr15_CpG_SNP_index <- image_res_chr15$Loc

chr15_CpG_image <- read.table("24May_IMAGE_CpGs_chr15.txt", header=T)
chr15_SNP_image <- read.table("24May_IMAGE_SNPs_chr15.txt", header=T)

image_res_chr15$CpG <- chr15_CpG_image[chr15_CpG_SNP_index,]
image_res_chr15$SNP <- chr15_SNP_image[chr15_CpG_SNP_index,]

image_res_chr16 <- read.table("24May_IMAGE_results_chr16.txt", header=TRUE)
image_res_chr16$chr <- "16"

chr16_CpG_SNP_index <- image_res_chr16$Loc

chr16_CpG_image <- read.table("24May_IMAGE_CpGs_chr16.txt", header=T)
chr16_SNP_image <- read.table("24May_IMAGE_SNPs_chr16.txt", header=T)

image_res_chr16$CpG <- chr16_CpG_image[chr16_CpG_SNP_index,]
image_res_chr16$SNP <- chr16_SNP_image[chr16_CpG_SNP_index,]

image_res_chr17 <- read.table("24May_IMAGE_results_chr17.txt", header=TRUE)
image_res_chr17$chr <- "17"

chr17_CpG_SNP_index <- image_res_chr17$Loc

chr17_CpG_image <- read.table("24May_IMAGE_CpGs_chr17.txt", header=T)
chr17_SNP_image <- read.table("24May_IMAGE_SNPs_chr17.txt", header=T)

image_res_chr17$CpG <- chr17_CpG_image[chr17_CpG_SNP_index,]
image_res_chr17$SNP <- chr17_SNP_image[chr17_CpG_SNP_index,]

image_res_chr18 <- read.table("10May_IMAGE_results_chr18.txt", header=TRUE)
image_res_chr18$chr <- "18"

chr18_CpG_SNP_index <- image_res_chr18$Loc

chr18_CpG_image <- read.table("10May_IMAGE_CpGs_chr18.txt", header=T)
chr18_SNP_image <- read.table("10May_IMAGE_SNPs_chr18.txt", header=T)

image_res_chr18$CpG <- chr18_CpG_image[chr18_CpG_SNP_index,]
image_res_chr18$SNP <- chr18_SNP_image[chr18_CpG_SNP_index,]

#chr19
image_res_chr19_aa <- read.table("24June_IMAGE_results_chr19_aa.txt", header = TRUE)

chr19_aa_CpG_SNP_index <- image_res_chr19_aa$Loc

chr19_aa_CpG <- read.table("24June_IMAGE_CpGs_chr19_aa.txt", header = TRUE)
chr19_aa_SNP <- read.table("24June_IMAGE_SNPs_chr19_aa.txt", header = TRUE)

image_res_chr19_aa$CpG <- chr19_aa_CpG[chr19_aa_CpG_SNP_index,]
image_res_chr19_aa$SNP <- chr19_aa_SNP[chr19_aa_CpG_SNP_index,]

image_res_chr19_ab <- read.table("24June_IMAGE_results_chr19_ab.txt", header = TRUE)

chr19_ab_CpG_SNP_index <- image_res_chr19_ab$Loc

chr19_ab_CpG <- read.table("24June_IMAGE_CpGs_chr19_ab.txt", header = TRUE)
chr19_ab_SNP <- read.table("24June_IMAGE_SNPs_chr19_ab.txt", header = TRUE)

image_res_chr19_ab$CpG <- chr19_ab_CpG[chr19_ab_CpG_SNP_index,]
image_res_chr19_ab$SNP <- chr19_ab_SNP[chr19_ab_CpG_SNP_index,]

image_res_chr19_ac <- read.table("24June_IMAGE_results_chr19_ac.txt", header = TRUE)

chr19_ac_CpG_SNP_index <- image_res_chr19_ac$Loc

chr19_ac_CpG <- read.table("24June_IMAGE_CpGs_chr19_ac.txt", header = TRUE)
chr19_ac_SNP <- read.table("24June_IMAGE_SNPs_chr19_ac.txt", header = TRUE)

image_res_chr19_ac$CpG <- chr19_ac_CpG[chr19_ac_CpG_SNP_index,]
image_res_chr19_ac$SNP <- chr19_ac_SNP[chr19_ac_CpG_SNP_index,]

image_res_chr19_ad <- read.table("24June_IMAGE_results_chr19_ad.txt", header = TRUE)

chr19_ad_CpG_SNP_index <- image_res_chr19_ad$Loc

chr19_ad_CpG <- read.table("24June_IMAGE_CpGs_chr19_ad.txt", header = TRUE)
chr19_ad_SNP <- read.table("24June_IMAGE_SNPs_chr19_ad.txt", header = TRUE)

image_res_chr19_ad$CpG <- chr19_ad_CpG[chr19_ad_CpG_SNP_index,]
image_res_chr19_ad$SNP <- chr19_ad_SNP[chr19_ad_CpG_SNP_index,]

image_res_chr19_ae <- read.table("24June_IMAGE_results_chr19_ae.txt", header = TRUE)

chr19_ae_CpG_SNP_index <- image_res_chr19_ae$Loc

chr19_ae_CpG <- read.table("24June_IMAGE_CpGs_chr19_ae.txt", header = TRUE)
chr19_ae_SNP <- read.table("24June_IMAGE_SNPs_chr19_ae.txt", header = TRUE)

image_res_chr19_ae$CpG <- chr19_ae_CpG[chr19_ae_CpG_SNP_index,]
image_res_chr19_ae$SNP <- chr19_ae_SNP[chr19_ae_CpG_SNP_index,]

image_res_chr19 <- rbind(image_res_chr19_aa, image_res_chr19_ab, image_res_chr19_ac, image_res_chr19_ad, image_res_chr19_ae)
image_res_chr19$chr <- "19"

#chr20

image_res_chr20 <- read.table("24May_IMAGE_results_chr20.txt", header=TRUE)
image_res_chr20$chr <- "20"

chr20_CpG_SNP_index <- image_res_chr20$Loc

chr20_CpG_image <- read.table("24May_IMAGE_CpGs_chr20.txt", header=T)
chr20_SNP_image <- read.table("24May_IMAGE_SNPs_chr20.txt", header=T)

image_res_chr20$CpG <- chr20_CpG_image[chr20_CpG_SNP_index,]
image_res_chr20$SNP <- chr20_SNP_image[chr20_CpG_SNP_index,]

#merge
input_df <- rbind(image_res_chr1, image_res_chr2, image_res_chr3, image_res_chr4, image_res_chr5, image_res_chr6, image_res_chr7,
image_res_chr8, image_res_chr9, image_res_chr10, image_res_chr11, image_res_chr12, image_res_chr13,
image_res_chr14, image_res_chr15, image_res_chr16, image_res_chr17, image_res_chr18, image_res_chr19, image_res_chr20)

input_df_converged <- subset(input_df,converged=='TRUE')

# Pvals_col_name = pvalue 

#the permutations
#set working directory to permutation files 
setwd(/scratch/ccosta8/RRBS/Output/ASM_all/permutations)

permutation_files <- list.files(pattern = "\\.txt$")
permutation_res <- lapply(permutation_files, function(x) {read.table(x, header=TRUE)["pvalue"]})
permutation_pvalues <- unlist(permutation_res)
perm_df <- as.data.frame(permutation_pvalues)
colnames(perm_df) <- "pvalue"


qqplot(-log10(perm_df$pvalue), -log10(image_res$pvalue), xlab = "-log10(Permuted p-values)", ylab= "-log10(Observed p-values IMAGE)")
hist(perm_df$pvalue, breaks=40, col="light blue", xlab="permuted p-value distribution IMAGE")
hist(image_res$pvalue, breaks=40, col="light blue", xlab="observes IMAGE p-value distribution")

Pvals_col_name <- "pvalue"

name <- "FDR_IMAGE_res_all"

#FDR correction
fdr_corrected_chr_res <- perm.fdr(input_df = input_df_converged, perm_df = perm_df, Pvals_col_name = Pvals_col_name, name = name)

#save 
write.table(fdr_corrected_chr_res, file ="fdr_corrected_chr_res_all.txt")
