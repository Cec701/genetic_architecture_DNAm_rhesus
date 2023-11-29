#!/usr/bin/Rscript
#install.packages("sva")
install.packages("sva", repos = c("https://bioconductor.org/packages/3.12/bioc")) 
library(sva)

# SVA analysis (surrogate variable)
# load normalized gene expression data 
# make a design matrix with sex, age, and hurricane exposure (mod)

sorted_ge_wgs_SVA_matrix <- read.table(sorted_GE_wgs_SVA_matrix.txt, header=TRUE)
#load normalised gene count
load("norm_gene_counts_filtered.RData")

resid_GE_wgs_overlap_ordered <- read.table("resid_GE_wgs_overlap_ordered.txt", header = TRUE)
ordered_GE_wgs_LID <- colnames(resid_GE_wgs_overlap_ordered)
norm_gene_counts_filtered_GE_wgs_overlap_order <- norm_gene_counts_filtered[,ordered_GE_wgs_LID]

mod = model.matrix(~Sex + Chronological.age + Hurricane.exposure, data=sorted_ge_wgs_SVA_matrix)
mod0 = model.matrix(~1, data = sorted_ge_wgs_SVA_matrix)


sva1 = sva(norm_gene_counts_filtered_GE_wgs_overlap_order,mod,mod0,n.sv=15)

#...
#lm(gene expression ~ SVA1 + SVA2 + SVA3 +SVA4 + SVA5)

# only the fam files need to change
# because the bed files subset by snps within 200kb of that gene won't change

# tranpose the GE into a dataframe
# and filter through each gene 

norm_gene_counts_filtered_GE_wgs_overlap_order_transpose <- t(norm_gene_counts_filtered_GE_wgs_overlap_order)

df_norm_gene_counts_GE_wgs_order <- as.data.frame(norm_gene_counts_filtered_GE_wgs_overlap_order_transpose)
df_SVA_GE <- as.data.frame(sva1$sv)

norm_gene_counts_GE_wgs_order_wSVA <- cbind(df_norm_gene_counts_GE_wgs_order,df_SVA_GE)
sva1$sv
library(dplyr)


# 5 SVA
n <- (ncol(norm_gene_counts_GE_wgs_order_wSVA)-15)
# run n regressions
model_lms_ge_sva <- lapply(1:n, function(x) lm(norm_gene_counts_GE_wgs_order_wSVA[,x]~ V1 + V2 + V3 + V4 + V5, data = norm_gene_counts_GE_wgs_order_wSVA))
summaries_modellm <- lapply(model_lms_ge_sva, summary)
X_model_res <- lapply(summaries_modellm, function(x) x$res)

X_model_res_5_sva_df <- as.data.frame(do.call(cbind, X_model_res))
gene_ID_GE_ordered <- colnames(norm_gene_counts_filtered_GE_wgs_overlap_order_transpose)
colnames(X_model_res_5_sva_df) <- gene_ID_GE_ordered

#1SVA 

model_lms_ge_1_sva <- lapply(1:n, function(x) lm(norm_gene_counts_GE_wgs_order_wSVA[,x]~ V1, data = norm_gene_counts_GE_wgs_order_wSVA))
summaries_modellm_1_sva <- lapply(model_lms_ge_1_sva, summary)
X_model_res_1_sva <- lapply(summaries_modellm_1_sva, function(x) x$res)

#9 SVA

model_lms_ge_9_sva <- lapply(1:n, function(x) lm(norm_gene_counts_GE_wgs_order_wSVA[,x]~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9, data = norm_gene_counts_GE_wgs_order_wSVA))
summaries_modellm_9_sva <- lapply(model_lms_ge_9_sva, summary)
X_model_res_9_sva <- lapply(summaries_modellm_9_sva, function(x) x$res)


#10 SVA

model_lms_ge_10_sva <- lapply(1:n, function(x) lm(norm_gene_counts_GE_wgs_order_wSVA[,x]~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10, data = norm_gene_counts_GE_wgs_order_wSVA))
summaries_modellm_10_sva <- lapply(model_lms_ge_10_sva, summary)
X_model_res_10_sva <- lapply(summaries_modellm_10_sva, function(x) x$res)

#11 SVA


model_lms_ge_11_sva <- lapply(1:n, function(x) lm(norm_gene_counts_GE_wgs_order_wSVA[,x]~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11, data = norm_gene_counts_GE_wgs_order_wSVA))
summaries_modellm_11_sva <- lapply(model_lms_ge_11_sva, summary)
X_model_res_11_sva <- lapply(summaries_modellm_11_sva, function(x) x$res)

#15 SVA

model_lms_ge_15_sva <- lapply(1:n, function(x) lm(norm_gene_counts_GE_wgs_order_wSVA[,x]~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15, data = norm_gene_counts_GE_wgs_order_wSVA))
summaries_modellm_15_sva <- lapply(model_lms_ge_15_sva, summary)
X_model_res_15_sva <- lapply(summaries_modellm_15_sva, function(x) x$res)

# subset the GE residual dataframes (with SVA 1 regressed out, 5, 10, or 15) by chr
#SVA = 5 
#chr1

macaque_GE_genes_chr1_list <- macaque_GE_genes_chr1$ID

t_resid_GE_wgs_overlap_ordered_df_chr1_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr1_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr1_SVA_5, file = "resid_GE_wgs_chr1_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA = 10
X_model_res_10_sva_df <- as.data.frame(do.call(cbind, X_model_res_10_sva))
colnames(X_model_res_10_sva_df) <- gene_ID_GE_ordered

t_resid_GE_wgs_overlap_ordered_df_chr1_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr1_list]
write.table(t_resid_GE_wgs_overlap_ordered_df_chr1_SVA_10, file = "resid_GE_wgs_chr1_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)

#SVA = 1
X_model_res_1_sva_df <-  as.data.frame(do.call(cbind, X_model_res_1_sva))
colnames(X_model_res_1_sva_df) <- gene_ID_GE_ordered
t_resid_GE_wgs_overlap_ordered_df_chr1_SVA_1 <- X_model_res_1_sva_df[,macaque_GE_genes_chr1_list]
write.table(t_resid_GE_wgs_overlap_ordered_df_chr1_SVA_1, file = "resid_GE_wgs_chr1_ordered_SVA_1.txt", row.names = FALSE, col.names = FALSE)

###
#SVA =  15
#X_model_res_15_sva

X_model_res_15_sva_df <-  as.data.frame(do.call(cbind, X_model_res_15_sva))
colnames(X_model_res_15_sva_df) <- gene_ID_GE_ordered
t_resid_GE_wgs_overlap_ordered_df_chr1_SVA_15 <- X_model_res_15_sva_df[,macaque_GE_genes_chr1_list]
write.table(t_resid_GE_wgs_overlap_ordered_df_chr1_SVA_15, file = "resid_GE_wgs_chr1_ordered_SVA_15.txt", row.names = FALSE, col.names = FALSE)

#SVA=11
X_model_res_11_sva_df <- as.data.frame(do.call(cbind, X_model_res_11_sva))
colnames(X_model_res_11_sva_df) <- gene_ID_GE_ordered

t_resid_GE_wgs_overlap_ordered_df_chr1_SVA_11 <- X_model_res_11_sva_df[,macaque_GE_genes_chr1_list]
write.table(t_resid_GE_wgs_overlap_ordered_df_chr1_SVA_11, file = "resid_GE_wgs_chr1_ordered_SVA_11.txt", row.names = FALSE, col.names = FALSE)

#SVA=9
X_model_res_9_sva_df <- as.data.frame(do.call(cbind, X_model_res_9_sva))
colnames(X_model_res_9_sva_df) <- gene_ID_GE_ordered

t_resid_GE_wgs_overlap_ordered_df_chr1_SVA_9 <- X_model_res_9_sva_df[,macaque_GE_genes_chr1_list]
write.table(t_resid_GE_wgs_overlap_ordered_df_chr1_SVA_9, file = "resid_GE_wgs_chr1_ordered_SVA_9.txt", row.names = FALSE, col.names = FALSE)


#chr 2

macaque_GE_genes_chr2_list <- macaque_GE_genes_chr2$ID

t_resid_GE_wgs_overlap_ordered_df_chr2_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr2_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr2_SVA_5, file = "resid_GE_wgs_chr2_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)


#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr2_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr2_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr2_SVA_10, file = "resid_GE_wgs_chr2_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)


#chr 3

macaque_GE_genes_chr3_list <- macaque_GE_genes_chr3$ID

t_resid_GE_wgs_overlap_ordered_df_chr3_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr3_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr3_SVA_5, file = "resid_GE_wgs_chr3_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr3_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr3_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr3_SVA_10, file = "resid_GE_wgs_chr3_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)

#chr 4

macaque_GE_genes_chr4_list <- macaque_GE_genes_chr4$ID

t_resid_GE_wgs_overlap_ordered_df_chr4_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr4_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr4_SVA_5, file = "resid_GE_wgs_chr4_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr4_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr4_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr4_SVA_10, file = "resid_GE_wgs_chr4_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)


#chr 5

macaque_GE_genes_chr5_list <- macaque_GE_genes_chr5$ID

t_resid_GE_wgs_overlap_ordered_df_chr5_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr5_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr5_SVA_5, file = "resid_GE_wgs_chr5_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr5_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr5_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr5_SVA_10, file = "resid_GE_wgs_chr5_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)


#chr 6

macaque_GE_genes_chr6_list <- macaque_GE_genes_chr6$ID

t_resid_GE_wgs_overlap_ordered_df_chr6_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr6_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr6_SVA_5, file = "resid_GE_wgs_chr6_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr6_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr6_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr6_SVA_10, file = "resid_GE_wgs_chr6_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)


#chr 7

macaque_GE_genes_chr7_list <- macaque_GE_genes_chr7$ID

t_resid_GE_wgs_overlap_ordered_df_chr7_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr7_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr7_SVA_5, file = "resid_GE_wgs_chr7_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr7_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr7_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr7_SVA_10, file = "resid_GE_wgs_chr7_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)

#chr 8

macaque_GE_genes_chr8_list <- macaque_GE_genes_chr8$ID

t_resid_GE_wgs_overlap_ordered_df_chr8_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr8_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr8_SVA_5, file = "resid_GE_wgs_chr8_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr8_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr8_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr8_SVA_10, file = "resid_GE_wgs_chr8_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)

#chr 9

macaque_GE_genes_chr9_list <- macaque_GE_genes_chr9$ID

t_resid_GE_wgs_overlap_ordered_df_chr9_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr9_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr9_SVA_5, file = "resid_GE_wgs_chr9_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr9_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr9_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr9_SVA_10, file = "resid_GE_wgs_chr9_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)


#chr 10

macaque_GE_genes_chr10_list <- macaque_GE_genes_chr10$ID

t_resid_GE_wgs_overlap_ordered_df_chr10_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr10_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr10_SVA_5, file = "resid_GE_wgs_chr10_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr10_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr10_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr10_SVA_10, file = "resid_GE_wgs_chr10_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)

#chr 11

macaque_GE_genes_chr11_list <- macaque_GE_genes_chr11$ID

t_resid_GE_wgs_overlap_ordered_df_chr11_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr11_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr11_SVA_5, file = "resid_GE_wgs_chr11_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr11_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr11_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr11_SVA_10, file = "resid_GE_wgs_chr11_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)


#chr 12

macaque_GE_genes_chr12_list <- macaque_GE_genes_chr12$ID

t_resid_GE_wgs_overlap_ordered_df_chr12_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr12_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr12_SVA_5, file = "resid_GE_wgs_chr12_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr12_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr12_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr12_SVA_10, file = "resid_GE_wgs_chr12_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)


#chr 13

macaque_GE_genes_chr13_list <- macaque_GE_genes_chr13$ID

t_resid_GE_wgs_overlap_ordered_df_chr13_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr13_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr13_SVA_5, file = "resid_GE_wgs_chr13_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr13_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr13_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr13_SVA_10, file = "resid_GE_wgs_chr13_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)

#chr14 

macaque_GE_genes_chr14_list <- macaque_GE_genes_chr14$ID

t_resid_GE_wgs_overlap_ordered_df_chr14_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr14_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr14_SVA_5, file = "resid_GE_wgs_chr14_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr14_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr14_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr14_SVA_10, file = "resid_GE_wgs_chr14_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)

#chr15

macaque_GE_genes_chr15_list <- macaque_GE_genes_chr15$ID

t_resid_GE_wgs_overlap_ordered_df_chr15_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr15_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr15_SVA_5, file = "resid_GE_wgs_chr15_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr15_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr15_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr15_SVA_10, file = "resid_GE_wgs_chr15_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)

#chr16

macaque_GE_genes_chr16_list <- macaque_GE_genes_chr16$ID

t_resid_GE_wgs_overlap_ordered_df_chr16_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr16_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr16_SVA_5, file = "resid_GE_wgs_chr16_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr16_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr16_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr16_SVA_10, file = "resid_GE_wgs_chr16_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)

#chr17
macaque_GE_genes_chr17 <- subset(macaque_GE_genes, chr == 17)
macaque_GE_genes_chr17_list <- macaque_GE_genes_chr17$ID

t_resid_GE_wgs_overlap_ordered_df_chr17_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr17_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr17_SVA_5, file = "resid_GE_wgs_chr17_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr17_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr17_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr17_SVA_10, file = "resid_GE_wgs_chr17_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)

#chr18

macaque_GE_genes_chr18_list <- macaque_GE_genes_chr18$ID

t_resid_GE_wgs_overlap_ordered_df_chr18_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr18_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr18_SVA_5, file = "resid_GE_wgs_chr18_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr18_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr18_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr18_SVA_10, file = "resid_GE_wgs_chr18_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)

#chr19

macaque_GE_genes_chr19_list <- macaque_GE_genes_chr19$ID

t_resid_GE_wgs_overlap_ordered_df_chr19_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr19_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr19_SVA_5, file = "resid_GE_wgs_chr19_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)
#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr19_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr19_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr19_SVA_10, file = "resid_GE_wgs_chr19_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)


#chr20
t_resid_GE_wgs_overlap_ordered_df_chr20_SVA_5[,1:5]


macaque_GE_genes_chr20_list <- macaque_GE_genes_chr20$ID

t_resid_GE_wgs_overlap_ordered_df_chr20_SVA_5 <- X_model_res_5_sva_df[,macaque_GE_genes_chr20_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr20_SVA_5, file = "resid_GE_wgs_chr20_ordered_SVA_5.txt", row.names = FALSE, col.names = FALSE)

#SVA10

t_resid_GE_wgs_overlap_ordered_df_chr20_SVA_10 <- X_model_res_10_sva_df[,macaque_GE_genes_chr20_list]

write.table(t_resid_GE_wgs_overlap_ordered_df_chr20_SVA_10, file = "resid_GE_wgs_chr20_ordered_SVA_10.txt", row.names = FALSE, col.names = FALSE)

