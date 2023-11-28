#!/usr/bin/Rscript
###---
#chr = ARRAY_ID 
# Macaque genes cpg sites within 50kb

library(dplyr)

read.table("Macaca_mulatta.Mmul_10.104.gtf", header=FALSE, sep='\t')

macaque_genes <- subset(macaque_gtf, V3 == "gene")

colnames(macaque_genes) <- c("chr","source","type","start","end","score","strand","phase","attributes")

positive_index <- which(macaque_genes$strand == "+")
negative_index <- which(macaque_genes$strand == "-")

macaque_genes$TSS <- 0

macaque_genes[positive_index,]$TSS <- macaque_genes[positive_index,]$start
macaque_genes[negative_index,]$TSS <- macaque_genes[negative_index,]$end

#Array ID value from the environment variable passed from sbatch (array for each chromosome: 1-20) 
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
n <- as.numeric(slurm_arrayid)

macaque_genes_chrn <- subset(macaque_genes, chr == n)

fdr_all_chr_res <- read.table("fdr_corrected_image_res_annotated_allchr_wchromhmm_FINAL.txt", header= TRUE)

#we need ALL sites to get the proportion 
fdr_all_chr_res_chrn <- subset(fdr_all_chr_res, chr == n)

gene_list <- c()
TSS_list <- c()
gene_associated_CpG <- c()
CpGs <- c()
TSSs <- c()
qvalue <- c()
qvalue_list <- c()
distance <- c()
distance_list <- c()

for (x in 1:nrow(macaque_genes_chrn)){
  for (i in 1:nrow(fdr_all_chr_res_chrn)){
    if ((abs(macaque_genes_chrn[x,]$TSS - fdr_all_chr_res_chrn[i,]$CpG)) < 50000) {
      distance <- abs(macaque_genes_chrn[x,]$TSS - fdr_all_chr_res_chrn[i,]$CpG)
      CpGs <- fdr_all_chr_res_chrn[i,]$CpG
      TSSs <- macaque_genes_chrn[x,]$TSS
      gene <- macaque_genes_chrn[x,]$attributes
      qvalue <- fdr_all_chr_res_chrn[i,]$fdr_FDR_IMAGE_res_all
      gene_associated_CpG <- c(gene_associated_CpG, CpGs)
      TSS_list <- c(TSS_list, TSSs)
      gene_list <- c(gene_list, gene) 
      qvalue_list <- c(qvalue_list, qvalue)
      distance_list <- c(distance_list, distance)
        }
  }
}

# Dataframe with the cpg sites and their associated gene  
# add a column/variable for distance 
chrn_cpg_tss_50kb <- as.data.frame(TSS_list)
chrn_cpg_tss_50kb$CpG <- gene_associated_CpG
chrn_cpg_tss_50kb$gene <- gene_list
chrn_cpg_tss_50kb$qvalue <- qvalue_list
chrn_cpg_tss_50kb$chr <- n
chrn_cpg_tss_50kb$distance <- distance_list

library(data.table)

#Filter by minimum distance to get nearest gene
filterbymindist_chrn <- setDT(chrn_cpg_tss_50kb)[chrn_cpg_tss_50kb[, .I[which.min(distance)], by=CpG]$V1]

#save 
write.table(filterbymindist_chrn, file= paste("chr", n, "cpg_TSS_50kb_nearestgene", ".txt", sep = "_")

#Get the empirical values in this script as an array
#then load in separate script
# and permute, 100 + times 

