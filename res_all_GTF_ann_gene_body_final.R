#!/usr/bin/Rscript
# Annotations - gene bodies
# Add to results

macaque_gtf <- read.table("Macaca_mulatta.Mmul_10.104.gtf", header=FALSE, sep='\t')

macaque_genes <- subset(macaque_gtf, V3 == "gene")

colnames(macaque_genes) <- c("chr","source","type","start","end","score","strand","phase","attributes")

positive_index <- which(macaque_genes$strand == "+")
negative_index <- which(macaque_genes$strand == "-")

macaque_genes$TSS <- 0

macaque_genes[positive_index,]$TSS <- macaque_genes[positive_index,]$start
macaque_genes[negative_index,]$TSS <- macaque_genes[negative_index,]$end

fdr_corrected_res_annotated <- read.table("fdr_corrected_res_all_gtf_chromhmm_annotations.txt")

#gene body
fdr_corrected_res_annotated$gene_body <- "NO"

#split by chr 
macaque_genes_chr1 <- subset(macaque_genes, chr == 1)
macaque_genes_chr2 <- subset(macaque_genes, chr == 2)
macaque_genes_chr3 <- subset(macaque_genes, chr == 3)
macaque_genes_chr4 <- subset(macaque_genes, chr == 4)
macaque_genes_chr5 <- subset(macaque_genes, chr == 5)
macaque_genes_chr6 <- subset(macaque_genes, chr == 6)
macaque_genes_chr7 <- subset(macaque_genes, chr == 7)
macaque_genes_chr8 <- subset(macaque_genes, chr == 8)
macaque_genes_chr9 <- subset(macaque_genes, chr == 9)
macaque_genes_chr10 <- subset(macaque_genes, chr == 10)
macaque_genes_chr11 <- subset(macaque_genes, chr == 11)
macaque_genes_chr12 <- subset(macaque_genes, chr == 12)
macaque_genes_chr13 <- subset(macaque_genes, chr == 13)
macaque_genes_chr14 <- subset(macaque_genes, chr == 14)
macaque_genes_chr15 <- subset(macaque_genes, chr == 15)
macaque_genes_chr16 <- subset(macaque_genes, chr == 16)
macaque_genes_chr17 <- subset(macaque_genes, chr == 17)
macaque_genes_chr18 <- subset(macaque_genes, chr == 18)
macaque_genes_chr19 <- subset(macaque_genes, chr == 19)
macaque_genes_chr20 <- subset(macaque_genes, chr == 20)


#chr1
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr1)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "1" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr1[i,]$start : macaque_genes_chr1[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr2
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr2)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "2" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr2[i,]$start : macaque_genes_chr2[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr3
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr3)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "3" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr3[i,]$start : macaque_genes_chr3[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr4
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr4)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "4" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr4[i,]$start : macaque_genes_chr4[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr5
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr5)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "5" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr5[i,]$start : macaque_genes_chr5[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr6
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr6)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "6" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr6[i,]$start : macaque_genes_chr6[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr7
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr7)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "7" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr6[i,]$start : macaque_genes_chr6[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr8
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr8)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "8" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr8[i,]$start : macaque_genes_chr8[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr9
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr9)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "9" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr9[i,]$start : macaque_genes_chr9[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr10
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr10)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "10" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr10[i,]$start : macaque_genes_chr10[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr11
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr11)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "11" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr11[i,]$start : macaque_genes_chr11[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr12
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr12)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "12" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr12[i,]$start : macaque_genes_chr12[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr13
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr13)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "13" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr13[i,]$start : macaque_genes_chr13[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr14
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr14)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "14" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr14[i,]$start : macaque_genes_chr14[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr15
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr15)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "15" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr15[i,]$start : macaque_genes_chr15[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr16
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr16)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "16" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr16[i,]$start : macaque_genes_chr16[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr17
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr17)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "17" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr17[i,]$start : macaque_genes_chr17[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr18
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr18)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "18" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr18[i,]$start : macaque_genes_chr18[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr19
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr19)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "19" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr19[i,]$start : macaque_genes_chr19[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

#chr20
gene_body_index <- c()
for (i in 1:nrow(macaque_genes_chr20)) {
  gene_body_index <- which(fdr_corrected_res_annotated$chr == "20" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_chr20[i,]$start : macaque_genes_chr20[i,]$end) == TRUE)
  if  (sum(gene_body_index) > 0){
    fdr_corrected_res_annotated[gene_body_index,]$gene_body <- "YES"
  }
}

write.table(fdr_corrected_res_annotated, file ="fdr_corrected_res_gtf_chromhmm_annotations_ALL.txt")