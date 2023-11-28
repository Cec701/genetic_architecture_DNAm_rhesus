#!/usr/bin/Rscript
# annotations - promoters, exons
 
macaque_gtf <- read.table("Macaca_mulatta.Mmul_10.104.gtf", header=FALSE, sep='\t')

macaque_genes <- subset(macaque_gtf, V3 == "gene")

colnames(macaque_genes) <- c("chr","source","type","start","end","score","strand","phase","attributes")

positive_index <- which(macaque_genes$strand == "+")
negative_index <- which(macaque_genes$strand == "-")

macaque_genes$TSS <- 0

macaque_genes[positive_index,]$TSS <- macaque_genes[positive_index,]$start
macaque_genes[negative_index,]$TSS <- macaque_genes[negative_index,]$end

macaque_genes$promoter_start <- 0

macaque_genes[positive_index,]$promoter_start <- macaque_genes[positive_index,]$TSS - 2000
macaque_genes[negative_index,]$promoter_start <- macaque_genes[negative_index,]$TSS + 2000

#actually need to subset positive and negative strands then
macaque_genes_positive <- macaque_genes[positive_index,]
macaque_genes_negative <- macaque_genes[negative_index,]

#split by chr positive 
macaque_genes_positive_chr1 <- subset(macaque_genes_positive, chr == 1)
macaque_genes_positive_chr2 <- subset(macaque_genes_positive, chr == 2)
macaque_genes_positive_chr3 <- subset(macaque_genes_positive, chr == 3)
macaque_genes_positive_chr4 <- subset(macaque_genes_positive, chr == 4)
macaque_genes_positive_chr5 <- subset(macaque_genes_positive, chr == 5)
macaque_genes_positive_chr6 <- subset(macaque_genes_positive, chr == 6)
macaque_genes_positive_chr7 <- subset(macaque_genes_positive, chr == 7)
macaque_genes_positive_chr8 <- subset(macaque_genes_positive, chr == 8)
macaque_genes_positive_chr9 <- subset(macaque_genes_positive, chr == 9)
macaque_genes_positive_chr10 <- subset(macaque_genes_positive, chr == 10)
macaque_genes_positive_chr11 <- subset(macaque_genes_positive, chr == 11)
macaque_genes_positive_chr12 <- subset(macaque_genes_positive, chr == 12)
macaque_genes_positive_chr13 <- subset(macaque_genes_positive, chr == 13)
macaque_genes_positive_chr14 <- subset(macaque_genes_positive, chr == 14)
macaque_genes_positive_chr15 <- subset(macaque_genes_positive, chr == 15)
macaque_genes_positive_chr16 <- subset(macaque_genes_positive, chr == 16)
macaque_genes_positive_chr17 <- subset(macaque_genes_positive, chr == 17)
macaque_genes_positive_chr18 <- subset(macaque_genes_positive, chr == 18)
macaque_genes_positive_chr19 <- subset(macaque_genes_positive, chr == 19)
macaque_genes_positive_chr20 <- subset(macaque_genes_positive, chr == 20)

#split by chr negative 
macaque_genes_negative_chr1 <- subset(macaque_genes_negative, chr == 1)
macaque_genes_negative_chr2 <- subset(macaque_genes_negative, chr == 2)
macaque_genes_negative_chr3 <- subset(macaque_genes_negative, chr == 3)
macaque_genes_negative_chr4 <- subset(macaque_genes_negative, chr == 4)
macaque_genes_negative_chr5 <- subset(macaque_genes_negative, chr == 5)
macaque_genes_negative_chr6 <- subset(macaque_genes_negative, chr == 6)
macaque_genes_negative_chr7 <- subset(macaque_genes_negative, chr == 7)
macaque_genes_negative_chr8 <- subset(macaque_genes_negative, chr == 8)
macaque_genes_negative_chr9 <- subset(macaque_genes_negative, chr == 9)
macaque_genes_negative_chr10 <- subset(macaque_genes_negative, chr == 10)
macaque_genes_negative_chr11 <- subset(macaque_genes_negative, chr == 11)
macaque_genes_negative_chr12 <- subset(macaque_genes_negative, chr == 12)
macaque_genes_negative_chr13 <- subset(macaque_genes_negative, chr == 13)
macaque_genes_negative_chr14 <- subset(macaque_genes_negative, chr == 14)
macaque_genes_negative_chr15 <- subset(macaque_genes_negative, chr == 15)
macaque_genes_negative_chr16 <- subset(macaque_genes_negative, chr == 16)
macaque_genes_negative_chr17 <- subset(macaque_genes_negative, chr == 17)
macaque_genes_negative_chr18 <- subset(macaque_genes_negative, chr == 18)
macaque_genes_negative_chr19 <- subset(macaque_genes_negative, chr == 19)
macaque_genes_negative_chr20 <- subset(macaque_genes_negative, chr == 20)

fdr_corrected_res_annotated <- read.table("fdr_corrected_image_res_annotated_allchr_wchromhmm.txt")

fdr_corrected_res_annotated <- fdr_corrected_res_annotated[,-c(1,3)]

colnames(fdr_corrected_res_annotated)[colnames(fdr_corrected_res_annotated) == 'V1.1'] <- 'Loc'

write.table(fdr_corrected_res_annotated, file= "fdr_corrected_image_res_annotated_allchr_wchromhmm_FINAL.txt")

fdr_corrected_res_annotated$promoter <- "NO"

#positive
promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr1)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "1" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr1[i,]$promoter_start : macaque_genes_positive_chr1[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr1)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "1" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr1[i,]$TSS : macaque_genes_negative_chr1[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}


#chr2

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr2)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "2" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr2[i,]$promoter_start : macaque_genes_positive_chr2[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr2)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "2" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr2[i,]$TSS : macaque_genes_negative_chr2[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr3

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr3)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "3" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr3[i,]$promoter_start : macaque_genes_positive_chr3[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr3)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "3" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr3[i,]$TSS : macaque_genes_negative_chr3[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr4

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr4)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "4" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr4[i,]$promoter_start : macaque_genes_positive_chr4[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr4)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "4" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr4[i,]$TSS : macaque_genes_negative_chr4[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr5

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr5)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "5" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr5[i,]$promoter_start : macaque_genes_positive_chr5[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr5)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "5" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr5[i,]$TSS : macaque_genes_negative_chr5[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr6

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr6)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "6" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr6[i,]$promoter_start : macaque_genes_positive_chr6[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr6)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "6" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr6[i,]$TSS : macaque_genes_negative_chr6[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr7

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr7)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "7" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr7[i,]$promoter_start : macaque_genes_positive_chr7[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr7)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "7" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr7[i,]$TSS : macaque_genes_negative_chr7[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr8

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr8)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "8" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr8[i,]$promoter_start : macaque_genes_positive_chr8[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr8)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "8" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr8[i,]$TSS : macaque_genes_negative_chr8[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr9

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr9)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "9" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr9[i,]$promoter_start : macaque_genes_positive_chr9[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr9)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "9" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr9[i,]$TSS : macaque_genes_negative_chr9[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr10

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr10)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "10" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr10[i,]$promoter_start : macaque_genes_positive_chr10[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr10)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "10" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr10[i,]$TSS : macaque_genes_negative_chr10[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr11

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr11)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "11" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr11[i,]$promoter_start : macaque_genes_positive_chr11[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr11)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "11" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr11[i,]$TSS : macaque_genes_negative_chr11[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr12

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr12)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "12" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr12[i,]$promoter_start : macaque_genes_positive_chr12[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr12)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "12" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr12[i,]$TSS : macaque_genes_negative_chr12[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr13

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr13)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "13" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr13[i,]$promoter_start : macaque_genes_positive_chr13[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr13)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "13" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr13[i,]$TSS : macaque_genes_negative_chr13[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr14

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr14)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "14" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr14[i,]$promoter_start : macaque_genes_positive_chr14[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr14)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "14" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr14[i,]$TSS : macaque_genes_negative_chr14[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr15

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr15)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "15" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr15[i,]$promoter_start : macaque_genes_positive_chr15[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr15)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "15" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr15[i,]$TSS : macaque_genes_negative_chr15[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr16

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr16)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "16" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr16[i,]$promoter_start : macaque_genes_positive_chr16[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr16)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "16" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr16[i,]$TSS : macaque_genes_negative_chr16[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr17

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr17)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "17" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr17[i,]$promoter_start : macaque_genes_positive_chr17[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr17)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "17" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr17[i,]$TSS : macaque_genes_negative_chr17[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr18

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr18)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "18" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr18[i,]$promoter_start : macaque_genes_positive_chr18[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr18)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "18" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr18[i,]$TSS : macaque_genes_negative_chr18[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr19

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr19)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "19" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr19[i,]$promoter_start : macaque_genes_positive_chr19[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr19)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "19" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr19[i,]$TSS : macaque_genes_negative_chr19[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#chr20

promoter_index_1 <- c()
for (i in 1:nrow(macaque_genes_positive_chr20)) {
  promoter_index_1 <- which(fdr_corrected_res_annotated$chr == "20" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_positive_chr20[i,]$promoter_start : macaque_genes_positive_chr20[i,]$TSS) == TRUE)
  if  (sum(promoter_index_1) > 0){
    fdr_corrected_res_annotated[promoter_index_1,]$promoter <- "YES"
  }
}


#negative 
promoter_index_2 <- c()
for (i in 1:nrow(macaque_genes_negative_chr20)) {
  promoter_index_2 <- which(fdr_corrected_res_annotated$chr == "20" & fdr_corrected_res_annotated$CpG %in% (macaque_genes_negative_chr20[i,]$TSS : macaque_genes_negative_chr20[i,]$promoter_start) == TRUE)
  if  (sum(promoter_index_2) > 0){
    fdr_corrected_res_annotated[promoter_index_2,]$promoter <- "YES"
  }
}

#exons

#let's look for sites in exons (features) and gene bodies (gene start:gene end)

macaque_exons <- subset(macaque_gtf, V3 == "exon")

colnames(macaque_exons) <- c("chr","source","type","start","end","score","strand","phase","attributes")

fdr_corrected_res_annotated$exon <- "NO"

#split by chr 
macaque_exons_chr1 <- subset(macaque_exons, chr == 1)
macaque_exons_chr2 <- subset(macaque_exons, chr == 2)
macaque_exons_chr3 <- subset(macaque_exons, chr == 3)
macaque_exons_chr4 <- subset(macaque_exons, chr == 4)
macaque_exons_chr5 <- subset(macaque_exons, chr == 5)
macaque_exons_chr6 <- subset(macaque_exons, chr == 6)
macaque_exons_chr7 <- subset(macaque_exons, chr == 7)
macaque_exons_chr8 <- subset(macaque_exons, chr == 8)
macaque_exons_chr9 <- subset(macaque_exons, chr == 9)
macaque_exons_chr10 <- subset(macaque_exons, chr == 10)
macaque_exons_chr11 <- subset(macaque_exons, chr == 11)
macaque_exons_chr12 <- subset(macaque_exons, chr == 12)
macaque_exons_chr13 <- subset(macaque_exons, chr == 13)
macaque_exons_chr14 <- subset(macaque_exons, chr == 14)
macaque_exons_chr15 <- subset(macaque_exons, chr == 15)
macaque_exons_chr16 <- subset(macaque_exons, chr == 16)
macaque_exons_chr17 <- subset(macaque_exons, chr == 17)
macaque_exons_chr18 <- subset(macaque_exons, chr == 18)
macaque_exons_chr19 <- subset(macaque_exons, chr == 19)
macaque_exons_chr20 <- subset(macaque_exons, chr == 20)

#chr1
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr1)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "1" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr1[i,]$start : macaque_exons_chr1[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr2
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr2)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "2" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr2[i,]$start : macaque_exons_chr2[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr3
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr3)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "3" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr3[i,]$start : macaque_exons_chr3[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr4
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr4)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "4" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr4[i,]$start : macaque_exons_chr4[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr5
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr5)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "5" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr5[i,]$start : macaque_exons_chr5[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr6
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr6)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "6" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr6[i,]$start : macaque_exons_chr6[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr7
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr7)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "7" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr7[i,]$start : macaque_exons_chr7[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr8
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr8)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "8" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr8[i,]$start : macaque_exons_chr8[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr9
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr9)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "9" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr9[i,]$start : macaque_exons_chr9[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr10
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr10)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "10" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr10[i,]$start : macaque_exons_chr10[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr11
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr11)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "11" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr11[i,]$start : macaque_exons_chr11[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr12
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr12)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "12" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr12[i,]$start : macaque_exons_chr12[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr13
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr13)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "13" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr13[i,]$start : macaque_exons_chr13[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr14
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr14)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "14" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr14[i,]$start : macaque_exons_chr14[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr15
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr15)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "15" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr15[i,]$start : macaque_exons_chr15[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr16
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr16)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "16" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr16[i,]$start : macaque_exons_chr16[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr17
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr17)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "17" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr17[i,]$start : macaque_exons_chr17[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr18
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr18)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "18" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr18[i,]$start : macaque_exons_chr18[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr19
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr19)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "19" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr19[i,]$start : macaque_exons_chr19[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}

#chr20
exon_index <- c()
for (i in 1:nrow(macaque_exons_chr20)) {
  exon_index <- which(fdr_corrected_res_annotated$chr == "20" & fdr_corrected_res_annotated$CpG %in% (macaque_exons_chr20[i,]$start : macaque_exons_chr20[i,]$end) == TRUE)
  if  (sum(exon_index) > 0){
    fdr_corrected_res_annotated[exon_index,]$exon <- "YES"
  }
}


write.table(fdr_corrected_res_annotated, file ="fdr_corrected_res_all_gtf_chromhmm_annotations.txt")