#!/usr/bin/Rscript
library(dplyr)
library(data.table)

# load empirical 
chr_cpg_TSS_50kb_sig_counts_greater10 <- read.table("allchr_cpg_TSS_50kb_sig_counts_empirical_proportions.txt", header = TRUE)
library(dplyr)
library(tidyr)

empirical_genes <- separate(data = chr_cpg_TSS_50kb_sig_counts_greater10, col = gene, into = c("ID", "info"), sep = "; gene_version")
empirical_genes$ID <- gsub("gene_id ", "", as.character(empirical_genes$ID))

empirical_genes_order <- empirical_genes[order(empirical_genes$proportion_sig, decreasing = TRUE),]
empirical_proportion_genes_protein_coding_order <- filter(empirical_genes_order, grepl("protein_coding",info))

# plots and stats
# Number of CpGs within 50kb of each gene
ggplot(empirical_proportion_genes_protein_coding_order,aes(x=total_sites)) + geom_histogram(bins=40, alpha=0.6, fill="dark blue") + theme_bw() + xlab("Number of CpGs within 50kb macaque protein coding genes (n=10,415 genes)")
mean(empirical_proportion_genes_protein_coding_order$total_sites)
#30.5494
median(empirical_proportion_genes_protein_coding_order$total_sites)
#23
sd(empirical_proportion_genes_protein_coding_order$total_sites)
#24.24092

# Proportion of significant CpGs within 50kb of each gene
ggplot(empirical_proportion_genes_protein_coding_order,aes(x=proportion_sig)) + geom_histogram(bins=40, alpha=0.6, fill="dark red") + theme_bw() + xlab("Proportion of CpGs with significant meQTL within 50kb macaque protein coding genes (n=10,415 genes)")
mean(empirical_proportion_genes_protein_coding_order$proportion_sig)
#0.7326375 sig/total CpGs
median(empirical_proportion_genes_protein_coding_order$proportion_sig)
#0.754717
sd(empirical_proportion_genes_protein_coding_order$proportion_sig)
#0.1840721
                                                                                     
empirical_proportion_genes_protein_coding_order_2 <- empirical_proportion_genes_protein_coding_order[,c("ID","proportion_sig")]
write.table(empirical_proportion_genes_protein_coding_order_2, file="sigmeqtl_proportion_50kb_proteingenes_sorted_forGSEA.txt", row.names = FALSE)

# Genes (any with 10 or more sites) = background
# All that fall above the 95%tile = the "significantly more meqtl than expected by chance" genes

# protein coding only
empirical_proportion_genes_protein_coding <- filter(empirical_genes, grepl("protein_coding",info))
# background genes
background_genes_50kbcpg_chr_morethan10site <- empirical_proportion_genes_protein_coding$ID

write.table(background_genes_50kbcpg_chr_morethan10site, file ="background_genes_50kbcpg_allchr_morethan10site_proteincoding.txt", row.names = FALSE)

# load all 100 permuted proportion distributions
temp = list.files(pattern="*_.txt")
permutationfiles = lapply(temp, read.table)
chr_100_permutations <- do.call(rbind, permutationfiles)
hist(chr_100_permutations$proportion_sig)

# Get the 95% interval/percentile 
# which proportions ares greater than 95% of the null/permuted distribution - 

# multiply k percent (0.95) by the total number of values n (27,800)

values <- dim(chr_100_permutations)[1]

percentile_95 <- 0.95*values

#0.95*27800
#26410

#sort(chr_100_permutations$proportion_sig)[26410]
perm_percentile_number <- sort(chr_100_permutations$proportion_sig)[percentile_95]

# 95th percentile for this null distribution, should then take any genes with actual proportions higher than this

chr_genes_50kbsigcpgsite_min10site_proporgreater_than_95percentile <- empirical_proportion_genes_protein_coding$ID[empirical_proportion_genes_protein_coding$proportion_sig > perm_percentile_number]

#save
write.table(chr_genes_50kbsigcpgsite_min10site_proporgreater_than_95percentile, file ="greaterthan95percentile_allchr_nearestgene_50kb_sigcpg_proteincoding.txt", row.names = FALSE)
