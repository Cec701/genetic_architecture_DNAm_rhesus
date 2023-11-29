#!/usr/bin/Rscript
#######  
# Get genotype and methylation info for IMAGE pairs to run in PQLseq
# UPDATED FINAL 

# enter chromosome number
n = 1

### load meta
meta <- read.table("REAL_FINAL_ORDERED_SUBSET_RRBSMETA_NOMIXUP_573")
###### get genotypes (from RRBS VCF converted to .traw with PLINK v1.9) 
image_snps_chrn <- read.table(paste("chr", n, "_image_snps.traw", sep=""), header=TRUE)
#### To get the cpgsnplist and count matrices
load(paste("merge_filter_573_chr",n,".RData", sep=""))
cpgslistdf <- as.data.frame(CpGlist)
cpgslistdf$chr <- n
cpgslistdf$ID <- paste(cpgslistdf$chr,cpgslistdf$CpGlist,sep="_")
cpgslistdf$SNP <- SNPlist

write.table(cpgslistdf, file=paste("cpg_snp_listdf_countorder_chr",n,".txt", sep=""))
row.names(rm) <- SNPlist
rownames(ym) <- SNPlist

write.table(rm, file=paste("chr",n,"_totalreadcounts.txt", sep=""))
write.table(ym, file=paste("chr",n,"_methreadcounts.txt", sep=""))

rm(ym, r1,r2,y1,y2,y2m,y1m,r1m,r2m)
rm(rm)
rm(SNP)

# snps

image_snps_chrn_edit <- image_snps_chrn[,-c(1:3, 5:6)]
rownames(image_snps_chrn_edit) <- image_snps_chrn_edit$POS
image_snps_chrn_edit_2 <- image_snps_chrn_edit[,-1]
colnames(image_snps_chrn_edit_2) <- gsub("\\..*", "", colnames(image_snps_chrn_edit_2))
colnames(image_snps_chrn_edit_2) <- gsub("X", "", colnames(image_snps_chrn_edit_2))
image_snps_chrn_edit_2_nomixup <- image_snps_chrn_edit_2[,-c(47,228,484,485,498,524,542,550)]
colnames_animalid <- meta$animal_id[match(names(image_snps_chrn_edit_2_nomixup[,43:572]), meta$lid)]
colnames(image_snps_chrn_edit_2_nomixup)[43:572] <- colnames_animalid
test_sortng_names <- meta$animal_id
image_snps_chrn_edit_2_nomixup_ORDERED <- image_snps_chrn_edit_2_nomixup[,test_sortng_names]
# save genotypes for input to PQLseq
write.table(image_snps_chrn_edit_2_nomixup_ORDERED, file =paste("image_genotypes_chr",n,"_573_ORDERED.txt", sep=""))
rm(image_snps_chr1, image_snps_chrn_edit, image_snps_chrn_edit_2, image_snps_chrn_edit_2_nomixup, image_snps_chrn_edit_2_nomixup_ORDERED)
