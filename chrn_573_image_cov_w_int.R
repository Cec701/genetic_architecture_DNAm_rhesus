#!/usr/bin/Rscript

#Enter chromosome number (1-20) or run as an array
CHROMNUM=8

#load filtered matrices 
load(file="merge_filter_573_chr8.RData")

#######################
# meQTL analysis with image
#######################

# run with all samples (N=573)
library(devtools)
install_github("fanyue322/IMAGE")
library(IMAGE)

rel=read.table('FINAL_REPEAT_PEDIGREE_RRBS_573_NOMIXUP', check.names = FALSE)
K_573 <- as.matrix(rel)

#add covariate matrix 
cov <- read.table('numeric_REAL_COVARIATE_MATRIX_ORDERED_FINAL_573',check.names = FALSE)
#add intercept
cov_int <- cbind(int=1, cov)
cov_573_int <- as.matrix(cov_int)

res_chr8_573_cov_int = IMAGE::image(geno,data,K_573,Covariates = cov_573_int, numCore=24)


write.table(res_chr8_573_cov_int,paste('20May_IMAGE_results_chr',CHROMNUM,'.txt',sep=''),row.names=F,sep='\t')
write.table(CpGlist,paste('20May_IMAGE_CpGs_chr',CHROMNUM,'.txt',sep=''),row.names=F,sep='\t')
write.table(SNPlist,paste('20May_IMAGE_SNPs_chr',CHROMNUM,'.txt',sep=''),row.names=F,sep='\t')

pdf(file='chr8_573_cov_int.pdf')
par(mfrow=c(2,2))
hist(res_chr8_573_cov_int$pvalue,breaks=100,main='')
hist(res_chr8_573_cov_int$h2,breaks=100,main='')
dev.off()

