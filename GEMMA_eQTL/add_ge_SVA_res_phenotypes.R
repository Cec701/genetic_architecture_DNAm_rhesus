#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
num = args[1]

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
n <- as.numeric(slurm_arrayid)
library("tidyverse", lib.loc="/packages/7x/R/4.0.2/lib64/R/library")

merged.chr.GP.HW.MAF.GE.with_repeats.sorted.gene_fam <- read.table(paste("merged.chr",num,".GP.HW.MAF.GE.with_repeats.sorted.gene",n,".fam", sep=""), header = FALSE)

# load output from SVA regression (first 5 SVs regressed out) 
resid_GE_wgs_chr_ordered <- read.table(paste("resid_GE_wgs_chr",num,"_ordered_SVA_5.txt",sep=""))

merged.chr.GP.HW.MAF.GE.with_repeats.sorted.gene_fam_update <- merged.chr.GP.HW.MAF.GE.with_repeats.sorted.gene_fam[,1:5]

fam_file_wgs_ge_bind_chr <- bind_cols(merged.chr.GP.HW.MAF.GE.with_repeats.sorted.gene_fam_update, resid_GE_wgs_chr_ordered)

setwd("SVA_5/")

write.table(fam_file_wgs_ge_bind_chr, file = paste("merged.chr",num,".GP.HW.MAF.GE.with_repeats.sorted.gene", n, ".fam", sep = ""), row.names = FALSE, col.names = FALSE)
