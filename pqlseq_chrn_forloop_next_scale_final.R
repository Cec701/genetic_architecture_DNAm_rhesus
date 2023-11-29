#!/usr/bin/Rscript
##PQLSEQ on IMAGE SNP-CpG PAIRS 
# TO CALCULATE PVE 

library(doParallel)
library(parallel)
library(Rcpp)
library(MuMIn)
library(Matrix)
library(dplyr, lib.loc = "/packages/7x/R/4.0.2/lib64/R/library")
source("PQLseq_MW.R") 
sourceCpp("PQLseq.cpp")

#Or array
args = commandArgs(trailingOnly=TRUE)
#args[1]
#SET CHR (run for 1-20)
CHR = args[1]

#load genotypes
image_snps_chr <- read.table(paste("image_genotypes_chr", CHR, "_573_ORDERED.txt", sep=""), header=TRUE)
#load covariates
covariates <- read.table("numeric_REAL_COVARIATE_MATRIX_ORDERED_FINAL_573")
#load total read counts
rm_total <- read.table(paste("chr", CHR, "_totalreadcounts.txt", sep=""), row.names=NULL)
rm_total <- rm_total[,-1]
#load methylated read counts
ym_total <- read.table(paste("chr", CHR, "_methreadcounts.txt", sep=""), row.names=NULL)
ym_total <- ym_total[,-1]
#load SNP-CpG list original 
cpgslistdf <- read.table(paste("cpg_snp_listdf_countorder_chr", CHR, ".txt", sep=""), header=TRUE)
# get df with snps with geno
has_geno_idx <- which(cpgslistdf$SNP %in% rownames(image_snps_chr))
cpgslistdf_withgeno <- cpgslistdf[has_geno_idx,]
#load relatedness matrix
pedigree_relatedmatrix <- read.table("FINAL_REPEAT_PEDIGREE_RRBS_573_NOMIXUP")

#scale variables
covariates$age <- scale(as.numeric(covariates$age))
covariates$sex_binary <- scale(as.numeric(covariates$sex_binary))
covariates$batch <- scale(as.numeric(covariates$batch))

#get distance
cpgslistdf_withgeno$distance <- abs(cpgslistdf_withgeno$SNP - cpgslistdf_withgeno$CpGlist)

#set names
rownames(pedigree_relatedmatrix) <- colnames(pedigree_relatedmatrix)
rownames(covariates) <- colnames(pedigree_relatedmatrix)
colnames(ym_total) <- colnames(pedigree_relatedmatrix)
colnames(rm_total) <- colnames(pedigree_relatedmatrix)

## each chunk will be separate df for the pqlseq run
## and will correspond to the corresponding cpg-snp pair in the cpglistdf
the_index_list <- which(!duplicated(cpgslistdf$SNP)== TRUE)
split_rm_total_counts_df <- split.data.frame(rm_total, findInterval(1:nrow(rm_total), the_index_list))
split_ym_meth_counts_df <- split.data.frame(ym_total, findInterval(1:nrow(ym_total), the_index_list))
count_rowname <- unique(cpgslistdf$SNP)
names(split_rm_total_counts_df) <- count_rowname
names(split_ym_meth_counts_df) <- count_rowname

# index count lists by snps with geno
count_matrix_idx <- which(names(split_rm_total_counts_df) %in% unique(cpgslistdf_withgeno$SNP))
split_rm_total_counts_withgeno <- split_rm_total_counts_df[count_matrix_idx]

count_matrix_idx <- which(names(split_ym_meth_counts_df) %in% unique(cpgslistdf_withgeno$SNP))
split_ym_meth_counts_withgeno <- split_ym_meth_counts_df[count_matrix_idx]

# correct snp-cpg df
the_index_list_withgeno <- which(!duplicated(cpgslistdf_withgeno$SNP)== TRUE)
split_cpglist_df <- split.data.frame(cpgslistdf_withgeno, findInterval(1:nrow(cpgslistdf_withgeno), the_index_list_withgeno))

# make df
pedigree_relatedmatrix_df <- as.data.frame(pedigree_relatedmatrix)

pqlseq_colnames <- c("numIDV", "beta_age", "beta_batch", "beta_sex_binary", "beta_genotype", "se_age", "se_batch", "se_sex_binary", "se_genotype", "pval_age",
 "pval_batch", "pval_sex_binary", "pval_genotype", "h2", "sigma2", "LL", "AIC", "converged")


#### RUN PQLSEQ
# for loop for model output

model_DNA_meth_list <- list()
for (i in 1:length(split_ym_meth_counts_withgeno)){
	genotype_i <- scale(as.numeric(image_snps_chr[i,]))
	model_DNA_meth_allres = try(pqlseq(RawCountDataSet= as.data.frame(split_ym_meth_counts_withgeno[i]), Phenotypes=genotype_i, Covariates=covariates,
                                 RelatednessMatrix=pedigree_relatedmatrix_df, LibSize=as.data.frame(split_rm_total_counts_withgeno[i]),
                                 fit.model="BMM"))
    
	if(class(model_DNA_meth_allres) == "try-error") {
		print("tryerror")
		next
 	} else {
	  print("else")
	  model_DNA_meth_allres <- as.data.frame(model_DNA_meth_allres)
	  colnames(model_DNA_meth_allres) <- pqlseq_colnames
	  split_cpglist_df_i <- as.data.frame(split_cpglist_df[[i]])
          model_DNA_meth_allres$CpG <- split_cpglist_df_i$CpGlist
          model_DNA_meth_allres$SNP <- split_cpglist_df_i$SNP
          model_DNA_meth_allres$dist <- split_cpglist_df_i$distance
    
##Get the variance for each model (each uses different subset)
##Don't need to scale these b/c the scaling was done outside of the model matrix 
##which animals have NA 
##animal names are colnames for the snps
    	  notmissinggeno_index <- which(!is.na(genotype_i))
    #and then filter the covariates file on that
    	  covariates_modeli <- covariates[notmissinggeno_index,]
    	  image_snps_chr_noNA <- genotype_i[notmissinggeno_index,]
    
    	  model_DNA_meth_allres$A_age <- ((model_DNA_meth_allres$beta_age)^2)*as.numeric(var(covariates_modeli$age))
    	  model_DNA_meth_allres$A_batch <- ((model_DNA_meth_allres$beta_batch)^2)*as.numeric(var(covariates_modeli$batch))
    	  model_DNA_meth_allres$A_sex_binary <- ((model_DNA_meth_allres$beta_sex_binary)^2)*as.numeric(var(covariates_modeli$sex_binary))
    	  model_DNA_meth_allres$A_genotype <- ((model_DNA_meth_allres$beta_genotype)^2)*as.numeric(var(image_snps_chr_noNA))
    
    #
    	  model_DNA_meth_list[[i]] <- model_DNA_meth_allres
	  
  	}
}                


final_model_DNA_meth = do.call(rbind, model_DNA_meth_list)      

write.table(final_model_DNA_meth, file=paste("chr", CHR, "_snps_pqlseqout_scaled_final.txt", sep=""))      
