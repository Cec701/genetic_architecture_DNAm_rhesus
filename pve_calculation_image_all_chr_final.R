#!/usr/bin/Rscript

#PVE calculation meQTL SNP-CpG pairs 
#load methylation levels (all chromosomes, rbind)

methylationlevel_chr1 <- read.table("image_res_meth_levels_chr1_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr2 <- read.table("image_res_meth_levels_chr2_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr3 <- read.table("image_res_meth_levels_chr3_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr4 <- read.table("image_res_meth_levels_chr4_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr5 <- read.table("image_res_meth_levels_chr5_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr6 <- read.table("image_res_meth_levels_chr6_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr7 <- read.table("image_res_meth_levels_chr7_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr8 <- read.table("image_res_meth_levels_chr8_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr9 <- read.table("image_res_meth_levels_chr9_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr10 <- read.table("image_res_meth_levels_chr10_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr11 <- read.table("image_res_meth_levels_chr11_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr12 <- read.table("image_res_meth_levels_chr12_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr13 <- read.table("image_res_meth_levels_chr13_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr14 <- read.table("image_res_meth_levels_chr14_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr15 <- read.table("image_res_meth_levels_chr15_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr16 <- read.table("image_res_meth_levels_chr16_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr17 <- read.table("image_res_meth_levels_chr17_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr18 <- read.table("image_res_meth_levels_chr18_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr19 <- read.table("image_res_meth_levels_chr19_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylationlevel_chr20 <- read.table("image_res_meth_levels_chr20_copy.txt", header= FALSE)
meth_header <- read.table("meth_level_header.txt", header=FALSE)

methylation_levels <- rbind(methylationlevel_chr1, methylationlevel_chr2, methylationlevel_chr3, methylationlevel_chr4, methylationlevel_chr5,
                            methylationlevel_chr6, methylationlevel_chr7, methylationlevel_chr8, methylationlevel_chr9, methylationlevel_chr10,
                            methylationlevel_chr11, methylationlevel_chr12, methylationlevel_chr13, methylationlevel_chr14, methylationlevel_chr15,
                            methylationlevel_chr16, methylationlevel_chr17, methylationlevel_chr18, methylationlevel_chr19, methylationlevel_chr20)


uniquecpg_methlevels <- methylation_levels[!duplicated(methylation_levels$`"cpg"`),]

rownames(uniquecpg_methlevels) <- uniquecpg_methlevels$`“cpg”`

uniquecpg_methlevels <- uniquecpg_methlevels[,-1]

t_uniquecpg_methlevels <- as.data.frame(t(uniquecpg_methlevels))

covariates <- read.table("numeric_REAL_COVARIATE_MATRIX_ORDERED_FINAL_573", check.names = FALSE)

PVE_dataframe <- cbind(t_uniquecpg_methlevels, covariates)

write.table(PVE_dataframe, file ="for_pvecalc_methw_covariate.txt")

n <- (ncol(PVE_dataframe)-3)

# run n regressions
model_lms <- lapply(1:n, function(x) lm(PVE_dataframe[,x]~ sex + age + batch, data = PVE_dataframe))
summaries_modellm <- lapply(model_lms, summary)

X_model_res <- lapply(summaries_modellm, function(x) x$res)

# The equation is:
#PVE<-(beta^2*var(X_res))/(beta^2*var(X_res)+sigma^2)
# where var(X_model_res_chr1[[1]]) = var(X_res)

# Take variance
var_x_res_model <- lapply(X_model_res, function(x) var(x))

cpg_IDs_pve <- colnames(PVE_dataframe[1:n])

names(var_x_res_model) <- cpg_IDs_pve

#load image results
res_macaque_meqtl <- read.table("fdr_corrected_image_res_annotated_allchr_wchromhmm_FINAL.txt", header=TRUE)

res_macaque_meqtl$ID <- paste(res_macaque_meqtl$chr, res_macaque_meqtl$CpG, sep="_")
var_x_res_model_dataframe <- as.data.frame(unlist(var_x_res_model))
var_x_res_model_dataframe$ID <- rownames(var_x_res_model_dataframe)
colnames(var_x_res_model_dataframe) <- c("var_x_res", "ID")

res_with_covariate_varx <- merge(var_x_res_model_dataframe, res_macaque_meqtl, by = "ID")

#PVE<-(beta^2*var(X_res))/(beta^2*var(X_res)+sigma^2)
#sigma^2 from model output

# calculate PVE
res_with_covariate_varx$PVE <- ((res_with_covariate_varx$beta)^2*res_with_covariate_varx$var_x_res)/((res_with_covariate_varx$beta)^2*res_with_covariate_varx$var_x_res + res_with_covariate_varx$sigma2)

write.table(res_with_covariate_varx, file = "fdr_corrected_image_res_final_wPVE.txt")

res_with_covariate_varx_meqtl <- subset(res_with_covariate_varx, fdr_FDR_IMAGE_res_all < 0.05)
res_with_covariate_varx_meqtl_filt <- subset(res_with_covariate_varx_meqtl, sigma2 > 0)
print(mean(res_with_covariate_varx_meqtl_filt$PVE))

library(ggplot2)
#pdf
pdf(file='image_meqtl_res_pve.pdf')
par(mfrow=c(1,2))
hist(res_with_covariate_varx_meqtl_filt$PVE, main = "meQTLs Percent Variance Explained (PVE)", xlab = "PVE", col = "sky blue")
ggplot(data = res_with_covariate_varx_meqtl_filt, aes(x=PVE, y=abs(beta))) + geom_point() + theme_bw() + xlab("PVE") + ylab("abs(beta)") +geom_smooth(method='lm')
