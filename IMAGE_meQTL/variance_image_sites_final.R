#!/usr/bin/Rscript
# Get variance for methylation level at each CpG 

# load methylated and total read counts and cpg snp df

# chr1 
#load total read counts
rm_total_1a <- read.table("chr1_totalreadcounts_1.txt", row.names=NULL)
rm_total_1a <- rm_total_1a[,-1]
#load meth read counts
ym_total_1a <- read.table("chr1_methreadcounts_1.txt", row.names=NULL)
ym_total_1a <- ym_total_1a[,-1]
#get meth as proportion
propor_meth_1a <- ym_total_1a/rm_total_1a
#load cpg_snp df
cpgslistdf_1a <- read.table("cpg_snp_listdf_countorder_chr1_1.txt", header=TRUE)
dim(cpgslistdf_1a)[1]
propr_meth_1a_1_vector <- as.vector(propor_meth_1a[1,])

#get variance for each row
var_cpg_1a <- apply(propor_meth_1a, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_1a$meth_var <- var_cpg_1a
hist(cpgslistdf_1a$meth_var)
plot(density(cpgslistdf_1a$meth_var))
#save cpglist with variance
write.table(cpgslistdf_1a, file="cpg_snp_listdf_countorder_chr1_1_wvar.txt")
rm(rm_total_1a, ym_total_1a)

#load total read counts
rm_total_1b <- read.table("chr1_totalreadcounts_2.txt", row.names=NULL)
rm_total_1b <- rm_total_1b[,-1]
#load meth read counts
ym_total_1b <- read.table("chr1_methreadcounts_2.txt", row.names=NULL)
ym_total_1b <- ym_total_1b[,-1]
#get meth as proportion
propor_meth_1b <- ym_total_1b/rm_total_1b
#load cpg_snp df
cpgslistdf_1b <- read.table("cpg_snp_listdf_countorder_chr1_2.txt", header=TRUE)
dim(cpgslistdf_1a)[1]
#get variance for each row
var_cpg_1b <- apply(propor_meth_1b, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_1b$meth_var <- var_cpg_1b
hist(cpgslistdf_1b$meth_var)
plot(density(cpgslistdf_1b$meth_var))
#save cpglist with variance
write.table(cpgslistdf_1b, file="cpg_snp_listdf_countorder_chr1_2_wvar.txt")
rm(rm_total_1b, ym_total_1b)

# load total read counts
rm_total_1c <- read.table("chr1_totalreadcounts_3.txt", row.names=NULL)
rm_total_1c <- rm_total_1c[,-1]
# load meth read counts
ym_total_1c <- read.table("chr1_methreadcounts_3.txt", row.names=NULL)
ym_total_1c <- ym_total_1c[,-1]
# get meth as proportion
propor_meth_1c <- ym_total_1c/rm_total_1c
# load cpg_snp df
cpgslistdf_1c <- read.table("cpg_snp_listdf_countorder_chr1_3.txt", header=TRUE)
# get variance for each row
var_cpg_1c <- apply(propor_meth_1c, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_1c$meth_var <- var_cpg_1c
hist(cpgslistdf_1c$meth_var)
plot(density(cpgslistdf_1c$meth_var))
# save cpglist with variance
write.table(cpgslistdf_1c, file="cpg_snp_listdf_countorder_chr1_3_wvar.txt")
rm(rm_total_1c, ym_total_1c)

# load total read counts
rm_total_1d <- read.table("chr1_totalreadcounts_4.txt", row.names=NULL)
rm_total_1d <- rm_total_1d[,-1]
#load meth read counts
ym_total_1d <- read.table("chr1_methreadcounts_4.txt", row.names=NULL)
ym_total_1d <- ym_total_1d[,-1]
#get meth as proportion
propor_meth_1d <- ym_total_1d/rm_total_1d
#load cpg_snp df
cpgslistdf_1d <- read.table("cpg_snp_listdf_countorder_chr1_4.txt", header=TRUE)
#get variance for each row
var_cpg_1d <- apply(propor_meth_1d, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_1d$meth_var <- var_cpg_1d
hist(cpgslistdf_1d$meth_var)
plot(density(cpgslistdf_1d$meth_var))
#save cpglist with variance
write.table(cpgslistdf_1d, file="cpg_snp_listdf_countorder_chr1_4_wvar.txt")
rm(rm_total_1d, ym_total_1d)

#load total read counts
rm_total_1e <- read.table("chr1_totalreadcounts_5.txt", row.names=NULL)
rm_total_1e <- rm_total_1e[,-1]
#load meth read counts
ym_total_1e <- read.table("chr1_methreadcounts_5.txt", row.names=NULL)
ym_total_1e <- ym_total_1e[,-1]
#get meth as proportion
propor_meth_1e <- ym_total_1e/rm_total_1e
#load cpg_snp df
cpgslistdf_1e <- read.table("cpg_snp_listdf_countorder_chr1_5.txt", header=TRUE)
#get variance for each row
var_cpg_1e <- apply(propor_meth_1e, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_1e$meth_var <- var_cpg_1e
hist(cpgslistdf_1e$meth_var)
plot(density(cpgslistdf_1e$meth_var))
#save cpglist with variance
write.table(cpgslistdf_1e, file="cpg_snp_listdf_countorder_chr1_5_wvar.txt")
rm(rm_total_1e, ym_total_1e)

#chr2 
#load total read counts
rm_total_2 <- read.table("chr2_totalreadcounts.txt", row.names=NULL)
rm_total_2 <- rm_total_2[,-1]
#load meth read counts
ym_total_2 <- read.table("chr2_methreadcounts.txt", row.names=NULL)
ym_total_2 <- ym_total_2[,-1]
#get meth as proportion
propor_meth_2 <- ym_total_2/rm_total_2
#load cpg_snp df
cpgslistdf_2 <- read.table("cpg_snp_listdf_countorder_chr2.txt", header=TRUE)
#get variance for each row
var_cpg_2 <- apply(propor_meth_2, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_2$meth_var <- var_cpg_2
hist(cpgslistdf_2$meth_var)
plot(density(cpgslistdf_2$meth_var))
#save cpglist with variance
write.table(cpgslistdf_2, file="cpg_snp_listdf_countorder_chr2_wvar.txt")
rm(rm_total_2, ym_total_2)

#chr3
#load total read counts
rm_total_3a <- read.table("chr3_totalreadcounts_1.txt", row.names=NULL)
rm_total_3a <- rm_total_3a[,-1]
#load meth read counts
ym_total_3a <- read.table("chr3_methreadcounts_1.txt", row.names=NULL)
ym_total_3a <- ym_total_3a[,-1]
#get meth as proportion
propor_meth_3a <- ym_total_3a/rm_total_3a
#load cpg_snp df
cpgslistdf_3a <- read.table("cpg_snp_listdf_countorder_chr3_1.txt", header=TRUE)
#get variance for each row
var_cpg_3a <- apply(propor_meth_3a, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_3a$meth_var <- var_cpg_3a
hist(cpgslistdf_3a$meth_var)
plot(density(cpgslistdf_3a$meth_var))
#save cpglist with variance
write.table(cpgslistdf_3a, file="cpg_snp_listdf_countorder_chr3_1_wvar.txt")
rm(rm_total_3a, ym_total_3a)

#load total read counts
rm_total_3b <- read.table("chr3_totalreadcounts_2.txt", row.names=NULL)
rm_total_3b <- rm_total_3b[,-1]
#load meth read counts
ym_total_3b <- read.table("chr3_methreadcounts_2.txt", row.names=NULL)
ym_total_3b <- ym_total_3b[,-1]
#get meth as proportion
propor_meth_3b <- ym_total_3b/rm_total_3b
#load cpg_snp df
cpgslistdf_3b <- read.table("cpg_snp_listdf_countorder_chr3_2.txt", header=TRUE)
#get variance for each row
var_cpg_3b <- apply(propor_meth_3b, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_3b$meth_var <- var_cpg_3b
hist(cpgslistdf_3b$meth_var)
plot(density(cpgslistdf_3b$meth_var))
#save cpglist with variance
write.table(cpgslistdf_3b, file="cpg_snp_listdf_countorder_chr3_2_wvar.txt")
rm(rm_total_3b, ym_total_3b)

#load total read counts
rm_total_3c <- read.table("chr3_totalreadcounts_3.txt", row.names=NULL)
rm_total_3c <- rm_total_3c[,-1]
#load meth read counts
ym_total_3c <- read.table("chr3_methreadcounts_3.txt", row.names=NULL)
ym_total_3c <- ym_total_3c[,-1]
#get meth as proportion
propor_meth_3c <- ym_total_3c/rm_total_3c
#load cpg_snp df
cpgslistdf_3c <- read.table("cpg_snp_listdf_countorder_chr3_3.txt", header=TRUE)
#get variance for each row
var_cpg_3c <- apply(propor_meth_3c, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_3c$meth_var <- var_cpg_3c
hist(cpgslistdf_3c$meth_var)
plot(density(cpgslistdf_3c$meth_var))
#save cpglist with variance
write.table(cpgslistdf_3c, file="cpg_snp_listdf_countorder_chr3_3_wvar.txt")
rm(rm_total_3c, ym_total_3c)

#load total read counts
rm_total_3d <- read.table("chr3_totalreadcounts_4.txt", row.names=NULL)
rm_total_3d <- rm_total_3d[,-1]
#load meth read counts
ym_total_3d <- read.table("chr3_methreadcounts_4.txt", row.names=NULL)
ym_total_3d <- ym_total_3d[,-1]
#get meth as proportion
propor_meth_3d <- ym_total_3d/rm_total_3d
#load cpg_snp df
cpgslistdf_3d <- read.table("cpg_snp_listdf_countorder_chr3_4.txt", header=TRUE)
#get variance for each row
var_cpg_3d <- apply(propor_meth_3d, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_3d$meth_var <- var_cpg_3d
hist(cpgslistdf_3d$meth_var)
plot(density(cpgslistdf_3d$meth_var))
#save cpglist with variance
write.table(cpgslistdf_3d, file="cpg_snp_listdf_countorder_chr3_4_wvar.txt")
rm(rm_total_3d, ym_total_3d)

#load total read counts
rm_total_3e <- read.table("chr3_totalreadcounts_5.txt", row.names=NULL)
rm_total_3e <- rm_total_3e[,-1]
#load meth read counts
ym_total_3e <- read.table("chr3_methreadcounts_5.txt", row.names=NULL)
ym_total_3e <- ym_total_3e[,-1]
#get meth as proportion
propor_meth_3e <- ym_total_3e/rm_total_3e
#load cpg_snp df
cpgslistdf_3e <- read.table("cpg_snp_listdf_countorder_chr3_5.txt", header=TRUE)
#get variance for each row
var_cpg_3e <- apply(propor_meth_3e, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_3e$meth_var <- var_cpg_3e
hist(cpgslistdf_3e$meth_var)
plot(density(cpgslistdf_3e$meth_var))
#it feels like the max variance has been like 0.25.. but the mean was 0.6 based on sigma2 estimate from image?
#save cpglist with variance
write.table(cpgslistdf_3e, file="cpg_snp_listdf_countorder_chr3_5_wvar.txt")
rm(rm_total_3e, ym_total_3e)

#load total read counts
rm_total_3f <- read.table("chr3_totalreadcounts_6.txt", row.names=NULL)
rm_total_3f <- rm_total_3f[,-1]
#load meth read counts
ym_total_3f <- read.table("chr3_methreadcounts_6.txt", row.names=NULL)
ym_total_3f <- ym_total_3f[,-1]
#get meth as proportion
propor_meth_3f <- ym_total_3f/rm_total_3f
#load cpg_snp df
cpgslistdf_3f <- read.table("cpg_snp_listdf_countorder_chr3_6.txt", header=TRUE)
#get variance for each row
var_cpg_3f <- apply(propor_meth_3f, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_3f$meth_var <- var_cpg_3f
hist(cpgslistdf_3f$meth_var)
plot(density(cpgslistdf_3f$meth_var))
#it feels like the max variance has been like 0.25.. but the mean was 0.6 based on sigma2 estimate from image?
#save cpglist with variance
write.table(cpgslistdf_3f, file="cpg_snp_listdf_countorder_chr3_6_wvar.txt")
rm(rm_total_3f, ym_total_3f)

#load total read counts
rm_total_3g <- read.table("chr3_totalreadcounts_7.txt", row.names=NULL)
rm_total_3g <- rm_total_3g[,-1]
#load meth read counts
ym_total_3g <- read.table("chr3_methreadcounts_7.txt", row.names=NULL)
ym_total_3g <- ym_total_3g[,-1]
#get meth as proportion
propor_meth_3g <- ym_total_3g/rm_total_3g
#load cpg_snp df
cpgslistdf_3g <- read.table("cpg_snp_listdf_countorder_chr3_7.txt", header=TRUE)
#get variance for each row
var_cpg_3g <- apply(propor_meth_3g, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_3g$meth_var <- var_cpg_3g
hist(cpgslistdf_3g$meth_var)
plot(density(cpgslistdf_3g$meth_var))
#it feels like the max variance has been like 0.25.. but the mean was 0.6 based on sigma2 estimate from image?
#save cpglist with variance
write.table(cpgslistdf_3g, file="cpg_snp_listdf_countorder_chr3_7_wvar.txt")
rm(rm_total_3g, ym_total_3g)

#load total read counts
rm_total_3h <- read.table("chr3_totalreadcounts_8.txt", row.names=NULL)
rm_total_3h <- rm_total_3h[,-1]
#load meth read counts
ym_total_3h <- read.table("chr3_methreadcounts_8.txt", row.names=NULL)
ym_total_3h <- ym_total_3h[,-1]
#get meth as proportion
propor_meth_3h <- ym_total_3h/rm_total_3h
#load cpg_snp df
cpgslistdf_3h <- read.table("cpg_snp_listdf_countorder_chr3_8.txt", header=TRUE)
#get variance for each row
var_cpg_3h <- apply(propor_meth_3h, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_3h$meth_var <- var_cpg_3h
hist(cpgslistdf_3h$meth_var)
plot(density(cpgslistdf_3h$meth_var))
#save cpglist with variance
write.table(cpgslistdf_3h, file="cpg_snp_listdf_countorder_chr3_8_wvar.txt")
rm(rm_total_3h, ym_total_3h)

#load total read counts
rm_total_3i <- read.table("chr3_totalreadcounts_9.txt", row.names=NULL)
rm_total_3i <- rm_total_3i[,-1]
#load meth read counts
ym_total_3i <- read.table("chr3_methreadcounts_9.txt", row.names=NULL)
ym_total_3i <- ym_total_3i[,-1]
#get meth as proportion
propor_meth_3i <- ym_total_3i/rm_total_3i
#load cpg_snp df
cpgslistdf_3i <- read.table("cpg_snp_listdf_countorder_chr3_9.txt", header=TRUE)
#get variance for each row
var_cpg_3i <- apply(propor_meth_3i, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_3i$meth_var <- var_cpg_3i
hist(cpgslistdf_3i$meth_var)
plot(density(cpgslistdf_3i$meth_var))
#save cpglist with variance
write.table(cpgslistdf_3i, file="cpg_snp_listdf_countorder_chr3_9_wvar.txt")
rm(rm_total_3i, ym_total_3i)

#load total read counts
rm_total_3j <- read.table("chr3_totalreadcounts_10.txt", row.names=NULL)
rm_total_3j <- rm_total_3j[,-1]
#load meth read counts
ym_total_3j <- read.table("chr3_methreadcounts_10.txt", row.names=NULL)
ym_total_3j <- ym_total_3j[,-1]
#get meth as proportion
propor_meth_3j <- ym_total_3j/rm_total_3j
#load cpg_snp df
cpgslistdf_3j <- read.table("cpg_snp_listdf_countorder_chr3_10.txt", header=TRUE)
#get variance for each row
var_cpg_3j <- apply(propor_meth_3j, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_3j$meth_var <- var_cpg_3j
hist(cpgslistdf_3j$meth_var)
plot(density(cpgslistdf_3j$meth_var))
#save cpglist with variance
write.table(cpgslistdf_3j, file="cpg_snp_listdf_countorder_chr3_10_wvar.txt")
rm(rm_total_3j, ym_total_3j)

#load total read counts
rm_total_4 <- read.table("chr4_totalreadcounts.txt", row.names=NULL)
rm_total_4 <- rm_total_4[,-1]
#load meth read counts
ym_total_4 <- read.table("chr4_methreadcounts.txt", row.names=NULL)
ym_total_4 <- ym_total_4[,-1]
#get meth as proportion
propor_meth_4 <- ym_total_4/rm_total_4
#load cpg_snp df
cpgslistdf_4 <- read.table("cpg_snp_listdf_countorder_chr4.txt", header=TRUE)
#get variance for each row
var_cpg_4 <- apply(propor_meth_4, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_4$meth_var <- var_cpg_4
hist(cpgslistdf_4$meth_var)
plot(density(cpgslistdf_4$meth_var))
#save cpglist with variance
write.table(cpgslistdf_4, file="cpg_snp_listdf_countorder_chr4_wvar.txt")
rm(rm_total_4, ym_total_4)

#load total read counts
rm_total_5 <- read.table("chr5_totalreadcounts.txt", row.names=NULL)
rm_total_5 <- rm_total_5[,-1]
#load meth read counts
ym_total_5 <- read.table("chr5_methreadcounts.txt", row.names=NULL)
ym_total_5 <- ym_total_5[,-1]
#get meth as proportion
propor_meth_5 <- ym_total_5/rm_total_5
#load cpg_snp df
cpgslistdf_5 <- read.table("cpg_snp_listdf_countorder_chr5.txt", header=TRUE)
#get variance for each row
var_cpg_5 <- apply(propor_meth_5, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_5$meth_var <- var_cpg_5
hist(cpgslistdf_5$meth_var)
plot(density(cpgslistdf_5$meth_var))
#save cpglist with variance
write.table(cpgslistdf_5, file="cpg_snp_listdf_countorder_chr5_wvar.txt")
rm(rm_total_5, ym_total_5)

#load total read counts
rm_total_6 <- read.table("chr6_totalreadcounts.txt", row.names=NULL)
rm_total_6 <- rm_total_6[,-1]
#load meth read counts
ym_total_6 <- read.table("chr6_methreadcounts.txt", row.names=NULL)
ym_total_6 <- ym_total_6[,-1]
#get meth as proportion
propor_meth_6 <- ym_total_6/rm_total_6
#load cpg_snp df
cpgslistdf_6 <- read.table("cpg_snp_listdf_countorder_chr6.txt", header=TRUE)
#get variance for each row
var_cpg_6 <- apply(propor_meth_6, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_6$meth_var <- var_cpg_6
hist(cpgslistdf_6$meth_var)
plot(density(cpgslistdf_6$meth_var))
#save cpglist with variance
write.table(cpgslistdf_6, file="cpg_snp_listdf_countorder_chr6_wvar.txt")
rm(rm_total_6, ym_total_6)

#load total read counts
rm_total_7 <- read.table("chr7_totalreadcounts.txt", row.names=NULL)
rm_total_7 <- rm_total_7[,-1]
#load meth read counts
ym_total_7 <- read.table("chr7_methreadcounts.txt", row.names=NULL)
ym_total_7 <- ym_total_7[,-1]
#get meth as proportion
propor_meth_7 <- ym_total_7/rm_total_7
#load cpg_snp df
cpgslistdf_7 <- read.table("cpg_snp_listdf_countorder_chr7.txt", header=TRUE)
#get variance for each row
var_cpg_7 <- apply(propor_meth_7, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_7$meth_var <- var_cpg_7
hist(cpgslistdf_7$meth_var)
plot(density(cpgslistdf_7$meth_var))
#save cpglist with variance
write.table(cpgslistdf_7, file="cpg_snp_listdf_countorder_chr7_wvar.txt")
rm(rm_total_7, ym_total_7)

#load total read counts
rm_total_8 <- read.table("chr8_totalreadcounts.txt", row.names=NULL)
rm_total_8 <- rm_total_8[,-1]
#load meth read counts
ym_total_8 <- read.table("chr8_methreadcounts.txt", row.names=NULL)
ym_total_8 <- ym_total_8[,-1]
#get meth as proportion
propor_meth_8 <- ym_total_8/rm_total_8
#load cpg_snp df
cpgslistdf_8 <- read.table("cpg_snp_listdf_countorder_chr8.txt", header=TRUE)
#get variance for each row
var_cpg_8 <- apply(propor_meth_8, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_8$meth_var <- var_cpg_8
hist(cpgslistdf_8$meth_var)
plot(density(cpgslistdf_8$meth_var))
#save cpglist with variance
vector_cpg1_chr8 <- as.vector(propor_meth_8[1,])
var(unlist(vector_cpg1_chr8), na.rm = TRUE)
write.table(cpgslistdf_8, file="cpg_snp_listdf_countorder_chr8_wvar.txt")
rm(rm_total_8, ym_total_8)

#load total read counts
rm_total_9 <- read.table("chr9_totalreadcounts.txt", row.names=NULL)
rm_total_9 <- rm_total_9[,-1]
#load meth read counts
ym_total_9 <- read.table("chr9_methreadcounts.txt", row.names=NULL)
ym_total_9 <- ym_total_9[,-1]
#get meth as proportion
propor_meth_9 <- ym_total_9/rm_total_9
#load cpg_snp df
cpgslistdf_9 <- read.table("cpg_snp_listdf_countorder_chr9.txt", header=TRUE)
#get variance for each row
var_cpg_9 <- apply(propor_meth_9, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_9$meth_var <- var_cpg_9
hist(cpgslistdf_9$meth_var)
plot(density(cpgslistdf_9$meth_var))
#save cpglist with variance
write.table(cpgslistdf_9, file="cpg_snp_listdf_countorder_chr9_wvar.txt")
rm(rm_total_9, ym_total_9)

#load total read counts
rm_total_10 <- read.table("chr10_totalreadcounts.txt", row.names=NULL)
rm_total_10 <- rm_total_10[,-1]
#load meth read counts
ym_total_10 <- read.table("chr10_methreadcounts.txt", row.names=NULL)
ym_total_10 <- ym_total_10[,-1]
#get meth as proportion
propor_meth_10 <- ym_total_10/rm_total_10
#load cpg_snp df
cpgslistdf_10 <- read.table("cpg_snp_listdf_countorder_chr10.txt", header=TRUE)
#get variance for each row
var_cpg_10 <- apply(propor_meth_10, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_10$meth_var <- var_cpg_10
hist(cpgslistdf_10$meth_var)
plot(density(cpgslistdf_10$meth_var))
#save cpglist with variance
write.table(cpgslistdf_10, file="cpg_snp_listdf_countorder_chr10_wvar.txt")
rm(rm_total_10, ym_total_10)

#load total read counts
rm_total_11 <- read.table("chr11_totalreadcounts.txt", row.names=NULL)
rm_total_11 <- rm_total_11[,-1]
#load meth read counts
ym_total_11 <- read.table("chr11_methreadcounts.txt", row.names=NULL)
ym_total_11 <- ym_total_11[,-1]
#get meth as proportion
propor_meth_11 <- ym_total_11/rm_total_11
#load cpg_snp df
cpgslistdf_11 <- read.table("cpg_snp_listdf_countorder_chr11.txt", header=TRUE)
#get variance for each row
var_cpg_11 <- apply(propor_meth_11, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_11$meth_var <- var_cpg_11
hist(cpgslistdf_11$meth_var)
plot(density(cpgslistdf_11$meth_var))
#save cpglist with variance
write.table(cpgslistdf_11, file="cpg_snp_listdf_countorder_chr11_wvar.txt")
rm(rm_total_11, ym_total_11)

#load total read counts
rm_total_12 <- read.table("chr12_totalreadcounts.txt", row.names=NULL)
rm_total_12 <- rm_total_12[,-1]
#load meth read counts
ym_total_12 <- read.table("chr12_methreadcounts.txt", row.names=NULL)
ym_total_12 <- ym_total_12[,-1]
#get meth as proportion
propor_meth_12 <- ym_total_12/rm_total_12
#load cpg_snp df
cpgslistdf_12 <- read.table("cpg_snp_listdf_countorder_chr12.txt", header=TRUE)
#get variance for each row
var_cpg_12 <- apply(propor_meth_12, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_12$meth_var <- var_cpg_12
hist(cpgslistdf_12$meth_var)
plot(density(cpgslistdf_12$meth_var))
#save cpglist with variance
write.table(cpgslistdf_12, file="cpg_snp_listdf_countorder_chr12_wvar.txt")
rm(rm_total_12, ym_total_12)

#load total read counts
rm_total_13 <- read.table("chr13_totalreadcounts.txt", row.names=NULL)
rm_total_13 <- rm_total_13[,-1]
#load meth read counts
ym_total_13 <- read.table("chr13_methreadcounts.txt", row.names=NULL)
ym_total_13 <- ym_total_13[,-1]
#get meth as proportion
propor_meth_13 <- ym_total_13/rm_total_13
#load cpg_snp df
cpgslistdf_13 <- read.table("cpg_snp_listdf_countorder_chr13.txt", header=TRUE)
#get variance for each row
var_cpg_13 <- apply(propor_meth_13, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_13$meth_var <- var_cpg_13
hist(cpgslistdf_13$meth_var)
plot(density(cpgslistdf_13$meth_var))
#save cpglist with variance
write.table(cpgslistdf_13, file="cpg_snp_listdf_countorder_chr13_wvar.txt")
rm(rm_total_13, ym_total_13)

#load total read counts
rm_total_14 <- read.table("chr14_totalreadcounts.txt", row.names=NULL)
rm_total_14 <- rm_total_14[,-1]
#load meth read counts
ym_total_14 <- read.table("chr14_methreadcounts.txt", row.names=NULL)
ym_total_14 <- ym_total_14[,-1]
#get meth as proportion
propor_meth_14 <- ym_total_14/rm_total_14
#load cpg_snp df
cpgslistdf_14 <- read.table("cpg_snp_listdf_countorder_chr14.txt", header=TRUE)
#get variance for each row
var_cpg_14 <- apply(propor_meth_14, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_14$meth_var <- var_cpg_14
hist(cpgslistdf_14$meth_var)
plot(density(cpgslistdf_14$meth_var))
#save cpglist with variance
write.table(cpgslistdf_14, file="cpg_snp_listdf_countorder_chr14_wvar.txt")
rm(rm_total_14, ym_total_14)

#load total read counts
rm_total_15 <- read.table("chr15_totalreadcounts.txt", row.names=NULL)
rm_total_15 <- rm_total_15[,-1]
#load meth read counts
ym_total_15 <- read.table("chr15_methreadcounts.txt", row.names=NULL)
ym_total_15 <- ym_total_15[,-1]
#get meth as proportion
propor_meth_15 <- ym_total_15/rm_total_15
#load cpg_snp df
cpgslistdf_15 <- read.table("cpg_snp_listdf_countorder_chr15.txt", header=TRUE)
#get variance for each row
var_cpg_15 <- apply(propor_meth_15, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_15$meth_var <- var_cpg_15
hist(cpgslistdf_15$meth_var)
plot(density(cpgslistdf_15$meth_var))
#save cpglist with variance
write.table(cpgslistdf_15, file="cpg_snp_listdf_countorder_chr15_wvar.txt")
rm(rm_total_15, ym_total_15)

#load total read counts
rm_total_16 <- read.table("chr16_totalreadcounts.txt", row.names=NULL)
rm_total_16 <- rm_total_16[,-1]
#load meth read counts
ym_total_16 <- read.table("chr16_methreadcounts.txt", row.names=NULL)
ym_total_16 <- ym_total_16[,-1]
#get meth as proportion
propor_meth_16 <- ym_total_16/rm_total_16
#load cpg_snp df
cpgslistdf_16 <- read.table("cpg_snp_listdf_countorder_chr16.txt", header=TRUE)
#get variance for each row
var_cpg_16 <- apply(propor_meth_16, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_16$meth_var <- var_cpg_16
hist(cpgslistdf_16$meth_var)
plot(density(cpgslistdf_16$meth_var))
#save cpglist with variance
write.table(cpgslistdf_16, file="cpg_snp_listdf_countorder_chr16_wvar.txt")
rm(rm_total_16, ym_total_16)

#load total read counts
rm_total_17 <- read.table("chr17_totalreadcounts.txt", row.names=NULL)
rm_total_17 <- rm_total_17[,-1]
#load meth read counts
ym_total_17 <- read.table("chr17_methreadcounts.txt", row.names=NULL)
ym_total_17 <- ym_total_17[,-1]
#get meth as proportion
propor_meth_17 <- ym_total_17/rm_total_17
#load cpg_snp df
cpgslistdf_17 <- read.table("cpg_snp_listdf_countorder_chr17.txt", header=TRUE)
#get variance for each row
var_cpg_17 <- apply(propor_meth_17, 1, function(x) var(x,na.rm=TRUE))
hist(var_cpg_17)
cpgslistdf_17$meth_var <- var_cpg_17
hist(cpgslistdf_17$meth_var)
plot(density(cpgslistdf_17$meth_var))
#save cpglist with variance
write.table(cpgslistdf_17, file="cpg_snp_listdf_countorder_chr17_wvar.txt")
rm(rm_total_17, ym_total_17)
var(as.numeric(propor_meth_17[1,]), na.rm = TRUE)

#load total read counts
rm_total_18 <- read.table("chr18_totalreadcounts.txt", row.names=NULL)
rm_total_18 <- rm_total_18[,-1]
#load meth read counts
ym_total_18 <- read.table("chr18_methreadcounts.txt", row.names=NULL)
ym_total_18 <- ym_total_18[,-1]
#get meth as proportion
propor_meth_18 <- ym_total_18/rm_total_18
#load cpg_snp df
cpgslistdf_18 <- read.table("cpg_snp_listdf_countorder_chr18.txt", header=TRUE)
#get variance for each row
var_cpg_18 <- apply(propor_meth_18, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_18$meth_var <- var_cpg_18
hist(cpgslistdf_18$meth_var)
plot(density(cpgslistdf_18$meth_var))
#save cpglist with variance
write.table(cpgslistdf_18, file="cpg_snp_listdf_countorder_chr18_wvar.txt")
rm(rm_total_18, ym_total_18)


#load total read counts
rm_total_19_a <- read.table("chr19_totalreadcounts_1.txt", row.names=NULL)
rm_total_19_a <- rm_total_19_a[,-1]
#load meth read counts
ym_total_19_a <- read.table("chr19_methreadcounts_1.txt", row.names=NULL)
ym_total_19_a <- ym_total_19_a[,-1]
#get meth as proportion
propor_meth_19_a <- ym_total_19_a/rm_total_19_a
#load cpg_snp df
cpgslistdf_19_a <- read.table("cpg_snp_listdf_countorder_chr19_1.txt", header=TRUE)
#get variance for each row
var_cpg_19_a <- apply(propor_meth_19_a, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_19_a$meth_var <- var_cpg_19_a
hist(cpgslistdf_19_a$meth_var)
plot(density(cpgslistdf_19_a$meth_var))
#save cpglist with variance
write.table(cpgslistdf_19_a, file="cpg_snp_listdf_countorder_chr19_1_wvar.txt")
rm(rm_total_19_a, ym_total_19_a)

#load total read counts
rm_total_19_b <- read.table("chr19_totalreadcounts_2.txt", row.names=NULL)
rm_total_19_b <- rm_total_19_b[,-1]
#load meth read counts
ym_total_19_b <- read.table("chr19_methreadcounts_2.txt", row.names=NULL)
ym_total_19_b <- ym_total_19_b[,-1]
#get meth as proportion
propor_meth_19_b <- ym_total_19_b/rm_total_19_b
#load cpg_snp df
cpgslistdf_19_b <- read.table("cpg_snp_listdf_countorder_chr19_2.txt", header=TRUE)
#get variance for each row
var_cpg_19_b <- apply(propor_meth_19_b, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_19_b$meth_var <- var_cpg_19_b
hist(cpgslistdf_19_b$meth_var)
plot(density(cpgslistdf_19_b$meth_var))
#save cpglist with variance
write.table(cpgslistdf_19_b, file="cpg_snp_listdf_countorder_chr19_2_wvar.txt")
rm(rm_total_19_b, ym_total_19_b)

#load total read counts
rm_total_19_c <- read.table("chr19_totalreadcounts_3.txt", row.names=NULL)
rm_total_19_c <- rm_total_19_c[,-1]
#load meth read counts
ym_total_19_c <- read.table("chr19_methreadcounts_3.txt", row.names=NULL)
ym_total_19_c <- ym_total_19_c[,-1]
#get meth as proportion
propor_meth_19_c <- ym_total_19_c/rm_total_19_c
#load cpg_snp df
cpgslistdf_19_c <- read.table("cpg_snp_listdf_countorder_chr19_3.txt", header=TRUE)
#get variance for each row
var_cpg_19_c <- apply(propor_meth_19_c, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_19_c$meth_var <- var_cpg_19_c
hist(cpgslistdf_19_c$meth_var)
plot(density(cpgslistdf_19_c$meth_var))
#save cpglist with variance
write.table(cpgslistdf_19_c, file="cpg_snp_listdf_countorder_chr19_3_wvar.txt")
rm(rm_total_19_c, ym_total_19_c)

#load total read counts
rm_total_19_d <- read.table("chr19_totalreadcounts_4.txt", row.names=NULL)
rm_total_19_d <- rm_total_19_d[,-1]
#load meth read counts
ym_total_19_d <- read.table("chr19_methreadcounts_4.txt", row.names=NULL)
ym_total_19_d <- ym_total_19_d[,-1]
#get meth as proportion
propor_meth_19_d <- ym_total_19_d/rm_total_19_d
#load cpg_snp df
cpgslistdf_19_d <- read.table("cpg_snp_listdf_countorder_chr19_4.txt", header=TRUE)
#get variance for each row
var_cpg_19_d <- apply(propor_meth_19_d, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_19_d$meth_var <- var_cpg_19_d
hist(cpgslistdf_19_d$meth_var)
plot(density(cpgslistdf_19_d$meth_var))
#save cpglist with variance
write.table(cpgslistdf_19_d, file="cpg_snp_listdf_countorder_chr19_4_wvar.txt")
rm(rm_total_19_d, ym_total_19_d)

#load total read counts
rm_total_19_e <- read.table("chr19_totalreadcounts_5.txt", row.names=NULL)
rm_total_19_e <- rm_total_19_e[,-1]
#load meth read counts
ym_total_19_e <- read.table("chr19_methreadcounts_5.txt", row.names=NULL)
ym_total_19_e <- ym_total_19_e[,-1]
#get meth as proportion
propor_meth_19_e <- ym_total_19_e/rm_total_19_e
#load cpg_snp df
cpgslistdf_19_e <- read.table("cpg_snp_listdf_countorder_chr19_5.txt", header=TRUE)
#get variance for each row
var_cpg_19_e <- apply(propor_meth_19_e, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_19_e$meth_var <- var_cpg_19_e
hist(cpgslistdf_19_e$meth_var)
plot(density(cpgslistdf_19_e$meth_var))
#save cpglist with variance
write.table(cpgslistdf_19_e, file="cpg_snp_listdf_countorder_chr19_5_wvar.txt")
rm(rm_total_19_e, ym_total_19_e)

#load total read counts
rm_total_20 <- read.table("chr20_totalreadcounts.txt", row.names=NULL)
rm_total_20 <- rm_total_20[,-1]
#load meth read counts
ym_total_20 <- read.table("chr20_methreadcounts.txt", row.names=NULL)
ym_total_20 <- ym_total_20[,-1]
#get meth as proportion
propor_meth_20 <- ym_total_20/rm_total_20
#load cpg_snp df
cpgslistdf_20 <- read.table("cpg_snp_listdf_countorder_chr20.txt", header=TRUE)
#get variance for each row
var_cpg_20 <- apply(propor_meth_20, 1, function(x) var(x,na.rm=TRUE))
cpgslistdf_20$meth_var <- var_cpg_20
hist(cpgslistdf_20$meth_var)
plot(density(cpgslistdf_20$meth_var))
#save cpglist with variance
write.table(cpgslistdf_20, file="cpg_snp_listdf_countorder_chr20_wvar.txt")
rm(rm_total_20, ym_total_20)

var_files <- list.files(pattern="*wvar.txt")
var_files_df <- lapply(var_files, read.table)
var_files_df_combine <- do.call(rbind, var_files_df)
hist(var_files_df_combine$meth_var)
plot(density(var_files_df_combine$meth_var))
mean(var_files_df_combine$meth_var)
#0.02922362
median(var_files_df_combine$meth_var)
#0.01970066
sd(var_files_df_combine$meth_var)
#0.02780466
