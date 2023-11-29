#GSEA(Gene Set Enrichment) Analysis of sig meQTL proportion genes 
#Genes that are enriched for nearby meqtl or depleted 
#Sorted by proportion of signifcant CpGs/total tested CpGs (control for number of tested cpgs near each gene)


BiocManager::install("fgsea")
library(fgsea)
install.packages("msigdbr")
library(msigdbr)
library(qvalue)
library(readr)
library(data.table)

BiocManager::install("umap")

empirical_proportion_genes_protein_coding_order_2 <- read.table("sigmeqtl_proportion_50kb_proteingenes_sorted_forGSEA.txt", header=TRUE)

##################use the ensembl IDs of hallmark dataset#################

hallmark.msigdb = msigdbr(species = "Macaca mulatta", category = "H")
hallmark.msigdb = split(x = hallmark.msigdb$ensembl_gene, f = hallmark.msigdb$gs_name)

GO.msigdb.BP = msigdbr(species = "Macaca mulatta", category = "C5", subcategory = "GO:BP")
GO.msigdb.BP = split(x = GO.msigdb.BP$ensembl_gene, f = GO.msigdb.BP$gs_name)

GO.msigdb.C5 = msigdbr(species = "Macaca mulatta", category = "C5")
GO.msigdb.C5 = split(x = GO.msigdb.C5$ensembl_gene, f = GO.msigdb.C5$gs_name)

GO.msigdb.MF = msigdbr(species = "Macaca mulatta", category = "C5", subcategory = "GO:MF")
GO.msigdb.MF = split(x = GO.msigdb.MF$ensembl_gene, f = GO.msigdb.MF$gs_name)

GO.msigdb.CC = msigdbr(species = "Macaca mulatta", category = "C5", subcategory = "GO:CC")
GO.msigdb.CC = split(x = GO.msigdb.CC$ensembl_gene, f = GO.msigdb.CC$gs_name)

GO.msigdb.HPO = msigdbr(species = "Macaca mulatta", category = "C5", subcategory = "HPO")
GO.msigdb.HPO = split(x = GO.msigdb.HPO$ensembl_gene, f = GO.msigdb.HPO$gs_name)

GO.c2 = msigdbr(species = "Macaca mulatta", category = "C2", subcategory = "CP:KEGG")
GO.c2 = split(x = GO.c2$ensembl_gene, f = GO.c2$gs_name)

GO.c7 = msigdbr(species = "Macaca mulatta", category = "C7")
GO.c7 = split(x = GO.c7$ensembl_gene, f = GO.c7$gs_name)

GO.c8 = msigdbr(species = "Macaca mulatta", category = "C8")
GO.c8 = split(x = GO.c8$ensembl_gene, f = GO.c8$gs_name)

#################meqtl genes####################
# rank genes all, using proportion sig meQTL
#vector of proportions 
meqtl_genes_proportionsig <- as.vector(empirical_proportion_genes_protein_coding_order_2$proportion_sig)
names(meqtl_genes_proportionsig) <- empirical_proportion_genes_protein_coding_order_2$ID

#GSEA, all genes and proportions
#library(fgsea)
set.seed(7)

fgsea.proportionssig.hallmark <- fgsea(pathways = hallmark.msigdb, #GO.c2, #hallmark.msigdb, GO.msigdb
                     stats    = meqtl_genes_proportionsig,
                     minSize  = 10,
                     maxSize  = 500, #2000 for the GO pathways
                     eps = 0.0)

library(dplyr)
cutoff=.2
meqtl_genes_proportionsig_sig_hallmark <- fgsea.proportionssig.hallmark %>% filter(padj<cutoff)
str_remove(meqtl_genes_proportionsig_sig$pathway, "HALLMARK_")

##
fgsea.proportionssig.GO.BP <- fgsea(pathways = GO.msigdb.BP, #GO.c2, #hallmark.msigdb, GO.msigdb
                              stats    = meqtl_genes_proportionsig,
                              minSize  = 10,
                              maxSize  = 2000, #2000 for the GO pathways
                              eps = 0.0)


cutoff=.2
meqtl_genes_proportionsig_sig_GO.BP_0.2 <- subset(fgsea.proportionssig.GO.BP, padj<cutoff)
write.table(meqtl_genes_proportionsig_sig_GO.BP_0.2[,1:7], file="gsea_go_bp_fdr0.2.txt", row.names = FALSE)

##
fgsea.proportionssig.GO.MF <- fgsea(pathways = GO.msigdb.MF, #GO.c2, #hallmark.msigdb, GO.msigdb
                                    stats    = meqtl_genes_proportionsig,
                                    minSize  = 10,
                                    maxSize  = 2000, #2000 for the GO pathways
                                    eps = 0.0)


cutoff=.2
meqtl_genes_proportionsig_sig_GO.MF <- subset(fgsea.proportionssig.GO.MF, padj<cutoff)
write.table(meqtl_genes_proportionsig_sig_GO.MF[,1:7], file="gsea_go_mf_fdr0.2.txt", row.names = FALSE)

##
fgsea.GO.c7.proportionssig <- fgsea(pathways = GO.c7, #GO.c2, #hallmark.msigdb, GO.msigdb
                               stats    = meqtl_genes_proportionsig,
                               minSize  = 10,
                               maxSize  = 2000, #2000 for the GO pathways
                               eps = 0.0)


cutoff=.2
meqtl_genes_proportionsig_sig_GO.c7_0.2 <- subset(fgsea.GO.c7.proportionssig, padj<cutoff)
write.table(meqtl_genes_proportionsig_sig_GO.c7_0.2[,1:7], file="gsea_go_C7_fdr0.2.txt", row.names = FALSE)

###
fgsea.GO.c2.proportionssig <- fgsea(pathways = GO.c2, #GO.c2, #hallmark.msigdb, GO.msigdb
                                 stats    = meqtl_genes_proportionsig,
                                 minSize  = 15,
                                 maxSize  = 2000, #2000 for the GO pathways
                                 eps = 0.0)

cutoff=.2
meqtl_genes_proportionsig_sig_GO.c2 <- fgsea.GO.c2.proportionssig %>% filter(padj<cutoff)



##GO.c2
fgsea.GO.c2.proportionssig <- fgsea(pathways = GO.c2, #GO.c2, #hallmark.msigdb, GO.msigdb
                                    stats    = meqtl_genes_proportionsig,
                                    minSize  = 15,
                                    maxSize  = 2000, #2000 for the GO pathways
                                    eps = 0.0)

cutoff=.2
meqtl_genes_proportionsig_sig_GO.c2 <- fgsea.GO.c2.proportionssig %>% filter(padj<cutoff)

library(ggplot2)
#plotEnrichment(GO.c7[["GSE17721_LPS_VS_CPG_1H_BMDC_UP"]], meqtl_genes_proportionsig) + labs(title="Up-regulated in dendritic cells (DC) stimulated with LPS vs CpG DNA (TLR9 agonist)")
#plotEnrichment(GO.msigdb.BP[["GOBP_REGULATION_OF_DEFENSE_RESPONSE_TO_VIRUS"]], meqtl_genes_proportionsig) + labs(title="Regulation of defense response to virus")
#plotEnrichment(GO.msigdb.BP[["GOBP_CYTOPLASMIC_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY"]], meqtl_genes_proportionsig) + labs(title="Cytoplasmic pattern recognition signaling")
#plotGseaTable(GO.msigdb[topPathwaysUp], meqtl_genes_proportionsig, fgsea.proportionssig.GO, 
#              gseaParam=0.5)
plotEnrichment(GO.msigdb.BP[["GOBP_REGULATION_OF_NEUROBLAST_PROLIFERATION"]], meqtl_genes_proportionsig) + labs(title="Regulation of neuroblast proliferation")
#
#(GO.msigdb[["GOBP_CYTOPLASMIC_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY"]], meqtl_genes_proportionsig)

##GO.c8
fgsea.GO.c8.proportionssig <- fgsea(pathways = GO.c8, #GO.c2, #hallmark.msigdb, GO.msigdb
                                    stats    = meqtl_genes_proportionsig,
                                    minSize  = 10,
                                    maxSize  = 2000, #2000 for the GO pathways
                                    eps = 0.0)


cutoff=.2
meqtl_genes_proportionsig_sig_GO.c8_0.2 <- fgsea.GO.c8.proportionssig %>% filter(padj<cutoff)
write.table(meqtl_genes_proportionsig_sig_GO.c8_0.2[,1:7], file="gsea_go_C8_fdr0.2.txt", row.names = FALSE)


meqtl_genes_proportionsig_sig_GO.MF$pathway <- gsub("GOMF_", "", meqtl_genes_proportionsig_sig_GO.MF$pathway)
meqtl_genes_proportionsig_sig_GO_0.2$pathway <- gsub("GOBP_", "", meqtl_genes_proportionsig_sig_GO_0.2$pathway)
meqtl_genes_proportionsig_sig_GO.BP_0.2$pathway <- gsub("GOBP_", "", meqtl_genes_proportionsig_sig_GO.BP_0.2$pathway)

#GO.msigdb.HPO

fgsea.GO.HPO.proportionssig <- fgsea(pathways = GO.msigdb.HPO, #GO.c2, #hallmark.msigdb, GO.msigdb
                                    stats    = meqtl_genes_proportionsig,
                                    minSize  = 15,
                                    maxSize  = 2000, #2000 for the GO pathways
                                    eps = 0.0)


meqtl_genes_proportionsig_sig_GO.HPO <- fgsea.GO.HPO.proportionssig  %>% filter(padj<cutoff)
#20%FDR

#GO.msigdb.CC

fgsea.GO.CC.proportionssig <- fgsea(pathways = GO.msigdb.CC, #GO.c2, #hallmark.msigdb, GO.msigdb
                                     stats    = meqtl_genes_proportionsig,
                                     minSize  = 15,
                                     maxSize  = 2000, #2000 for the GO pathways
                                     eps = 0.0)

cutoff=.2
meqtl_genes_proportionsig_sig_GO.CC_0.2 <- fgsea.GO.CC.proportionssig %>% filter(padj<cutoff)


#GO.msigdb.C5 
fgsea.GO.C5.proportionssig <- fgsea(pathways = GO.msigdb.C5, #GO.c2, #hallmark.msigdb, GO.msigdb
                                    stats    = meqtl_genes_proportionsig,
                                    minSize  = 15,
                                    maxSize  = 2000, #2000 for the GO pathways
                                    eps = 0.0)

cutoff=.2
meqtl_genes_proportionsig_sig_GO.C5_0.2 <- fgsea.GO.C5.proportionssig %>% filter(padj<cutoff)


library(ggplot2)
library(dplyr)
meqtl_genes_proportionsig_sig_GO.c7_0.2$pathway <- gsub("GSE[0-9]+_","", meqtl_genes_proportionsig_sig_GO.c7_0.2$pathway)
meqtl_genes_proportionsig_sig_GO.BP_0.2$pathway <- gsub("GOBP_","", meqtl_genes_proportionsig_sig_GO.BP_0.2$pathway)

fgseaResTidy <- meqtl_genes_proportionsig_sig_GO.c7_0.2 %>%  as_tibble() %>%  arrange(desc(NES))
fgseaResTidy_immunesignature_filt <- fgseaResTidy[1:15,]
ggplot(fgseaResTidy_immunesignature_filt, aes(reorder(pathway, NES), NES, fill=padj)) +
  geom_col() +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Top Immunologic Signatures pathways GSEA") + 
  theme_minimal() + scale_fill_gradient(low="orange", high="navy blue") + theme(text = element_text(size=15))

fgseaResTidy_GOBP <- meqtl_genes_proportionsig_sig_GO.BP_0.2 %>%  as_tibble() %>%  arrange(desc(NES))
fgseaResTidy_GOBP_filt <- fgseaResTidy_GOBP[1:15,]
ggplot(fgseaResTidy_GOBP_filt, aes(reorder(pathway, NES), NES, fill=padj)) +
  geom_col() +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=" Top GO:BP pathways GSEA") + 
  theme_minimal() + scale_fill_gradient(low="orange", high="navy blue") + theme(text = element_text(size=15))

fgseaResTidy_GOMF <- meqtl_genes_proportionsig_sig_GO.MF %>%  as_tibble() %>%  arrange(desc(NES))
fgseaResTidy_GOMF_filt <- fgseaResTidy_GOMF[1:15,]
ggplot(fgseaResTidy_GOMF_filt, aes(reorder(pathway, NES), NES, fill=padj)) +
  geom_col() +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Top GO:MF pathways GSEA") + 
  theme_minimal() + scale_fill_gradient(low="orange", high="navy blue") + theme(text = element_text(size=15))

fgseaResTidy_celltype <- meqtl_genes_proportionsig_sig_GO.c8_0.2 %>%  as_tibble() %>%  arrange(desc(NES))
fgseaResTidy_celltype_filt <- fgseaResTidy_celltype[1:15,]
ggplot(fgseaResTidy_celltype_filt, aes(reorder(pathway, NES), NES, fill=padj)) +
  geom_col() +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Cell type pathways GSEA") + 
  theme_minimal() + scale_fill_gradient(low="orange", high="navy blue") + theme(text = element_text(size=15))

fgseaResTidy_CC <- meqtl_genes_proportionsig_sig_GO.CC %>%  as_tibble() %>%  arrange(desc(NES))
ggplot(fgseaResTidy_CC, aes(reorder(pathway, NES), NES)) +
  geom_col() +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO:CC pathways NES from GSEA") + 
  theme_minimal()


#need to 
fgseaResTidy_GOBP_0.2 <- meqtl_genes_proportionsig_sig_GO_0.2 %>%  as_tibble() %>%  arrange(desc(NES))
fgseaResTidy_GOBP_0.2_filt <- fgseaResTidy_GOBP_0.2[1:50,]

ggplot(fgseaResTidy_GOBP_0.2_filt, aes(reorder(pathway, NES), NES)) +
  geom_col() +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Top GO:BP pathways NES from GSEA (FDR 20%)") + 
  theme_minimal()


