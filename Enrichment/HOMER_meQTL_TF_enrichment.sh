#!/bin/sh


#HOMER TF Enrichment Analysis 
# v4.11

findMotifsGenome.pl IMAGE_SNP_sig_0.05_200bp_region_unique.bed Genome/Macaca_mulatta.Mmul_10.dna.toplevel.fa homer_unique_cpgflag_meqtl_FINAL/ 
-size 200 -cpg -bg image_all_res_unique_SNPs_homer.bed -mset vertebrates
