#!/bin/sh
## If a site is less than genotype probability 90% for an individual, set it to '.' (missing) for that individual

module load bcftools/1.10.2

## Filter for sites with high genotype probability from imputation and annotate ID column-need for plink LD filtering
bcftools filter -i 'FMT/GP>0.9' --set-GTs . /scratch/mwatowic/gxe_project/gxe_modeling/cayo_vcfs_impute/merged_imputed/merged.chr${SLURM_ARRAY_TASK_ID}.vcf.gz | \
    bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o merged.chr${SLURM_ARRAY_TASK_ID}.GP90.vcf.gz

