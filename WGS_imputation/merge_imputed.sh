#!/bin/bash

module load bcftools/1.10.2
module load tabix/0.2.6

chr=`sed -n ${SLURM_ARRAY_TASK_ID}p chroms` 

bcftools merge imputed_filtered/*_${chr}.vcf.gz | bgzip -c > imputed_merge/all_cayo_wgs_${chr}.vcf.gz
tabix -p vcf imputed_merge/all_cayo_wgs_${chr}.vcf.gz

