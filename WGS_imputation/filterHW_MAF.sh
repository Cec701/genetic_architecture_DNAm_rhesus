#!/bin/sh
module load plink/1.9.0
module load bedtools2/2.24.0
module load tabix/0.2.6

## Filter out sites based on MAF, HWE
plink --vcf merged.chr${SLURM_ARRAY_TASK_ID}.GP90.vcf.gz --maf 0.05 --hwe 10e-8 --out merged.chr${SLURM_ARRAY_TASK_ID}.GP.HW.MAF --recode vcf

bgzip merged.chr${SLURM_ARRAY_TASK_ID}.GP.HW.MAF.vcf

