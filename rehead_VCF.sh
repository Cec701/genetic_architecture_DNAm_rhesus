#!/bin/sh

#SBATCH -A ccosta8
#SBATCH -o RRBS_Head_Job.%j.out
#SBATCH -e RRBS_Head_Job.%j.err 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cec701@nyu.edu


# Reheader each vcf file 
# replace NA001 with Sample ID 

module load bcftools/1.10.2

for indv in `cat sample_names`; do rm temp_${indv}.txt; touch temp_${indv}.txt; echo $indv >> temp_${indv}.txt;
bcftools reheader -s temp_${indv}.txt /scratch/ccosta8/RRBS/Output/${indv}.R1_CGmap.PASS.DP5.vcf.gz > /scratch/ccosta8/RRBS/Output/${indv}.R1_CGmap.PASS2.DP5.vcf.gz;
tabix -p vcf /scratch/ccosta8/RRBS/Output/${indv}.R1_CGmap.PASS2.DP5.vcf.gz; rm temp_${indv}.txt; done
