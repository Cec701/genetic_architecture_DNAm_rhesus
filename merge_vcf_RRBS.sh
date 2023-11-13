#!/bin/sh

#SBATCH -A ccosta8
#SBATCH -o RRBS_Merge_Job.%j.out
#SBATCH -e RRBS_Merge_Job.%j.err 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cec701@nyu.edu

# Merge all VCF files from RRBS genotyping 

module load bcftools/1.10.2

bcftools merge /scratch/ccosta8/RRBS/Output/*.R1_CGmap.PASS2.DP5.vcf.gz --output /scratch/ccosta8/RRBS/Output/all_RRBS_samples.R1_CGmap.PASS2.DP5.vcf.gz --output-type z 

