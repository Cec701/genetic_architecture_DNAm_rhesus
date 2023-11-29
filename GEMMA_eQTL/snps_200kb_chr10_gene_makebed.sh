#!/bin/sh

#SBATCH -A ccosta8
#SBATCH -o snps_200kb_chr10_gene.%j.out
#SBATCH -e snps_200kb_chr10_gene.%j.err 
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cec701@nyu.edu

# Run separate for each chromosome 1-20 (Example chr10)
# Array = 1-number of genes for that chromosome 
module load plink/1.9.0 

# create file with genomic range of 200kb from start and end of macaque gene
sed -n ${SLURM_ARRAY_TASK_ID}p GE_genes_range_file_200kb_chr10.txt > GE_genes_range_file_200kb_chr10_gene${SLURM_ARRAY_TASK_ID}.txt

# Extract all SNPs that fall within 200kb of that gene 
plink --bfile merged.chr10.GP.HW.MAF.GE.with_repeats.sorted --extract range GE_genes_range_file_200kb_chr10_gene${SLURM_ARRAY_TASK_ID}.txt --make-bed --out chr10/merged.chr10.GP.HW.MAF.GE.with_repeats.sorted.gene${SLURM_ARRAY_TASK_ID}

