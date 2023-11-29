#!/bin/sh

#SBATCH -A ccosta8
#SBATCH -o ge_pheno.%j.out
#SBATCH -e ge_pheno.%j.err 
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cec701@nyu.edu

module load r/4.0.2

#first argument is chr number 

Rscript add_ge_SVA_res_phenotypes.R 1
