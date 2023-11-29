#!/bin/sh

#SBATCH -A ccosta8
#SBATCH -o chr1_perm_gemma.%j.out
#SBATCH -e chr1_perm_gemma.%j.err 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cec701@nyu.edu

# run eQTL GEMMA permutation chromsome 1

FILE=merged.chr1.GP.HW.MAF.GE.with_repeats.sorted.gene${SLURM_ARRAY_TASK_ID}.bed

for i in {1..10}; do 
	cd perm$i
	if [ -f "$FILE" ]; then
		../../../../../../Programs/./gemma-0.98.5-linux-static-AMD64 -bfile merged.chr1.GP.HW.MAF.GE.with_repeats.sorted.gene${SLURM_ARRAY_TASK_ID} -k ../../merged.GP.HW.MAF.GE.with_repeats.sorted.all_chr.GRM.cXX.txt -lmm 1 -n ${SLURM_ARRAY_TASK_ID} -c ../../sorted_ge_wgs_GEMMA_covariates.txt -o gemma_eqtl_chr1_gene${SLURM_ARRAY_TASK_ID}_perm$i
	else
		echo "$FILE does not exist."
	fi
	cd ..
done
