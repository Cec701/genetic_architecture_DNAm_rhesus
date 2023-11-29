#!/bin/sh

#SBATCH -A ccosta8
#SBATCH -o chr1_permutations_makefiles.%j.out
#SBATCH -e chr1_permutations_makefiles.%j.err 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cec701@nyu.edu

# eQTL GEMMA permutations: make permuted files 

NUM=$(( ($SLURM_ARRAY_TASK_ID + 5) )) 
echo $NUM

for i in {1..10}; do 
	cd perm$i
	awk -v var="$NUM" '{print $var}' ../merged.chr1.GP.HW.MAF.GE.with_repeats.sorted.gene${SLURM_ARRAY_TASK_ID}.fam | shuf > merged.chr1.GP.HW.MAF.GE.with_repeats.sorted.gene${SLURM_ARRAY_TASK_ID}_perm$i.txt
	awk -v var="$NUM" 'NR==FNR{a[NR]=$0;next}{$var=a[FNR]}1' merged.chr1.GP.HW.MAF.GE.with_repeats.sorted.gene${SLURM_ARRAY_TASK_ID}_perm$i.txt ../merged.chr1.GP.HW.MAF.GE.with_repeats.sorted.gene${SLURM_ARRAY_TASK_ID}.fam > merged.chr1.GP.HW.MAF.GE.with_repeats.sorted.gene${SLURM_ARRAY_TASK_ID}.fam
	cd ..
done
