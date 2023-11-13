#!/bin/sh

#SBATCH -A ccosta8
#SBATCH -o RRBS_ASM_Job.%j.out
#SBATCH -e RRBS_ASM_Job.%j.err 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cec701@nyu.edu

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p sample_names`




# Array job : Allele-specific methylation


path_cgmap=/scratch/ccosta8/RRBS/Programs/xz_mod_cgmaptools/cgmaptools
path_ref=/scratch/ccosta8/RRBS/Genomes/Macaca_mulatta.Mmul_10.dna.toplevel.fa
path_bam2=/scratch/nsnyderm/meQTL/bams/${sample}.R1_val_1_bismark_bt2_pe.bam

# extract sites to keep
module load vcftools (#add version)
module load samtools/1.9

samtools index $path_bam2

gunzip /scratch/ccosta8/RRBS/Output/${indv}.R1_CGmap.PASS2.DP5.vcf.gz

$path_cgmap asm -r $path_ref -b $path_bam2 -l /scratch/ccosta8/RRBS/Output/${indv}.R1_CGmap.PASS2.DP5.vcf -m ass > /scratch/ccosta8/RRBS/Output/ASM/${indv}.R1_CGmap.PASS2.DP5.ASM
