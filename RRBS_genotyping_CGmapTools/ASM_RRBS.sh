#!/bin/sh

#SBATCH -A ccosta8
#SBATCH -o RRBS_ASM_Job.%j.out
#SBATCH -e RRBS_ASM_Job.%j.err 
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cec701@nyu.edu

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p sample_names`


# Array job : Allele-specific methylation

path_cgmap=/scratch/ccosta8/RRBS/Programs/xz_mod_cgmaptools/cgmaptools
path_ref=/scratch/ccosta8/RRBS/Genomes/Macaca_mulatta.Mmul_10.dna.toplevel.fa
path_bam2=/scratch/ccosta8/RRBS/Intermediate/${sample}.R1_val_1_bismark_bt2_pe.sort.bam

# extract sites to keep
module load vcftools/0.1.12b
module load samtools/1.9

samtools index $path_bam2

gunzip /scratch/ccosta8/RRBS/Output/${sample}.R1_CGmap.PASS2.DP5.vcf.gz

# call allele-specific methylation/create asm files
$path_cgmap asm -r $path_ref -b $path_bam2 -l /scratch/ccosta8/RRBS/Output/${sample}.R1_CGmap.PASS2.DP5.vcf -m ass > /scratch/ccosta8/RRBS/Output/ASM/${sample}.R1_CGmap.PASS2.DP5.ASM
