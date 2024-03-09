#!/bin/bash
## written for use on the ASU HPC (Agave)

#SBATCH -o mpileup.%j.out
#SBATCH -e mpileup.%j.err
#SBATCH -A ccosta8
#SBATCH --mem=100gb                             # Job memory request
#SBATCH --time=0-2                              # Time limit days-hrs
#SBATCH --cpus-per-task=21                      # cpus needed (so how many parallel threads you want), 1 per chrom

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/ccosta8/RRBS/Output/GxE/WGS_pipeline/aln/rerun_sample_names.txt`
bam=/scratch/ccosta8/RRBS/Output/GxE/WGS_pipeline/aln/${sample}.mmul10.bam

module load speedseq/0.1.2
module load samtools/1.9

cat chroms | parallel --verbose -j 20 "samtools mpileup -r {1} -l /scratch/mwatowic/impute/mgap_new/mgap.chr{1}.phased.AF0.vcf.gz \
    ${bam} | gzip -c > pileups/${sample}.{1}.gz"

