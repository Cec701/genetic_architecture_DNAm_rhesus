#!/bin/bash
## written for use on the ASU HPC (Agave) 

##SBATCH --mail-type=END,FAIL                  # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=150gb                            # Job memory request
#SBATCH --time=0-4                             # Time limit days-hrs
#SBATCH --cpus-per-task=1                      # cpus needed (so how many parallel threads you want), 1 per chrom

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p /data/CEM/smacklab/data/mapped_reads/cayo_WGS/samples_cayo`
chrom=$1

module load speedseq/0.1.2
module load samtools/1.9
module load singularity/3.6.3
module load bcftools/1.10.2
module load tabix/0.2.6

singularity exec -B /scratch/mwatowic/cayo_WGS_impute loimpute_latest.sif loimpute \
    -ne 1000 \
    -id ${sample} \
    -i pileups/${sample}.${chrom}.gz \
    -h mgap_new/mgap.chr${chrom}.phased.AF0.vcf.gz \
    -o imputed/${sample}_${chrom}

tabix -p vcf imputed/${sample}_${chrom}.vcf.gz
