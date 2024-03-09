#!/bin/bash
#SBATCH -G a100:2
#SBATCH -t 0-4
#SBATCH --mem=100G
#SBATCH --cpus-per-task=12

genome=mmul10 #shorthand name for genome (for bam name)
genome_path=/data/CEM/smacklab/genomes/bwa/Mmul_10 #path to bwa indexed genome
genome_prefix=Macaca_mulatta.Mmul_10.dna.toplevel.fa #indexed genome prefix
output_path=aln #path to directory for alignments
fastq_path=${1}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${fastq_path}/samples`

echo $sample
echo $fastq_path

#cat ${fastq_path}/${sample}*R1*gz > merged_fqs/${sample}.R1.fastq.gz
#cat ${fastq_path}/${sample}*R2*gz > merged_fqs/${sample}.R2.fastq.gz

r1=`ls ${fastq_path}/${sample}*R1*.fastq.gz`
r2=`ls ${fastq_path}/${sample}*R2*.fastq.gz`

apptainer exec --nv -B ${genome_path},${fastq_path},${PWD},/scratch/nsnyderm/cayo_wgs /packages/apps/simg/parabricks-4.0.sif\
	pbrun fq2bam \
	--ref ${genome_path}/${genome_prefix} \
	--in-fq ${r1} ${r2} \
	--read-group-sm ${sample} \
	--out-bam ${output_path}/${sample}.${genome}.bam

module load samtools-1.16-gcc-11.2.0

samtools coverage -r 1  \
	-o ${sample}.cov \
	${output_path}/${sample}.${genome}.bam
