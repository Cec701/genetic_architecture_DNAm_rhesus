#!/bin/bash
module load anaconda/py3
source activate bismark
module load trim_galore/0.4.0

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p libraries` #libraries is a file with one library per line

genome_path=/data/CEM/smacklab/genomes/bismark/mmul10/
fastq_path=/scratch/nsnyderm/sge_rrbs/fastqs
out_path=/scratch/mwatowic/dnam/RRBS/bams
cov_path=/scratch/mwatowic/dnam/RRBS/cov_files

## trim adaptors
trim_galore --rrbs --non_directional -o ${fastq_path} --gzip ${fastq_path}/${sampleID}*R1*fastq.gz

## map the reads
bismark --genome_folder ${genome_path} --non_directional -o ${out_path} --score_min L,0,-0.6 -R 12 --parallel 4 ${fastq_path}/${sampleID}*trimmed.fq.gz

## extract methylation data
bismark_methylation_extractor -p -o ${cov_path} --no_overlap --bedGraph --comprehensive --parallel 24 --merge_non_CpG --genome_folder ${genome_path} ${out_path}/${sampleID}*.bam

conda deactivate
## sbatch --exclusive -p mrline -q wildfire -t 24:00:00 --array=1-581%50 trim_map.sh