# Individual Variant Calling (SNPs from Cayo RRBS BAM Files) 

#!/bin/sh 

#SBATCH -A ccosta8
#SBATCH -o RRBS_SNP_Job.%j.out
#SBATCH -e RRBS_SNP_Job.%j.err 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cec701@nyu.edu

#run below for each sampleID name in sample_names list (replace sampleID in each instance of ${sample})
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p sample_names`

module load samtools/1.9

path_cgmap=/scratch/ccosta8/RRBS/Programs/xz_mod_cgmaptools/cgmaptools
path_ref=/scratch/ccosta8/RRBS/Genomes/Macaca_mulatta.Mmul_10.dna.toplevel.fa

path_bam=/scratch/nsnyderm/meQTL/bams/${sample}.R1_val_1_bismark_bt2_pe.bam #(cycle through all sample IDs)
path_bam_sort=/scratch/ccosta8/RRBS/Intermediate/${sample}.R1_val_1_bismark_bt2_pe.sort.bam #(what the name of the file output should be once you do the bam file sort step)

path_out=/scratch/ccosta8/RRBS/Output/${sample}.R1_CGmap #(for all the outputs) 

#sort BAM files

samtools sort -o $path_bam_sort $path_bam
samtools view $path_bam | wc -l

#convert BAM to CGmap format 
$path_cgmap convert bam2cgmap -b $path_bam_sort -g $path_ref -o $path_out

#call SNPs from CGmap file
$path_cgmap snv -i $path_out.ATCGmap.gz -m bayes -v $path_out.vcf --bayes-dynamicP -o $path_out.out --bayes-e=0.01 -a

#filter for DP>4 (read depth), remove "VAGUE" calls, and only include variants also present in Cayo genetic dataset (cayo.PASS.vcf.gz)
module bcftools/1.10.2
module load tabix/0.2.6

bgzip -c $path_out.vcf > $path_out.vcf.gz
tabix -p vcf $path_out.vcf.gz
bcftools view -f 'PASS,.' $path_out.vcf.gz --output-type z > $path_out.PASS.vcf.gz
tabix -p vcf $path_out.PASS.vcf.gz
bcftools filter -i 'FORMAT/DP>4' $path_out.PASS.vcf.gz --output $path_out.PASS.DP5.vcf.gz --output-type z --regions-file /scratch/nsnyderm/meQTL/cayo.PASS.vcf.gz
tabix -p vcf $path_out.PASS.DP5.vcf.gz

rm $path_out.PASS.vcf.gz
rm $path_out.out
rm $path_out.ATCGmap.gz
