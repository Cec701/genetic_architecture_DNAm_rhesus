#!/bin/sh

#SBATCH -A ccosta8
#SBATCH -o RRBS_Head_Job.%j.out
#SBATCH -e RRBS_Head_Job.%j.err 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cec701@nyu.edu


# Reheader each vcf file 
# replace NA001 with Sample ID 

module load bcftools/1.10.2

for indv in `cat sample_names`; do rm temp_${indv}.txt; touch temp_${indv}.txt; echo $indv >> temp_${indv}.txt;
bcftools reheader -s temp_${indv}.txt /scratch/ccosta8/RRBS/Output/${indv}.R1_CGmap.PASS.DP5.vcf.gz > /scratch/ccosta8/RRBS/Output/${indv}.R1_CGmap.PASS2.DP5.vcf.gz;
tabix -p vcf /scratch/ccosta8/RRBS/Output/${indv}.R1_CGmap.PASS2.DP5.vcf.gz; rm temp_${indv}.txt; done



for f in *.g.vcf.gz; do d=$(echo "$f" | sed -e 's/\.g.vcf.gz$//');echo "filename: $f : $d";done





for indv in `cat wgs_animal_ids.txt`; do rm temp_${indv}.txt; touch temp_${indv}.txt; echo $indv >> temp_${indv}.txt; bcftools reheader -s temp_${indv}.txt ${indv}.g.vcf.gz > ${indv}.rehead.vcf.gz; tabix -p vcf ${indv}.rehead.vcf.gz; rm temp_${indv}.txt; done

for indv in `cat new_wgs_animal_ID.txt`; do rm temp_${indv}.txt; touch temp_${indv}.txt; echo $indv >> temp_${indv}.txt; bcftools reheader -s temp_${indv}.txt ${indv}.g.vcf.gz > ${indv}.rehead.vcf.gz; tabix -p vcf ${indv}.rehead.vcf.gz; rm temp_${indv}.txt; done
for indv in `cat new_wgs_animal_ID.txt`; do touch temp_${indv}.txt; echo $indv >> temp_${indv}.txt; bcftools reheader -s temp_${indv}.txt ${indv}.g.vcf.gz > ${indv}.rehead.vcf.gz; tabix -p vcf ${indv}.rehead.vcf.gz; rm temp_${indv}.txt; done

#why is 44T
#and 4C1 in there twice? 
#doesn't explain the error for 22R though which is only there once.. 
#let's just rehead everything.. 

for indv in `cat all_WGS_animalID_feb2023.txt`; do touch temp_${indv}.txt; echo $indv >> temp_${indv}.txt; bcftools reheader -s temp_${indv}.txt ${indv}.g.vcf.gz > ${indv}.rehead.vcf.gz; tabix -p vcf ${indv}.rehead.vcf.gz; rm temp_${indv}.txt; done

sbatch -p general -t 1-00:00:00 --mem=20G