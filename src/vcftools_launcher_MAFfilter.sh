#!/bin/bash
#SBATCH -p short
#SBATCH --nodes 1
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --mem=10G
#SBATCH -o vcftools.%N.%j.out
#SBATCH -e vcftools.%N.%j.err
#SBATCH --job-name="vcftools_MAFsubset"
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=schaal.s@northeastern.edu

vcftools --vcf 3383643_InversionVCF.vcf --maf 0.01 --out 3383643_InversionVCF_MAF01 --recode --recode-INFO-all