#!/bin/bash
#SBATCH --job-name=3384626_slimInv.txt
#SBATCH --mem=10Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --output=../figures/test/clustOut/3384626.%j.out
#SBATCH --error=../figures/test/clustOut/3384626.%j.err
module load lotterhos/2019-11-15
srun slim -d seed=3384626 -d muBase=1e-06 -d n=1000 -d sigmaK=3 -d alpha=0.002 -d rep=3 -d muInv=0 -d enVar=0 -d mig1=0.001 -d mig2=0.001 -d burnin=10000 -d chromNum=21 -d r=1e-06 -d theta1=1 -d theta2=-1 -d dom=2 -m -t 20210429_InversionModel.slim
echo "finished selection sim"
srun slim -d seed=3384626 -d muBase=1e-06 -d n=1000 -d sigmaK=3 -d alpha=0.002 -d rep=3 -d muInv=0 -d enVar=0 -d mig1=0.001 -d mig2=0.001 -d burnin=10000 -d chromNum=21 -d r=1e-06 -d theta1=1 -d theta2=-1 -d dom=2 -m -t 20210429_InversionModelnoSel.slim
wait
echo "finished no selection sim"
vcftools --vcf ../results/test/3384626_InversionVCF.vcf --maf 0.01 --out ../results/test/3384626_InversionVCF_MAF01 --recode --recode-INFO-all
vcftools --vcf ../results/test/3384626noSel_InversionVCF.vcf --maf 0.01 --out ../results/test/3384626noSel_InversionVCF_MAF01 --recode --recode-INFO-all
Rscript --vanilla ProcessSingleSimFile.R "../results/test/" "../figures/test/" 3384626
