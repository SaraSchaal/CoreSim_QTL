#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --job-name=NA_slimInv.txt
#SBATCH --mem=100Mb
#SBATCH --partition=general
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=NA.output
#SBATCH --error=NA.error
module load lotterhos lotterhos/2019-11-15
srun slim -d seed=NA -d muBase=NA -d n=NA -d sigmaK=NA -d alpha=NA -d rep=NA -d muInv=NA -d enVar=NA -d mig1=NA -d mig2=NA -d burnin=NA -d chromNum=NA -d r=NA -d theta1=NA -d theta2=NA -d dom=NA -m -t InversionModel.slim
