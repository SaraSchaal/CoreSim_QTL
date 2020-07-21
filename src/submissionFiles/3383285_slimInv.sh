#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --job-name=3383285_slimInv.txt
#SBATCH --mem=100Mb
#SBATCH --partition=general
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=outFiles/3383285.out
#SBATCH --error=outFiles3383285.err
module load lotterhos/2019-11-15
srun slim -d seed=3383285 -d muBase=0.0000001 -d n=1000 -d sigmaK=0.9 -d alpha=0.002 -d rep=1 -d muInv=0 -d enVar=0 -d mig1=0.5 -d mig2=0.5 -d burnin=10000 -d chromNum=21 -d r=0.000001 -d theta1=1 -d theta2=-1 -d dom=FALSE -m -t InversionModel.slim
