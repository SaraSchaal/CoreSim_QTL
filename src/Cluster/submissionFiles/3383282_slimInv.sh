#!/bin/bash
#SBATCH --job-name=3383282_slimInv.txt
#SBATCH --mem=2Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --output=outFiles/3383282.%j.out
#SBATCH --error=outFiles/3383282.%j.err
module load lotterhos/2019-11-15
srun slim -d seed=3383282 -d muBase=1e-07 -d n=1000 -d sigmaK=0.75 -d alpha=0.002 -d rep=1 -d muInv=0 -d enVar=0 -d mig1=0.001 -d mig2=0.001 -d burnin=10000 -d chromNum=21 -d r=1e-06 -d theta1=1 -d theta2=-1 -d dom=2 -m -t 20201013_InversionModel.slim
