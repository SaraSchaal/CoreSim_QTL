#!/bin/bash
#SBATCH --job-name=3384695_slimInv.txt
#SBATCH --mem=1Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --output=../figures/test/clustOut/20210719/3384695.%j.out
#SBATCH --error=../figures/test/clustOut/20210719/3384695.%j.err
module load lotterhos/2019-11-15
Rscript --vanilla ProcessSingleSimFile_20210706.R "../results/test/" "../figures/test/" 3384695