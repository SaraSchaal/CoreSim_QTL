#!/bin/bash
#SBATCH --job-name=3383979_slimInvGenome
#SBATCH --mem=10Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --output=../figures/test/clustOut/invGenome/3383979.%j.out
#SBATCH --error=../figures/test/clustOut/invGenome/3383979.%j.err
module load lotterhos/2019-11-15
Rscript --vanilla invGenomeCalc.R "../results/test/" "../figures/test/" 3383979