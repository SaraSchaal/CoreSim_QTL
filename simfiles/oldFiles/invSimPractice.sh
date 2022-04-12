#!/bin/bash
#SBATCH --job-name=slimInvTest                    # sets the job name
#SBATCH -n 1                                 	  # reserves 1 machine
#SBATCH -N 1
#SBATCH --mem=50Gb                                # reserves 50 GB memory
#SBATCH --partition=lotterhos                     # requests that the job is executed in partition my partition
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --time=24:00:00                            # reserves machines/cores for 24 hours.
#SBATCH --output=slimInvTest.%j.out                # sets the standard output to be stored in file my_nice_job.%j.out, where %j is the job id)
#SBATCH --error=slimInvTest.%j.err                 # sets the standard error to be stored in file my_nice_job.%j.err, where %j is the job id)

module load lotterhos/2019-11-15

srun slim -d muBase=1e-03 -d n=1000 -d mig1=0.01 -d mig2=0.01 -d sigmaK=0.9 -d alpha=0.002 -d seed=3383283 -d enVar=0.01 -d muInv=1e-3 -d dom=F -m -t InversionModel.slim 