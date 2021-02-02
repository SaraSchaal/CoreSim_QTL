#!/bin/bash
#SBATCH --job-name=simRun3384725 	        	  # sets the job name
#SBATCH -n 1                                 	  # reserves 1 machine
#SBATCH -N 1
#SBATCH --mem=10Gb                                # reserves 50 GB memory
#SBATCH --partition=lotterhos                     # requests that the job is executed in partition my partition
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --time=100:00:00                            # reserves machines/cores for 72 hours.
#SBATCH --output=slimInvMovie.%j.out                # sets the standard output to be stored in file my_nice_job.%j.out, where %j is the job id)
#SBATCH --error=slimInvMovie.%j.err                 # sets the standard error to be stored in file my_nice_job.%j.err, where %j is the job id)

module load lotterhos/2020-08-24

srun slim -d chromNum=21 -d r=1e-06 -d theta1=1 -d theta2=-1 -d burnin=10000 -d muBase=1e-06 -d n=1000 -d mig1=0.25 -d mig2=0.25 -d sigmaK=0.75 -d alpha=0.002 -d seed=3384725 -d enVar=0 -d muInv=1e-3 -d dom=2 -m -t 20210106_InversionModel.slim 
