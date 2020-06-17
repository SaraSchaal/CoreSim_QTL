#!/bin/bash
#SBATCH --job-name=stacksTest                     # sets the job name
#SBATCH -n 1                                 	  # reserves 1 machine
#SBATCH -N 1
#SBATCH --mem=100Gb                               # reserves 100 GB memory
#SBATCH --partition=lotterhos                     # requests that the job is executed in partition my partition
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --time=3:00:00                            # reserves machines/cores for 4 hours.
#SBATCH --output=stacksTest.%j.out                # sets the standard output to be stored in file my_nice_job.%j.out, where %j is the job id)
#SBATCH --error=stacksTest.%j.err                 # sets the standard error to be stored in file my_nice_job.%j.err, where %j is the job id)

module load lotterhos/2019-11-15

srun slim 