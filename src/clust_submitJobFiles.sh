#!/bin/bash
#SBATCH --job-name=pracSubmitJob			      
#SBATCH --mem=100Mb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1                        
#SBATCH --output=slimInvTest.%j.out                
#SBATCH --error=slimInvTest.%j.err                

for file in submissionFiles/*.sh 
do 
	echo $file
	sbatch $file
done

