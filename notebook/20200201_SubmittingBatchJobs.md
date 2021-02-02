# Description of files to run

To run a loop to submit a batch of simulation to the cluster you will need two files:
invSim_jobSubmit.R
invSimParams.txt
submitBatchScript.sh
clust_submitJobFiles.sh

The first file, invSim_jobSubmit.R, creates submission files and requires you have a folder within your current working directory called "submissionFiles". 
This script will write a new submission file for every row of parameters in the other required file which is invSimParams.txt. This txt file needs to be in
the current working directory. To create the files for submission to the cluster you then run the submitBatchScript.sh file using the following line of code:
```sbatch submitBatchScript.sh```
This creates the first 1000 job submissions to your submissionFiles folder and names them with the seed value for that simulation. 
The cluster can only handle 1000 jobs at a time. Now you should have 1000 jobs in the submission folder and you can run the next bash script:
```clust_submitJobFiles.sh```
This then submits all 1000 jobs to the cluster. If you have more than 1000 jobs then once these 1000 complete, you either delete or move the submission files and
change the for loop in your invSim_jobSubmit.R file to be for(i in 1001:2000) or less than 2000 if you have fewer than 2000 jobs. Keep doing this until you have run 
through all rows in your parameters file.
 
