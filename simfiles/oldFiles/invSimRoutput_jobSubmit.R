#######################################
#### Submit jobs for Inversion Sim ####
#######################################

#params <- read.table("src/rerun3_invSimParams.txt", header = TRUE)
params <- removeLowInvMu_df
for(i in 1:nrow(params)){
  
  filename <- paste(params$seed[i], "_slimInv.sh",sep="")
  fileConn <- file(print(paste("src/submissionFiles/20220220_finalrerun/R/", params$seed[i],"_slimInv.sh",sep="")))
  if(i <= 500){
    partition <- "short"
  } else if(i > 500 & i <= 700){
    partition <- "long"
  } else {
    partition <- "lotterhos"
  }
  writeLines(c("#!/bin/bash",
               paste0("#SBATCH --job-name=", params$seed[i],"_slimInv.txt",sep=""),
               "#SBATCH --mem=1Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               paste0("#SBATCH --partition=", partition),
               "#SBATCH --time=2:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --tasks-per-node=1",
               paste0("#SBATCH --output=/scratch/schaal.s/InversionSimulations/figures/20220220/clustOut_R/",params$seed[i],"NewRun.%j.out"),
               paste0("#SBATCH --error=/scratch/schaal.s/InversionSimulations/figures/20220220/clustOut_R/",params$seed[i],"NewRun.%j.err"),
               "module load lotterhos/2020-08-24",
        
               paste0("Rscript --vanilla ProcessSingleSimFile_20220205.R \"/scratch/schaal.s/InversionSimulations/results/20220220/10gen/\" \"/scratch/schaal.s/InversionSimulations/figures/20220220/\" ", params$seed[i])
               
  ), fileConn)
  
}
