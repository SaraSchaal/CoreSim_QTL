#######################################
#### Submit jobs for Inversion Sim ####
#######################################

params <- read.table("src/rerun3_invSimParams.txt", header = TRUE)
params <- df.params.2
for(i in 1:nrow(params)){
  
  filename <- paste(params$seed[i], "_slimInv.sh",sep="")
  fileConn <- file(print(paste("src/submissionFiles/20210621_additionalParams/", params$seed[i],"_slimInv.sh",sep="")))
  
  writeLines(c("#!/bin/bash",
               paste0("#SBATCH --job-name=", params$seed[i],"_slimInv.txt",sep=""),
               "#SBATCH --mem=1Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=short",
               "#SBATCH --time=1:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --tasks-per-node=1",
               paste0("#SBATCH --output=../figures/test/clustOut/20210719/",params$seed[i],".%j.out"),
               paste0("#SBATCH --error=../figures/test/clustOut/20210719/",params$seed[i],".%j.err"),
               "module load lotterhos/2019-11-15",
        
               paste0("Rscript --vanilla ProcessSingleSimFile_20210706.R \"../results/test/\" \"../figures/test/\" ", params$seed[i])
               
  ), fileConn)
  
}
