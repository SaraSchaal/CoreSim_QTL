#######################################
#### Submit jobs for Inversion Sim ####
#######################################

params <- read.table("src/rerun3_invSimParams.txt", header = TRUE)
params <- as.data.frame(df.params )
for(i in 1:nrow(params)){
  
  filename <- paste(params$seed[i], "_slimInv.sh",sep="")
  fileConn <- file(print(paste("src/submissionFiles/20210930_fixOrigin/R_script/", params$seed[i],"_slimInv.sh",sep="")))
  
  writeLines(c("#!/bin/bash",
               paste0("#SBATCH --job-name=", params$seed[i],"_slimInv.txt",sep=""),
               "#SBATCH --mem=1Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=lotterhos",
               "#SBATCH --time=1:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --tasks-per-node=1",
               paste0("#SBATCH --output=../figures/20210930/clustOut_R/",params$seed[i],".%j.out"),
               paste0("#SBATCH --error=../figures/20210930/clustOut_R/",params$seed[i],".%j.err"),
               "module load lotterhos/2019-11-15",
        
               paste0("Rscript --vanilla ProcessSingleSimFile_20210803.R \"../results/20210930/\" \"../figures/20210930/\" ", params$seed[i])
               
  ), fileConn)
  
}
