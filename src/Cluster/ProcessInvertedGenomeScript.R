#######################################
#### Submit jobs for Inversion Sim ####
#######################################

params <- read.table("src/invSimParams.txt", header = TRUE)
params <- df.params.2
for(i in 1:nrow(params)){
  
  filename <- paste(params$seed[i], "_slimInvGenome.sh",sep="")
  fileConn <- file(print(paste("src/submissionFiles/20210614_invGenome/Already_run/", params$seed[i],"_slimInvGenome.sh",sep="")))
  
  writeLines(c("#!/bin/bash",
               paste0("#SBATCH --job-name=", params$seed[i],"_slimInvGenome",sep=""),
               "#SBATCH --mem=10Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=short",
               "#SBATCH --time=4:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --tasks-per-node=1",
               paste0("#SBATCH --output=../figures/test/clustOut/invGenome/",params$seed[i],".%j.out"),
               paste0("#SBATCH --error=../figures/test/clustOut/invGenome/",params$seed[i],".%j.err"),
               "module load lotterhos/2019-11-15",
                
               paste0("Rscript --vanilla invGenomeCalc.R \"../results/test/\" \"../figures/test/\" ", params$seed[i])

  ), fileConn)
  
}


