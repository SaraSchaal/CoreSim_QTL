#######################################
#### Submit jobs for Inversion Sim ####
#######################################

params =read.table("./invSimParams.txt", header = TRUE)

for(i in 1:1000){

  filename <- paste(params$seed[i], "_slimInv.sh",sep="")
  fileConn<-file(print(paste("./submissionFiles/", params$seed[i],"_slimInv.sh",sep="")))
  
  # if(params$mu_base[i]>0.000001){
  #  time <- "#SBATCH --time=72:00:00"
  # } else {
  #  time <- "#SBATCH --time=120:00:00"
  # }
  
  writeLines(c("#!/bin/bash",
               paste0("#SBATCH --job-name=", params$seed[i],"_slimInv.txt",sep=""),
               "#SBATCH --mem=2Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=lotterhos",
               "#SBATCH --time=48:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --tasks-per-node=1",
               paste0("#SBATCH --output=outFiles/",params$seed[i],".%j.out"),
               paste0("#SBATCH --error=outFiles/",params$seed[i],".%j.err"),
               "module load lotterhos/2019-11-15",
               paste0("srun slim -d seed=", params$seed[i], " -d muBase=", params$muBase[i], " -d n=",
                      params$n[i], " -d sigmaK=", params$sigmaK[i], " -d alpha=", params$alpha[i], " -d rep=",
                      params$rep[i], " -d muInv=", params$muInv[i]," -d enVar=", params$enVar[i], " -d mig1=",
                      params$mig1[i]," -d mig2=", params$mig2[i]," -d burnin=", params$burnin[i], " -d chromNum=",
                      params$chromNum[i]," -d r=", params$r[i], " -d theta1=", params$theta1[i], " -d theta2=", 
                      params$theta2[i], " -d dom=", params$dom[i], " -m -t 20201013_InversionModel.slim")

  ), fileConn)
  
}
