#######################################
#### Submit jobs for Inversion Sim ####
#######################################

#params <- read.table("src/rerun3_invSimParams.txt", header = TRUE)
params <- df.params.3

for(i in 1:nrow(params)){
  
  filename <- paste(params$seed[i], "_slimInv.sh",sep="")
  fileConn <- file(print(paste("src/submissionFiles/20210621_additionalParams/rerun/", params$seed[i],"_slimInv.sh",sep="")))
  
  writeLines(c("#!/bin/bash",
               paste0("#SBATCH --job-name=", params$seed[i],"_slimInv.txt",sep=""),
               "#SBATCH --mem=10Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=short",
               "#SBATCH --time=24:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --tasks-per-node=1",
               paste0("#SBATCH --output=../figures/test/clustOut/",params$seed[i],".%j.out"),
               paste0("#SBATCH --error=../figures/test/clustOut/",params$seed[i],".%j.err"),
               "module load lotterhos/2019-11-15",
               paste0("srun slim -d seed=", params$seed[i], " -d muBase=", params$muBase[i], " -d n=",
                      params$n[i], " -d sigmaK=", params$sigmaK[i], " -d alpha=", params$alpha[i], " -d rep=",
                      params$rep[i], " -d muInv=", params$muInv[i]," -d enVar=", params$enVar[i], " -d mig1=",
                      params$mig1[i]," -d mig2=", params$mig2[i]," -d burnin=", params$burnin[i], " -d chromNum=",
                      params$chromNum[i]," -d r=", params$r[i], " -d theta1=", params$theta1[i], " -d theta2=", 
                      params$theta2[i], " -d dom=", params$dom[i], " -m -t 20210429_InversionModel.slim"),
               
               paste0("echo ", "\"finished selection sim\""),
               
               paste0("srun slim -d seed=", params$seed[i], " -d muBase=", params$muBase[i], " -d n=",
                      params$n[i], " -d sigmaK=", params$sigmaK[i], " -d alpha=", params$alpha[i], " -d rep=",
                      params$rep[i], " -d muInv=", params$muInv[i]," -d enVar=", params$enVar[i], " -d mig1=",
                      params$mig1[i]," -d mig2=", params$mig2[i]," -d burnin=", params$burnin[i], " -d chromNum=",
                      params$chromNum[i]," -d r=", params$r[i], " -d theta1=", params$theta1[i], " -d theta2=", 
                      params$theta2[i], " -d dom=", params$dom[i], " -m -t 20210429_InversionModelnoSel.slim"),
               "wait",
               "echo \"finished no selection sim\"",
               
               paste0("vcftools --vcf ../results/test/", params$seed[i], "_InversionVCF.vcf --maf 0.01 --out ../results/test/",
                      params$seed[i], "_InversionVCF_MAF01 --recode --recode-INFO-all"),
               
               paste0("vcftools --vcf ../results/test/", params$seed[i], "noSel_InversionVCF.vcf --maf 0.01 --out ../results/test/",
                      params$seed[i], "noSel_InversionVCF_MAF01 --recode --recode-INFO-all"),
               
               paste0("Rscript --vanilla ProcessSingleSimFile_20210706.R \"../results/test/\" \"../figures/test/\" ", params$seed[i])

  ), fileConn)
  
}
