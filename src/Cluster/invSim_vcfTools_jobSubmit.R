#######################################
#### Submit jobs for Inversion Sim ####
#######################################

#params <- read.table("src/rerun3_invSimParams.txt", header = TRUE)
params <- 3384725
params <- as.data.frame(params)
colnames(params) <- "seed"
for(i in 1:nrow(params)){
  
  filename <- paste(params$seed[i], "_slimInv.sh",sep="")
  fileConn <- file(print(paste("src/submissionFiles/20211203_errorInOrigin/", params$seed[i],"_slimInv.sh",sep="")))
  
  writeLines(c("#!/bin/bash",
               paste0("#SBATCH --job-name=", params$seed[i],"_slimInv.txt",sep=""),
               "#SBATCH --mem=10Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=lotterhos",
               "#SBATCH --time=6:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --tasks-per-node=1",
               paste0("#SBATCH --output=../figures/20211203/clustOut/",params$seed[i],".%j.out"),
               paste0("#SBATCH --error=../figures/20211203/clustOut/",params$seed[i],".%j.err"),
               "module load lotterhos/2019-11-15",
               
               paste0("vcftools --vcf ../results/20211203/", params$seed[i], "_InversionVCF.vcf --maf 0.01 --out ../results/20211203/",
                      params$seed[i], "_InversionVCF_MAF01 --recode --recode-INFO-all"),
               
               paste0("vcftools --vcf ../results/20211203/", params$seed[i], "noSel_InversionVCF.vcf --maf 0.01 --out ../results/20211203/",
                      params$seed[i], "noSel_InversionVCF_MAF01 --recode --recode-INFO-all")
               
               #paste0("Rscript --vanilla ProcessSingleSimFile_20210803.R \"../results/20210930/\" \"../figures/20210930/\" ", params$seed[i])

  ), fileConn)
  
}
