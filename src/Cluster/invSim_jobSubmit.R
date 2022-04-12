#######################################
#### Submit jobs for Inversion Sim ####
#######################################

#params <- read.table("src/rerun3_invSimParams.txt", header = TRUE)
params <- removeLowInvMu_df[removeLowInvMu_df$seed >= 3384722,] # NEED MUTS
params <- removeLowInvMu_df[removeLowInvMu_df$seed >= 3384779 & removeLowInvMu_df$seed <= 3384896,] # muts DONE
params <- removeLowInvMu_df[removeLowInvMu_df$seed >= 3384722 & removeLowInvMu_df$seed <= 3384778,] # missed NOT RUN YET
params <- removeLowInvMu_df[removeLowInvMu_df$seed >= 3385005,] # running on lotterhos
params <- removeLowInvMu_df[removeLowInvMu_df$seed <= 3385004 & removeLowInvMu_df$seed >= 3384897 & removeLowInvMu_df$seed != 3384909 &
                            removeLowInvMu_df$seed != 3384905,]
params <- removeLowInvMu_df[removeLowInvMu_df$seed >= 3383900 & removeLowInvMu_df$seed <= 3384721,]
params <- removeLowInvMu_df[removeLowInvMu_df$seed <= 3384721 & removeLowInvMu_df$seed >=3384611,]
params <- removeLowInvMu_df[removeLowInvMu_df$seed == 3384905,]
for(i in 1:nrow(params)){
  
  filename <- paste(params$seed[i], "_slimInv.sh",sep="")
  fileConn <- file(print(paste("src/submissionFiles/20220220_finalrerun/noSel/", params$seed[i],"_noSelslimInv.sh",sep="")))
  if(i <= 25){
    partition <- "lotterhos"
  } else {
    partition <- "lotterhos"
  }
  writeLines(c("#!/bin/bash",
               paste0("#SBATCH --job-name=", params$seed[i],"_slimInv",sep=""),
               "#SBATCH --mem=5Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               paste0("#SBATCH --partition=", partition),
               "#SBATCH --time=48:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --tasks-per-node=1",
               paste0("#SBATCH --output=/scratch/schaal.s/InversionSimulations/figures/20220220/clustOut/",params$seed[i],".%j.out"),
               paste0("#SBATCH --error=/scratch/schaal.s/InversionSimulations/figures/20220220/clustOut/",params$seed[i],".%j.err"),
               "module load lotterhos/2020-08-24",
               # paste0("srun slim -d seed=", params$seed[i], " -d muProp=", params$muProp[i], " -d n=",
               #         params$n[i], " -d sigmaK=", params$sigmaK[i], " -d alpha=", params$alpha[i], " -d rep=",
               #         params$rep[i], " -d muInv=", params$muInv[i]," -d enVar=", params$enVar[i], " -d mig1=",
               #         params$mig1[i]," -d mig2=", params$mig2[i]," -d burnin=", params$burnin[i], " -d chromNum=",
               #         params$chromNum[i]," -d r=", params$r[i], " -d theta1=", params$theta1[i], " -d theta2=",
               #         params$theta2[i], " -d dom=", params$dom[i], " -m -t 20220324_InversionModel10gen.slim"),
               # 
               # paste0("echo ", "\"finished selection sim\""),

               paste0("srun slim -d seed=", params$seed[i], " -d muProp=", params$muProp[i], " -d n=",
                      params$n[i], " -d sigmaK=", params$sigmaK[i], " -d alpha=", params$alpha[i], " -d rep=",
                      params$rep[i], " -d muInv=", params$muInv[i]," -d enVar=", params$enVar[i], " -d mig1=",
                      params$mig1[i]," -d mig2=", params$mig2[i]," -d burnin=", params$burnin[i], " -d chromNum=",
                      params$chromNum[i]," -d r=", params$r[i], " -d theta1=", params$theta1[i], " -d theta2=",
                      params$theta2[i], " -d dom=", params$dom[i], " -m -t 20220215_InversionModelnoSel10gen.slim"),
               "wait",
               "echo \"finished no selection sim\"",
               "module unload lotterhos/2020-08-24",
               "module load lotterhos/2019-11-15",
               # paste0("vcftools --vcf /scratch/schaal.s/InversionSimulations/results/20220220/10gen/", params$seed[i], "_InversionVCF.vcf --maf 0.01 --out /scratch/schaal.s/InversionSimulations/results/20220220/10gen/",
               #         params$seed[i], "_InversionVCF_MAF01 --recode --recode-INFO-all"),
               
               paste0("vcftools --vcf /scratch/schaal.s/InversionSimulations/results/20220220/10gen/", params$seed[i], "noSel_InversionVCF.vcf --maf 0.01 --out /scratch/schaal.s/InversionSimulations/results/20220220/10gen/",
                       params$seed[i], "noSel_InversionVCF_MAF01 --recode --recode-INFO-all")

               #paste0("Rscript --vanilla ProcessSingleSimFile_20210803.R \"../results/20210930/\" \"../figures/20210930/\" ", params$seed[i])

  ), fileConn)
  
}
