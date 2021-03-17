#### Process a single script for simulations ####

### Load Packages and Download Data Files ###
## List Packages Needed 
packages_needed <- c("IntegratedMRF", "vcfR", "distances","ggplot2", "metR", "fields",
                     "MultivariateRandomForest", "gridExtra", "akima",
                     "MLmetrics", "ash", "plotly", "stringr", "tidyverse",
                     "bigsnpr", "bigstatsr", "ggpubr", "purrr", "dplyr")

## install packages that aren't installed already
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

## load each library
for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}

# download data
folder <- "results/Inversion/20210220_inOutInvFST/"
#folder <- "results/20210129_Movie/"
seed <- "3384725"
#df.params <- read.table("src/InvSimParams.txt", header = TRUE)
df.invTime <- read.table(paste0(folder, seed, "_outputInvTime.txt", sep = ""), header = TRUE)
df.invData <- read.table(paste0(folder, seed, "_outputInvSumInfo.txt", sep = ""), header = TRUE)
df.muts <- read.table(paste0(folder, seed, "_outputMutations.txt", sep = ""), header = TRUE)
df.popDyn <- read.table(paste0(folder, seed, "_outputPopDynam.txt", sep = ""), header = TRUE)
df.indPheno <- read.table(paste0(folder, seed, "_outputIndPheno.txt", sep = ""), header = TRUE)
df.invQTNData <- read.table(paste0(folder, seed, "_outputInvQtnSumInfo.txt", sep = ""), header = TRUE)
df.invQTNTime <- read.table(paste0(folder, seed, "_outputInvQtn.txt", sep = ""), header = TRUE)
df.params <- read.table(paste0(folder, seed, "_outputSimStats.txt", sep = ""), header = FALSE)
colnames(df.params) <- c("seed", "mig1", "mig2", "N1", "N2", "r", "muInv", 
                         "muBase", "alpha", "sigmaK", "burnin", "rep", "enVar")

df.invTime.NS <- read.table(paste0(folder, seed, "noSel_outputInvTime.txt", sep = ""), header = TRUE)
df.invData.NS <- read.table(paste0(folder, seed, "noSel_outputInvSumInfo.txt", sep = ""), header = TRUE)
df.muts.NS <- read.table(paste0(folder, seed, "noSel_outputMutations.txt", sep = ""), header = TRUE)
df.popDyn.NS <- read.table(paste0(folder, seed, "noSel_outputPopDynam.txt", sep = ""), header = TRUE)


### Manipulate data frames for merging
### SELECTION ###
df.invAllDatatemp <- merge(df.invData, df.invTime, by.x = c("inv_id"), 
                           by.y = c("inv_id"), all.y = TRUE)
colnames(df.invAllDatatemp)[9] <- "sim_gen"
df.invAllData <- df.invAllDatatemp[,-2]
df.invAllData$inv_age <- df.invAllData$sim_gen - df.invAllData$inv_originGen

### NO SELECTION ###
df.invAllDatatemp.NS <- merge(df.invData.NS, df.invTime.NS, by.x = c("inv_id"), 
                              by.y = c("inv_id"), all.y = TRUE)
colnames(df.invAllDatatemp.NS)[9] <- "sim_gen"
df.invAllData.NS <- df.invAllDatatemp.NS[,-2]
df.invAllData.NS$inv_age <- df.invAllData.NS$sim_gen - df.invAllData.NS$inv_originGen


## COPY AND PASTE WHERE NEEDED
pdf(paste0("figures/", seed, "XXXX.pdf"), height = 15, width = 15)

dev.off()

png(paste0("figures/", seed, "XXXX.png"), width = 480, height = 480, units = "px")

dev.off()