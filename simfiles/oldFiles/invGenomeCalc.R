#### Process a single script for simulations
### Sara M. Schaal

######################################################################################################
### Load Packages and Download Data Files ###
## List Packages Needed 
packages_needed <- c("scales","IntegratedMRF", "vcfR", "distances","ggplot2", "metR", "RColorBrewer",
                     "MultivariateRandomForest", "gridExtra", "akima", "fields", "ggnewscale",
                     "MLmetrics", "ash", "plotly", "stringr", "tidyverse", "viridisLite", "remotes",
                     "bigsnpr", "bigstatsr", "ggpubr", "purrr", "dplyr", "lfmm", "pcadapt", "BiocManager" )
bioc_packages <- c("LEA", "qvalue")


## Install packages that aren't installed already and load them
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library( packages_needed[i], character.only = TRUE)
}

for(i in 1:length(bioc_packages)){
  if(!(bioc_packages[i] %in% installed.packages())){BiocManager::install(bioc_packages[i])}
  library(bioc_packages[i], character.only = TRUE)
}

if(!("OutFLANK" %in% installed.packages())){remotes::install_github("whitlock/OutFLANK")}
library(OutFLANK)

### Download Data
args = commandArgs(trailingOnly=TRUE)
folderIn <- args[1] # "results/Inversion/20210525_fullData/"
folderOut <- args[2] #"figures/20210523_fullData/" 
seed <- args[3] #"3383648"  

df.invTime <- read.table(paste0(folderIn, seed, "_outputInvTime.txt", sep = ""), header = TRUE)
df.invData <- read.table(paste0(folderIn, seed, "_outputInvSumInfo.txt", sep = ""), header = TRUE)
df.muts <- read.table(paste0(folderIn, seed, "_outputMutations.txt", sep = ""), header = TRUE)
df.popDyn <- read.table(paste0(folderIn, seed, "_outputPopDynam.txt", sep = ""), header = TRUE)
df.indPheno <- read.table(paste0(folderIn, seed, "_outputIndPheno.txt", sep = ""), header = TRUE)
df.invQTNData <- read.table(paste0(folderIn, seed, "_outputInvQtnSumInfo.txt", sep = ""), header = TRUE)
df.invQTNTime <- read.table(paste0(folderIn, seed, "_outputInvQtn.txt", sep = ""), header = TRUE)
df.params <- read.table(paste0(folderIn, seed, "_outputSimStats.txt", sep = ""), header = FALSE)
colnames(df.params) <- c("seed", "mig1", "mig2", "N1", "N2", "muBase", "muInv", 
                         "r", "alpha", "sigmaK", "burnin", "rep", "enVar")

df.invTime.NS <- read.table(paste0(folderIn, seed, "noSel_outputInvTime.txt", sep = ""), header = TRUE)
df.invData.NS <- read.table(paste0(folderIn, seed, "noSel_outputInvSumInfo.txt", sep = ""), header = TRUE)
df.muts.NS <- read.table(paste0(folderIn, seed, "noSel_outputMutations.txt", sep = ""), header = TRUE)
df.popDyn.NS <- read.table(paste0(folderIn, seed, "noSel_outputPopDynam.txt", sep = ""), header = TRUE)
df.indPheno.NS <- read.table(paste0(folderIn, seed, "noSel_outputIndPheno.txt", sep = ""), header = TRUE)

## Functions 
# This function saves the legend into an object to use for plotting as an element in ggarrange
g_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
######################################################################################################

# If there are no inversions in the simulation we do not want to evaluate the majority of the following
# code. The only part we need is the amount of local adaptation that the simulation reached. 
if(df.params$muInv != 0){
  ######################################################################################################
  #### Manipulate data frames for merging ####
  ### SELECTION ###
  df.invAllDatatemp <- merge(df.invData, df.invTime, by.x = c("inv_id"), 
                             by.y = c("inv_id"), all.y = TRUE)
  colnames(df.invAllDatatemp)[9] <- "sim_gen"
  df.invAllData <- df.invAllDatatemp[,-2]
  df.invAllData$inv_age <- df.invAllData$sim_gen - df.invAllData$inv_originGen
  df.invAllData$num_qtns_Lscaled <- df.invAllData$num_qtns/df.invAllData$inv_length
  
  
  ### NO SELECTION ###
  df.invAllDatatemp.NS <- merge(df.invData.NS, df.invTime.NS, by.x = c("inv_id"), 
                                by.y = c("inv_id"), all.y = TRUE)
  colnames(df.invAllDatatemp.NS)[9] <- "sim_gen"
  df.invAllData.NS <- df.invAllDatatemp.NS[,-2]
  df.invAllData.NS$inv_age <- df.invAllData.NS$sim_gen - df.invAllData.NS$inv_originGen
  df.invAllData.NS$num_qtns_Lscaled <- df.invAllData.NS$num_qtns/df.invAllData.NS$inv_length
  
  ######################################################################################################
  
  
  
  ######################################################################################################
  #### Calculate Inversion Windows ####
  
  # filter the inversion data to get the final generation
  df.invDataFinalGen <- filter(df.invAllData, sim_gen == 50000)
  df.invDataFinalGen.NS <- filter(df.invAllData.NS, sim_gen == 50000)
  
  ## CHECK TO SEE IF THERE ARE ANY INVERSIONS IN THE FINAL GENERATION
  if(nrow(df.invDataFinalGen) > 0 & nrow(df.invDataFinalGen.NS) > 0){
    
    # identify all the positions in the genome that are within the inversion windows
    invWind <- function(x){
      invWindowBases <- NULL
      numSegInv <- 0
      inversion.ID <- NULL
      for(i in 1:nrow(x)) {
        bases <- seq(from = x$inv_pos[i],to = x$inv_end[i], by = 1)
        invWindowBases <- c(invWindowBases, bases)
        numSegInv <- numSegInv + 1
        inversion.ID <- c(inversion.ID, rep(x$inv_id[i], length(bases)))
      }
      
      # identify all the positions in the genome that are within the inversion windows
      # filtered for MAF of 0.01
      invWindowBasesMAF <- NULL
      numSegInvMAF <- 0
      inversion.ID.MAF <- NULL
      for(i in 1:nrow(x)) {
        Inv_wtype_freq <-  1 - x$freq[i]					
        MAF <-  min(x$freq[i], Inv_wtype_freq)				
        if(MAF > 0.01){
          basesMAF <- seq(from = x$inv_pos[i], to = x$inv_end[i], by = 1)
          invWindowBasesMAF <- c(invWindowBasesMAF, basesMAF)
          numSegInvMAF <- numSegInvMAF + 1
          inversion.ID.MAF <- c(inversion.ID.MAF, rep(x$inv_id[i], length(basesMAF)))
        }
      }
      
      df.invWindow <- data.frame("invWindBases" = invWindowBases, "invID" = inversion.ID)
      df.invWindowMAF <- data.frame("invWindBasesMAF" = invWindowBasesMAF, "invIDMAF" = inversion.ID.MAF)
      obj_list <- list("df.invWind" = df.invWindow,
                       "df.invWindMAF" = df.invWindowMAF ,
                       "numSeqInv" = numSegInv, "numSeqInvMAF" = numSegInvMAF)
      return(obj_list)
    }
    
    df.invWind <- invWind(df.invDataFinalGen)
    invWindBases.MAF <- invWind(df.invDataFinalGen)[[2]]
    df.invWindNS <- invWind(df.invDataFinalGen.NS) 
    invWindBasesNS.MAF <- invWind(df.invDataFinalGen.NS)[[2]]
    
    uniqueBases <- length(unique(invWindBases.MAF$invWindBasesMAF))
    num.overlap <- length(invWindBases.MAF$invWindBasesMAF) - uniqueBases
    percGenome <- (uniqueBases/2000000)*100
    df.invGenome <- c(uniqueBases, num.overlap, round(percGenome, digits = 1))
    
    write.table(t(c(seed, df.invGenome)), paste0(folderOut, "outputInvGenome.txt"), col.names = F, row.names = F, append = T)
  }
    
} else {
  write.table(t(c(seed, 0,0,0)), paste0(folderOut, "outputInvGenome.txt"), col.names = F, row.names = F, append = T)
  
}
    
### end output inverted genome size
######################################################################################################  

    