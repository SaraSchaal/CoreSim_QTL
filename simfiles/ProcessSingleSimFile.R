#### Process a single script for simulations
### Sara M. Schaal

######################################################################################################
### Load Packages and Download Data Files ###
## List Packages Needed 
packages_needed <- c("IntegratedMRF", "vcfR", "distances","ggplot2", "metR", 
                     "MultivariateRandomForest", "gridExtra", "akima", "fields",
                     "MLmetrics", "ash", "plotly", "stringr", "tidyverse",
                     "bigsnpr", "bigstatsr", "ggpubr", "purrr", "dplyr", "lfmm", "pcadapt" )

## Install packages that aren't installed already
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

## Load each library
for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LEA")
library(LEA)
BiocManager::install("qvalue")
library(qvalue)
install.packages("remotes")
remotes::install_github("whitlock/OutFLANK")
library(OutFLANK)

### Download Data
folderIn <- "results/Inversion/20210321_runLowMig/" #args[1]
folderOut <- "figures/20210321_lowVhighMig/High_Mig/" #args[2]
seed <- "3384725" #args[3]

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

## Functions 
# This function saves the legend into an object to use for plotting as an element in ggarrange
g_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
######################################################################################################



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

### end inversion window calculation
######################################################################################################  



######################################################################################################  
#### Identify which QTNs fall within inversion windows ####

## filter for MAF 
# SELECTION
str(df.muts)
df.muts$FST <- as.numeric(as.character(df.muts$FST))
freq <- df.muts$freq
altFreq <- 1-freq
df.calc <- cbind(freq, altFreq)
df.muts$MAF <- apply(df.calc, 1, min)
df.muts.MAF <- subset(df.muts, subset = MAF >= 0.01)

# NO SELECTION
df.muts.NS$FST <- as.numeric(as.character(df.muts.NS$FST))
freq <- df.muts.NS$freq
altFreq <- 1-freq
df.calc <- cbind(freq,altFreq)
df.muts.NS$MAF <- apply(df.calc, 1, min)
df.muts.NS.MAF <- subset(df.muts.NS, subset = MAF >= 0.01)

# identify which qtns overlap with inv window locations
invWinQTNrows <- which(df.muts.MAF$position %in% invWindBases.MAF$invWindBasesMAF)
invWinQTNrows.NS <- which(df.muts.NS.MAF$position %in% invWindBasesNS.MAF$invWindBasesMAF)

# use this to identify how many positions have multiple qtns 
#mult.qtns <- length(which(qtnMuts$position %in% invWindBases)) - length(intersect(invWindBases, qtnMuts$position))
#mult.qtns.NS <- length(which(qtnMuts.NS$position %in% invWindBasesNS)) - length(intersect(invWindBasesNS, qtnMuts.NS$position))

## Selection
for(i in 1:nrow(df.muts.MAF)){
  if(df.muts.MAF$type[i] == "m2"){
    if(i %in% invWinQTNrows){
      df.muts.MAF$inOut[i] <- "in"
    } else {
      df.muts.MAF$inOut[i] <- "out"
    }
  } else if(df.muts.MAF$type[i] == "m1"){
    df.muts.MAF$inOut[i] <- "neut"
  } else {
    df.muts.MAF$inOut[i] <- "inv"
  }
}
df.muts.MAF$inOut <- as.factor(df.muts.MAF$inOut)

## No selection
for(i in 1:nrow(df.muts.NS.MAF)){
  if(df.muts.NS.MAF$type[i] == "m2"){
    if(i %in% invWinQTNrows.NS){
      df.muts.NS.MAF$inOut[i] <- "in"
    } else {
      df.muts.NS.MAF$inOut[i] <- "out"
    }
  } else if(df.muts.NS.MAF$type[i] == "m1"){
    df.muts.NS.MAF$inOut[i] <- "neut"
  } else {
    df.muts.NS.MAF$inOut[i] <- "inv"
  }
}
df.muts.NS.MAF$inOut <- as.factor(df.muts.NS.MAF$inOut)


### end identifying inversion window QTNs
######################################################################################################  




######################################################################################################  
### Identify adaptive inversions ####

# visualize the FST values inside inversions for each sim
par(mfrow= c(1,2), mar = c(2,2,3,2), oma = c(1.5,1.5,3,0))
df.invQTNs <- subset(df.muts.MAF, subset = df.muts.MAF$inOut == "in")
df.invQTNs$FST <- as.numeric(as.character(df.invQTNs$FST))
hist(df.invQTNs$FST, main = "Selection",
     xlab = "FST") 

df.invQTNs.NS <- subset(df.muts.NS.MAF, subset = df.muts.NS.MAF$inOut == "in")
df.invQTNs.NS$FST <- as.numeric(as.character(df.invQTNs.NS$FST))
hist(df.invQTNs.NS$FST, main = "No Selection",
     xlab = "FST") 
mtext(expression(bold("Inversion Window QTN FST values")), side = 3, outer = TRUE, cex = 1.5)
mtext(expression(bold("FST")), side = 1, outer = TRUE)
mtext(expression(bold("Frequency")), side = 2, outer = TRUE)

# create null distributions for outlier criteria
null <- df.muts.NS.MAF$FST[!is.nan(df.muts.NS.MAF$FST) & df.muts.NS.MAF$type == "m2"] # criteria 1: compare to a null distribution of all QTNs in no-selec sim
#null <- df.muts.NS.MAF$FST[df.muts.NS.MAF$inOut == "in"]
#null_crit1_2 <- df.muts.NS.MAF$FST[df.muts.NS.MAF$inOut == "in"]
null_neut <- df.muts.MAF$FST[df.muts.NS.MAF$type == "m1"] # criteria 2: compare to a null distribution of neutral QTNs in selec sim

# subset for just qtns
df.qtnMuts.MAF <- df.muts.MAF[df.muts.MAF$type == "m2" | df.muts.MAF$type == "m1",]
df.qtnMuts.NS.MAF <- df.muts.MAF[df.muts.MAF$type == "m2" | df.muts.MAF$type == "m1",]

# criteria calculation loop 
df.qtnMuts.MAF$crit1_p.value <- NULL
#df.qtnMuts.MAF$crit1_p.value.2 <- NULL
df.qtnMuts.MAF$crit2_p.value <- NULL
df.qtnMuts.MAF$crit3_Va <- NULL
for(i in 1:nrow(df.qtnMuts.MAF)){
  if(!is.nan(df.qtnMuts.MAF$FST[i])){
    obs <- df.qtnMuts.MAF$FST[i]
    df.qtnMuts.MAF$crit1_p.value[i] <- 1-rank(c(null, obs))[length(null)+1]/(length(null)+1)
    #df.qtnMuts.MAF$crit1_p.value.2[i] <- 1-rank(c(null_crit1_2, obs))[length(null_crit1_2)+1]/(length(null_crit1_2)+1)
    df.qtnMuts.MAF$crit2_p.value[i] <- 1-rank(c(null_neut, obs))[length(null_neut)+1]/(length(null_neut)+1)
    df.qtnMuts.MAF$crit3_Va[i] <- (df.qtnMuts.MAF$selCoef[i]^2)*df.qtnMuts.MAF$freq[i]*(1-df.qtnMuts.MAF$freq[i])
  } else {
    df.qtnMuts.MAF$crit1_p.value[i] <- NA
    #df.qtnMuts.MAF$crit1_p.value.2[i] <- NA
    df.qtnMuts.MAF$crit2_p.value[i] <- NA
    df.qtnMuts.MAF$crit3_Va[i] <- NA
  }
}

# 3rd criteria -- added genetic variation for each mutation and the proportion of added genetic variation
df.qtnMuts.MAF$crit3_Va_prop <- df.qtnMuts.MAF$crit3_Va/sum(df.qtnMuts.MAF$crit3_Va, na.rm = TRUE)
df.qtnMuts.MAF$crit3_Va_perc <- df.qtnMuts.MAF$crit3_Va_prop * 100
# side note: try to get some simulations to show drift play with parameters

# plot for criteria 1 
df.qtnMuts.MAF$inOut <- recode_factor(df.qtnMuts.MAF$inOut, 'in' = 'Inside Inversion', 'out' = 'Outside Inversion')
crit1.plot <- ggplot(df.qtnMuts.MAF[df.qtnMuts.MAF$type == "m2",], aes(x = crit1_p.value, fill = inOut)) +
                geom_histogram(position = "identity", alpha = 0.8) + 
                scale_fill_manual(values = c("red", "blue")) +
                labs(title =  "Criteria 1 - compare to null of all \nQTNs no-selection simulation",
                     x = "empirical p-value") +
                theme(legend.position = "none")
  

# plot for criteria 2
crit2.plot <- ggplot(df.qtnMuts.MAF[df.qtnMuts.MAF$type == "m2",], aes(x = crit2_p.value, fill = inOut)) +
                geom_histogram(position = "identity", alpha = 0.8) + 
                scale_fill_manual(values = c("red", "blue")) +
                labs(title =  "Criteria 2 - compare to null of neutral \nQTNs selection simulation",
                     x = "empirical p-value") +
                theme(legend.position = "none") 

# plot Va percent
crit3.plot.leg <- ggplot(df.qtnMuts.MAF[df.qtnMuts.MAF$type == "m2",], aes(x = crit3_Va_perc, fill = inOut)) +
                     geom_histogram(position = "identity", alpha = 0.8) + 
                     scale_fill_manual(name = "QTN Location", 
                                       values = c("red", "blue")) +
                    labs(title = "", x = "")  + 
                    theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank())
                     

# plot Va percent - subset for above 0.01
crit3.subplot <- ggplot(df.qtnMuts.MAF[df.qtnMuts.MAF$type == "m2" & df.qtnMuts.MAF$crit3_Va_perc > 0.01, ], 
                        aes(x = crit3_Va_perc, fill = inOut)) +
                  geom_histogram(position = "identity", alpha = 0.8) + 
                  scale_fill_manual(values = c("red", "blue")) +
                  labs(title = "Criteria 3 - compare % Va within \nselection simulation subset for > 0.01",
                       x = "percent of Va explained") +
                  theme(legend.position = "none")
  
  
 crit.leg <- g_legend(crit3.plot.leg)

  pdf(paste0(folderOut, seed, "_adaptInvCriteria.pdf"), height = 5, width = 15)

    ggarrange(crit1.plot, crit2.plot, crit3.subplot, crit3.plot.leg, ncol = 4, widths = c(2.3,2.3,2.3,0.8))

  dev.off()


 ## summarize criteria for each inversion
# first find average for each criteria for non-inverted segments 
outQTNs <- df.qtnMuts.MAF[df.qtnMuts.MAF$inOut == "Outside Inversion",]
genome.sizeNoNeut <- 2000000
inverted.genome.size <- nrow(invWindBases.MAF)
non.inverted.genome.size <- genome.sizeNoNeut - inverted.genome.size

non.inv.crit1 <- sum(outQTNs$crit1_p.value == 0)/non.inverted.genome.size
non.inv.crit2 <- sum(outQTNs$crit2_p.value == 0)/non.inverted.genome.size
Va_perc_Out <- sum(outQTNs$crit3_Va_perc)
non.inv.crit3 <- Va_perc_Out/non.inverted.genome.size

## loop over each inversion that is MAF filtered
inv.IDs <- unique(invWindBases.MAF$invIDMAF)
crit.output <- matrix(NA, nrow = length(inv.IDs), ncol = 12)
colnames(crit.output) <- c("invWindID", "length", "first.base", "final.base", "num.QTNs", "crit1", "crit2", 
                           "inv_Va_perc", "crit3", "crit1_YN", "crit2_YN", "crit3_YN")
adapt.inv <- NULL

for(i in 1:length(inv.IDs)){
  focalWindow <- invWindBases.MAF[invWindBases.MAF$invIDMAF == inv.IDs[i],]
  wind.length <- nrow(focalWindow)
  first.wind.base <- focalWindow$invWindBasesMAF[1]
  final.wind.base <- focalWindow$invWindBasesMAF[wind.length]
  focalWindQTNs <- subset(df.qtnMuts.MAF, subset = position >= first.wind.base & position <= final.wind.base)
  focalWind.crit1 <- sum(focalWindQTNs$crit1_p.value == 0)/wind.length
  focalWind.crit2 <- sum(focalWindQTNs$crit2_p.value == 0)/wind.length
  Va_focInv <- sum(focalWindQTNs$crit3_Va_perc)
  focalWind.crit3 <- sum(focalWindQTNs$crit3_Va_perc)/wind.length
  crit1_YN <- ifelse(focalWind.crit1 > non.inv.crit1, 1, 0)
  crit2_YN <- ifelse(focalWind.crit2 > non.inv.crit2, 1, 0)
  crit3_YN <- ifelse(focalWind.crit3 > non.inv.crit3, 1, 0)
  crit.output[i,]<- c(inv.IDs[i], wind.length, first.wind.base, final.wind.base, nrow(focalWindQTNs), 
                      focalWind.crit1, focalWind.crit2, Va_focInv, focalWind.crit3, crit1_YN, crit2_YN, crit3_YN)
  if(sum(crit1_YN, crit2_YN, crit3_YN) == 3){
    adapt.inv <- c(adapt.inv, inv.IDs[i])
  }
}
df.crit.output <- as.data.frame(crit.output)
df.crit.output$crit1_thres <- non.inv.crit1
df.crit.output$crit2_thres <- non.inv.crit2
df.crit.output$crit3_thres <- non.inv.crit3 

df.crit.output
inQTNs <- df.qtnMuts.MAF[df.qtnMuts.MAF$inOut == "Inside Inversion",]
Va_perc_In <- sum(inQTNs$crit3_Va_perc)
#Va_perc_In2 <- sum(crit.output[,8]) # Extra ten percent because of overlap between some inversions so some qtns get counted more than once. 
#Need to calulate this with the first way. 

total.perc <- sum(Va_perc_In,Va_perc_Out) 
total.perc == 100
#sanity check???? why doesn't this work

### end identifying inversion window QTNs
######################################################################################################  

######################################################################################################
#### Subset Inversions ####

### SELECTION SIMS ###
## Adaptive Inversions
## Summary stats for each gen for line graph
adapt.inv.data <- df.invAllData %>%
  filter(inv_id %in% adapt.inv) %>%
  group_by(sim_gen) %>%
  summarise_at(c("inv_age", "mean_qtnSelCoef", "num_qtns", "inv_length", "num_qtns_Lscaled"), 
               mean, .groups = "keep") %>%
  rename(inv_ageAdapt = inv_age, mean_qtnSelCoefAdapt = mean_qtnSelCoef, 
         num_qtnsAdapt = num_qtns, inv_lengthAdapt = inv_length, 
         num_qtns_LscaledAdapt = num_qtns_Lscaled)

# no averaging for boxplots
adapt.inv.data.nosum <- df.invAllData %>%
  group_by(sim_gen) %>%
  filter(inv_id %in% adapt.inv)

adapt.inv.data.nosum$adaptInv <- "Adaptive"

## Non adaptive inversions
## Summary stats for each gen for line graph
non.adapt.inv.data <- df.invAllData %>%
  filter(!inv_id %in% adapt.inv) %>%
  group_by(sim_gen) %>%
  summarise_at(c("inv_age", "mean_qtnSelCoef", "num_qtns", "inv_length", "num_qtns_Lscaled"), 
               mean, .groups = "keep") %>%
  rename(inv_ageNonAdapt = inv_age, mean_qtnSelCoefNonAdapt = mean_qtnSelCoef, 
         numqtnsNonAdapt = num_qtns, inv_lengthNonAdapt = inv_length,
         num_qtns_LscaledNonAdapt = num_qtns_Lscaled) 

# no averaging for boxplots
non.adapt.inv.data.nosum <- df.invAllData %>%
  group_by(sim_gen) %>%
  filter(!inv_id %in% adapt.inv)
non.adapt.inv.data.nosum$adaptInv <- "Nonadaptive"

# bind together adaptive and non adaptive column together 
#adaptInvCol <- as.data.frame(c(adapt.inv.data.nosum$adaptInv, non.adapt.inv.data.nosum$adaptInv))

# bind together datasets for boxplots 
df.adaptSplitbox <- rbind(as.data.frame(adapt.inv.data.nosum), as.data.frame(non.adapt.inv.data.nosum))
#df.adaptSplitbox <- cbind(df.adaptSplitbox, adaptInvCol)
#colnames(df.adaptSplitbox)[ncol(df.adaptSplitbox)] <- "adaptInv"

# Standard deviation
sd.Adapt <- aggregate(cbind(inv_age, mean_qtnSelCoef, num_qtns, inv_length, num_qtns_Lscaled)~sim_gen, 
                      data = adapt.inv.data.nosum, FUN = sd)
colnames(sd.Adapt)[2:6] <- c("sd_inv_ageAdapt", "sd_qtnSelCoefAdapt", "sd_num_qtnsAdapt", 
                             "sd_inv_lengthAdapt", "sd_num_qtns_LscaledAdapt")
sd.NonAdapt <- aggregate(cbind(inv_age, mean_qtnSelCoef, num_qtns, inv_length, 
                               num_qtns_Lscaled)~sim_gen, data = non.adapt.inv.data.nosum, FUN = sd)
colnames(sd.NonAdapt)[2:6] <- c("sd_inv_ageNonAdapt", "sd_qtnSelCoefNonAdapt", "sd_num_qtnsNonAdapt", 
                                "sd_inv_lengthNonAdapt", "sd_num_qtns_LscaledNonAdapt")              

## Join dataframes with parameters
df.AdaptSplitTb <- full_join(adapt.inv.data, non.adapt.inv.data, by = "sim_gen")

## convert to data frame and factor parameter columns
df.AdaptSplitTemp <- as.data.frame(df.AdaptSplitTb)  
df.AdaptSplitTemp2 <- left_join(df.AdaptSplitTemp, sd.Adapt, by = "sim_gen")
df.AdaptSplit <- left_join(df.AdaptSplitTemp2, sd.NonAdapt, by = "sim_gen")

### NO SELECTION SIMS ###

## Summary stats for each gen for line graph
inv.data.NS <- df.invAllData.NS %>%
  group_by(sim_gen) %>%
  summarise_at(c("inv_age", "mean_qtnSelCoef", "num_qtns", "inv_length", "num_qtns_Lscaled"), 
                                mean, .groups = "keep")

sd.NS <- aggregate(cbind(inv_age, mean_qtnSelCoef, num_qtns, inv_length, num_qtns_Lscaled)~sim_gen, 
                      data = df.invAllData.NS, FUN = sd)
colnames(sd.NS)[2:6] <- c("sd_inv_age_NS", "sd_qtnSelCoef_NS", "sd_num_qtns_NS", 
                             "sd_inv_length_NS", "sd_num_qtns_Lscaled_NS")

inv.data.NS <- as.data.frame(inv.data.NS)
colnames(inv.data.NS)[2:6] <- c("inv_age_NS", "mean_qtnSelCoef_NS", "num_qtns_NS", "inv_length_NS", "num_qtns_Lscaled_NS")

df.inv.data.NS <- full_join(inv.data.NS, sd.NS, by = "sim_gen")

df.AdaptNSsplit <- full_join(df.AdaptSplit, df.inv.data.NS, by = "sim_gen")

# no averaging for boxplots
df.invAllData.NS$adaptInv <- "No selection"

# bind all data for boxplotting together
df.adaptSplitboxNS <- rbind(df.adaptSplitbox, df.invAllData.NS)
df.adaptSplitboxNS$adaptInv <- as.factor(df.adaptSplitboxNS$adaptInv)


## end subset inversions
######################################################################################################


#OUTPUT: Va_perc_In, LA final Gen, 

######################################################################################################
#### Inital full data plotting #####

## Local adaptation
  LA.title <- bquote(atop(paste(bold("Local Adaptation: "), italic(m), " = ", .(df.params$mig1), 
                                " ", mu[inv], " = ",  .(df.params$muInv)), 
                          atop(textstyle(paste(" ", mu, " = ", .(format(df.params$muBase, scientific = T, digits = 2)), 
                                              " ", alpha, " = ", .(df.params$alpha), " ", sigma[K], 
                                              " = ", .(df.params$sigmaK),
                                              " ", sigma[env], " = ", .(df.params$enVar))))))

  LA.plot <- ggplot(data = df.popDyn, 
                    aes(x = sim_gen, y = localAdaptSA)) + 
    geom_line(color = "cadetblue3", size = 0.75) + 
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = LA.title,
        y = "Local Adaptation",
        x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(-0.1, 1))

  par(mar = c(1,1,1,2))
  
  pdf(paste0(folderOut, seed, "_LA.pdf"), height = 5, width = 7)

    LA.plot + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

  dev.off()


## Phenotypes ##
  df.pheno <- pivot_longer(df.popDyn[, c(1,10,14)], cols = c(meanPhenoP1, meanPhenoP2),
                            names_to = "pop", values_to = "meanPheno")
  pheno.plot.leg <- ggplot(data = df.pheno, 
                        aes(x = sim_gen, y = meanPheno, group = pop)) + 
    geom_line(aes(color = pop), size = 0.75) + 
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = expression(bold("Selection")), y = "Phenotype", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(min(df.popDyn$meanPhenoP2 - df.popDyn$sdPhenoP2- 0.1),
                                                    max(df.popDyn$meanPhenoP1 + df.popDyn$sdPhenoP1 + 0.1))) + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
    scale_color_manual(name = "Population",
                       values = c( "cadetblue3", "navy"),
                       labels = c("Pop 1", "Pop 2"))

  pheno.plot <- ggplot(data = df.popDyn, 
                         aes(x = sim_gen, y = meanPhenoP1)) + 
    geom_line(color = "cadetblue3", size = 0.75) + 
    geom_line(aes(y = meanPhenoP2, x = sim_gen), color = "navy") + 
    geom_ribbon(aes(ymin= meanPhenoP1 - sdPhenoP1, ymax= meanPhenoP1 + sdPhenoP1), fill = "cadetblue3", alpha=0.2) +
    geom_ribbon(aes(ymin= meanPhenoP2 - sdPhenoP2, ymax= meanPhenoP2 + sdPhenoP2), fill = "navy", alpha=0.2) +
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = expression(bold("Selection")), y = "Phenotype", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(min(df.popDyn$meanPhenoP2 - df.popDyn$sdPhenoP2- 0.1),
                                                    max(df.popDyn$meanPhenoP1 + df.popDyn$sdPhenoP1 + 0.1))) + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
   
 

  pheno.plot.NS <- ggplot(data = df.popDyn.NS, 
                            aes(x = sim_gen, y = meanPhenoP1)) + 
    geom_line(size = 0.75, color = "cadetblue3") + 
    geom_ribbon(aes(ymin= meanPhenoP1 - sdPhenoP1, ymax= meanPhenoP1 + sdPhenoP1), fill = "cadetblue3", alpha=0.2) +
    geom_line(aes(y = meanPhenoP2, x = sim_gen), color = "navy") + 
    geom_ribbon(aes(ymin= meanPhenoP2 - sdPhenoP2, ymax= meanPhenoP2 + sdPhenoP2), fill = "navy", alpha=0.2) +
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = expression(bold("No Selection")), y = " ", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(min(df.popDyn$meanPhenoP2 - df.popDyn$sdPhenoP2- 0.1),
                                                  max(df.popDyn$meanPhenoP1 + df.popDyn$sdPhenoP1 + 0.1))) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  

  pheno.leg <- g_legend(pheno.plot.leg)

  pdf(paste0(folderOut, seed, "_pheno.pdf"), height = 5, width = 7)

    ggarrange(pheno.plot, pheno.plot.NS, pheno.leg, ncol = 3, widths = c(2.3,2.3,0.8))

  dev.off()

## Fitnesses ##
  fitP1.plot <- ggplot(data = df.popDyn, 
                       aes(x = sim_gen, y = meanFitP1)) + 
    geom_line(color = "cadetblue3", size = 0.75) + 
    geom_ribbon(aes(ymin=  meanFitP1 - sdFitP1, ymax= meanFitP1 + sdFitP1), fill = "cadetblue3", alpha=0.2) +
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = expression(bold("Population 1")), y = "Fitness", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits =c((min(df.popDyn$meanFitP1 - df.popDyn$sdFitP1, 
                                                        df.popDyn$meanFitP2 - df.popDyn$sdFitP2) - 0.1),
                                                  ((max(df.popDyn$meanFitP1 + df.popDyn$sdFitP1, 
                                                        df.popDyn$sdFitP2 + df.popDyn$sdFitP2) + 0.1)))) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 


  fitP2.plot <- ggplot(data = df.popDyn, 
                       aes(x = sim_gen, y = meanFitP2)) + 
    geom_line(color = "navy", size = 0.75) + 
    geom_ribbon(aes(ymin=  meanFitP2 - sdFitP2, ymax= meanFitP2 + sdFitP2), fill = "navy", alpha=0.2) +
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = expression(bold("Population 2")), y = " ", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c((min(df.popDyn$meanFitP1 - df.popDyn$sdFitP1, 
                                                        df.popDyn$meanFitP2 - df.popDyn$sdFitP2) - 0.1),
                                                    ((max(df.popDyn$meanFitP1 + df.popDyn$sdFitP1, 
                                                          df.popDyn$sdFitP2 + df.popDyn$sdFitP2) + 0.1)))) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 

  pdf(paste0(folderOut, seed, "_fitness.pdf"), height = 5, width = 7)

    ggarrange(fitP1.plot, fitP2.plot)

  dev.off()

## Average Inversion Age ##
  library("viridisLite")
  # SELECTION #
  df.invage.temp <- pivot_longer(df.AdaptNSsplit[,c(1,2,7,22)], cols = c(inv_ageAdapt, inv_ageNonAdapt, inv_age_NS),
                            names_to = "Adaptsplit", values_to = "inv_age")
  df.invage <- df.invage.temp %>%
    mutate(Adaptsplit = fct_relevel(Adaptsplit, "inv_ageAdapt","inv_ageNonAdapt","inv_age_NS"))
  df.AdaptSplit$inv_ageAdaptLower <- df.AdaptSplit$inv_ageAdapt - df.AdaptSplit$sd_inv_ageAdapt
  df.AdaptSplit$inv_ageNonAdaptLower <- df.AdaptSplit$inv_ageNonAdapt - df.AdaptSplit$sd_inv_ageNonAdapt
  df.inv.data.NS$inv_ageNSLower <- df.inv.data.NS$inv_age_NS - df.inv.data.NS$sd_inv_age_NS
  df.AdaptSplit$inv_ageAdaptLower[df.AdaptSplit$inv_ageAdaptLower<0] <- 0
  df.AdaptSplit$inv_ageNonAdaptLower[df.AdaptSplit$inv_ageNonAdaptLower<0] <- 0
  df.inv.data.NS$inv_ageNSLower[df.inv.data.NS$inv_ageNSLower<0] <- 0
  
  inv.age.plot <- ggplot(data = df.invage, 
                         aes(x = sim_gen, y = inv_age, group = Adaptsplit)) + 
    geom_line(data = df.invage, aes(color = Adaptsplit), size = 0.75) + 
    geom_ribbon(data = df.AdaptSplit, aes(x = sim_gen, ymin=  inv_ageAdaptLower, 
                                           ymax= inv_ageAdapt + sd_inv_ageAdapt), 
               fill = "#3E4A89FF", alpha=0.2, inherit.aes = FALSE) +
    geom_ribbon(data = df.AdaptSplit, aes(x = sim_gen, ymin = inv_ageNonAdaptLower,
                                           ymax = inv_ageNonAdapt + sd_inv_ageNonAdapt),
                fill = "firebrick", alpha=0.2, inherit.aes = FALSE) +
    geom_ribbon(data = df.inv.data.NS, aes(x = sim_gen, ymin =  inv_ageNSLower,
                                          ymax = inv_age_NS + sd_inv_age_NS),
               fill = "black", alpha=0.2, inherit.aes = FALSE) +
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = " ", y = "Average Inversion Age", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_color_manual(name = "", labels = c("Adaptive", "Nonadaptive", "No Selection"), 
                       values=c( "#3E4A89FF", "red", "black")) +
    #theme(legend.position = "none") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.invage$inv_age))) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  
  
  # NO SELECTION #
  #df.invage.NS <- pivot_longer(df.AdaptSplit.NS[, c(1,2,7)], cols = c(inv_ageAdapt, inv_ageNonAdapt),
                              # names_to = "AdaptSplit", values_to = "inv_age")
  # 
  # inv.age.plot.NS <- ggplot(data = inv.data.NS, 
  #                           aes(x = sim_gen, y = inv_age)) + 
  #   geom_line(aes(color = FSTsplit), size = 0.75) + 
  #   geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
  #   labs(title = " ", y = "", x = "Generation") +
  #   theme_classic() +
  #   theme(panel.background = element_blank(), 
  #         strip.background = element_rect(colour = "white", fill = "grey92")) +
  #   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  #   scale_y_continuous(expand = c(0, 0), limits = c(0, max(inv.data.NS$inv_age))) +
  #   theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  # 
  # leg <- g_legend(inv.age.plot.NS)
  # 
  # inv.age.plot.NS.noleg <- inv.age.plot.NS + theme(legend.position = "none")
  # 
  pdf(paste0(folderOut, seed, "_invAge3.pdf"), height = 5, width = 7)
   inv.age.plot
  dev.off()

  
  # Subset for final generation
  final.inv <- df.adaptSplitboxNS[df.adaptSplitboxNS$sim_gen == 50000,]
  final.inv$adaptInv <- factor(final.inv$adaptInv, levels = c("Adaptive", "Nonadaptive", "No selection"))
  
  pdf(paste0(folderOut, seed, "_invAgebox.pdf"), height = 5, width = 7)
    ggplot(data = final.inv, aes(x = adaptInv, y= inv_age, fill = adaptInv)) +
      geom_boxplot() + 
      scale_fill_manual(values = c( "#3E4A89FF", "firebrick", "black")) + 
      theme_classic() +
      theme(panel.background = element_blank(), 
            strip.background = element_rect(colour = "white", fill = "grey92")) +
      labs(title = " ", y = "Average Inversion Age", x = "Inversion Status") 
  dev.off()
  
## Average Inversion Length ##
  # SELECTION #  
  df.invlength.temp <- pivot_longer(df.AdaptNSsplit[, c(1,5,10,25)], cols = c(inv_lengthAdapt, inv_lengthNonAdapt, inv_length_NS),
                               names_to = "Adaptsplit", values_to = "inv_length")
  df.invlength <- df.invlength.temp %>%
    mutate(Adaptsplit = fct_relevel(Adaptsplit, "inv_lengthAdapt","inv_lengthNonAdapt","inv_length_NS"))
  df.AdaptSplit$inv_lengthAdaptLower <- df.AdaptSplit$inv_lengthAdapt - df.AdaptSplit$sd_inv_lengthAdapt
  df.AdaptSplit$inv_lengthNonAdaptLower <- df.AdaptSplit$inv_lengthNonAdapt - df.AdaptSplit$sd_inv_lengthNonAdapt
  df.inv.data.NS$inv_lengthNSLower <- df.inv.data.NS$inv_length_NS - df.inv.data.NS$sd_inv_length_NS
  df.AdaptSplit$inv_lengthAdaptLower[df.AdaptSplit$inv_lengthAdaptLower<0] <- 0
  df.AdaptSplit$inv_lengthNonAdaptLower[df.AdaptSplit$inv_lengthNonAdaptLower<0] <- 0
  df.inv.data.NS$inv_lengthNSLower[df.inv.data.NS$inv_lengthNSLower<0] <- 0
  #"paleturquoise2", "skyblue", "skyblue4"
  inv.length.plot <- ggplot(data = df.invlength, 
                            aes(x = sim_gen, y = inv_length, group = Adaptsplit)) + 
    geom_line(aes(color = Adaptsplit), size = 0.75, alpha = 0.9) + 
    geom_ribbon(data = df.AdaptSplit, aes(x = sim_gen, ymin=  inv_lengthAdaptLower, 
                                          ymax= inv_lengthAdapt + sd_inv_lengthAdapt), 
                fill = "#3E4A89FF", alpha=0.2, inherit.aes = FALSE) +
    geom_ribbon(data = df.AdaptSplit, aes(x = sim_gen, ymin = inv_lengthNonAdaptLower,
                                          ymax = inv_lengthNonAdapt + sd_inv_lengthNonAdapt),
                fill = "firebrick", alpha=0.2, inherit.aes = FALSE) +
    geom_ribbon(data = df.inv.data.NS, aes(x = sim_gen, ymin =  inv_lengthNSLower,
                                           ymax = inv_length_NS + sd_inv_length_NS),
                fill = "black", alpha=0.2, inherit.aes = FALSE) +
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = " ", y = "Average Inversion Length", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_color_manual(name = "", labels = c( "Adaptive", "Nonadaptive","No Selection"),
                       values=c("#3E4A89FF", "firebrick", "black")) +
   # theme(legend.position = "none") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 50000)) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  
  
  # NO SELECTION #
  # df.invlength.NS <- pivot_longer(df.FSTsplit.NS[, c(1,5,10)], cols = c(inv_lengthT10, inv_lengthB90),
  #                                 names_to = "FSTsplit", values_to = "inv_length")
  # 
  # inv.length.plot.NS <- ggplot(data = df.invlength.NS, 
  #                              aes(x = sim_gen, y = inv_length, group = FSTsplit)) + 
  #   geom_line(aes(color = FSTsplit), size = 0.75, alpha= 0.9) + 
  #   geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
  #   labs(title = " ", y = "", x = "Generation") +
  #   theme_classic() +
  #   theme(panel.background = element_blank(), 
  #         strip.background = element_rect(colour = "white", fill = "grey92")) +
  #   scale_color_manual(labels = c("Top 10%", "Bot 90%"),
  #                      values=c("lightsalmon1", "lightsalmon3")) +
  #   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  #   scale_y_continuous(expand = c(0, 0), limits = c(0, 50000)) +
  #   theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  # 
  # legLeng <- g_legend(inv.length.plot.NS)
  # 
  # inv.length.plot.NS.noLeg <- inv.length.plot.NS + theme(legend.position = "none")
  # 
  pdf(paste0(folderOut, seed, "_invLength.pdf"), height = 5, width = 7)
  # ggarrange(inv.length.plot, inv.length.plot.NS.noLeg, legLeng, labels = c("Selection", "No Selection"),
  #           ncol = 3, widths = c(2.3,2.3,0.8))
    inv.length.plot
  dev.off()

  pdf(paste0(folderOut, seed, "_invLengthbox.pdf"), height = 5, width = 7)
  ggplot(data = final.inv, aes(x = adaptInv, y= inv_length, fill = adaptInv)) +
    geom_boxplot() + 
    scale_fill_manual(values = c("#3E4A89FF", "firebrick", "black")) + 
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    labs(title = " ", y = "Average Inversion Length", x = "Inversion Status") 
  dev.off()
## Inversion QTNs Length Scaled ##
  # Selection
  df.invQTNsLscaled.temp <- pivot_longer(df.AdaptNSsplit[, c(1,6,11,26)], cols = c(num_qtns_LscaledAdapt, num_qtns_LscaledNonAdapt, num_qtns_Lscaled_NS),
                                    names_to = "Adaptsplit", values_to = "inv_numQTNs")
  df.invQTNsLscaled <- df.invQTNsLscaled.temp %>%
    mutate(Adaptsplit = fct_relevel(Adaptsplit, "num_qtns_LscaledAdapt", "num_qtns_LscaledNonAdapt", "num_qtns_Lscaled_NS"))
  
  df.AdaptSplit$inv_numQTNsLscaledAdaptLower <- df.AdaptSplit$num_qtns_LscaledAdapt - df.AdaptSplit$sd_num_qtns_LscaledAdapt
  df.AdaptSplit$inv_numQTNsLscaledNonAdaptLower <- df.AdaptSplit$num_qtns_LscaledNonAdapt - df.AdaptSplit$sd_num_qtns_LscaledNonAdapt
  df.inv.data.NS$inv_numQTNsLscaledNSLower <- df.inv.data.NS$num_qtns_Lscaled_NS - df.inv.data.NS$sd_num_qtns_Lscaled_NS
  df.AdaptSplit$inv_numQTNsLscaledAdaptLower[df.AdaptSplit$inv_numQTNsLscaledAdaptLower<0] <- 0
  df.AdaptSplit$inv_numQTNsLscaledNonAdaptLower[df.AdaptSplit$inv_numQTNsLscaledNonAdaptLower<0] <- 0
  df.inv.data.NS$inv_numQTNsLscaledNSLower[df.inv.data.NS$inv_numQTNsLscaledNSLower<0] <- 0
  # thistle, plum4, thistle4
  inv.qtns.Lscaled.plot <- ggplot(data = df.invQTNsLscaled, 
                                  aes(x = sim_gen, y = inv_numQTNs, group = Adaptsplit)) + 
    geom_line(aes(color = Adaptsplit), size = 0.75) + 
    geom_ribbon(data = df.AdaptSplit, aes(x = sim_gen, ymin= inv_numQTNsLscaledAdaptLower, 
                                          ymax= num_qtns_LscaledAdapt + sd_num_qtns_LscaledAdapt), 
                fill = "#3E4A89FF", alpha=0.2, inherit.aes = FALSE) +
    geom_ribbon(data = df.AdaptSplit, aes(x = sim_gen, ymin = inv_numQTNsLscaledNonAdaptLower,
                                          ymax = num_qtns_LscaledNonAdapt + sd_num_qtns_LscaledNonAdapt),
                fill = "firebrick", alpha=0.2, inherit.aes = FALSE) +
    geom_ribbon(data = df.inv.data.NS, aes(x = sim_gen, ymin =  inv_numQTNsLscaledNSLower,
                                           ymax = num_qtns_Lscaled_NS + sd_num_qtns_Lscaled_NS),
                fill = "black", alpha=0.2, inherit.aes = FALSE) +
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = "",
         y = "Average number of inversion QTNs \nscaled by inversion length",
         x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_color_manual(name = "", labels = c( "Adaptive", "Nonadaptive","No Selection"),
                       values=c("#3E4A89FF", "firebrick", "black")) +
    #theme(legend.position = "none") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.invQTNsLscaled$inv_numQTNs))) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  
  # No Selection    
  # df.invQTNs.Lscaled.NS <- pivot_longer(df.FSTsplit.NS[, c(1,6,11)], cols = c(num_qtns_LscaledT10, num_qtns_LscaledB90),
  #                                       names_to = "FSTsplit", values_to = "inv_qtnNum")
  # 
  # inv.qtns.Lscaled.plot.NS <- ggplot(data = df.invQTNs.Lscaled.NS, 
  #                                    aes(x = sim_gen, y = inv_qtnNum, group = FSTsplit)) + 
  #   geom_line(aes(color = FSTsplit), size = 0.75) + 
  #   geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
  #   labs(title = " ", y = "", x = "Generation") +
  #   theme_classic() +
  #   theme(panel.background = element_blank(), 
  #         strip.background = element_rect(colour = "white", fill = "grey92")) +
  #   scale_color_manual(labels = c("Top 10%", "Bot 90%"), 
  #                      values=c("plum4", "thistle")) +
  #   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  #   scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.invQTNsLscaled$inv_qtnNum))) +
  #   theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  # 
  # legQTNsLscaled <- g_legend(inv.qtns.Lscaled.plot.NS)
  # 
  # inv.qtns.Lscaled.plot.NS.noleg <- inv.qtns.Lscaled.plot.NS + theme(legend.position = "none")
  # 
  pdf(paste0(folderOut, seed, "_invQTNsLscaled2.pdf"), height = 5, width = 7)
  inv.qtns.Lscaled.plot
  # ggarrange(inv.qtns.Lscaled.plot, inv.qtns.Lscaled.plot.NS.noleg, legQTNsLscaled, labels = c("Selection", "No Selection"),
  #           ncol = 3, widths = c(2.3,2.3,0.8))
  dev.off()
  
  pdf(paste0(folderOut, seed, "_invQTNsLscaledbox.pdf"), height = 5, width = 7)
    ggplot(data = final.inv, aes(x = adaptInv, y= num_qtns_Lscaled, fill = adaptInv)) +
      geom_boxplot() + 
      scale_fill_manual(values = c("#3E4A89FF", "firebrick", "black")) + 
      theme_classic() +
      theme(panel.background = element_blank(), 
            strip.background = element_rect(colour = "white", fill = "grey92")) +
      labs(title = " ", y = "Average number of inversion QTNs \nscaled by inversion length", x = "Inversion Status") 
  dev.off()
### end initial full data plotting
######################################################################################################
  
  
  
 
######################################################################################################  
#### Manhattan plot ####
  # SELECTION #
  df.neutQTNmuts <- df.muts.MAF[df.muts.MAF$inOut != "inv",]
  df.neutQTNmuts$inOut <- factor(df.neutQTNmuts$inOut)
  manh.plot <- ggplot(df.neutQTNmuts, aes(x = position, y = FST, group = inOut)) + 
    geom_point(aes(color = inOut)) + 
    scale_color_manual(values=c( "red", "goldenrod", "navy")) +
    xlim(0, 2100000) + 
    ylim(0, max(df.neutQTNmuts$FST)) +
    labs(title = expression(bold("Selection"))) + 
    theme(legend.position = "none") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92"),
          text = element_text(size = 11)) 
  
  # NO SELECTION #
  df.neutQTNmuts.NS <- df.muts.NS.MAF[df.muts.NS.MAF$inOut != "inv",]
  df.neutQTNmuts.NS$inOut <- factor(df.neutQTNmuts.NS$inOut)
  manh.plot.NS <- ggplot(df.neutQTNmuts.NS, aes(x = position, y = FST, group = inOut)) + 
    geom_point(aes(color = inOut)) + 
    scale_color_manual(values=c( "red", "goldenrod", "navy")) +
    xlim(0, 2100000) + 
    ylim(0, max(df.neutQTNmuts$FST)) +
    labs(title = expression(bold("No Selection")))  + 
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92"),
          text = element_text(size = 11)) 
  
  # No legend
  manh.plot.NS.noleg <- manh.plot.NS + theme(legend.position = "none")
  manh.plot.noleg <-  manh.plot+ theme(legend.position = "none")
  legManh <- g_legend(manh.plot.NS)
  
  ## TO DO: ADD CHROMOSOME COLORS
  png(paste0(folderOut, seed, "_manh.png"), width = 700, height = 400, units = "px")
  ggarrange(manh.plot.noleg, manh.plot.NS.noleg, legManh, ncol = 3, widths = c(2.3,2.3,0.8))
  dev.off()
  
  pdf(paste0(folderOut, seed, "_manh.pdf"), height = 5, width = 7)
  ggarrange(manh.plot.noleg, manh.plot.NS.noleg, legManh, ncol = 3, widths = c(2.3,2.3,0.8))
  dev.off()
  
#### end manhattan plot
######################################################################################################  

  
  
######################################################################################################   
#### plot inversion origin dynamics ####
  # Subset for final generation
  df.invFinalGen <- subset(df.invTime, subset = sim_gen == 50000)
  df.InvDataOrigin <- left_join(df.invFinalGen, df.invData, by = "inv_id")
  df.invFinalGen.NS <- subset(df.invTime.NS, subset = sim_gen == 50000)
  df.InvDataOrigin.NS <- left_join(df.invFinalGen.NS, df.invData.NS, by = "inv_id")
  
  # subset for MAF > 0.01
  df.InvDataOriginMAF <- df.InvDataOrigin[df.InvDataOrigin$freq > 0.01,]
  for(i in 1:nrow(df.InvDataOriginMAF)){
    if(df.InvDataOriginMAF$freq_p1[i] > df.InvDataOriginMAF$freq_p2[i]){
      df.InvDataOriginMAF$pop[i] <- "pop1"
    } else {
      df.InvDataOriginMAF$pop[i] <- "pop2"
      
    }
  }
  
  df.InvDataOriginMAF.NS <- df.InvDataOrigin.NS[df.InvDataOrigin.NS$freq > 0.01,]
  for(i in 1:nrow(df.InvDataOriginMAF.NS)){
    if(df.InvDataOriginMAF.NS$freq_p1[i] > df.InvDataOriginMAF.NS$freq_p2[i]){
      df.InvDataOriginMAF.NS$pop[i] <- "pop1"
    } else {
      df.InvDataOriginMAF.NS$pop[i] <- "pop2"
      
    }
  }
  
  ## Subset dataframe to get how the MAF filtered inversions change through time
  inv.IDs <- as.vector(df.InvDataOriginMAF$inv_id)
  df.invFinalAllData <- df.invTime[df.invTime$inv_id %in% inv.IDs, ]
  df.invFinalAllData$qtnSelCoefsum <- df.invFinalAllData$mean_qtnSelCoef*df.invFinalAllData$num_qtns
  df.invFinalAllDataPop <- left_join(df.invFinalAllData,
                                     df.InvDataOriginMAF[c(2,13:19)], 
                                     by = "inv_id")
  inv.IDs.NS <- as.vector(df.InvDataOriginMAF.NS$inv_id)
  df.invFinalAllData.NS <- df.invTime.NS[df.invTime.NS$inv_id %in% inv.IDs.NS, ]
  df.invFinalAllData.NS$qtnSelCoefsum <- df.invFinalAllData.NS$mean_qtnSelCoef*df.invFinalAllData.NS$num_qtns
  df.invFinalAllDataPop.NS <- left_join(df.invFinalAllData.NS,
                                     df.InvDataOriginMAF.NS[c(2,13:19)], 
                                     by = "inv_id")
  
  df.invFinalsubset <- df.invFinalAllDataPop %>% filter(sim_gen %in% seq(0, 50000, by = 1000)) 
  df.invFinalsubset.NS <- df.invFinalAllDataPop.NS %>% filter(sim_gen %in% seq(0, 50000, by = 1000)) 
  
  plot.inv.orig <- ggplot(df.invFinalsubset, aes(x = sim_gen, y = qtnSelCoefsum)) + 
    geom_point(aes(color = pop, size = inv_FST, alpha = inv_id)) + 
    geom_line(aes(color = pop, group = inv_id, alpha = inv_id)) + 
    scale_color_manual(values=c("navy", "red")) + 
    scale_size(range = c(0.5, 4), breaks = c(0.00001, 0.05, 0.15, 0.2)) +   
    #scale_alpha(range = c( 1, 0.2)) +
    #scale_alpha_discrete(range = c(0.35, 0.9)) +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92"),
          text = element_text(size = 15)) +
    labs(title = expression(bold("Selection")),
         y = "sum of each Inversion QTNs \neffects on phenotype",
         x = "Generation") +
    ylim(-0.2, 0.2) +
    xlim(0,50000) +
    #theme(legend.position = "none")
    guides(color = guide_legend(title = "Pop with Highest\nFrequency of Inv")) +
    guides(size = guide_legend(title = "Inversion FST")) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  plot.inv.orig.NS <- ggplot(df.invFinalsubset.NS, aes(x = sim_gen, y = qtnSelCoefsum)) + 
    geom_point(aes(color = pop, size = inv_FST), alpha = 0.8) + 
    geom_line(aes(color = pop, group = inv_id), alpha = 0.8) + 
    scale_color_manual(values=c("navy", "red")) + 
    #scale_color_gradient(low = "gainsboro", high = "darkorange", 
                         #name = "Days Since \nStart of \nExperiment", 
                         #breaks = c(4, 8, 12, 16), 
                         #labels = c("Day 6 PM", "Day 8 AM", "Day 10 AM", "Day 13 PM")) + 
    scale_size(range = c(0.5, 4), breaks = c(0.00001, 0.05, 0.15, 0.2)) + 
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92"),
          text = element_text(size = 15)) +
    labs(title = expression(bold("No Selection")),
         y = "",
         x = "Generation") +
    ylim(-0.2, 0.2) +
    xlim(0,50000) +
    theme(legend.position = "none") +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    #guides(color = guide_legend(title = "Pop with Highest\nFrequency of Inv")) +
    #guides(size = guide_legend(title = "Inversion FST")) 
  
  plot.inv.orig.noleg <- plot.inv.orig + theme(legend.position = "none")
  
  legInvOrig <- g_legend(plot.inv.orig)
  pdf(paste0(folderOut, seed, "_invOrigin.pdf"), height = 5, width = 12)
  ggarrange( plot.inv.orig.noleg, plot.inv.orig.NS, legInvOrig, ncol = 3, widths = c(2.3,2.3,0.8))
  dev.off()
  
#### end plot origin dynamics
######################################################################################################  

######################################################################################################  
#### process VCF ####
  
  # Inspect Individual Data 
  dim(df.indPheno)
  head(df.indPheno)
  tail(df.indPheno)
  
  # Mutation stats at end of sim
  # for all muts
  dim(df.muts)
  head(df.muts)
  table(df.muts$type) 
  hist(as.numeric(as.character(df.muts$FST)))
  hist(df.muts$freq)
  
  # inversion summary
  head(df.invData)
  
  # VCF file
  #vcffile <- list.files(path=path, pattern=paste0(".vcf"))
  #vcf <- read.vcfR(paste0(folder,seed, "_InversionVCF.vcf"))
  vcf.MAF <- read.vcfR(paste0(folderIn,seed, "_InversionVCF_MAF01.recode.vcf"))
  
  head(vcf.MAF)
  head(vcf.MAF@fix, 50)
  dim(vcf.MAF@fix)
  
  # subset mutations file for MAF > 0.01
  freq <- df.muts$freq
  altFreq <- 1 - df.muts$freq
  df.calc <- cbind(freq, altFreq)
  df.muts$MAF <- apply(df.calc, 1, min)
  df.mutsMAF <- subset(df.muts, subset = MAF >= 0.01)
  dim(df.mutsMAF)
  df.mutsMAF$FST <- as.numeric(as.character(df.mutsMAF$FST))
  colnames(df.mutsMAF)[2] <- "LocusName"
  
  # example of how to find a specific mutation in the vcf file
  df.mutsMAF[2,]
  vcf.MAF@fix[grep(df.mutsMAF$LocusName[1], vcf.MAF@fix[,"INFO"]),]
  
  geno <- vcf.MAF@gt[,-1] # this gets individual ids and genotypes
  position <- getPOS(vcf.MAF) # this gets position of mutations
  chromosome <- getCHROM(vcf.MAF) # Identify chromosome ID for each mutation
  
  if (sum(duplicated(position)) != 0){
    print("This simulation needs to be checked for duplicated locus positions")
  }
  
  # convert the genome data to only have three options 0 (homoz), 1 (heteroz), or 2 (homoz alt)
  G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
  G[geno %in% c("0/0", "0|0")] <- 0
  G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
  G[geno %in% c("1/1", "1|1")] <- 2
  
  # calculate allele frequencies
  # sum across rows to get how many copies of the allele then divide by 2 times the number of columns (ind x diploid)
  a_freq <- rowSums(G)/(2*ncol(G))
  hist(a_freq) 
  
  # get individual name information from vcf file (e.g., i999) 
  vcf_ind <- data.frame(vcf_ind=colnames(vcf.MAF@gt)[-1])
  # store mutation metadata which includes mut ID (MID), selection coef (S)
  # dominance (DOM), population origin (PO), ? (GO), mutation type (MT), ? (AC), ? (DP)
  meta <- vcf.MAF@fix[,"INFO"]
  head(meta)
  length(meta)
  
  # identify mutation ids and make sure its the number of muts you expect
  LocusNameOGorder <- regmatches(meta, regexpr("[0-9]+[0-9]", meta))
  length(LocusNameOGorder)
  
  ## make a data frame where we have the locus and position in the same dataframe (useful later)
  df.OGlocusInfo <- data.frame(LocusName = LocusNameOGorder, position = position, chromosome = chromosome)
  df.ord <- df.OGlocusInfo[order(as.numeric(df.OGlocusInfo$position)),]
  
  # sort mutations and make sure vcf matches slim output
  # slim output is one base off so add a base to match vcf position
  df.mutsMAF$position_vcf <- df.mutsMAF$position + 1
  df.mutsMAFord <- df.mutsMAF[order(df.mutsMAF$position_vcf),]
  head(df.mutsMAFord)
  vcf_pos <- as.numeric(vcf.MAF@fix[,"POS"])
  head(vcf_pos)
 
  # 
  vcf_muts <- data.frame(vcf_muts=vcf.MAF@fix[,8])
  colnames(G) <- vcf_ind$vcf_ind # adds individual ids as column names
  rownames(G) <- regmatches(meta, regexpr("[0-9]+[0-9]", meta)) #ADD MUTATION NAMES TO G
  head(G[,1:5])
  dim(G)
  
  dim(vcf_ind)
  head(vcf_ind)
  head(df.indPheno)
  df.indPheno$vcf_ind <- paste0("i",0:1999) # hard coding
  # The individual IDs in Slim do not match the IDs in the VCF file. 
  # I will assume they are in the same order
  tail(df.indPheno)
  
  # Add vcf individual IDs with slim individual IDs and pop numbers 
  indPhen_df_vcf <- merge(vcf_ind, df.indPheno, by="vcf_ind")
  dim(df.indPheno)
  dim(indPhen_df_vcf)
  
  # reorder the dataframe so that it is by subpop first then ind ID
  indPhen_df_vcf <- indPhen_df_vcf[order(indPhen_df_vcf$subpop, indPhen_df_vcf$id),]
  head(indPhen_df_vcf)
  tail(indPhen_df_vcf)
  
  # split individuals by subpop
  pop1_ids <- which(indPhen_df_vcf$subpop==1)
  pop2_ids <- which(indPhen_df_vcf$subpop==2)
  
  # split G by subpop
  G_pop1 <- G[, pop1_ids]
  G_pop2 <- G[, pop2_ids]
  dim(G_pop1)
  head(G_pop1[,1:5])
  head(G_pop2[,1:5])
  
  ## clustering ##
  # first transform matrix so that columns are mutations
  fordistp1 <- as.data.frame(t(G_pop1))
  fordistp2 <- as.data.frame(t(G_pop2))
  # this calculates the distance matrix for individuals based on similar
  # mutation combinations and euclidean is the square root of sum of squares 
  # of mutational dissimilarities
  dist_matp1 <- dist(fordistp1, method="euclidean")
  dist_matp2 <- dist(fordistp2, method="euclidean")
  # This then clusters individuals into distinct groups based on genetic 
  # distances calculated previously 
  pop1_clust <- hclust(dist_matp1, method = "ward.D")
  pop2_clust <- hclust(dist_matp2, method = "ward.D")
  str(pop1_clust)
  str(pop2_clust)
  pop1_order <- pop1_clust$order
  pop2_order <- pop2_clust$order

  # Why are we doing this for all allele affect sizes? shouldn't we just be looking at m2 mutations?
  # create variables that identify which mutations are which
  whichinversionmuts <- grep("MT=3", vcf.MAF@fix[,"INFO"]) #inversions
  whichqtnmuts <- grep("MT=2", vcf.MAF@fix[,"INFO"]) #qtns
  whichneutmuts <- grep("MT=1", vcf.MAF@fix[,"INFO"]) #neut
  
  #vcf.MAF@fix[whichinversionmuts,"INFO"] # if you want info for specific mutation types
  # store info for mutation sin another variable and split data into columns
  info <- str_split(vcf.MAF@fix[,"INFO"], pattern =";", simplify=TRUE)
  head(info)
  dim(info)
  
  # find allele effect sizes
  a <- as.numeric(substring(info[,2], first=3)) 
  head(a)
  hist(a, breaks=seq(-0.01, 0.01, length.out=101))
  summary(a)
  length(a)
  dim(G)

  #G * a gives the overall effect size of the mutations on the phenotype
  G1_alpha <- G_pop1*a 
  head(G1_alpha[,1:10])
  hist(G1_alpha, breaks=seq(-0.02, 0.02, length.out=101))
  
  G2_alpha <- G_pop2*a 
  head(G2_alpha[,1:10])
  hist(G2_alpha, breaks=seq(-0.02, 0.02, length.out=101))
  
  # this gives the distribution of phenotypes in each population 
  # population 1 is evolving to an optimum of 1
  # population 2 is evolving to an optimum of -1
  hist(colSums(G1_alpha))
  hist(colSums(G2_alpha))
  
  # Sanity check - mutations in rows
  head(G[1:100,1:10])
  t(G1_alpha[1:100,1:2])
  
  hist(G1_alpha)
  hist(G2_alpha)
  dim(G1_alpha)
  
  # get position of all mutations and plot a hist of mutation locations
  # the regions with higher frequencies of mutations are potentially inverted regions
  hist(vcf_pos, breaks=seq(0,2100000, length.out=100))
  hist(df.mutsMAFord$position, breaks=seq(0,2100000, length.out=100))
  

  df.mutsMAFord$is_vcf <- NA
  
  # Check to make sure FST values match between files
  # this is slow, but correct
  G_FST <- rep(NA, nrow(G)) 
  count <- 0
  for (i in 1:nrow(df.mutsMAFord)){
    x <- grep(df.mutsMAFord$LocusName[i], vcf.MAF@fix[,"INFO"])
    if (length(x)==1){
      G_FST[x] <- df.mutsMAFord$FST[i]
    }
    count <- count+1
    print(count)
  }
  
  hist(df.mutsMAFord$FST, breaks = 50)
  head(G_FST)
  hist(G_FST, breaks = 50)
  
  sum(is.na(G_FST))
  # mutations that don't match up
  
  length(a)
  length(G_FST)
  
  head(df.mutsMAFord)
  #df.mutsMAFord<- df.mutsMAF[order(df.mutsMAF$position),]
  #G.ord <- G[order(G)]
  
#### end process VCF
######################################################################################################    

  
  
  
######################################################################################################      
#### plot heatmaps
  hist(a)
  a2 <- a
  # make an arbitary cutoff to visualize loci effect on phenotype
  a2[a>0.001] <- 1
  a2[a<0.001] <- -1
  
  G1_alpha <- G_pop1*a2*G_FST # make sure G and a line up
  G2_alpha <- G_pop2*a2*G_FST # make sure G and a line up
  
  hist(G_pop1*a2)

  pdf(paste0(folderOut, seed, "_heatmapPop1alphaFST.pdf"), height = 5, width = 7)
  heatmap(t(G1_alpha[, pop1_order]),   
          main="Pop1 G*a+-*FST",cexCol = 0.3,
          Colv = NA, useRaster=TRUE,
          scale="none",
          # Rowv = NA, 
          col=two.colors(100, start = "blue", end="red", middle="white"))
  # ADDING BREAKS SCREWS UP EVERYTHING
  #breaks=seq(-0.005, 0.005, length.out = 101))
  dev.off()
  
  pdf(paste0(folderOut, seed, "_heatmapPop2alphaFST.pdf"), height = 5, width = 7)
  
  heatmap(t(G2_alpha[, pop2_order]),   
          main="Pop2 G*a+-*FST",cexCol = 0.3,
          Colv = NA, useRaster=TRUE,
          scale="none",
          #  Rowv = NA, 
          col=two.colors(100, start = "blue", end="red", middle="white") )
  
  dev.off()
  
  G_ref1 <- G_pop1
  G_ref2 <- G_pop2
  
  af_pop1 <- rowSums(G_pop1)/(2*ncol(G_pop1))
  af_pop2 <- rowSums(G_pop2)/(2*ncol(G_pop2))
  hist(af_pop1)
  hist(af_pop2)
  todo <- which(af_pop1>af_pop2)
  
  todo_homo <- which(G[,1]==2)
  for (i in todo_homo){
    G_ref1[i,] <- abs(G_pop1[i,]-2)
    G_ref2[i,] <- abs(G_pop2[i,]-2)
  }
  
  todo_hetero <- which(G[,1]==1 & af_pop1>af_pop2)
  for (i in todo_hetero){
    G_ref1[i,] <- abs(G_pop1[i,]-2)
    G_ref2[i,] <- abs(G_pop2[i,]-2)
  }
  
  table(G_ref1[,1])
  table(G_ref1)
  table(G_ref2)

  pdf(paste0(folderOut, seed, "_heatmapPop1geno.pdf"), height = 5, width = 7)
  
  heatmap(t(G_ref1[,pop1_order]), Rowv = NA,  main="Pop1 genotypes",cexCol = 0.3,
          Colv = NA, useRaster=TRUE,
          scale="none")
  dev.off()
  
  pdf(paste0(folderOut, seed, "_heatmapPop2geno.pdf"), height = 5, width = 7)
  heatmap(t(G_ref2[,pop2_order]), Rowv = NA,  main="Pop2 genotypes",cexCol = 0.3,
          Colv = NA, useRaster=TRUE,
          scale="none")
  dev.off()
#### end plot heatmanps
######################################################################################################    

######################################################################################################
#### BigSnpr ####
## This chunk of code finds a quasi-independent set of SNPs
  ## the g matrix needs to be order by position in the genome can not reorder by any other value

# list the G matrix with the position of each mutation and the chromosome 
# it is on. This is used for
  training <- list(G = G, 
                  position = df.ord$position,
                  chromosome = df.ord$chromosome)
# confirm that these are all in the proper order - now yes!
  
# puts it in the raw format and stores likelihood genotype probability
  G_coded <- add_code256(big_copy(t(training$G),
                                  type = "raw"), 
                        code=bigsnpr:::CODE_012)

# this is doing SNP pruning - removing correlated SNPs
  newpc <- snp_autoSVD(G=G_coded, infos.chr = as.integer(training$chromosome),
                       infos.pos = training$position)
  # take snps with highest MAF and correlate snps around it
  # Snps with R^2 > 0.2 are removed
  # the subset is the indexes of the remaining SNPs
  str(newpc) #2760

# These are the indexes of the quasi-independent 
# set of loci that are kept after pruning for LD
  which_pruned = attr(newpc, which = "subset")
  length(which_pruned)

  training$G_coded <- G_coded
  training$G_pruned <- training$G[which_pruned,]
  training$which_pruned <- which_pruned

  df.mutsMAF$quasi_indep <- FALSE
  df.mutsMAF$quasi_indep[training$which_pruned] <- TRUE

#### end bigsnpr
######################################################################################################
  
######################################################################################################
#### PCADAPT
  
  # check structure of data
  training$G[1:6, 1:6]
  head(df.mutsMAFord)
  head(df.mutsMAFord[order(df.mutsMAFord$position_vcf), 1:ncol(df.mutsMAFord)])
  
  # calculate pca loadings from genotype matrix 
  # first you need to convert to a lfmm file
  gename <- paste0(seed, "_genotypes.lfmm")
  write.lfmm(t(training$G), gename)
  pcafile <- read.pcadapt(gename, type="lfmm")
  pca_all <- pcadapt(pcafile,K=2)
  head(pca_all$loadings)
  str(pca_all)
  
  plot(pca_all$scores[,1], pca_all$scores[,2], col = c(rep("red", 1000), rep("blue", 1000)))
  
  head(df.mutsMAF)
  head(pca_all)
  
  # add pca loadings to muts dataframe 
  df.mutsMAFord$pca_ALL_PC1_loadings <- pca_all$loadings[,1]
  df.mutsMAFord$pca_ALL_PC2_loadings <- pca_all$loadings[,2]
  head(df.mutsMAFord)

  plot(df.mutsMAFord$position_vcf, df.mutsMAFord$pca_ALL_PC1_loadings)

### PCA loadings for pruned data ####
  gename2 <- paste0(seed, "_genotypes_pruned.lfmm")
  write.lfmm(t(training$G_pruned), gename2)
  pcafile2 <- read.pcadapt(gename2, type="lfmm")
  pca_pruned <- pcadapt(pcafile2,K=2)
  
  plot(pca_pruned$scores[,1], pca_pruned$scores[,2], col = c(rep("red", 1000), rep("blue", 1000)))
  
  # add pca loadings to muts dataframe
  df.mutsMAFord$pca_PRUNED_PC1_loadings <- NA 
  df.mutsMAFord$pca_PRUNED_PC2_loadings <- NA
  df.mutsMAFord$pca_PRUNED_PC1_loadings[training$which_pruned] <- pca_pruned$loadings[,1]
  df.mutsMAFord$pca_PRUNED_PC2_loadings[training$which_pruned] <- pca_pruned$loadings[,2] 
  
  # look at the scree plot to see if I need to add more pops because of all the inversions
  plot(df.mutsMAFord$position_vcf, df.mutsMAFord$pca_PRUNED_PC1_loadings)
  plot(df.mutsMAFord$position_vcf, df.mutsMAFord$pca_PRUNED_PC2_loadings)
  
  # check to see if pruned dataset correlates with pca results from full dataset
  cor.test(df.mutsMAFord$pca_ALL_PC1_loadings, df.mutsMAFord$pca_PRUNED_PC1_loadings) # yes
  cor.test(df.mutsMAFord$pca_ALL_PC2_loadings, df.mutsMAFord$pca_PRUNED_PC2_loadings) # no

  ### outlier detection ###
  ## all data
  df.mutsMAFord$pcadapt_4.3.3_ALL_chisq <- as.numeric(pca_all$chi2.stat)
  df.mutsMAFord$pcadapt_4.3.3_ALL_log10p <- -log10(pca_all$pvalues)
  
  # plot to look at results
  plot(df.mutsMAFord$position_vcf, df.mutsMAFord$pcadapt_4.3.3_ALL_chisq)
  plot(df.mutsMAFord$position_vcf, df.mutsMAFord$pcadapt_4.3.3_ALL_log10p)

  ### outlier detection ### 
  ## pruned data
  test <- snp_gc(snp_pcadapt(training$G_coded, U.row = newpc$u[,1]))
  df.mutsMAFord$pcadapt_4.3.3_PRUNED_log10p <- -predict(test,log10=T)
  df.mutsMAFord$pcadapt_4.3.3_PRUNED_pvalue <- 10^(-df.mutsMAFord$pcadapt_4.3.3_PRUNED_log10p)
  df.mutsMAFord$qvalue <- qvalue(df.mutsMAFord$pcadapt_4.3.3_PRUNED_pvalue)$qvalues
  df.mutsMAFord$pcadapt_outlier <- ifelse(df.mutsMAFord$qvalue > 0.01, FALSE, TRUE)
  df.mutsMAFord$pcadapt_outlier <- as.factor(df.mutsMAFord$pcadapt_outlier)
  
  plot(newpc$u[,1], newpc$u[,2], col = c(rep("red", 1000), rep("blue", 1000)))
  
  pcadapt.log10p <- ggplot(df.mutsMAFord, aes(x = position_vcf, y = pcadapt_4.3.3_PRUNED_log10p)) +
    geom_point(aes(color = pcadapt_outlier)) +
    scale_color_manual(values = c("black",  "red")) +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92"),
          text = element_text(size = 11)) +
    labs(title = "PCAdapt",
         y = "-log10(p-value)",
         x = "Genome Position") + 
    ylim(c(0,30))+
    theme(legend.position = "none")
 
#### end PCADAPT
######################################################################################################



######################################################################################################
#### OutFLANK

FstDataFrame <- MakeDiploidFSTMat(t(training$G),rownames(training$G),
                                  c(rep("Pop1", 1000), rep("Pop2", 1000)))
colnames(FstDataFrame)[3] <- "FST_outflank"

# ask about qthreshold and trim fraction
out_pruned <- OutFLANK(FstDataFrame[training$which_pruned,], NumberOfSamples=2, 
                       LeftTrimFraction=0.05, RightTrimFraction=0.05,
                        Hmin=0.1, qthreshold=0.01)     
str(out_pruned)
head(df.mutsMAFord)
colnames(df.mutsMAFord)[10] <- "FST_slim"
df.FST.temp <- merge(df.mutsMAFord, FstDataFrame, by = "LocusName")
df.FST <- df.FST.temp[order(df.FST.temp$position_vcf),]
head(df.FST)

df.out <- pOutlierFinderChiSqNoCorr(df.FST, 
                                Fstbar = out_pruned$FSTNoCorrbar, 
                                dfInferred = out_pruned$dfInferred, Hmin=0.1)

df.out$OutFLANK_0.2_PRUNED_log10p <- -log10(df.out$pvaluesRightTail)
df.out$OutFLANK_0.2_PRUNED_log10p_add <- -log10(df.out$pvaluesRightTail + 1/10000000000000000000)

OutFLANKResultsPlotter(out_pruned, Zoom = T)

outflank.fst <- ggplot(df.out, aes(x = position_vcf, y = FST_outflank)) +
  geom_point(aes(color = OutlierFlag))+
  scale_color_manual(values = c("black", "red")) +
  theme_classic() +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92"),
        text = element_text(size = 11)) +
  labs(title ="OutFLANK",
       y = "FST",
       x = "Genome Position") + 
  ylim(c(0, max(c(df.out$FST_outflank, df.neutQTNmuts$FST)))) +
  theme(legend.position = "none")


outflank.log10p <- ggplot(df.out, aes(x = position_vcf, y = OutFLANK_0.2_PRUNED_log10p_add)) + 
  geom_point(aes(color = OutlierFlag))+
  scale_color_manual(values = c("black", "red")) +
  theme_classic() +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92"),
        text = element_text(size = 11)) +
  labs(title = "OutFLANK",
       y = "-log10(p-values)",
       x = "Genome Position") + 
  ylim(c(0, 30)) + 
  theme(legend.position = "none")

manh.plot.slim <- manh.plot + labs(title = "Simulation Output") +
  geom_point(aes(color = inOut), alpha = 0.8) + 
  scale_color_manual(values=c( "firebrick", "goldenrod", "navy")) +
  theme(legend.position = "none") + ylim(c(0, max(c(df.out$FST_outflank, df.neutQTNmuts$FST))))
pdf(paste0(folderOut, seed, "_outlierTests.pdf"), height = 7, width = 5)
ggarrange(manh.plot.slim, outflank.fst, outflank.log10p, pcadapt.log10p, nrow = 4, ncol = 1)
dev.off()

#### end OutFLANK
######################################################################################################

######################################################################################################    
## Outlier identification from genome scans

colnames(df.qtnMuts.MAF)[2] <- "LocusName"
df.outlierComp <- merge(df.qtnMuts.MAF[c(2,12:ncol(df.qtnMuts.MAF))], df.out, by = "LocusName")
head(df.outlierComp)
df.outlierComp$testNum <- NA
df.outlierComp$whichTests <- NA
for(i in 1:nrow(df.outlierComp)){
  if(!is.na(df.outlierComp$OutlierFlag[i]) & !is.na(df.outlierComp$pcadapt_outlier[i])){
    if(df.outlierComp$crit1_p.value[i] != 0 | df.outlierComp$crit2_p.value[i] != 0 & df.outlierComp$OutlierFlag[i] == FALSE & df.outlierComp$pcadapt_outlier[i] == FALSE){
      df.outlierComp$testNum[i] <- 0
      df.outlierComp$whichTests[i] <- "none"
    } else if(df.outlierComp$crit1_p.value[i] == 0 & df.outlierComp$crit2_p.value[i] == 0 & df.outlierComp$OutlierFlag[i] == TRUE & df.outlierComp$pcadapt_outlier[i] == TRUE){
      df.outlierComp$testNum[i] <- 3
      df.outlierComp$whichTests[i] <- "crit_OutFLANK_PCAdapt"
    } else if(df.outlierComp$crit1_p.value[i] == 0 & df.outlierComp$crit2_p.value[i] == 0 & df.outlierComp$OutlierFlag[i] == TRUE & df.outlierComp$pcadapt_outlier[i] == FALSE){
      df.outlierComp$testNum[i] <- 2
      df.outlierComp$whichTests[i] <- "crit_OutFLANK"
    } else if(df.outlierComp$crit1_p.value[i] == 0 & df.outlierComp$crit2_p.value[i] == 0 & df.outlierComp$OutlierFlag[i] == FALSE & df.outlierComp$pcadapt_outlier[i] == TRUE){
      df.outlierComp$testNum[i] <- 2
      df.outlierComp$whichTests[i] <- "crit_PCAdapt"
    } else if(df.outlierComp$crit1_p.value[i] != 0 | df.outlierComp$crit2_p.value[i] != 0 & df.outlierComp$OutlierFlag[i] == TRUE & df.outlierComp$pcadapt_outlier[i] == FALSE){
      df.outlierComp$testNum[i] <- 1
      df.outlierComp$whichTests[i] <- "OutFLANK"
    } else if(df.outlierComp$crit1_p.value[i] != 0 | df.outlierComp$crit2_p.value[i] != 0 & df.outlierComp$OutlierFlag[i] == FALSE & df.outlierComp$pcadapt_outlier[i] == TRUE){ 
      df.outlierComp$testNum[i] <- 1
      df.outlierComp$whichTests[i] <- "PCAdapt"
    }
  }
}
df.outlierComp$whichTests <- as.factor(df.outlierComp$whichTests)
df.outlierComp$whichTests <- relevel(df.outlierComp$whichTests, "none")
levels(df.outlierComp$inOut)
overall.outliers <- ggplot(df.outlierComp[!is.na(df.outlierComp$whichTests),], aes(x = position_vcf, y = FST_slim)) + 
  geom_point(aes(color = whichTests), alpha = 0.8) +
  facet_wrap(~inOut, nrow = 4, ncol = 1) + 
  scale_color_manual(values = c("black", "blue", "red", "goldenrod")) +
  theme_classic() +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92"),
        text = element_text(size = 11)) +
  labs(title = "Outlier Comparison",
       y = "FST",
       x = "Genome Position") + 
  ylim(c(0, 0.3))


#### end Outlier identification
######################################################################################################

######################################################################################################    
## Date to output to table
## Q1&2: VA in inversions; final LA
## Q3: average difference in # QTNs start vs. end for adaptive inversions; average difference in FST start vs end
## Q4: average difference in inversion characterists (age, length, num qtns/length) between adaptive, non-adaptive and no selection
## Q5:

#c(Va_perc_In, LA_final, )


#### end Outlier identification
######################################################################################################


######################################################################################################    
## COPY AND PASTE WHERE NEEDED
pdf(paste0(outFolder, seed, "_heatmapPop1geno.pdf"), height = 5, width = 7)

dev.off()

png(paste0("figures/", seed, "XXXX.png"), width = 480, height = 480, units = "px")

dev.off()