######################################################################################################  
#### The following code will process all simulation files and output results figures
#### Sara M. Schaal

######################################################################################################
### Load Packages and Download Data Files ###
## List Packages Needed 
packages_needed <- c("scales", "vcfR", "ggplot2", "RColorBrewer",
                      "gridExtra", "akima", "fields", "ggnewscale",
                      "ash", "plotly", "stringr", "tidyverse",  
                      "viridisLite", "ggpubr", "purrr", "dplyr")



## Install packages that aren't installed already and load them
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library( packages_needed[i], character.only = TRUE)
}

g_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
######################################################################################################  
#### Run the following code chunk once to get full data files to do further analyses on ####
df <- read.table("figures/20210514_noAdaptInv/outputAdaptInvCrit.txt")
## Inversion Through Time
df.invTime <- NULL
count <- 0
for(i in 1:nrow(df.params)){
  seed <- df.params$seed[i]
  invTimeNewFile <- read.table(paste(folder, seed, "_outputInvTime.txt", sep=""), header = TRUE, stringsAsFactors = FALSE)
  if(nrow(invTimeNewFile) > 0){
    invTimeNewFile$seed <- seed
    df.invTime <-  rbind(df.invTime, invTimeNewFile)
  }
  count <- count + 1
  print(count)
}
write.table(df.invTime, "FullSet_invTime.txt", row.names = FALSE)

## Inversion Summary Data
df.invData <- NULL
count <- 0
for(i in 1:nrow(df.params)){
  seed <- df.params$seed[i]
  invData <- read.table(paste(folder, seed, "_outputInvSumInfo.txt", sep=""), header = TRUE,
                        stringsAsFactors = FALSE)
  if(nrow(invData) > 0){
    invData$seed <- seed
    df.invData <- rbind(df.invData, invData)
  }
  count <- count + 1
  print(count)
}
write.table(df.invData, "FullSet_invData.txt", row.names = FALSE)

## QTN Inversion Through Time
df.invQTNTime <- NULL
count <- 0
for(i in 1:nrow(df.params)){
  seed <- df.params$seed[i]
  invQTNTimeNewFile <- read.table(paste(folder, seed, "_outputInvQtn.txt", sep=""), header = TRUE, stringsAsFactors = FALSE)
  if(nrow(invQTNTimeNewFile) > 0){
    invQTNTimeNewFile$seed <- seed
    df.invQTNTime <-  rbind(df.invQTNTime, invQTNTimeNewFile)
  }
  count <- count + 1
  print(count)
}
write.table(df.invQTNTime, "FullSet_invQTNTime.txt", row.names = FALSE)

## QTN Inversion Summary Data
df.invQTNData <- NULL
count <- 0
for(i in 1:nrow(df.params)){
  seed <- df.params$seed[i]
  invQTNDataNewFile <- read.table(paste(folder, seed, "_outputInvQtnSumInfo.txt", sep=""), header = TRUE, stringsAsFactors = FALSE)
  if(nrow(invQTNDataNewFile) > 0){
    invQTNDataNewFile$seed <- seed
    df.invQTNData <-  rbind(df.invQTNData, invQTNDataNewFile)
  }
  count <- count + 1
  print(count)
}
write.table(df.invQTNData, "FullSet_invQTNSumInfo.txt", row.names = FALSE)

## Pop Dynamics
df.popDyn <- NULL
count <- 0
for(i in 1:nrow(df.params)){
  seed <- df.params$seed[i]
  popDynNewFile <- read.table(paste(folder, seed, "_outputPopDynam.txt", sep=""), header = TRUE, stringsAsFactors = FALSE)
  if(nrow(popDynNewFile) > 0){
    popDynNewFile$seed <- seed
    df.popDyn <-  rbind(df.popDyn, popDynNewFile)
  }
  count <- count + 1
  print(count)
}
write.table(df.popDyn, "FullSet_popDyn.txt", row.names = FALSE)

## SimStats
df.simStats <- NULL
count <- 0
for(i in 1:nrow(df.params)){
  seed <- df.params$seed[i]
  simStatsNewFile <- read.table(paste(folder, seed, "_outputSimStats.txt", sep=""), stringsAsFactors = FALSE)
  simStatsNewFile$seed <- seed
  df.simStats <-  rbind(df.simStats, simStatsNewFile)
  count <- count + 1
  print(count)
}
df.simStats <- df.simStats[,2:ncol(df.simStats)]
colnames(df.simStats) <- c("mig1", "mig2", "pop1N", "pop2N", "mu_base", "mu_inv", "r", "alpha", "sigmaK", "burnin", "dom", "enVar", "Seed")
write.table(df.simStats, "FullSet_simStats.txt", row.names = FALSE)

## Mutations File
df.finalMuts <- NULL
count <- 0
no.Data <- NULL
for(i in 1:nrow(df.params)){
  seed <- df.params$seed[i]
  if(file.exists(paste(folder, seed, "_outputMutations.txt", sep=""))){
    finalMutsNewFile <- read.table(paste(folder, seed, "_outputMutations.txt", sep=""), stringsAsFactors = FALSE, header = TRUE)
    finalMutsNewFile$seed <- seed
    df.finalMuts <-  rbind(df.finalMuts, finalMutsNewFile)
  } else {
    no.Data <- c(no.Data, seed)
  }
  count <- count + 1
  print(count)
}
colnames(df.finalMuts)[2] <- "mut_id" 
write.table(df.finalMuts, "FullSet_finalMuts.txt", row.names = FALSE)

folderIn <- "./results/Inversion/20210719_fullSummaryData/"
df.summary <- read.table(paste0(folderIn, "outputSumData.txt"))
colnames(df.summary) <- c("seed", "Va_perc_In", "LA_final", "num_inv", "num_inv_NA", "num_inv_NS",
                          "capture_gain_p", "capture_no_gain_p", "neutral_gain_p", "neutral_no_gain_p",
                          "capture_gain_p_NA", "capture_no_gain_p_NA", "neutral_gain_p_NA", "neutral_no_gain_p_NA",
                          "capture_gain_p_NS", "capture_no_gain_p_NS", "neutral_gain_p_NS", "neutral_no_gain_p_NS",
                          "ave_start_QTNs", "ave_start_QTNs_NA", "ave_start_QTNs_NS",
                          "ave_end_QTNs", "ave_end_QTNs_NA", "ave_end_QTNs_NS", 
                          "ave_start_FST", "ave_start_FST_NA", "ave_start_FST_NS", 
                          "ave_end_FST", "ave_end_FST_NA", "ave_end_FST_NS",  
                          "ave_abV_start_qtnSelCoef", "ave_abV_start_qtnSelCoef_NA", "ave_abV_start_qtnSelCoef_NS",
                          "ave_abV_end_qtnSelCoef", "ave_abV_end_qtnSelCoef_NA", "ave_abV_end_qtnSelCoef_NS", 
                          "true_pos_pcadapt", "true_pos_outflank", "false_neg_pcadapt", "false_neg_outflank", 
                          "true_neg_pcadapt", "true_neg_outflank", "false_pos_pcadapt", "false_pos_outflank",
                          "true_neg_pcadapt_NS",  "false_pos_pcadapt_NS", 
                          "true_neg_outflank_NS", "false_pos_outflank_NS")

#df.simStats_2 <- read.table("src/invSimParams_2.txt", header = TRUE)
#df.simStats_1 <- read.table("src/invSimParams.txt", header = TRUE)
#df.simStats <- rbind(df.simStats_1, df.simStats_2)
df.simStats <- read.table("src/invSimParams.txt", header = TRUE)
df.invCharFinalGen <- read.table(paste0(folderIn, "outputInvChar_finalGen.txt"))
colnames(df.invCharFinalGen) <- c("seed", "inv_id", "inv_age", "inv_length", "num_qtns_Lscaled", "adaptInv")
df.invCharAllData <- read.table(paste0(folderIn, "outputInvChar_allData.txt"), fill = TRUE)

#### END Run the following code chunk once to get full data files to do further analyses on ####
######################################################################################################  

## Sanity check this should be 0
dim(df.simStats[!df.simStats$seed %in% df.summary$seed,])

######################################################################################################  
#### Unique Parameters ####
# Subset original parameters dataframe to identify all unique parameter combinations. 
# Then using those parameter values identify the simulation seeds of replicate simulations 
# for each set of parameters.
unique.params <- unique(df.simStats[c("chromNum", "mig1", "mig2", "n", "muBase", "muInv", 
                                      "r", "alpha", "sigmaK", "theta1", "theta2", "burnin", "dom", "enVar")])

unique.params$params <- NULL
for(i in 1:nrow(unique.params)){
  unique.params$params[i] <- paste(unique.params[i,1:14], collapse = " ")
}

df.simStats$params <- NULL
for(i in 1:nrow(df.simStats)){
  df.simStats$params[i] <- paste(df.simStats[i,c(14, 9, 10, 7, 2, 3, 15, 5, 4, 11, 12, 13, 16, 8)], 
                                 collapse = " ")
}

reps <- 5
for(i in 1:reps){
  seedColName <- paste("Seed", i, sep = "")
  unique.params[,seedColName] <- NA
}

for(i in 1:nrow(unique.params)){
 seeds <- df.simStats[unique.params$params[i] == df.simStats[,17], ]$seed
 unique.params[i, 16:20] <- seeds
}

######################################################################################################    



######################################################################################################    
## Calculate average LA across rep sims 
### Average Replicate Simulations - Pop dynamics file ONLY NEED THIS CHUNK ONCE
# Create a dataframe that has the average of all the replicate parameter values. 
# Then bind the corresponding parameters to the dataframe for ease of plotting and analysis later. 


########################
### average for loop ###
###################################################################

# df.LA.average <- NULL
# df.LA.stderr <- NULL
# df.LA.all <- NULL
# count <- 0
# # step through each line of the unique parameters dataframe
# for(i in 1:nrow(unique.params)){
#   ptm <- proc.time()
#   # create empty variable for putting each data set that should be average
#   df.average <- NULL
#   av.seeds <- NULL
#   # step through the population dynamics dataframe
#   for(j in 1:nrow(df.popDyn)){
#     # step through the different seeds that have the unique parameters
#     for(k in 1:reps){
#       seedCol <- paste("Seed", k, sep="")
#       # if that seed is not an NA then do the next step
#       if(!is.na(unique.params[i, seedCol])){
#         # if the seed from unique parameters matches the seed of the pop dynam dataset store it
#         if(df.popDyn$seed[j] == unique.params[i, seedCol]){
#           # bind data from same parameter combination together for averaging later
#           df.average <- rbind(df.average, df.popDyn[j,])
#           # store unique seeds to this dataset
#           av.seeds <- unique(c(av.seeds, unique.params[i, seedCol]))
#         }
#       }
#     }
#   }
#   
#   
#   ## average columns of interest
#   # df.averageSubset <- df.average[,vect.aggCols]
#   new.average <- aggregate(.~sim_gen, data = df.average, FUN = mean)
#   SE <- function(x){
#     sd(x)/sqrt(length(x))
#   }
#   new.se <- aggregate(.~sim_gen, data = df.average, FUN = SE)
#   df.simParams <- data.frame(mu_inv = rep(unique.params$mu_inv[i], nrow(new.average)), 
#                              mig = rep(unique.params$mig1[i], nrow(new.average)),
#                              alpha = rep(unique.params$alpha[i], nrow(new.average)), 
#                              sigmaK = rep(unique.params$sigmaK[i], nrow(new.average)), 
#                              enVar = rep(unique.params$enVar[i], nrow(new.average)), 
#                              mu_base = rep(unique.params$mu_base[i], nrow(new.average)))
#   
#   df.simParamsAll <- data.frame(mu_inv = rep(unique.params$mu_inv[i], nrow(df.average)), 
#                                 mig = rep(unique.params$mig1[i], nrow(df.average)),
#                                 alpha = rep(unique.params$alpha[i], nrow(df.average)), 
#                                 sigmaK = rep(unique.params$sigmaK[i], nrow(df.average)), 
#                                 enVar = rep(unique.params$enVar[i], nrow(df.average)), 
#                                 mu_base = rep(unique.params$mu_base[i], nrow(df.average)))
#   
#   # create seed columns so we can keep track of relevant seed names
#   df.Seedcolumns <- NULL
#   df.SeedcolumnsAll <- NULL
#   vect.colNames <- NULL
#   for(m in 1:reps){
#     seedCol <- paste("Seed", m, sep = "")
#     if(!is.na(av.seeds[m])){
#       df.Seedcolumns <- cbind(df.Seedcolumns, rep(av.seeds[m], nrow(new.average)))
#       df.SeedcolumnsAll <- cbind(df.SeedcolumnsAll, rep(av.seeds[m], nrow(df.average)))
#     } else {
#       df.Seedcolumns <- cbind(df.Seedcolumns, rep(NA, nrow(new.average)))
#       df.SeedcolumnsAll <- cbind(df.SeedcolumnsAll, rep(NA, nrow(df.average)))
#     }
#     vect.colNames <- c(vect.colNames, seedCol)
#   }
#   colnames(df.Seedcolumns) <- vect.colNames
#   colnames(df.SeedcolumnsAll) <- vect.colNames
#   
#   # finally bind together all the data in the final dataframe
#   df.LA.average <- rbind(df.LA.average, cbind(new.average, df.simParams, df.Seedcolumns))
#   df.LA.stderr <- rbind(df.LA.stderr, cbind(new.se, df.simParams, df.Seedcolumns))
#   df.LA.all <- rbind(df.LA.all, cbind(df.average, df.simParamsAll, df.SeedcolumnsAll))
#   count <- count + 1
#   print(count)
#   print(proc.time() - ptm)
# } # close average for loop
# ###################################################################
# df.LA.average <- subset(df.LA.average, select = -seed)
# df.LA.stderr <- subset(df.LA.stderr, select = -seed)
# df.LA.all <- subset(df.LA.all, select = -seed)
# 
# #df.LA.average <- read.table("FullSet_LAaverage.txt", header = TRUE)
# write.table(df.LA.average, "FullSet_LAaverage.txt", row.names = FALSE)
# write.table(df.LA.stderr, "FullSet_LAstderr.txt", row.names = FALSE)
# write.table(df.LA.all, "FullSet_LAallData.txt", row.names = FALSE)
# 
# # this should work I don't understand why it doesn't
# #df.LA.average[,df.LA.average] <- lapply(df.LA.average[,LA.cols], function(x) as.factor(as.character(x)))
# df.LA.average <- subset(df.LA.average, select = -seed)
# 
# # convert every column of parameters data to factor 
# for(i in 16:ncol(df.LA.average)){
#   df.LA.average[,i] <- as.factor(as.character(df.LA.average[,i]))
# }
# 
# # Run this chunk if average file has already been made
# df.LA.average <- read.table("src/results/FullSet_LAaverage.txt", header = TRUE, stringsAsFactors = FALSE)
# 
# for(i in 16:ncol(df.LA.average)){
#   df.LA.average[,i] <- as.factor(as.character(df.LA.average[,i]))
# }

#df.sumInfo <- read.table("./results/Inversion/20210525_fulldata/finalData/outputSumData.txt")

#df.sumInfoP <- left_join(df.sumInfo, df.params, by = "seed")

#### end LA analysis
######################################################################################################


######################################################################################################
#### Q1 & 2: Plot Percent Additive Genetic Variance in Adaptive Inversions

## Join summary data together with parameters
df.sumData <- left_join(df.summary, df.simStats[,1:(ncol(df.simStats)-1)], by = "seed")

## Separate data into three different inversion mutation rates
df.muInv0 <- df.sumData[df.sumData$muInv == 0,]
df.muInv3 <- df.sumData[df.sumData$muInv == 1e-03,]
df.muInv6 <- df.sumData[df.sumData$muInv == 1e-06,]

## Join Inv mu 0 and with both of the other Inv mutation rates to calculate local adaptation difference
df.muInv6_0 <- full_join(df.muInv6[,c(2:3,49:53, 55:57)], df.muInv0[,c(3,49:53, 55:57)], by = c("muBase", "sigmaK", "alpha", "enVar", "mig1", "mig2", "rep"))
colnames(df.muInv6_0)[c(1,2,4,11:12)] <- c("VA_perc_In_6", "LA_final_6", "muInv_6", "LA_final_0", "muInv_0")
df.muInv3_0 <- full_join(df.muInv3[,c(2:3,49:53, 55:57)], df.muInv0[,c(3,49:53, 55:57)], by = c("muBase", "sigmaK", "alpha", "enVar", "mig1", "mig2", "rep"))
colnames(df.muInv3_0)[c(1,2,4,11:12)] <- c("VA_perc_In_3", "LA_final_3", "muInv_3", "LA_final_0", "muInv_0")
df.muInv3_0$LA_diff <- df.muInv3_0$LA_final_3 - df.muInv3_0$LA_final_0
df.muInv6_0$LA_diff <- df.muInv6_0$LA_final_6 - df.muInv6_0$LA_final_0

## Summarize the additive genetic variation and Local adaptation difference across replicates 
# muInv = 1e-03
df.muInv3_0_av <- aggregate(cbind(VA_perc_In_3, LA_diff, LA_final_3, LA_final_0)~muBase + sigmaK + alpha + enVar + mig1 + mig2, 
                            FUN=mean, data = df.muInv3_0)
df.muInv3_0_sd <- aggregate(cbind(VA_perc_In_3, LA_diff, LA_final_3, LA_final_0)~muBase + sigmaK + alpha + enVar + mig1 + mig2, 
                            FUN=sd, data = df.muInv3_0)
df.muInv3_0_av$LA_diff_sd <-  df.muInv3_0_sd$LA_diff
df.muInv3_0_av$LA_diff_upSD <- df.muInv3_0_av$LA_diff + df.muInv3_0_av$LA_diff_sd
df.muInv3_0_av$LA_diff_lowSD <- df.muInv3_0_av$LA_diff - df.muInv3_0_av$LA_diff_sd
df.muInv3_0_av$LA_final_3_sd <-  df.muInv3_0_sd$LA_final_3
df.muInv3_0_av$LA_final_3_upSD <- df.muInv3_0_av$LA_final_3 + df.muInv3_0_av$LA_final_3_sd
df.muInv3_0_av$LA_final_3_lowSD <- df.muInv3_0_av$LA_final_3 - df.muInv3_0_av$LA_final_3_sd
df.muInv3_0_av$LA_final_0_sd <-  df.muInv3_0_sd$LA_final_0
df.muInv3_0_av$LA_final_0_upSD <- df.muInv3_0_av$LA_final_0 + df.muInv3_0_av$LA_final_0_sd
df.muInv3_0_av$LA_final_0_lowSD <- df.muInv3_0_av$LA_final_0 - df.muInv3_0_av$LA_final_0_sd
df.muInv3_0_av$VA_perc_sd <-  df.muInv3_0_sd$VA_perc_In_3
df.muInv3_0_av$VA_perc_upSD <- df.muInv3_0_av$VA_perc_In_3 + df.muInv3_0_av$VA_perc_sd
df.muInv3_0_av$VA_perc_lowSD <- df.muInv3_0_av$VA_perc_In_3 - df.muInv3_0_av$VA_perc_sd

## Plot local adaptation difference
df.muInv3_0$muBase <- as.factor(df.muInv3_0$muBase)
df.muInv3_0$mig1 <- as.factor(df.muInv3_0$mig1)
df.muInv3_0$sigmaK <- as.factor(df.muInv3_0$sigmaK)
df.muInv6_0$muBase <- as.factor(df.muInv6_0$muBase)
df.muInv6_0$mig1 <- as.factor(df.muInv6_0$mig1)
df.muInv6_0$sigmaK <- as.factor(df.muInv6_0$sigmaK)
sigmaK.labels <- c("0.75" = "Strong Selection", "1.5" = "Moderate Selection", "3" = "Weak Selection")

# make all subsetting columns factors
for(i in 1:6){
  df.muInv3_0_av[,i] <- as.factor(df.muInv3_0_av[,i])
}

# recode necessary factor levels for plotting 
df.muInv3_0_av$muBase <- recode_factor(df.muInv3_0_av$muBase,'0.0000001' = '0.0002', '0.000001' = '0.002')
df.muInv3_0_av$muBase <- recode_factor(df.muInv3_0_av$muBase,'1e-07' = '0.0002', '1e-06' = '0.002')
df.muInv3_0_av$sigmaK <- factor(df.muInv3_0_av$sigmaK, c(3, 1.5, 0.75))
  
#### Local Adaptation Plots ####
# Difference in local adaptation
# ploygenic architecture
plot.LA_diff_inv3_pgen <- ggplot(df.muInv3_0_av[df.muInv3_0_av$enVar == 0 & 
                                              df.muInv3_0_av$alpha == 0.002 &
                                              df.muInv3_0_av$muBase == 0.002,], 
                            aes(x = mig1, y = LA_diff, group = sigmaK)) + 
  geom_errorbar(aes(ymin=LA_diff_lowSD, ymax=LA_diff_upSD), size =  0.2, width=.5) +
  geom_point(aes(color = viridis(4)[3]), shape = 19, size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  labs(y = expression("LA"[Inv]*" - LA"[noInv]),
       x = " ") +
  guides(color = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold"))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  scale_color_manual(values=viridis(4)[2]) +
  ylim(-0.04, 0.6)

# oligogenic architecture
plot.LA_diff_inv3_ogen <- ggplot(df.muInv3_0_av[df.muInv3_0_av$enVar == 0 & 
                                                df.muInv3_0_av$alpha == 0.2 &
                                                df.muInv3_0_av$muBase == "0.0002",], 
                               aes(x = mig1, y = LA_diff, group = sigmaK)) + 
  geom_errorbar(aes(ymin=LA_diff_lowSD, ymax=LA_diff_upSD), size =  0.2, width=.5) +
  geom_point(aes(color = viridis(4)[4]), shape = 17, size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  guides(color = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)")) +
  labs(x = " ", y = " ") + 
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold"))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  scale_color_manual(values=viridis(4)[4]) +
  ylim(-0.04, 0.6)



#### Amount of Additive Genetic Variance Plots ####
df.invGenome <- read.table(paste0(folderIn, "outputInvGenome_allData.txt"))
head(df.invGenome)
tail(df.invGenome)
colnames(df.invGenome) <- c("seed", "uniqueBases", "numOverlap", "percGenome")
df.invGenomeParam <- full_join(df.invGenome, df.simStats, by = "seed")
df.invGenome_av <- aggregate(percGenome~muBase + sigmaK + muInv + alpha + enVar + mig1 + mig2, 
                             FUN=mean, data = df.invGenomeParam)
df.invGenome_sd <- aggregate(percGenome~muBase + sigmaK + muInv + alpha + enVar + mig1 + mig2, 
                             FUN=sd, data = df.invGenomeParam)
df.invGenome_av$percGenome_sd <- df.invGenome_sd$percGenome
df.invGenome_av$percGenome_lowSD <- df.invGenome_av$percGenome - df.invGenome_av$percGenome_sd
df.invGenome_av$percGenome_upSD <- df.invGenome_av$percGenome + df.invGenome_av$percGenome_sd

for(i in 1:7 ){
  df.invGenome_av[,i] <- as.factor(df.invGenome_av[,i])
}

#df.invGenome_av$muBase <- recode_factor(df.invGenome_av$muBase, '0.0000001' = '0.0002', '0.000001' = '0.002')
df.invGenome_av$muBase <- recode_factor(df.invGenome_av$muBase, '1e-07' = '0.0002', '1e-06' = '0.002')

df.invGenome_av$sigmaK <- factor(df.invGenome_av$sigmaK, c(3, 1.5, 0.75))
df.VA <- left_join(df.muInv3_0_av, df.invGenome_av[df.invGenome_av$muInv == 0.001,], by = c("muBase", "sigmaK", "alpha", 
                                                  "enVar", "mig1", "mig2"))
df.VA$mig2 <- as.numeric(as.character(df.VA$mig2))

# polygenic architecture
plot.VA_in_inv3_pgen <- ggplot(df.VA[df.VA$enVar == 0 & 
                                df.VA$alpha == 0.002 &
                                df.VA$muBase == 0.002,], 
                          aes(x = mig1, y = VA_perc_In_3, group = sigmaK)) + 
  geom_errorbar(aes(ymin=VA_perc_lowSD, ymax=VA_perc_upSD), size =  0.2, width=.5) +
  geom_ribbon(aes(ymin= percGenome_lowSD, ymax= percGenome_upSD), fill = "navy", alpha=0.2) +
  geom_point(aes(color = muBase), shape = 19, size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) )+ 
  labs(y = expression("%VA"[inv]),
       x = " ") +
  guides(color = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold"))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  scale_color_manual(values=viridis(4)[2]) + 
  ylim(c(0, 75))

# oligogenic architecture
plot.VA_in_inv3_ogen <- ggplot(df.VA[df.VA$enVar == 0 & 
                                             df.VA$alpha == 0.2 &
                                             df.VA$muBase == "0.0002",], 
                          aes(x = mig1, y = VA_perc_In_3, group = sigmaK)) + 
  geom_errorbar(aes(ymin=VA_perc_lowSD, ymax=VA_perc_upSD), size = 0.2, width=.5) +
  geom_ribbon(aes(ymin= percGenome_lowSD, ymax= percGenome_upSD), fill = "goldenrod", alpha=0.2) +
  geom_point(aes(color = muBase, shape = muBase), size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) )+ 
  labs(y = " ",
       x = " ") +
  guides(color = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold"))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  scale_color_manual(values=viridis(4)[4]) + 
  ylim(c(0, 75))
         

## total LA
# polygenic architecture
plot.totalLA_pgen <- ggplot(df.muInv3_0_av[df.muInv3_0_av$enVar == 0 & 
                                           df.muInv3_0_av$alpha == 0.002 &
                                           df.muInv3_0_av$muBase == 0.002,], 
                            aes(x = mig1, y = LA_final_3, group = sigmaK)) + 
  geom_errorbar(aes(ymin=LA_final_3_lowSD, ymax=LA_final_3_upSD), size = 0.2, width=.5) +
  geom_ribbon(aes(ymin= LA_final_0_lowSD, ymax= LA_final_0_upSD), fill = "navy", alpha=0.2) +
  geom_point(aes(color = muBase, shape = muBase), size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) )+ 
  labs(title = "Polygenic Architecture",
       y = expression("LA"[inv]),
       x = " ") +
  guides(color = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold")) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  scale_color_manual(values=viridis(4)[2]) + 
  ylim(c(-0.01,1))


plot.totalLA_ogen <- ggplot(df.muInv3_0_av[df.muInv3_0_av$enVar == 0 & 
                                           df.muInv3_0_av$alpha == 0.2 &
                                           df.muInv3_0_av$muBase == "0.0002",], 
                            aes(x = mig1, y = LA_final_3, group = sigmaK)) + 
  geom_errorbar(aes(ymin=LA_final_3_lowSD, ymax=LA_final_3_upSD), size = 0.2, width=.5) +
  geom_ribbon(aes(ymin= LA_final_0_lowSD, ymax= LA_final_0_upSD), fill = "goldenrod", alpha=0.2) +
  geom_point(aes(color = muBase, shape = muBase), size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) )+ 
  labs(title = "Oligogenic Architecture",
       y = " ",
       x = " ") +
  guides(color = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=18,)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  scale_color_manual(values=viridis(4)[4]) + 
  ylim(c(-0.01,1))

#### Number of Inversions Plots ####
numInv_av <- aggregate(num_inv~muBase + sigmaK + alpha + enVar + mig1 + mig2, 
                       FUN=mean, data = df.muInv3)
numInv_sd <- aggregate(num_inv~muBase + sigmaK + alpha + enVar + mig1 + mig2, 
                       FUN=sd, data = df.muInv3)
#numInv6_av <- aggregate(num_inv~muBase + sigmaK + alpha + enVar + mig1 + mig2, 
#                       FUN=mean, data = df.muInv6)
#numInv6_sd <- aggregate(num_inv~muBase + sigmaK + alpha + enVar + mig1 + mig2, 
#                       FUN=sd, data = df.muInv6)

numInv_av$numInv_sd <- numInv_sd$num_inv
numInv_av$numInv_lowSD <- numInv_av$num_inv - numInv_av$numInv_sd
numInv_av$numInv_upSD <- numInv_av$num_inv + numInv_av$numInv_sd

for(i in 1:6 ){
  numInv_av[,i] <- as.factor(numInv_av[,i])
}

#numInv_av$muBase <- recode_factor(numInv_av$muBase, '0.0000001' = '0.0002', '0.000001' = '0.002')
numInv_av$muBase <- recode_factor(numInv_av$muBase,  '1e-07' = '0.0002', '1e-06' = '0.002')

numInv_av$sigmaK <- factor(numInv_av$sigmaK, c(3, 1.5, 0.75))
plot.numInv.pgen <- ggplot(numInv_av[numInv_av$enVar == 0 & numInv_av$alpha == 0.002 & 
                                  numInv_av$muBase == 0.002 & numInv_av$num_inv != 0,], 
                      aes(x = mig1, y = num_inv, group = sigmaK)) + 
  geom_errorbar(aes(ymin=numInv_lowSD, ymax=numInv_upSD), width=.5) +
  geom_point(color = viridis(4)[2], size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) )+ 
  labs(y = expression(bar(N)[inv]),
       x = "Migration Rate") +
  guides(color = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold"))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  scale_color_manual(values=viridis(4)[2]) + ylim(-0.5, 15)


plot.numInv.ogen <- ggplot(numInv_av[numInv_av$enVar == 0 & numInv_av$alpha == 0.2
                                           & numInv_av$muBase == "0.0002" & numInv_av$num_inv != 0,], 
                                 aes(x = mig1, y = num_inv, group = sigmaK)) + 
  geom_errorbar(aes(ymin=numInv_lowSD, ymax=numInv_upSD), width=.5) +
  geom_point(color = viridis(4)[4], shape = 17, size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) ) + 
  labs(y = " ",
       x = "Migration Rate") +
  guides(color = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold"))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  scale_color_manual(values=viridis(4)[4]) + ylim(-0.5, 15)


folderOutFig <- "./figures/ManuscriptFigs/"
pdf(file = paste0(folderOutFig, "fig1_LAplots.pdf"), width = 15, height = 15)
ggarrange(plot.totalLA_pgen, plot.totalLA_ogen, plot.LA_diff_inv3_pgen, plot.LA_diff_inv3_ogen, 
          plot.VA_in_inv3_pgen, plot.VA_in_inv3_ogen, plot.numInv.pgen, plot.numInv.ogen,
            nrow = 4, ncol = 2, labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
dev.off()


#ggarrange(plot.numInv, plot.numInv_largeAlpha, labels = c("A",  "B"))


# # I could think of main figure that is a 5 row plot, with migration rate on the x-axis 
# (and other subsets of sims that we decide on), showing on the y axis in each row: 
# (i) the amount of local adaptation relative to no inversions, 
# (ii) the percent of additive genetic variance in adaptive inversions, 
# (iii) and the number of adaptive inversions, 
# (iv) the average age of inversions in different categories, 
# (v) the length of inversions in different categories.
  

## Summarize the additive genetic variation and Local adaptation difference across replicates 
# muInv = 1e-06
 
#### end Q1 & 2 plotting
######################################################################################################

######################################################################################################
#### Q3: Plot Proportion of Adaptive inversions in each evo history category

head(df.muInv3)
dim(df.muInv3)
df.muInv3$capture_gain_count <- df.muInv3$num_inv*df.muInv3$capture_gain_p
df.muInv3$capture_no_gain_count <- df.muInv3$num_inv*df.muInv3$capture_no_gain_p
df.muInv3$neutral_gain_count <- df.muInv3$num_inv*df.muInv3$neutral_gain_p
df.muInv3$neutral_no_gain_count <- df.muInv3$num_inv*df.muInv3$neutral_no_gain_p

df.ev.hist.av <- aggregate(cbind(capture_gain_count, capture_no_gain_count, 
                neutral_gain_count, neutral_no_gain_count)~muBase + sigmaK + alpha + enVar + mig1 + mig2, 
          FUN=mean, data = df.muInv3)
df.ev.hist.long <- as.data.frame(pivot_longer(df.ev.hist.av, cols = c(capture_gain_count, capture_no_gain_count, 
                                     neutral_gain_count, neutral_no_gain_count),
             names_to = "evHist", values_to = "count"))
for(i in 1:7){
  df.ev.hist.long[,i] <- as.factor(df.ev.hist.long[,i])
}

df.ev.hist.long$evHist <- recode_factor(df.ev.hist.long$evHist, 'capture_gain_p' = 'Capture & Gain', 
                                        'capture_no_gain_p' = 'Capture & No Gain', 
                                        'neutral_gain_p' = 'Neutral & Gain', 
                                        'neutral_no_gain_p'  = 'Neutral & No Gain')
plot.evHist_highQTN <- ggplot(df.ev.hist.long[df.ev.hist.long$enVar == 0 & df.ev.hist.long$alpha == 0.002 & df.ev.hist.long$muBase == 1e-06,], 
                      aes(x = mig1, y = prop, fill = evHist, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(5)[4:1]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "High QTN Mutation Rate", y = "Proportion", x = "Migration Rate")

plot.evoHist_lowQTN <- ggplot(df.ev.hist.long[df.ev.hist.long$enVar == 0 & df.ev.hist.long$alpha == 0.002 & df.ev.hist.long$muBase == 1e-07,], 
                      aes(x = mig1, y = prop, fill = evHist, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(5)[4:1]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "Low QTN Mutation Rate", y = " ", x = "Migration Rate")

plot.evHist_highQTN_largeAlpha <- ggplot(df.ev.hist.long[df.ev.hist.long$enVar == 0 & df.ev.hist.long$alpha == 0.2 & df.ev.hist.long$muBase == 1e-06,], 
                              aes(x = mig1, y = prop, fill = evHist, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(5)[4:1]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "High QTN Mutation Rate", y = "Proportion", x = "Migration Rate")

plot.evoHist_lowQTN_largeAlpha <- ggplot(df.ev.hist.long[df.ev.hist.long$enVar == 0 & df.ev.hist.long$alpha == 0.2 & df.ev.hist.long$muBase == 1e-07,], 
                              aes(x = mig1, y = prop, fill = evHist, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(5)[4:1]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "Low QTN Mutation Rate", y = " ", x = "Migration Rate")


evohist.leg <- g_legend(plot.evoHist_lowQTN)
plot.evoHist_lowQTN.noLeg <- plot.evoHist_lowQTN + theme(legend.position = "none")
ggarrange(plot.evHist_highQTN, plot.evoHist_lowQTN.noLeg, evohist.leg, nrow = 1, ncol = 3,
          widths = c(2.3,2.3,0.8))

evohist_largeAlpha.leg <- g_legend(plot.evoHist_lowQTN_largeAlpha)
plot.evoHist_lowQTN_largeAlpha.noLeg <- plot.evoHist_lowQTN_largeAlpha + theme(legend.position = "none")
ggarrange(plot.evHist_highQTN_largeAlpha, plot.evoHist_lowQTN_largeAlpha.noLeg, evohist_largeAlpha.leg, nrow = 1, ncol = 3,
          widths = c(2.3,2.3,0.8))

#### end Q3 plotting
######################################################################################################

######################################################################################################
#### Q4: Plot Characteristics

head(df.invCharFinalGen)
for(i in 2:ncol(df.simStats)){
  df.simStats[,i] <- as.factor(as.character(df.simStats[,i]))
}
df.invCharParam <- left_join(df.invCharFinalGen, df.simStats, by = "seed")
df.invChar.muInv3 <- df.invCharParam[df.invCharParam$muInv == 0.001,]
options(scipen = 999)
df.invChar.muInv3$muBase <- recode_factor(df.invChar.muInv3$muBase, '1e-09' = '0.000002', '1e-08' = '0.00002', '1e-07' = '0.0002', '1e-06' = '0.002')
df.invChar.muInv3$sigmaK <- factor(df.invChar.muInv3$sigmaK, c(3, 1.5, 0.75))

df.invChar.muInv3[df.invChar.muInv3$enVar == 0 & df.invChar.muInv3$muBase == 0.002 & df.invChar.muInv3$alpha == 0.002 & df.invChar.muInv3$mig1 == 0.001 
                  & df.invChar.muInv3$sigmaK == 0.75,]


color_scale <- c(inferno(4)[3], "darkgrey", inferno(4)[1])
df.invChar.muInv3$adaptInv <- factor(df.invChar.muInv3$adaptInv, levels = c("Adaptive", "Nonadaptive", "No selection"))

## subset to get data for only parameter sets that contain all 3 factor levels
plot.output <- NULL
for(i in 1:length(unique(df.invChar.muInv3$params))){
  df <- df.invChar.muInv3[df.invChar.muInv3$params==unique(df.invChar.muInv3$params)[i],]
  if(length(unique(df$adaptInv)) == 3){
    plot.output <- rbind(plot.output, df)
  }
}

plot.age.pgen <- ggplot(data = plot.output[plot.output$enVar == 0 & 
                                                      plot.output$muBase == 0.002 & 
                                                      plot.output$alpha == 0.002,], 
       aes(x = mig1, y= inv_age, fill = adaptInv, color = adaptInv)) +
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values = c("orange", "black", "grey21"))+
  scale_fill_manual(values = color_scale) + 
  guides(fill = guide_legend(title = "Inversion Status"),
         color = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "Polygenic Architecture", y = expression("Average Age"[inv]*" (gen)"), x = " ")

plot.age.ogen <- ggplot(data = plot.output[plot.output$enVar == 0 &
                                                   plot.output$muBase == "0.0002" &
                                                   plot.output$alpha == 0.2,], 
       aes(x = mig1, y= inv_age, fill = adaptInv, color = adaptInv)) +
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values = c("orange", "black", "grey21"))+
  scale_fill_manual(values = color_scale) + 
  guides(fill = guide_legend(title = "Inversion Status"),
         color = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "Oligogenic Architecture", y = " ", x = " ") +
  scale_x_discrete(" ", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE)


plot.length.pgen <- ggplot(data = plot.output[plot.output$enVar == 0 & 
                                                         plot.output$muBase == 0.002 & 
                                                         plot.output$alpha == 0.002,], 
       aes(x = mig1, y= inv_length, fill = adaptInv, color = adaptInv)) +
  facet_wrap(~sigmaK,  labeller = labeller(sigmaK = sigmaK.labels)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = color_scale) + 
  scale_color_manual(values = c("orange", "black", "grey21"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = guide_legend(title = "Inversion Status"),
         color = guide_legend(title = "Inversion Status")) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = expression("Average Length"[inv]*" (bp)"), x = " ")

plot.length.ogen <- ggplot(data = plot.output[plot.output$enVar == 0 & 
                                                      plot.output$muBase == "0.0002" & 
                                                      plot.output$alpha == 0.2,], 
       aes(x = mig1, y= inv_length, fill = adaptInv, color = adaptInv)) +
  facet_wrap(~sigmaK,  labeller = labeller(sigmaK = sigmaK.labels)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = color_scale) + 
  scale_color_manual(values = c("orange", "black", "grey21"))+
  theme_classic() +
  theme(legend.position = "none") +
  guides(fill = guide_legend(title = "Inversion Status"),
         color = guide_legend(title = "Inversion Status")) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = " ", x = " ") +
  scale_x_discrete("", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE)


plot.numQTN.pgen <- ggplot(data = plot.output[plot.output$enVar == 0 & 
                                                         plot.output$muBase == 0.002 & 
                                                         plot.output$alpha == 0.002,], 
       aes(x = mig1, y= num_qtns_Lscaled, fill = adaptInv, color = adaptInv)) +
  facet_wrap(~sigmaK,  labeller = labeller(sigmaK = sigmaK.labels)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = color_scale) + 
  scale_color_manual(values = c("orange", "black", "grey21")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = guide_legend(title = "Inversion Status"),
         color = guide_legend(title = "Inversion Status")) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = expression(bar(N)[QTNs] / "(Length"[inv]*")"), x = "Migration Rate") + 
  ylim(c(0,0.012))


plot.numQTN.ogen <- ggplot(data = plot.output[plot.output$enVar == 0  &
                                                      plot.output$muBase == "0.0002" &
                                                      plot.output$alpha == 0.2,], 
                              aes(x = mig1, y= num_qtns_Lscaled, fill = adaptInv, color = adaptInv)) +
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = color_scale) + 
  scale_color_manual(values = c("orange", "black", "grey21")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = guide_legend(title = "Inversion Status"),
         color = guide_legend(title = "Inversion Status")) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = " ") + 
  ylim(c(0,0.012)) +
  scale_x_discrete("Migration Rate", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE)


## blank graph
invChar.legend <- g_legend(plot.length.pgen)
plot.length.pgen.noLeg <- plot.length.pgen + theme(legend.position = "none")
blank <- ggplot() + theme_void()
ggarrange(plot.age.pgen, plot.age.ogen, blank, plot.length.pgen.noLeg, plot.length.ogen, 
          invChar.legend, plot.numQTN.pgen, plot.numQTN.ogen, blank, ncol = 3, nrow = 3,
          widths = c(2.3, 2.3, 0.8, 2.3, 2.3, 0.8, 2.3, 2.3, 0.8), 
          labels = c("A", "B", "", "C", "D", "", "E", "F"))
# pdf 15 wide by 12 high 

### end Q4 plotting
######################################################################################################

######################################################################################################
#### characteristics a different way

color_scale <- inferno(4)[3:1]
df.invChar.muInv3$adaptInv <- factor(df.invChar.muInv3$adaptInv, levels = c("Adaptive", "Nonadaptive", "No selection"))

## Get averages and SD for all characteristics across all replicates: mean(c(rep1, rep2, .. rep5))
df.invChar.muInv3.av <- aggregate(cbind(inv_age, inv_length, num_qtns_Lscaled)~adaptInv + muBase + muInv + sigmaK + alpha + enVar + mig1 + mig2, data = df.invChar.muInv3, FUN = mean)
colnames(df.invChar.muInv3.av)[9:11] <- c("inv_age_av", "inv_length_av", "num_qtns_Lscaled_av")
df.invChar.muInv3.sd <- aggregate(cbind(inv_age, inv_length, num_qtns_Lscaled)~adaptInv + muBase + muInv + sigmaK + alpha + enVar + mig1 + mig2, data = df.invChar.muInv3, FUN = sd)
colnames(df.invChar.muInv3.sd)[9:11] <- c("inv_age_sd", "inv_length_sd", "num_qtns_Lscaled_sd")
df.invChar.av.sd <- cbind(df.invChar.muInv3.av, df.invChar.muInv3.sd[9:11])
head(df.invChar.av.sd)

## Convert the adaptInv column for all the averages and SD to wide format for plotting
df.nonadapt <- df.invChar.av.sd[df.invChar.av.sd$adaptInv == "Nonadaptive",]
df.noselect <- df.invChar.av.sd[df.invChar.av.sd$adaptInv == "No selection",]
df.adaptive <- df.invChar.av.sd[df.invChar.av.sd$adaptInv == "Adaptive",]
colnames(df.adaptive)[9:14] <- c("inv_age_av_A", "inv_length_av_A", "num_qtns_Lscaled_av_A", "inv_age_sd_A", "inv_length_sd_A", "num_qtns_Lscaled_sd_A")
colnames(df.nonadapt)[9:14] <- c("inv_age_av_NA", "inv_length_av_NA", "num_qtns_Lscaled_av_NA", "inv_age_sd_NA", "inv_length_sd_NA", "num_qtns_Lscaled_sd_NA")
colnames(df.noselect)[9:14] <- c("inv_age_av_NS", "inv_length_av_NS", "num_qtns_Lscaled_av_NS", "inv_age_sd_NS", "inv_length_sd_NS", "num_qtns_Lscaled_sd_NS")

df.invChar.plot.temp <- left_join(df.adaptive, df.nonadapt, by = c("muBase", "muInv", "sigmaK", "alpha", "enVar", "mig1", "mig2"))
df.invChar.plot <- left_join(df.invChar.plot.temp, df.noselect, by = c("muBase", "muInv", "sigmaK", "alpha", "enVar", "mig1", "mig2"))

## START PLOTS
plot.age.adapt.pgen <- ggplot(data = df.invChar.plot[df.invChar.plot$enVar == 0 & 
                                                        df.invChar.plot$muBase == 0.002 & 
                                                        df.invChar.plot$alpha == 0.002,], 
                           aes(x = mig1, y= inv_age_av_A, group = sigmaK)) +
  geom_errorbar(aes(ymin=inv_age_av_A - inv_age_sd_A, 
                    ymax=inv_age_av_A + inv_age_sd_A), size = 0.4, width=0.5) +
  geom_ribbon(aes(ymin = inv_age_av_NA - inv_age_sd_NA, 
                  ymax = inv_age_av_NA + inv_age_sd_NA), fill = inferno(4)[2],  alpha = 0.5) + 
  geom_line(aes(y = inv_age_av_NA, x = mig1), linetype = "solid", color = inferno(4)[2]) + 
  geom_ribbon(aes(ymin = inv_age_av_NS - inv_age_sd_NS, 
                  ymax = inv_age_av_NS + inv_age_sd_NS), fill = inferno(4)[1],  alpha = 0.4) + 
  geom_line(aes(y = inv_age_av_NS, x = mig1), linetype = "dashed", color = inferno(4)[1]) + 
  geom_point(color = inferno(4)[3], size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "Polygenic Architecture", y = expression("Average Age"[inv]*" (Gen)"), x = " ") +
  ylim(-10000, 50000)

plot.age.adapt.ogen <- ggplot(data = df.invChar.plot[df.invChar.plot$enVar == 0 & 
                                                       df.invChar.plot$muBase == "0.0002" & 
                                                       df.invChar.plot$alpha == 0.2,], 
                              aes(x = mig1, y= inv_age_av_A, group = sigmaK)) +
  geom_errorbar(aes(ymin=inv_age_av_A - inv_age_sd_A, 
                    ymax=inv_age_av_A + inv_age_sd_A), size = 0.4, width=0.5) +
  geom_ribbon(aes(ymin = inv_age_av_NA - inv_age_sd_NA, 
                  ymax = inv_age_av_NA + inv_age_sd_NA), fill = inferno(4)[2],  alpha = 0.5) + 
  geom_line(aes(y = inv_age_av_NA, x = mig1), linetype = "solid", color = inferno(4)[2]) + 
  geom_ribbon(aes(ymin = inv_age_av_NS - inv_age_sd_NS, 
                  ymax = inv_age_av_NS + inv_age_sd_NS), fill = inferno(4)[1],  alpha = 0.4) + 
  geom_line(aes(y = inv_age_av_NS, x = mig1), linetype = "dashed", color = inferno(4)[1]) + 
  geom_point(color = inferno(4)[3], size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "Oligogenic Architecture", y = " ", x = " ") +
  ylim(-10000, 50000)

plot.length.adapt.pgen <- ggplot(data = df.invChar.plot[df.invChar.plot$enVar == 0 & 
                                                       df.invChar.plot$muBase == 0.002 & 
                                                       df.invChar.plot$alpha == 0.002,], 
                              aes(x = mig1, y= inv_length_av_A, group = sigmaK)) +
  geom_ribbon(aes(ymin = inv_length_av_NA - inv_length_sd_NA, 
                  ymax = inv_length_av_NA + inv_length_sd_NA), fill = inferno(4)[2],  alpha = 0.5) + 
  geom_line(aes(y = inv_length_av_NA, x = mig1), linetype = "solid", color = inferno(4)[2]) + 
  geom_ribbon(aes(ymin = inv_length_av_NS - inv_length_sd_NS, 
                  ymax = inv_length_av_NS + inv_length_sd_NS), fill = inferno(4)[1],  alpha = 0.4) + 
  geom_line(aes(y = inv_length_av_NS, x = mig1), linetype = "dashed", color = inferno(4)[1]) + 
  geom_errorbar(aes(ymin=inv_length_av_A - inv_length_sd_A, 
                    ymax=inv_length_av_A + inv_length_sd_A), size = 0.4, width=0.5) +
  geom_point(color = inferno(4)[3], size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = expression("Average Length"[inv]*" (bp)"), x = " ") +
  ylim(-5000, 60000)

plot.length.adapt.ogen <- ggplot(data = df.invChar.plot[df.invChar.plot$enVar == 0 & 
                                                       df.invChar.plot$muBase == "0.0002" & 
                                                       df.invChar.plot$alpha == 0.2,], 
                              aes(x = mig1, y= inv_length_av_A, group = sigmaK)) +
  geom_ribbon(aes(ymin = inv_length_av_NA - inv_length_sd_NA, 
                  ymax = inv_length_av_NA + inv_length_sd_NA), fill = inferno(4)[2],  alpha = 0.5) + 
  geom_line(aes(y = inv_length_av_NA, x = mig1), linetype = "solid", color = inferno(4)[2]) + 
  geom_ribbon(aes(ymin = inv_length_av_NS - inv_length_sd_NS, 
                  ymax = inv_length_av_NS + inv_length_sd_NS), fill = inferno(4)[1],  alpha = 0.4) + 
  geom_line(aes(y = inv_length_av_NS, x = mig1), linetype = "dashed", color = inferno(4)[1]) + 
  geom_errorbar(aes(ymin=inv_length_av_A - inv_length_sd_A, 
                    ymax=inv_length_av_A + inv_length_sd_A), size = 0.4, width=0.5) +
  geom_point(color = inferno(4)[3], size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = " ", x = " ") +
  ylim(-10000, 60000)

plot.QTN.adapt.pgen <- ggplot(data = df.invChar.plot[df.invChar.plot$enVar == 0 & 
                                                          df.invChar.plot$muBase == 0.002 & 
                                                          df.invChar.plot$alpha == 0.002,], 
                                 aes(x = mig1, y= num_qtns_Lscaled_av_A, group = sigmaK)) +
  geom_ribbon(aes(ymin = num_qtns_Lscaled_av_NA - num_qtns_Lscaled_sd_NA, 
                  ymax = num_qtns_Lscaled_av_NA + num_qtns_Lscaled_sd_NA), fill = inferno(4)[2],  alpha = 0.5) + 
  geom_line(aes(y = num_qtns_Lscaled_av_NA, x = mig1), linetype = "solid", color = inferno(4)[2]) + 
  geom_ribbon(aes(ymin = num_qtns_Lscaled_av_NS - num_qtns_Lscaled_sd_NS, 
                  ymax = num_qtns_Lscaled_av_NS + num_qtns_Lscaled_sd_NS), fill = inferno(4)[1],  alpha = 0.4) + 
  geom_line(aes(y = num_qtns_Lscaled_av_NS, x = mig1), linetype = "dashed", color = inferno(4)[1]) + 
  geom_errorbar(aes(ymin=num_qtns_Lscaled_av_A - num_qtns_Lscaled_sd_A, 
                    ymax=num_qtns_Lscaled_av_A + num_qtns_Lscaled_sd_A), size = 0.4, width=0.5) +
  geom_point(color = inferno(4)[3], size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = expression(bar(N)[QTNs] / "(Length"[inv]*" (bp))"), x = "Migration Rate") +
  ylim(min(df.invChar.plot$num_qtns_Lscaled_av_A - df.invChar.plot$num_qtns_Lscaled_sd_A), 
       max(df.invChar.plot$num_qtns_Lscaled_av_A - df.invChar.plot$num_qtns_Lscaled_sd_A))

plot.QTN.adapt.ogen <- ggplot(data = df.invChar.plot[df.invChar.plot$enVar == 0 & 
                                                          df.invChar.plot$muBase == "0.0002" & 
                                                          df.invChar.plot$alpha == 0.2,], 
                                 aes(x = mig1, y= num_qtns_Lscaled_av_A, group = sigmaK)) +
  geom_ribbon(aes(ymin = num_qtns_Lscaled_av_NA - num_qtns_Lscaled_sd_NA, 
                  ymax = num_qtns_Lscaled_av_NA + num_qtns_Lscaled_sd_NA), fill = inferno(4)[2],  alpha = 0.5) + 
  geom_line(aes(y = num_qtns_Lscaled_av_NA, x = mig1), linetype = "solid", color = inferno(4)[2]) + 
  geom_ribbon(aes(ymin = num_qtns_Lscaled_av_NS - num_qtns_Lscaled_sd_NS, 
                  ymax = num_qtns_Lscaled_av_NS + num_qtns_Lscaled_sd_NS), fill = inferno(4)[1],  alpha = 0.4) + 
  geom_line(aes(y = num_qtns_Lscaled_av_NS, x = mig1), linetype = "dashed", color = inferno(4)[1]) + 
  geom_errorbar(aes(ymin=num_qtns_Lscaled_av_A - num_qtns_Lscaled_sd_A, 
                    ymax=num_qtns_Lscaled_av_A + num_qtns_Lscaled_sd_A), size = 0.4, width=0.5) +
  geom_point(color = inferno(4)[3], size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = " ", x = "Migration Rate") +
  ylim(min(df.invChar.plot$num_qtns_Lscaled_av_A - df.invChar.plot$num_qtns_Lscaled_sd_A), 
       max(df.invChar.plot$num_qtns_Lscaled_av_A + df.invChar.plot$num_qtns_Lscaled_sd_A))

ggarrange(plot.age.adapt.pgen, plot.age.adapt.ogen, plot.length.adapt.pgen, plot.length.adapt.ogen, 
          plot.QTN.adapt.pgen, plot.QTN.adapt.ogen, ncol = 2, nrow = 3,
          widths = c(2.3, 2.3, 2.3, 2.3, 2.3, 2.3), labels = c("A", "B", "C", "D", "E", "F"))


#### characteristics 
######################################################################################################


######################################################################################################
#### Q5: Plot Genome scan results

df.outliers <- read.table("results/Inversion/20210525_fulldata/outputSumData.txt", header = FALSE,
                                 stringsAsFactors = FALSE)
head(df.muInv3)
colnames(df.outliers) <- colnames(df.summary)
df.outlierSumData <- left_join(df.outliers, df.simStats[,1:(ncol(df.simStats)-1)], by = "seed")
df.muInv3_outliers <- df.outlierSumData[df.outlierSumData$muInv == 1e-03,]

df.pcadapt.av <- aggregate(cbind(true_pos_pcadapt, false_neg_pcadapt, 
                                 true_neg_pcadapt, false_pos_pcadapt,
                                 true_neg_pcadapt_NS, false_pos_pcadapt_NS)~muBase + sigmaK + alpha + enVar + mig1 + mig2, 
                           FUN=mean, data = df.muInv3_outliers)
df.outflank.av <- aggregate(cbind(true_pos_outflank, false_neg_outflank, 
                                 true_neg_outflank, false_pos_outflank,
                                 true_neg_outflank_NS, false_pos_outflank_NS)~muBase + sigmaK + alpha + enVar + mig1 + mig2, 
                           FUN=mean, data = df.muInv3_outliers)

df.pcadapt.long <- as.data.frame(pivot_longer(df.pcadapt.av, cols = c(true_pos_pcadapt, false_neg_pcadapt, 
                                                                      true_neg_pcadapt, false_pos_pcadapt),
                                              names_to = "outcome", values_to = "count"))
for(i in 1:6 ){
  df.pcadapt.long[,i] <- as.factor(df.pcadapt.long[,i])
}

df.pcadapt.long$outcome <- as.factor(df.pcadapt.long$outcome)
df.pcadapt.sub$outcome <- recode_factor(df.pcadapt.sub$outcome, 'true_pos_pcadapt' = 'Adaptive Outlier', 'true_neg_pcadapt' = 'Nonadaptive Nonoutlier', 
                                         'false_pos_pcadapt' = 'Nonadaptive Outlier', 'false_neg_pcadapt' = 'Adaptive Nonoutlier')

df.pcadapt.sub$outcome <- recode_factor(df.pcadapt.sub$outcome, 'True Positive' = 'Adaptive Outlier', 'True Negative' = 'Nonadaptive Nonoutlier', 
                                        'False Positive' = 'Nonadaptive Outlier', 'False Negative' = 'Adaptive Nonoutlier')
#df.pcadapt.long$outcome <- factor(df.pcadapt.long$outcome, levels = c('True Positive', 'True Negative', 'False Positive', 'False Negative'))
df.pcadapt.subhigh <- df.pcadapt.long[df.pcadapt.long$enVar == 0 & df.pcadapt.long$alpha == "0.002" & df.pcadapt.long$muBase == "1e-06",]
 write.table(df.pcadapt.sub, "outlierHigh.txt")
plot.pcadapt.pgen <- ggplot(df.pcadapt.sub[order(df.pcadapt.sub$outcome),], 
                              aes(x = mig1, y = count, fill = outcome, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(4)[4:1]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "PCAdapt", y = "Polygenic Architecture\ncount", x = " ") + 
  ylim(c(0, 40))


df.pcadapt.sublow <- df.pcadapt.long[df.pcadapt.long$enVar == 0 & df.pcadapt.long$alpha == "0.2" & df.pcadapt.long$muBase == "1e-07",]
plot.pcadapt_lowQTN <- ggplot(df.pcadapt.sublow[order(df.pcadapt.sublow$outcome),], 
                               aes(x = mig1, y = count, fill = outcome, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(4)[4:1]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = "Oligogenic Architecture\ncount", x = "Migration Rate") + 
  ylim(c(0, 40)) +
  scale_x_discrete("Migration Rate", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE)



df.outflank.long <- as.data.frame(pivot_longer(df.outflank.av, cols = c(true_pos_outflank, false_neg_outflank, 
                                                                      true_neg_outflank, false_pos_outflank),
                                              names_to = "outcome", values_to = "count"))
for(i in 1:6 ){
  df.outflank.long[,i] <- as.factor(df.outflank.long[,i])
}

df.outflank.long$outcome <- recode_factor(df.outflank.long$outcome, 'true_pos_outflank' = 'True Positive', 'true_neg_outflank' = 'True Negative', 
                                         'false_pos_outflank' = 'False Positive', 'false_neg_outflank' = 'False Negative')
df.outflank.long$outcome <- factor(df.outflank.long$outcome, levels = c('True Positive', 'True Negative', 'False Positive', 'False Negative'))
df.outflank.subhigh <- df.outflank.long[df.outflank.long$enVar == 0 & df.outflank.long$alpha == 0.002 & df.outflank.long$muBase == 1e-06,]
plot.outflank_highQTN <- ggplot(df.outflank.subhigh[order(df.outflank.subhigh$outcome),], 
                               aes(x = mig1, y = count, fill = outcome, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(4)[4:1]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "OutFLANK", y = " ", x = " ") + 
  ylim(c(0, 40))

df.outflank.sublow <- df.outflank.long[df.outflank.long$enVar == 0 & df.outflank.long$alpha == 0.2 & df.outflank.long$muBase == 1e-07,]
plot.outflank_lowQTN <- ggplot(df.outflank.sublow[order(df.outflank.sublow$outcome),], 
                              aes(x = mig1, y = count, fill = outcome, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(4)[4:1]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = " ", x = "Migration Rate") + 
  ylim(c(0, 40)) +
  scale_x_discrete("Migration Rate", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE)


eval.legend <- g_legend(plot.pcadapt.pgen)
plot.pcadapt_highQTN.noLeg <- plot.pcadapt.pgen + theme(legend.position = "none")

ggarrange(plot.pcadapt_highQTN.noLeg, plot.outflank_highQTN, blank, 
          plot.pcadapt_lowQTN, plot.outflank_lowQTN, 
          eval.legend, ncol = 3, nrow = 2, widths = c(2.3,2.3,0.8,2.3,2.3,0.8),
          labels = c("A", "B", "", "C", "D"))

## width 15 height 12

## Large Alpha ##
df.pcadapt.subhigh_largeAlpha <- df.pcadapt.long[df.pcadapt.long$enVar == 0 & df.pcadapt.long$alpha == 0.2 & df.pcadapt.long$muBase == 1e-06,]
plot.pcadapt_highQTN_largeAlpha <- ggplot(df.pcadapt.subhigh_largeAlpha[order(df.pcadapt.subhigh_largeAlpha$outcome),], 
                               aes(x = mig1, y = count, fill = outcome, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(5)[4:1]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "PCAdapt", y = "High QTN Mutation Rate\ncount", x = " ") + 
  ylim(c(0, 40))


df.pcadapt.sublow_largeAlpha <- df.pcadapt.long[df.pcadapt.long$enVar == 0 & df.pcadapt.long$alpha == 0.2 & df.pcadapt.long$muBase == 1e-07,]
plot.pcadapt_lowQTN_largeAlpha <- ggplot(df.pcadapt.sublow_largeAlpha[order(df.pcadapt.sublow_largeAlpha$outcome),], 
                              aes(x = mig1, y = count, fill = outcome, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(5)[4:1]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = "Low QTN Mutation Rate\ncount", x = "Migration Rate") + 
  ylim(c(0, 40))

df.outflank.subhigh_largeAlpha <- df.outflank.long[df.outflank.long$enVar == 0 & df.outflank.long$alpha == 0.2 & df.outflank.long$muBase == 1e-06,]
plot.outflank_highQTN_largeAlpha <- ggplot(df.outflank.subhigh_largeAlpha[order(df.outflank.subhigh_largeAlpha$outcome),], 
                                aes(x = mig1, y = count, fill = outcome, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(5)[4:1]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "OutFLANK", y = " ", x = " ") + 
  ylim(c(0, 40))

df.outflank.sublow_largeAlpha <- df.outflank.long[df.outflank.long$enVar == 0 & df.outflank.long$alpha == 0.2 & df.outflank.long$muBase == 1e-07,]
plot.outflank_lowQTN_largeAlpha <- ggplot(df.outflank.sublow_largeAlpha[order(df.outflank.sublow_largeAlpha$outcome),], 
                               aes(x = mig1, y = count, fill = outcome, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(5)[4:1]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = " ", x = "Migration Rate") + 
  ylim(c(0, 40))

eval.legend <- g_legend(plot.outflank_lowQTN_largeAlpha)
plot.outflank_lowQTN_largeAlpha.noLeg <- plot.outflank_lowQTN_largeAlpha + theme(legend.position = "none")

ggarrange(plot.pcadapt_highQTN_largeAlpha, plot.outflank_highQTN_largeAlpha, blank, 
          plot.pcadapt_lowQTN_largeAlpha, plot.outflank_lowQTN_largeAlpha.noLeg, 
          eval.legend, ncol = 3, nrow = 2, widths = c(2.3,2.3,0.8,2.3,2.3,0.8), 
          labels = c("A", "B", "","C", "D"))


## No Selection Sims

df.pcadapt.NS.long <- as.data.frame(pivot_longer(df.pcadapt.av, cols = c(true_neg_pcadapt_NS, false_pos_pcadapt_NS),
                                              names_to = "outcome", values_to = "count"))
for(i in 1:6 ){
  df.pcadapt.NS.long[,i] <- as.factor(df.pcadapt.NS.long[,i])
}

df.pcadapt.NS.long$outcome <- recode_factor(df.pcadapt.NS.long$outcome, 'True Negative' = 'Nonadaptive Nonoutlier', 
                                         'False Positive' = 'Nonadaptive Outlier')
df.pcadapt.long$outcome <- factor(df.pcadapt.long$outcome, levels = c('Nonadaptive Nonoutlier', 'Nonadaptive Outlier'))
df.pcadapt.NS.subhigh <- df.pcadapt.NS.long[df.pcadapt.NS.long$enVar == 0 & df.pcadapt.NS.long$alpha == 0.002 & df.pcadapt.NS.long$muBase == 1e-06,]
plot.pcadapt_NS_highQTN <- ggplot(df.pcadapt.NS.subhigh[order(df.pcadapt.NS.subhigh$outcome),], 
                               aes(x = mig1, y = count, fill = outcome, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(4)[c(3,2)]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "PCAdapt", y = "Polygenic Architecture\ncount", x = " ")

df.pcadapt.NS.sublow <- df.pcadapt.NS.long[df.pcadapt.NS.long$enVar == 0 & df.pcadapt.NS.long$alpha == 0.2 & df.pcadapt.NS.long$muBase == 1e-07,]
plot.pcadapt_NS_lowQTN <- ggplot(df.pcadapt.NS.sublow[order(df.pcadapt.NS.sublow$outcome),], 
                              aes(x = mig1, y = count, fill = outcome, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(4)[c(3,2)]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = "Oligogenic Architecture\ncount", x = "Migration Rate") +
  scale_x_discrete("Migration Rate", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE)

df.outflank.NS.long <- as.data.frame(pivot_longer(df.outflank.av, cols = c(true_neg_outflank_NS, false_pos_outflank_NS),
                                                 names_to = "outcome", values_to = "count"))
for(i in 1:6){
  df.outflank.NS.long[,i] <- as.factor(df.outflank.NS.long[,i])
}

df.outflank.NS.long$outcome <- recode_factor(df.outflank.NS.long$outcome, 'True Negative' = 'Nonadaptive Nonoutlier', 
                                            'False Positive' = 'Nonadaptive Outlier')
df.outflank.NS.long$outcome <- factor(df.outflank.NS.long$outcome, levels = c('Nonadaptive Nonoutlier', 'Nonadaptive Outlier'))

df.outflank.NS.subhigh <- df.outflank.NS.long[df.outflank.NS.long$enVar == 0 & df.outflank.NS.long$alpha == 0.002 & df.outflank.NS.long$muBase == 1e-06,]
plot.outflank_NS_highQTN <- ggplot(df.outflank.NS.subhigh[order(df.outflank.NS.subhigh$outcome),], 
                                aes(x = mig1, y = count, fill = outcome, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(4)[c(3,2)]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "OutFLANK", y = " ", x = " ")

df.outflank.NS.sublow <- df.outflank.NS.long[df.outflank.NS.long$enVar == 0 & df.outflank.NS.long$alpha == 0.2 & df.outflank.NS.long$muBase == 1e-07,]
plot.outflank_NS_lowQTN <- ggplot(df.outflank.NS.sublow[order(df.outflank.NS.sublow$outcome),], 
                               aes(x = mig1, y = count, fill = outcome, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(4)[c(3,2)]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = " ", x = "Migration Rate") +
  scale_x_discrete("Migration Rate", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE)

genomeScans.legend <- g_legend(plot.pcadapt_NS_highQTN)
plot.pcadapt_NS_highQTN.noLeg <- plot.pcadapt_NS_highQTN + theme(legend.position = "none")

ggarrange(plot.pcadapt_NS_highQTN.noLeg, plot.outflank_NS_highQTN, blank, plot.pcadapt_NS_lowQTN, plot.outflank_NS_lowQTN, 
          genomeScans.legend, ncol = 3, nrow = 2, widths = c(2.3,2.3,0.8,2.3,2.3,0.8), 
          labels = c("A", "B", "", "C", "D"))

## Large Alpha ##
df.pcadapt.NS.subhigh_largeAlpha <- df.pcadapt.NS.long[df.pcadapt.NS.long$enVar == 0 & df.pcadapt.NS.long$alpha == 0.2 & df.pcadapt.NS.long$muBase == 1e-06,]
plot.pcadapt_NS_highQTN_largeAlpha <- ggplot(df.pcadapt.NS.subhigh_largeAlpha[order(df.pcadapt.NS.subhigh_largeAlpha$outcome),], 
                                  aes(x = mig1, y = count, fill = outcome, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(5)[c(3,2)]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "PCAdapt", y = "High QTN Mutation Rate\ncount", x = " ")

df.pcadapt.NS.sublow_largeAlpha <- df.pcadapt.NS.long[df.pcadapt.NS.long$enVar == 0 & df.pcadapt.NS.long$alpha == 0.2 & df.pcadapt.NS.long$muBase == 1e-07,]
plot.pcadapt_NS_lowQTN_largeAlpha <- ggplot(df.pcadapt.NS.sublow_largeAlpha[order(df.pcadapt.NS.sublow_largeAlpha$outcome),], 
                                 aes(x = mig1, y = count, fill = outcome, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(5)[c(3,2)]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = "Low QTN Mutation Rate\ncount", x = "Migration Rate")

df.outflank.NS.subhigh_largeAlpha <- df.outflank.NS.long[df.outflank.NS.long$enVar == 0 & df.outflank.NS.long$alpha == 0.2 & df.outflank.NS.long$muBase == 1e-06,]
plot.outflank_NS_highQTN_largeAlpha <- ggplot(df.outflank.NS.subhigh_largeAlpha[order(df.outflank.NS.subhigh_largeAlpha$outcome),], 
                                   aes(x = mig1, y = count, fill = outcome, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(5)[c(3,2)]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "OutFLANK", y = " ", x = " ")

df.outflank.NS.sublow_largeAlpha <- df.outflank.NS.long[df.outflank.NS.long$enVar == 0 & df.outflank.NS.long$alpha == 0.2 & df.outflank.NS.long$muBase == 1e-07,]
plot.outflank_NS_lowQTN_largeAlpha <- ggplot(df.outflank.NS.sublow_largeAlpha[order(df.outflank.NS.sublow_largeAlpha$outcome),], 
                                  aes(x = mig1, y = count, fill = outcome, group = sigmaK)) + 
  geom_bar(position = "stack", stat="identity", size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = viridis(5)[c(3,2)]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = " ", x = "Migration Rate")

genomeScans.legend_largeAlpha <- g_legend(plot.pcadapt_NS_highQTN_largeAlpha)
plot.pcadapt_NS_highQTN_largeAlpha.noLeg <- plot.pcadapt_NS_highQTN_largeAlpha + theme(legend.position = "none")

ggarrange(plot.pcadapt_NS_highQTN_largeAlpha.noLeg, plot.outflank_NS_highQTN_largeAlpha, blank, plot.pcadapt_NS_lowQTN_largeAlpha, plot.outflank_NS_lowQTN_largeAlpha, 
          genomeScans.legend_largeAlpha, ncol = 3, nrow = 2, widths = c(2.3,2.3,0.8,2.3,2.3,0.8))

#### end Q4 plotting
######################################################################################################


######################################################################################################
#### deleted code

# ## mu inversion 6
# for(i in 1:6){
#   df.muInv6_0_av[,i] <- as.factor(df.muInv6_0_av[,i])
# }
#df.muInv6_0_av$muBase <- recode_factor(df.muInv6_0_av$muBase, "0.000000001" = '0.000002', "0.00000001" = '0.00002','0.0000001' = '0.0002', '0.000001' = '0.002')

# plot.LA_diff_inv6_av <- ggplot(df.muInv6_0_av[df.muInv6_0_av$enVar == 0 & df.muInv6_0_av$alpha == 0.002,], 
#                                aes(x = mig1, y = LA_diff, group = interaction(muBase, sigmaK))) + 
#   geom_errorbar(aes(ymin=LA_diff_lowSD, ymax=LA_diff_upSD), width=.2) +
#   geom_point(aes(color = muBase, shape = muBase), size = 3 ) + 
#   facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
#   labs(title = "Difference in Amount of Local Adaptation",
#        y = expression("LA"[Inv]*" - LA"[noInv]),
#        x = "Migration Rate") +
#   guides(color = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)"),
#          shape = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)")) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90))  +
#   theme(panel.background = element_blank(), 
#         strip.background = element_rect(colour = "white", fill = "grey92")) +
#   scale_color_manual(values=viridis(4)[c(1:4)]) +
#   ylim(c(min(df.muInv6_0_av$LA_diff_lowSD[!is.na(df.muInv6_0_av$LA_diff_lowSD)]), 
#          max(df.muInv3_0_av$LA_diff_upSD[!is.na(df.muInv3_0_av$LA_diff_upSD)])))

## large alpha
# plot.LA_diff_inv3_largeAlpha_av <- ggplot(df.muInv3_0_av[df.muInv3_0_av$enVar == 0 & df.muInv3_0_av$alpha == 0.2,], 
#                                aes(x = mig1, y = LA_diff, group = interaction(muBase, sigmaK))) + 
#   geom_errorbar(aes(ymin=LA_diff_lowSD, ymax=LA_diff_upSD), width=.2) +
#   geom_point(aes(color = muBase, shape = muBase), size = 3 ) + 
#   facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
#   labs(title = "Difference in Amount of Local Adaptation",
#        y = expression("LA"[Inv]*" - LA"[noInv]),
#        x = "Migration Rate") +
#   guides(color = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)"),
#          shape = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)")) +
#   theme_classic() +
#   theme(legend.position = "none") + 
#   theme(axis.text.x = element_text(angle = 90))  +
#   theme(panel.background = element_blank(), 
#         strip.background = element_rect(colour = "white", fill = "grey92")) +
#   scale_color_manual(values=viridis(4)[c(2,4)]) +
#   ylim(min(df.muInv6_0_av$LA_diff_lowSD[!is.na(df.muInv6_0_av$LA_diff_lowSD)]), 
#        max(df.muInv3_0_av$LA_diff_upSD[!is.na(df.muInv3_0_av$LA_diff_upSD)]))
# 
# plot.LA_diff_inv6_largeAlpha_av <- ggplot(df.muInv6_0_av[df.muInv6_0_av$enVar == 0 & df.muInv6_0_av$alpha == 0.2,], 
#                                aes(x = mig1, y = LA_diff, group = interaction(muBase, sigmaK))) + 
#   geom_errorbar(aes(ymin=LA_diff_lowSD, ymax=LA_diff_upSD), width=.2) +
#   geom_point(aes(color = muBase, shape = muBase), size = 3 ) + 
#   facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
#   labs(title = "Difference in Amount of Local Adaptation",
#        y = expression("LA"[Inv]*" - LA"[noInv]),
#        x = "Migration Rate") +
#   guides(color = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)"),
#          shape = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)")) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90))  +
#   theme(panel.background = element_blank(), 
#         strip.background = element_rect(colour = "white", fill = "grey92")) +
#   scale_color_manual(values=viridis(4)[c(2,4)]) +
#   ylim(c(min(df.muInv6_0_av$LA_diff_lowSD[!is.na(df.muInv6_0_av$LA_diff_lowSD)]), 
#          max(df.muInv3_0_av$LA_diff_upSD[!is.na(df.muInv3_0_av$LA_diff_upSD)])))

# plot.LA_diff_inv3 <- ggplot(df.muInv3_0[df.muInv3_0$enVar == 0 & df.muInv3_0$alpha == 0.002,], 
#        aes(x = mig1, y = LA_diff, group = interaction(muBase, sigmaK))) + 
#   geom_point(aes(color = muBase, shape = muBase)) + 
#   facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) )+ 
#   labs(title = "High Inversion Mutation Rate (1e-3)",
#        y = "Difference in LA Between No Inv and Inv Sims",
#        x = "Migration Rate") +
#   guides(color = guide_legend(title = "QTN Mutation Rate"),
#          shape = guide_legend(title = "QTN Mutation Rate")) +
#   theme_classic() +
#   theme(panel.background = element_blank(), 
#         strip.background = element_rect(colour = "white", fill = "grey92")) +
#   scale_color_manual(values=viridis(4)[c(2,4)]) +
#   ylim(c(min(df.muInv6_0$LA_diff[!is.na(df.muInv6_0$LA_diff)]), 
#          max(df.muInv3_0$LA_diff[!is.na(df.muInv3_0$LA_diff)])))


# plot.LA_diff_inv6 <- ggplot(df.muInv6_0[df.muInv6_0$enVar == 0 & df.muInv6_0$alpha == 0.002,], 
#        aes(x = mig1, y = LA_diff, group = interaction(muBase, sigmaK))) + 
#   geom_point(aes(color = muBase, shape = muBase)) + 
#   facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) )+ 
#   labs(title = "Moderate Inversion Mutation Rate (1e-6)",
#        y = "Difference in LA Between No Inv and Inv Sims",
#        x = "Migration Rate") +
#   guides(color = guide_legend(title = "QTN Mutation Rate"),
#          shape = guide_legend(title = "QTN Mutation Rate")) +
#   theme_classic() +
#   theme(panel.background = element_blank(), 
#         strip.background = element_rect(colour = "white", fill = "grey92")) +
#   scale_color_manual(values=viridis(4)[c(2,4)]) + 
#   ylim(c(min(df.muInv6_0$LA_diff[!is.na(df.muInv6_0$LA_diff)]), 
#          max(df.muInv3_0$LA_diff[!is.na(df.muInv3_0$LA_diff)])))

# plot.VA_in_inv3 <- ggplot(df.muInv3_0[df.muInv3_0$enVar == 0 & df.muInv3_0$alpha == 0.002,], 
#        aes(x = mig1, y = VA_perc_In_3, group = interaction(muBase, sigmaK))) + 
#   geom_point(aes(color = muBase, shape = muBase)) + 
#   facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) )+ 
#   labs(title = "High Inversion Mutation Rate (1e-3)",
#        y = "Amount of Additive Genetic Variance Inside Inversions",
#        x = "Migration Rate") +
#   guides(color = guide_legend(title = "QTN Mutation Rate"),
#          shape = guide_legend(title = "QTN Mutation Rate")) +
#   theme_classic() +
#   theme(panel.background = element_blank(), 
#         strip.background = element_rect(colour = "white", fill = "grey92")) +
#   scale_color_manual(values=viridis(4)[c(2,4)]) + 
#   ylim(c(0, 
#          max(df.muInv3_0$VA_perc_In_3[!is.na(df.muInv3_0$VA_perc_In_3)])))

# muInv = 1e-06
# df.muInv6_0_av <- aggregate(LA_diff~muBase + sigmaK + alpha + enVar + mig1 + mig2, 
#                             FUN=mean, data = df.muInv6_0)
# df.muInv6_0_sd <- aggregate(LA_diff~muBase + sigmaK + alpha + enVar + mig1 + mig2, 
#                             FUN=sd, data = df.muInv6_0)
# df.muInv6_0_av$LA_diff_sd <-  df.muInv6_0_sd$LA_diff
# df.muInv6_0_av$LA_diff_upSD <- df.muInv6_0_av$LA_diff + df.muInv6_0_av$LA_diff_sd
# df.muInv6_0_av$LA_diff_lowSD <- df.muInv6_0_av$LA_diff - df.muInv6_0_av$LA_diff_sd

#### Amount of inverted regions Plots ####
# df.invGenome <- read.table(paste0(folderIn, "outputInvGenome_allData.txt"))
# head(df.invGenome)
# colnames(df.invGenome) <- c("seed", "uniqueBases", "numOverlap", "percGenome")
# df.invGenomeParam <- full_join(df.invGenome, df.simStats, by = "seed")
# df.invGenome_av <- aggregate(percGenome~muBase + sigmaK + muInv + alpha + enVar + mig1 + mig2, 
#                             FUN=mean, data = df.invGenomeParam)
# df.invGenome_sd <- aggregate(percGenome~muBase + sigmaK + muInv + alpha + enVar + mig1 + mig2, 
#                              FUN=sd, data = df.invGenomeParam)
# df.invGenome_av$percGenome_sd <- df.invGenome_sd$percGenome
# df.invGenome_av$percGenome_lowSD <- df.invGenome_av$percGenome - df.invGenome_av$percGenome_sd
# df.invGenome_av$percGenome_upSD <- df.invGenome_av$percGenome + df.invGenome_av$percGenome_sd
# 
# for(i in 1:6 ){
#   df.invGenome_av[,i] <- as.factor(df.invGenome_av[,i])
# }
# df.invGenome_av$muBase <- recode_factor(df.invGenome_av$muBase, "0.000000001" = '0.000002', "0.00000001" = '0.00002', '0.0000001' = '0.0002', '0.000001' = '0.002')
# plot.percInvGenome <- ggplot(df.invGenome_av[df.invGenome_av$enVar == 0 & df.invGenome_av$alpha == 0.002 & df.invGenome_av$muInv == 0.001,], 
#                                  aes(x = mig1, y = percGenome, group = interaction(muBase, sigmaK))) + 
#   geom_errorbar(aes(ymin=percGenome_lowSD, ymax=percGenome_upSD), width=.2) +
#   geom_point(aes(color = muBase, shape = muBase), size = 3) + 
#   facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) )+ 
#   labs(title = "Percent of Inverted Genome",
#        y = "Percent of Inverted Genome",
#        x = "Migration Rate") +
#   guides(color = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)"),
#          shape = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)")) +
#   theme_classic() +
#   theme(legend.position = "none") +
#   theme(axis.text.x = element_text(angle = 90))  +
#   theme(panel.background = element_blank(), 
#         strip.background = element_rect(colour = "white", fill = "grey92")) +
#   scale_color_manual(values=viridis(4)[c(2,4)]) + 
#   ylim(c(0, 
#          max(df.muInv3_0_av$VA_perc_In_3[!is.na(df.muInv3_0_av$VA_perc_upSD)]) + 10))
# 
# plot.percInvGenome_largeAlpha <- ggplot(df.invGenome_av[df.invGenome_av$enVar == 0 & df.invGenome_av$alpha == 0.2 & df.invGenome_av$muInv == 0.001,], 
#                                         aes(x = mig1, y = percGenome, group = interaction(muBase, sigmaK))) + 
#   geom_errorbar(aes(ymin=percGenome_lowSD, ymax=percGenome_upSD), width=.2) +
#   geom_point(aes(color = muBase, shape = muBase), size = 3) + 
#   facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) )+ 
#   labs(title = "Percent of Inverted Genome",
#        y = "Percent of Inverted Genome",
#        x = "Migration Rate") +
#   guides(color = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)"),
#          shape = guide_legend(title = "QTN Mutation Rate\n(Ne*mu)")) +
#   theme_classic() +
#   theme(legend.position = "none") +
#   theme(axis.text.x = element_text(angle = 90))  +
#   theme(panel.background = element_blank(), 
#         strip.background = element_rect(colour = "white", fill = "grey92")) +
#   scale_color_manual(values=viridis(4)[c(2,4)]) + 
#   ylim(c(0, 
#          max(df.muInv3_0_av$VA_perc_In_3[!is.na(df.muInv3_0_av$VA_perc_upSD)]) + 10))
# 
# plot.numInv.leg <- g_legend(plot.numInv)
# plot.numInv.noLeg <- plot.numInv + theme(legend.position = "none")
# 
# plot.numInv_largeAlpha.noLeg <- plot.numInv_largeAlpha + theme(legend.position = "none")
# plot.numInv_largeAlpha.leg <- g_legend(plot.numInv_largeAlpha)
# 
# ggarrange(plot.LA_diff_inv3_av, plot.VA_in_inv3, plot.percInvGenome, plot.numInv.noLeg, plot.numInv.leg, 
#           nrow = 1, ncol = 5, widths = c(2.3,2.3,2.3,2.3,0.8))
# ggarrange(plot.LA_diff_inv3_largeAlpha_av, plot.VA_in_inv3_largeAlpha, plot.percInvGenome_largeAlpha, plot.numInv_largeAlpha.noLeg, 
#           plot.numInv_largeAlpha.leg, 
#           nrow = 1, ncol = 5, widths = c(2.3,2.3,2.3,2.3,0.8))
# 
# ggarrange(plot.LA_diff_inv3_pgen, plot.VA_in_inv3, blank, plot.LA_diff_inv3_mgen, 
#           plot.VA_in_inv3_largeAlpha, plot.numInv.leg, widths = c(2.3,2.3,0.8,2.3,2.3,0.8), ncol = 3, nrow = 2)
# 

######################################################################################################    
## COPY AND PASTE WHERE NEEDED
pdf(paste0("figures/", seed, "_XXX.pdf"), height = 5, width = 7)

dev.off()

png(paste0("figures/", seed, "XXXX.png"), width = 480, height = 480, units = "px")

dev.off()