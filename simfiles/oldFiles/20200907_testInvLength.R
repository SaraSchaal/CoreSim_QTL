#### Test for inversion size issue ####
#df.params <- as.data.frame(unique(substr(list.files("./results/Inversion/20200907_testForInvLength", 
                                                    #pattern = "\\d{13}")[1], 1, 13)))

#colnames(df.params) <- "seed"

df.params <- data.frame(seed = 1)

#folder <- "results/Inversion/20200907_testForInvLength/"
folder <- "results/Inversion/"


## upload datasets
###################################################################################################
## Simulation parameters
df.simStats <- NULL
for(i in 1:nrow(df.params)){
  seed <- df.params$seed[i]
  if(seed %in% not.done.seeds){
  } else {
    simStatsNewFile <- read.table(paste(folder, seed, "_outputSimStats.txt", sep=""), 
                                  stringsAsFactors = FALSE)
    simStatsNewFile$seed <- seed
    df.simStats <-  rbind(df.simStats, simStatsNewFile)
  }
}
df.simStats <- df.simStats[,2:ncol(df.simStats)]
colnames(df.simStats) <- c("mig1", "mig2", "pop1N", "pop2N", "mu_base", "mu_inv", "r", 
                           "alpha", "sigmaK", "burnin", "dom", "enVar", "Seed")

unique.params <- df.simStats

## Inversion Through Time
df.invTime <- NULL
for(i in 1:nrow(df.params)){
  seed <- df.params$seed[i]
  invTimeNewFile <- read.table(paste(folder, seed, "_outputInvTime.txt", sep=""), header = TRUE, stringsAsFactors = FALSE)
  if(nrow(invTimeNewFile) > 0){
    invTimeNewFile$seed <- seed
    df.invTime <-  rbind(df.invTime, invTimeNewFile)
  }
}

df.invData <- NULL
for(i in 1:nrow(df.params)){
  seed <- df.params$seed[i]
  invData <- read.table(paste(folder, seed, "_outputInvSumInfo.txt", sep=""), header = TRUE,
                        stringsAsFactors = FALSE)
  if(nrow(invData) > 0){
    invData$seed <- seed
    df.invData <- rbind(df.invData, invData)
  }
}
#############################################################################

## inversion calculations
#############################################################################
#df.invTimeFreq <- subset(df.invTime, subset = df.invTime$freq > 0.05)
df.invData$inv_id <- as.numeric(df.invData$inv_id)

df.invAllData <- merge(df.invData, df.invTime, by.x = c("seed", "inv_id"), by.y = c("seed", "inv_id"), all.y = TRUE)
colnames(df.invAllData)[8] <- "sim_gen"

df.invLength.average <- NULL
inv.sims <- subset(unique.params, subset = unique.params$mu_inv > 0)

  
## average columns of interest
# df.averageSubset <- df.average[,vect.aggCols]
new.length.average <- aggregate(inv_length~sim_gen + seed, data = df.invAllData, FUN = mean)
new.length.SD <- aggregate(inv_length~sim_gen + seed, data = df.invAllData, FUN = sd)
new.num.inv <- aggregate(inv_length~sim_gen, data = df.length.average, FUN = length)
new.numQTNs.average <- aggregate(num_qtns~sim_gen, data = df.length.average, FUN = mean)
df.average.inv <- cbind(new.length.average, new.length.SD[,2],
                        new.numQTNs.average[,2], new.num.inv[,2])
df.simParams <- data.frame(mu_inv = rep(inv.sims$mu_inv[i], nrow(df.average.inv)), 
                           mig = rep(inv.sims$mig1[i], nrow(df.average.inv)),
                           alpha = rep(inv.sims$alpha[i], nrow(df.average.inv)), 
                           sigmaK = rep(inv.sims$sigmaK[i], nrow(df.average.inv)), 
                           enVar = rep(inv.sims$enVar[i], nrow(df.average.inv)), 
                           mu_base = rep(inv.sims$mu_base[i], nrow(df.average.inv)))
  

