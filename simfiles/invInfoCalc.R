df.invTime <- read.table("FullSet_invTime.txt", header = TRUE)
df.invData <- read.table("FullSet_invData.txt", header = TRUE)

df.invTimeFreq <- subset(df.invTime, subset = df.invTime$freq > 0.05)
df.invData$inv_id <- as.numeric(df.invData$inv_id)

df.invAllData <- merge( df.invData, df.invTimeFreq, by.x = c("seed", "inv_id"), by.y = c("seed", "inv_id"), all.y = TRUE)
colnames(df.invAllData)[8] <- "sim_gen"
df.invAllData$inv_age <- df.invAllData$sim_gen - df.invAllData$inv_originGen

df.invLength.average <- NULL
df.invLength.all <- NULL
inv.sims <- subset(unique.params, subset = unique.params$mu_inv > 0)
for(i in 1:nrow(inv.sims)){
  # create empty variable for putting each data set that should be average
  df.length.average <- NULL
  av.length.seeds <- NULL
  # step through the population dynamics dataframe
  for(j in 1:nrow(df.invAllData)){
    # step through the different seeds that have the unique parameters
    for(k in 1:reps){
      seedCol <- paste("Seed", k, sep="")
      # if that seed is not an NA then do the next step
      if(!is.na(inv.sims[i, seedCol])){
        # if the seed from unique parameters matches the seed of the invTime dataset store it
        if(df.invAllData$seed[j] == inv.sims[i, seedCol]){
          df.length.average <- rbind(df.length.average, df.invAllData[j,])
          av.length.seeds <- unique(c(av.length.seeds, inv.sims[i, seedCol]))
        }
      }
    }
  }
  
  
  ## average columns of interest
  # df.averageSubset <- df.average[,vect.aggCols]
  SE <- function(x){
    sd(x)/sqrt(length(x))
  }
  length.average <- aggregate(inv_length~sim_gen, data = df.length.average, FUN = mean)
  length.se <- aggregate(inv_length~sim_gen, data = df.length.average, FUN = SE)
  
  num.inv <- aggregate(inv_length~sim_gen, data = df.length.average, FUN = length)
  numQTNs.average <- aggregate(num_qtns~sim_gen, data = df.length.average, FUN = mean)
  av.inv.age <- aggregate(inv_age~sim_gen, data = df.length.average, FUN = mean)
  se.inv.age <- aggregate(inv_age~sim_gen, data = df.length.average, FUN = se)
  df.average.inv <- cbind(length.average, length.SD[,2], numQTNs.average[,2], 
                          num.inv[,2], av.inv.age[,2], se.inv.age[,2])
  df.simParams <- data.frame(mu_inv = rep(inv.sims$mu_inv[i], nrow(df.average.inv)), 
                             mig = rep(inv.sims$mig1[i], nrow(df.average.inv)),
                             alpha = rep(inv.sims$alpha[i], nrow(df.average.inv)), 
                             sigmaK = rep(inv.sims$sigmaK[i], nrow(df.average.inv)), 
                             enVar = rep(inv.sims$enVar[i], nrow(df.average.inv)), 
                             mu_base = rep(inv.sims$mu_base[i], nrow(df.average.inv)))
  
  # create seed columns so we can keep track of relevant seed names
  df.Seedcolumns <- NULL
  vect.colNames <- NULL
  #vect.NAseeds <- NULL
  for(m in 1:reps){
    seedCol <- paste("Seed", m, sep = "")
    if(!is.na(av.seeds[m])){
      df.Seedcolumns <- cbind(df.Seedcolumns, rep(av.seeds[m], nrow(df.average.inv)))
    } else {
      df.Seedcolumns <- cbind(df.Seedcolumns, rep(NA, nrow(df.average.inv)))
    }
    vect.colNames <- c(vect.colNames, seedCol)
  }
  colnames(df.Seedcolumns) <- vect.colNames
  
  # finally bind together all the data in the final dataframe
  df.invLength.average <- rbind(df.invLength.average, cbind(df.average.inv, df.simParams, df.Seedcolumns))
  df.invLength.all <- 
  df.LA.all <- rbind(df.LA.all, cbind(df.average, df.simParamsAll, df.SeedcolumnsAll))
} # close average for loop
###################################################################
colnames(df.invLength.average)[3:7] <- c("SD_length", "ave_num_qtns", "num_inv", "ave_inv_age", "SD_inv_age")

write.table(df.invLength.average, "FullSet_invInfo.txt")