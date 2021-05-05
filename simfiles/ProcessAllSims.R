######################################################################################################  
#### The following code will process all simulation files and output results figures
#### Sara M. Schaal
######################################################################################################  
#### Run the following code chunk once to get full data files to do further analyses on ####

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

#### END Run the following code chunk once to get full data files to do further analyses on ####
######################################################################################################  


######################################################################################################  
#### Unique Parameters ####
# Subset original parameters dataframe to identify all unique parameter combinations. 
# Then using those parameter values identify the simulation seeds of replicate simulations 
# for each set of parameters.


unique.params <- unique(df.simStats[c("mig1","mig2", "pop1N", "pop2N", "mu_base", "mu_inv", "r", "alpha", "sigmaK", "burnin", "dom", "enVar")])
reps <- 2
for(i in 1:reps){
  seedColName <- paste("Seed", i, sep = "")
  unique.params[,seedColName] <- NA
}

for(i in 1:nrow(unique.params)){
  for(j in 1:nrow(df.simStats)){
    if(apply(unique.params[i,1:12], 1, paste, collapse = " ") == 
       apply(df.simStats[j,c("mig1","mig2", "pop1N", "pop2N", "mu_base", "mu_inv", "r", "alpha", 
                             "sigmaK", "burnin", "dom", "enVar")], 1, paste, collapse = " ")){
      for(k in 1:reps){
        seedCol <- paste("Seed", k, sep="")
        if(is.na(unique.params[i, seedCol])){
          unique.params[i, seedCol] <- df.simStats$Seed[j]
          break
        } 
      }
    }
  }
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

df.LA.average <- NULL
df.LA.stderr <- NULL
df.LA.all <- NULL
count <- 0
# step through each line of the unique parameters dataframe
for(i in 1:nrow(unique.params)){
  ptm <- proc.time()
  # create empty variable for putting each data set that should be average
  df.average <- NULL
  av.seeds <- NULL
  # step through the population dynamics dataframe
  for(j in 1:nrow(df.popDyn)){
    # step through the different seeds that have the unique parameters
    for(k in 1:reps){
      seedCol <- paste("Seed", k, sep="")
      # if that seed is not an NA then do the next step
      if(!is.na(unique.params[i, seedCol])){
        # if the seed from unique parameters matches the seed of the pop dynam dataset store it
        if(df.popDyn$seed[j] == unique.params[i, seedCol]){
          # bind data from same parameter combination together for averaging later
          df.average <- rbind(df.average, df.popDyn[j,])
          # store unique seeds to this dataset
          av.seeds <- unique(c(av.seeds, unique.params[i, seedCol]))
        }
      }
    }
  }
  
  
  ## average columns of interest
  # df.averageSubset <- df.average[,vect.aggCols]
  new.average <- aggregate(.~sim_gen, data = df.average, FUN = mean)
  SE <- function(x){
    sd(x)/sqrt(length(x))
  }
  new.se <- aggregate(.~sim_gen, data = df.average, FUN = SE)
  df.simParams <- data.frame(mu_inv = rep(unique.params$mu_inv[i], nrow(new.average)), 
                             mig = rep(unique.params$mig1[i], nrow(new.average)),
                             alpha = rep(unique.params$alpha[i], nrow(new.average)), 
                             sigmaK = rep(unique.params$sigmaK[i], nrow(new.average)), 
                             enVar = rep(unique.params$enVar[i], nrow(new.average)), 
                             mu_base = rep(unique.params$mu_base[i], nrow(new.average)))
  
  df.simParamsAll <- data.frame(mu_inv = rep(unique.params$mu_inv[i], nrow(df.average)), 
                                mig = rep(unique.params$mig1[i], nrow(df.average)),
                                alpha = rep(unique.params$alpha[i], nrow(df.average)), 
                                sigmaK = rep(unique.params$sigmaK[i], nrow(df.average)), 
                                enVar = rep(unique.params$enVar[i], nrow(df.average)), 
                                mu_base = rep(unique.params$mu_base[i], nrow(df.average)))
  
  # create seed columns so we can keep track of relevant seed names
  df.Seedcolumns <- NULL
  df.SeedcolumnsAll <- NULL
  vect.colNames <- NULL
  for(m in 1:reps){
    seedCol <- paste("Seed", m, sep = "")
    if(!is.na(av.seeds[m])){
      df.Seedcolumns <- cbind(df.Seedcolumns, rep(av.seeds[m], nrow(new.average)))
      df.SeedcolumnsAll <- cbind(df.SeedcolumnsAll, rep(av.seeds[m], nrow(df.average)))
    } else {
      df.Seedcolumns <- cbind(df.Seedcolumns, rep(NA, nrow(new.average)))
      df.SeedcolumnsAll <- cbind(df.SeedcolumnsAll, rep(NA, nrow(df.average)))
    }
    vect.colNames <- c(vect.colNames, seedCol)
  }
  colnames(df.Seedcolumns) <- vect.colNames
  colnames(df.SeedcolumnsAll) <- vect.colNames
  
  # finally bind together all the data in the final dataframe
  df.LA.average <- rbind(df.LA.average, cbind(new.average, df.simParams, df.Seedcolumns))
  df.LA.stderr <- rbind(df.LA.stderr, cbind(new.se, df.simParams, df.Seedcolumns))
  df.LA.all <- rbind(df.LA.all, cbind(df.average, df.simParamsAll, df.SeedcolumnsAll))
  count <- count + 1
  print(count)
  print(proc.time() - ptm)
} # close average for loop
###################################################################
df.LA.average <- subset(df.LA.average, select = -seed)
df.LA.stderr <- subset(df.LA.stderr, select = -seed)
df.LA.all <- subset(df.LA.all, select = -seed)

#df.LA.average <- read.table("FullSet_LAaverage.txt", header = TRUE)
write.table(df.LA.average, "FullSet_LAaverage.txt", row.names = FALSE)
write.table(df.LA.stderr, "FullSet_LAstderr.txt", row.names = FALSE)
write.table(df.LA.all, "FullSet_LAallData.txt", row.names = FALSE)

# this should work I don't understand why it doesn't
#df.LA.average[,df.LA.average] <- lapply(df.LA.average[,LA.cols], function(x) as.factor(as.character(x)))
df.LA.average <- subset(df.LA.average, select = -seed)

# convert every column of parameters data to factor 
for(i in 16:ncol(df.LA.average)){
  df.LA.average[,i] <- as.factor(as.character(df.LA.average[,i]))
}


# Run this chunk if average file has already been made
df.LA.average <- read.table("src/results/FullSet_LAaverage.txt", header = TRUE, stringsAsFactors = FALSE)

for(i in 16:ncol(df.LA.average)){
  df.LA.average[,i] <- as.factor(as.character(df.LA.average[,i]))
}

#### end average LA calculation 
######################################################################################################


######################################################################################################
#### Q1 & 2: Plot Percent Additive Genetic Variance in Adaptive Inversions
read.table(paste(folder, seed, "_allSimSummaryCalcs.txt", sep=""), header = TRUE, stringsAsFactors = FALSE)







#### end Q1 & 2 plotting
######################################################################################################

######################################################################################################
#### Q3: Plot Proportion of Adaptive inversions in each evo history category








#### end Q3 plotting
######################################################################################################

######################################################################################################
#### Q4: Plot Characteristics








#### end Q4 plotting
######################################################################################################




######################################################################################################    
## COPY AND PASTE WHERE NEEDED
pdf(paste0("figures/", seed, "_XXX.pdf"), height = 5, width = 7)

dev.off()

png(paste0("figures/", seed, "XXXX.png"), width = 480, height = 480, units = "px")

dev.off()