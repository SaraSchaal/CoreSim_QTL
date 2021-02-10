
unique.params <- read.table("FullSet_uniqueParams.txt", header = TRUE)
df.popDyn <- read.table("FullSet_popDyn.txt", header = TRUE)

########################
### average for loop ###
###################################################################

df.LA.average <- NULL
# step through each line of the unique parameters dataframe
for(i in 1:nrow(unique.params)){
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
          df.average <- rbind(df.average, df.popDyn[j,])
          av.seeds <- unique(c(av.seeds, unique.params[i, seedCol]))
        }
      }
    }
  }
  
  ## average columns of interest
  # df.averageSubset <- df.average[,vect.aggCols]
  new.average <- aggregate(.~sim_gen, data = df.average, FUN = mean)
  df.simParams <- data.frame(mu_inv = rep(unique.params$mu_inv[i], nrow(new.average)), 
                             mig = rep(unique.params$mig1[i], nrow(new.average)),
                             alpha = rep(unique.params$alpha[i], nrow(new.average)), 
                             sigmaK = rep(unique.params$sigmaK[i], nrow(new.average)), 
                             enVar = rep(unique.params$enVar[i], nrow(new.average)), 
                             mu_base = rep(unique.params$mu_base[i], nrow(new.average)))
  
  
  # create seed columns so we can keep track of relevant seed names
  df.Seedcolumns <- NULL
  vect.colNames <- NULL
  vect.NAseeds <- NULL
  for(m in 1:reps){
    seedCol <- paste("Seed", m, sep = "")
    if(!is.na(av.seeds[m])){
      df.Seedcolumns <- cbind(df.Seedcolumns, rep(av.seeds[m], nrow(new.average)))
    } else {
      df.Seedcolumns <- cbind(df.Seedcolumns, rep(NA, nrow(new.average)))
    }
    vect.colNames <- c(vect.colNames, seedCol)
  }
  colnames(df.Seedcolumns) <- vect.colNames
  
  # finally bind together all the data in the final dataframe
  df.LA.average <- rbind(df.LA.average, cbind(new.average, df.simParams, df.Seedcolumns))
  
} # close average for loop
###################################################################
df.LA.average <- subset(df.LA.average, select = -seed)

write.table(df.LA.average, "FullSet_localAdaptation.txt")