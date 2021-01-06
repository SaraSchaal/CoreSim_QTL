# install packages
  install.packages("purrr")
  library(purrr)
  library(ggplot2)
  folder <- "results/Inversion/20201115_FullSet/"

# download data
  df.invTime <- read.table(paste(folder, "FullSet_invTime.txt", sep = ""), header = TRUE)
  df.invData <- read.table(paste(folder, "FullSet_invData.txt", sep = ""), header = TRUE)
  df.params <- read.table(paste(folder, "invSimParams.txt", sep = ""), header = TRUE)
  unique.params <- read.table(paste(folder, "FullSet_uniqueParams.txt", sep = ""), header = TRUE)
  inv.sims <- subset(unique.params, subset = unique.params$mu_inv > 0)

# equations needed
  SE <- function(x){
    sd(x)/sqrt(length(x))
  }

# manipulate dataframes
  df.invTimeFreq <- subset(df.invTime, subset = df.invTime$freq > 0.05)
  df.invData$inv_id <- as.numeric(df.invData$inv_id)
  df.invAllData <- merge(df.invData, df.invTimeFreq, by.x = c("seed", "inv_id"), 
                         by.y = c("seed", "inv_id"), all.y = TRUE)
  colnames(df.invAllData)[8] <- "sim_gen"
  df.invAllData$inv_age <- df.invAllData$sim_gen - df.invAllData$inv_originGen

# calculate summary stats on all data (not averaged across reps)
  ave.length.all <- aggregate(inv_length~sim_gen + seed, data = df.invAllData, FUN = mean)
  colnames(ave.length.all)[3] <- "aveLength"
  se.length.all <- aggregate(inv_length~sim_gen + seed, data = df.invAllData, FUN = SE)
  colnames(se.length.all)[3] <- "seLength"
  ave.numQTNs.all <- aggregate(num_qtns~sim_gen + seed, data = df.invAllData, FUN = mean)
  colnames(ave.numQTNs.all)[3] <- "aveQTNs"
  se.numQTNs.all <- aggregate(num_qtns~sim_gen + seed, data = df.invAllData, FUN = SE)
  colnames(se.numQTNs.all)[3] <- "seQTNs"
  ave.age.all <- aggregate(inv_age~sim_gen + seed, data = df.invAllData, FUN = mean)
  colnames(ave.age.all)[3] <- "aveAge"
  se.age.all <- aggregate(inv_age~sim_gen + seed, data = df.invAllData, FUN = SE)
  colnames(se.age.all)[3] <- "seAge"
  
  df.merged <- Reduce(function(x,y) {merge(x = x, y = y, by = c("seed", "sim_gen"))},
                      list(ave.length.all, se.length.all, ave.numQTNs.all, 
                           se.numQTNs.all, ave.age.all, se.age.all))

  # code to spot check merged worked properly
  ## mean(subset(df.invAllData, select = inv_length, subset = seed == 3383642 & sim_gen == 32400)[,1])
  
# merge parameters and seed 
  df.all.data <- merge(df.merged, df.params, by = "seed", all.x = TRUE)
  df.ave.data <- aggregate(cbind(aveLength, aveAge, aveQTNs)~muBase+sigmaK+muInv+alpha+mig1+enVar+sim_gen, data = df.all.data, FUN = mean)
 # df.all.ave.data <- merge(df.params, df.ave.data, by = "", all.y = TRUE)
  
 # more spot checking code
 ## mean(subset(df.all.ave.data, select = aveLength, subset = seed == 3385440 & sim_gen == 10800)[,1])
  
#########################################################################################################
### CODE FOR AVERAGING RELICATES ###
# 
#   df.invLength.average <- NULL
#   reps <- 5
#   for(i in 1:nrow(inv.sims)){
#     # create empty variable for putting each data set that should be average
#     df.length.average <- NULL
#     av.length.seeds <- NULL
#     # step through the inversion information dataframe
#     for(j in 1:nrow(df.invAllData)){
#       # step through the different seeds that have the unique parameters
#       for(k in 1:reps){
#         seedCol <- paste("Seed", k, sep="")
#         # if that seed is not an NA then do the next step
#         if(!is.na(inv.sims[i, seedCol])){
#           # if the seed from unique parameters matches the seed of the invTime dataset store it
#           if(df.invAllData$seed[j] == inv.sims[i, seedCol]){
#             df.length.average <- rbind(df.length.average, df.invAllData[j,])
#             av.length.seeds <- unique(c(av.length.seeds, inv.sims[i, seedCol]))
#           }
#         }
#       }
#     }
#     
#     
#     ## average columns of interest across reps
#     length.average <- aggregate(inv_length~sim_gen, data = df.length.average, FUN = mean)
#     length.se <- aggregate(inv_length~sim_gen, data = df.length.average, FUN = SE)
#     num.inv <- aggregate(inv_length~sim_gen, data = df.length.average, FUN = length)
#     numQTNs.average <- aggregate(num_qtns~sim_gen, data = df.length.average, FUN = mean)
#     av.inv.age <- aggregate(inv_age~sim_gen, data = df.length.average, FUN = mean)
#     se.inv.age <- aggregate(inv_age~sim_gen, data = df.length.average, FUN = SE)
#     df.average.inv <- cbind(length.average, length.se[,2], numQTNs.average[,2], 
#                             num.inv[,2], av.inv.age[,2], se.inv.age[,2])
#     df.simParams <- data.frame(mu_inv = rep(inv.sims$mu_inv[i], nrow(df.average.inv)), 
#                                mig = rep(inv.sims$mig1[i], nrow(df.average.inv)),
#                                alpha = rep(inv.sims$alpha[i], nrow(df.average.inv)), 
#                                sigmaK = rep(inv.sims$sigmaK[i], nrow(df.average.inv)), 
#                                enVar = rep(inv.sims$enVar[i], nrow(df.average.inv)), 
#                                mu_base = rep(inv.sims$mu_base[i], nrow(df.average.inv)))
#     
#     # create seed columns so we can keep track of relevant seed names
#     df.Seedcolumns <- NULL
#     vect.colNames <- NULL
#     for(m in 1:reps){
#       seedCol <- paste("Seed", m, sep = "")
#       if(!is.na(av.seeds[m])){
#         df.Seedcolumns <- cbind(df.Seedcolumns, rep(av.seeds[m], nrow(df.average.inv)))
#       } else {
#         df.Seedcolumns <- cbind(df.Seedcolumns, rep(NA, nrow(df.average.inv)))
#       }
#       vect.colNames <- c(vect.colNames, seedCol)
#     }
#     colnames(df.Seedcolumns) <- vect.colNames
#     
#     # finally bind together all the data in the final dataframe
#     df.invLength.average <- rbind(df.invLength.average, cbind(df.average.inv, df.simParams, df.Seedcolumns))
#   
#   } # close average for loop
#   ###################################################################
#   colnames(df.invLength.average)[3:7] <- c("SD_length", "ave_num_qtns", "num_inv", "ave_inv_age", "SD_inv_age")
#   
#   write.table(df.invLength.average, "FullSet_invInfo.txt")

### CLOSE CODE FOR AVERAGING REPS ###
#############################################################################################################
  
  
  
  # convert every column of parameters data to factor 
  for(i in 9:ncol(df.all.data)){
    df.all.data[,i] <- as.factor(as.character(df.all.data[,i]))
  }  
  
  for(i in 1:(ncol(df.ave.data)-4)){
    df.ave.data[,i] <- as.factor(as.character(df.ave.data[,i]))
  }
  
#############################################################################################################
### #PLOTTING ####
 # hmm <-df.ave.data[df.ave.data$sigmaK == 3 & df.ave.data$alpha == 0.002 & df.ave.data$mig1 == 0.001 & df.ave.data$muBase == 1e-07 & df.ave.data$sim_gen == 38000,]
  ## make labels for facet wrapping
    alpha.labels <- c("0.2" = "sigmaMu = 0.2", "0.002" = "sigmaMu = 0.002")
    mig.labels <- c("0.001" = "mig = 0.001", "0.01" = "mig = 0.01", "0.1" = "mig = 0.1", 
                  "0.25" = "mig = 0.25", "0.4" = "mig = 0.4", "0.5" = "mig = 0.5")
  # ggplot(data = df.ave.data[df.ave.data$sigmaK == 3 & df.ave.data$alpha == 0.002 & df.ave.data$mig1 == 0.001,],
  #        aes(x = sim_gen, y = aveAge, group = muBase)) + 
  #       geom_line(aes(color = muBase), size = 0.75)
  # 
  # ## plot inversion age  -- facets are alpha and migration
  #   inv.Age <-  ggplot(data = df.ave.data[df.ave.data$sigmaK == 0.75 & df.ave.data$enVar == 0 & df.ave.data$muInv == 0.001,], 
  #                      aes(x = sim_gen, y = aveAge, group = muBase)) + 
  #               geom_line(aes(color = muBase), size = 0.75) + 
  #               facet_wrap(~ alpha + mig1, labeller = labeller(alpha = alpha.labels, mig1 = mig.labels),
  #                          ncol = 6, nrow = 2) + 
  #               labs(title = "Average Inversions Age Through Time",
  #                    y = "Average Inversion Age",
  #                    x = "Generation") +
  #               guides(color = guide_legend(title = "QTN Mutation Rate"), 
  #                      linetype = guide_legend(title = "Migration Rate")) +
  #               theme_classic() +
  #               theme(panel.background = element_blank(), 
  #                     strip.background = element_rect(colour = "white", fill = "grey92")) +
  #               scale_color_manual(values=c( "cadetblue1","cadetblue3", "cornflowerblue", "navy", "black")) +
  #               scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  #               scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.ave.data$aveAge)))

  ## plot inversion age -- facets are migration and sigma K lines are QTN mutation rate and alpha
    inv.Age <-  ggplot(data = df.ave.data[df.ave.data$enVar == 0 & df.ave.data$muInv == 0.001,], 
                       aes(x = sim_gen, y = aveAge, group = interaction(muBase, alpha))) + 
                geom_line(aes(color = muBase, linetype = alpha), size = 0.75) + 
                facet_wrap(~ sigmaK + mig1 , labeller = labeller(mig1 = mig.labels),
                           ncol = 6, nrow = 3) + 
                labs(title = "Average Inversions Age Through Time",
                     y = "Average Inversion Age",
                     x = "Generation") +
                guides(color = guide_legend(title = "QTN Mutation Rate"), 
                       linetype = guide_legend(title = "Sigma Mu (Effect Size)")) +
                theme_classic() +
                theme(panel.background = element_blank(), 
                      strip.background = element_rect(colour = "white", fill = "grey92")) +
                scale_color_manual(values=c( "cadetblue3", "navy")) +
                scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
                scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.ave.data$aveAge)))
     
    ## plot all data         
    inv.Age <-  ggplot(data = df.all.data[df.all.data$enVar == 0 & df.all.data$muInv == 0.001,], 
                       aes(x = sim_gen, y = aveAge, group = interaction(muBase, alpha, rep))) + 
      geom_line(aes(color = muBase, linetype = alpha), size = 0.75, alpha = 0.3) + 
      geom_line(data = df.ave.data[df.ave.data$enVar == 0 & df.ave.data$muInv == 0.001,], 
                aes(x = sim_gen, y = aveAge, group = interaction(muBase, alpha), linetype = alpha, color = muBase), size = 1.2) + 
      facet_wrap(~ sigmaK + mig1 , labeller = labeller(mig1 = mig.labels),
                 ncol = 6, nrow = 3) + 
      labs(title = "Average Inversions Age Through Time",
           y = "Average Inversion Age",
           x = "Generation") +
      guides(color = guide_legend(title = "QTN Mutation Rate"), 
             linetype = guide_legend(title = "Sigma Mu (Effect Size)")) +
      theme_classic() +
      theme(panel.background = element_blank(), 
            strip.background = element_rect(colour = "white", fill = "grey92")) +
      scale_color_manual(values=c( "cadetblue3", "navy")) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.ave.data$aveAge)))
    
  ## violin plot inversion age
    # df.all.data$sim_gen_fact <- as.factor(as.character(df.all.data$sim_gen))
    # inv.Age.viol <- ggplot(data = df.all.data[df.all.data$sigmaK == 0.75 & df.all.data$alpha == 0.002 & df.all.data$mig1 == 0.001,], 
    #                        aes(x = sim_gen_fact, y = aveAge, fill = muBase)) +
    #                # facet_wrap(~ alpha + mig1, labeller = labeller(alpha = alpha.labels, mig1 = mig.labels)) + 
    #                 geom_violin() + 
    #                 stat_summary(fun=mean, geom="point", shape=23, size=2)
    # 
    # 
    # 
    
    inv.Length <-  ggplot(data = df.all.data[df.all.data$enVar == 0 & df.all.data$muInv == 0.001,], 
                       aes(x = sim_gen, y = aveLength, group = interaction(muBase, alpha, rep))) + 
      geom_line(aes(color = muBase, linetype = alpha), size = 0.75, alpha = 0.3) + 
      geom_line(data = df.ave.data[df.ave.data$enVar == 0 & df.ave.data$muInv == 0.001,], 
                aes(x = sim_gen, y = aveLength, group = interaction(muBase, alpha), linetype = alpha, color = muBase), size = 1.2) + 
      facet_wrap(~ sigmaK + mig1 , labeller = labeller(mig1 = mig.labels),
                 ncol = 6, nrow = 3) + 
      labs(title = "Average Inversion Length Through Time",
           y = "Average Inversion Length",
           x = "Generation") +
      guides(color = guide_legend(title = "QTN Mutation Rate"), 
             linetype = guide_legend(title = "Sigma Mu (Effect Size)")) +
      theme_classic() +
      theme(panel.background = element_blank(), 
            strip.background = element_rect(colour = "white", fill = "grey92")) +
      scale_color_manual(values=c( "cadetblue3", "navy")) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.ave.data$aveLength)))
    
    
    inv.numQTNs <-  ggplot(data = df.all.data[df.all.data$enVar == 0 & df.all.data$muInv == 0.001,], 
                          aes(x = sim_gen, y = aveQTNs, group = interaction(muBase, alpha, rep))) + 
      geom_line(aes(color = muBase, linetype = alpha), size = 0.75, alpha = 0.3) + 
      geom_line(data = df.ave.data[df.ave.data$enVar == 0 & df.ave.data$muInv == 0.001,], 
                aes(x = sim_gen, y = aveQTNs, group = interaction(muBase, alpha), linetype = alpha, color = muBase), size = 1.2) + 
      facet_wrap(~ sigmaK + mig1 , labeller = labeller(mig1 = mig.labels),
                 ncol = 6, nrow = 3) + 
      labs(title = "Average Number of Inversion QTNs Through Time",
           y = "Average Number of Inversion QTNs",
           x = "Generation") +
      guides(color = guide_legend(title = "QTN Mutation Rate"), 
             linetype = guide_legend(title = "Sigma Mu (Effect Size)")) +
      theme_classic() +
      theme(panel.background = element_blank(), 
            strip.background = element_rect(colour = "white", fill = "grey92")) +
      scale_color_manual(values=c( "cadetblue3", "navy")) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.ave.data$aveQTNs)))
    