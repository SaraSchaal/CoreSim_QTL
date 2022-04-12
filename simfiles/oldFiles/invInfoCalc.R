# install packages
  install.packages("purrr")
  library(purrr)
  library(ggplot2)
  folder <- "results/Inversion/20201115_FullSet/"

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
  
  # convert every column of parameters data to factor 
  for(i in 9:ncol(df.all.data)){
    df.all.data[,i] <- as.factor(as.character(df.all.data[,i]))
  }  
  
  for(i in 1:(ncol(df.ave.data)-4)){
    df.ave.data[,i] <- as.factor(as.character(df.ave.data[,i]))
  }
  
#############################################################################################################
#### Subset Inversions ####
  
  ## Top 10 percent of FST values
 top10.data <- df.invAllData %>%
    group_by(seed, sim_gen) %>%
    filter(inv_FST>=quantile(inv_FST, 0.9)) %>%
    summarise_at(c("inv_age", "mean_qtnSelCoef", "num_qtns", "inv_length"), 
                 mean, .groups = "keep") %>%
    rename(inv_ageT10 = inv_age, mean_qtnSelCoefT10 = mean_qtnSelCoef, 
           num_qtnsT10 = num_qtns, inv_lengthT10 = inv_length)
  
   ## Bottom 90 percent of FST values
  bottom90.data <- df.invAllData %>%
    group_by(seed, sim_gen) %>%
    filter(inv_FST<quantile(inv_FST, 0.9)) %>%
    summarise_at(c("inv_age", "mean_qtnSelCoef", "num_qtns", "inv_length"), 
                 mean, .groups = "keep") %>%
    rename(inv_ageB90 = inv_age, mean_qtnSelCoefB90 = mean_qtnSelCoef, 
           numqtnsB90 = num_qtns, inv_lengthB90 = inv_length)
  
   ## Join dataframes with parameters
  df.FSTsplitTb <- full_join(top10.data, bottom90.data, by = c("seed", "sim_gen"))
  df.FSTsplitTb <- left_join(df.FSTsplitTb, df.params, by = "seed")  
  
   ## convert to data frame and factor parameter columns
  df.FSTsplit <- as.data.frame(df.FSTsplitTb)
  for(i in 11:(ncol(df.FSTsplit))){
    df.FSTsplit[,i] <- as.factor(as.character(df.FSTsplit[,i]))
  }
   ## get average values across reps
  df.aveFSTsplit <- aggregate(cbind(inv_ageT10, mean_qtnSelCoefT10, num_qtnsT10, 
                                    inv_lengthT10, inv_ageB90, mean_qtnSelCoefB90,
                                    numqtnsB90, inv_lengthB90)~
                                muBase+sigmaK+muInv+alpha+mig1+enVar+sim_gen, data = df.FSTsplit, FUN = mean)

  ## check size is correct
  length(unique(paste(df.invAllData$seed, df.invAllData$sim_gen)))
      
#############################################################################################################
#### PLOTTING ####
#############################################################################################################
  
  ## make labels for facet wrapping
    alpha.labels <- c("0.2" = "sigmaMu = 0.2", "0.002" = "sigmaMu = 0.002")
    mig.labels <- c("0.001" = "mig = 0.001", "0.01" = "mig = 0.01", "0.1" = "mig = 0.1", 
                  "0.25" = "mig = 0.25", "0.4" = "mig = 0.4", "0.5" = "mig = 0.5")


  ## plot inversion age --  split by inversiion FST with reps
    ggplot(data = df.FSTsplit[df.FSTsplit$enVar == 0 & df.FSTsplit$muInv == 0.001 & 
                                df.FSTsplit$sigmaK == 0.75 & df.FSTsplit$alpha == 0.002,], 
           aes(x = sim_gen, y = inv_ageT10, group = interaction(muBase, rep))) + 
           geom_line(aes(color = muBase), size = 0.75, alpha = 0.3) + 
           geom_line(data = df.aveFSTsplit[df.aveFSTsplit$enVar == 0 & df.aveFSTsplit$muInv == 0.001 & 
                                        df.aveFSTsplit$sigmaK == 0.75 & df.aveFSTsplit$alpha == 0.002,], 
                     aes(x = sim_gen, y = inv_ageT10, group = muBase, color = muBase), size = 1.2) + 
           facet_wrap(~ mig1 , labeller = labeller(mig1 = mig.labels), ncol = 6, nrow = 1) + 
           labs(title = "Average Inversions Age Through Time",
                 y = "Average Inversion Age",
                 x = "Generation") +
           guides(color = guide_legend(title = "QTN Mutation Rate")) +
           theme_classic() +
           theme(panel.background = element_blank(), 
                 strip.background = element_rect(colour = "white", fill = "grey92")) +
           scale_color_manual(values=c( "cadetblue3", "navy")) +
           scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
           scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.FSTsplit$inv_ageT10)))
    
  ## plot inversion age -- split by inversion FST just average
    ggplot(data = df.aveFSTsplit[df.aveFSTsplit$enVar == 0 & df.aveFSTsplit$muInv == 0.001 & 
                              df.aveFSTsplit$sigmaK == 0.75 & df.aveFSTsplit$alpha == 0.002,], 
           aes(x = sim_gen, y = inv_ageT10)) + 
      geom_line(aes(color = mig1), size = 0.75, alpha = 0.95) + 
      geom_line(data = df.aveFSTsplit[df.aveFSTsplit$enVar == 0 & df.aveFSTsplit$muInv == 0.001 & 
                                      df.aveFSTsplit$sigmaK == 0.75 & df.aveFSTsplit$alpha == 0.002,], 
                 aes(x = sim_gen, y = inv_ageB90, color = mig1), 
                 size = 0.75, alpha = 0.95, linetype = 3) + 
      facet_wrap(~ muBase + mig1 , labeller = labeller(mig1 = mig.labels), ncol = 6, nrow = 2) + 
      labs(title = "Average Inversions Age Through Time",
           y = "Average Inversion Age",
           x = "Generation") +
      guides(color = guide_legend(title = "Migration Rate")) +
      scale_color_manual(values=c( "darkgrey","cadetblue1","cadetblue3", 
                                   "cornflowerblue", "navy", "black")) +
      theme_classic() +
      theme(panel.background = element_blank(), 
            strip.background = element_rect(colour = "white", fill = "grey92")) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.FSTsplit$inv_ageT10)))
    
    ggplot(data = df.aveFSTsplit[df.aveFSTsplit$enVar == 0 & df.aveFSTsplit$muInv == 0.001 & 
                                   df.aveFSTsplit$sigmaK == 0.75 & df.aveFSTsplit$alpha == 0.002,], 
           aes(x = sim_gen, y = inv_ageT10)) + 
      geom_line(aes(color = muBase), size = 0.75, alpha = 0.95) + 
      geom_line(data = df.aveFSTsplit[df.aveFSTsplit$enVar == 0 & df.aveFSTsplit$muInv == 0.001 & 
                                        df.aveFSTsplit$sigmaK == 0.75 & df.aveFSTsplit$alpha == 0.002,], 
                aes(x = sim_gen, y = inv_ageB90, color = muBase), 
                size = 0.75, alpha = 0.95, linetype = 3) + 
      facet_wrap(~ muBase + mig1 , labeller = labeller(mig1 = mig.labels), ncol = 6, nrow = 2) + 
      labs(title = "Average Inversions Age Through Time",
           y = "Average Inversion Age",
           x = "Generation") +
      guides(color= guide_legend(title = "QTN Mutation Rate")) +
      scale_color_manual(values=c( "cadetblue3", "dodgerblue3")) +
      theme_classic() +
      theme(panel.background = element_blank(), 
            strip.background = element_rect(colour = "white", fill = "grey92")) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.aveFSTsplit$inv_ageT10)))
   
  ## plot all data         
    inv.Age <-  ggplot(data = df.all.data[df.all.data$enVar == 0 & df.all.data$muInv == 0.001 &
                                          df.all.data$sigmaK == 0.75 & df.all.data$alpha == 0.002,], 
                       aes(x = sim_gen, y = aveAge, group = interaction(muBase, rep))) + 
      geom_line(aes(color = muBase), size = 0.75, alpha = 0.3) + 
      geom_line(data = df.ave.data[df.ave.data$enVar == 0 & df.ave.data$muInv == 0.001 & 
                                   df.ave.data$sigmaK == 0.75 & df.ave.data$alpha == 0.002,], 
                aes(x = sim_gen, y = aveAge, group = muBase, color = muBase), 
              size = 1.2) + 
      facet_wrap(~ mig1 , labeller = labeller(mig1 = mig.labels),
                 ncol = 6, nrow = 1) + 
      labs(title = "Average Inversions Age Through Time",
           y = "Average Inversion Age",
           x = "Generation") +
      guides(color = guide_legend(title = "QTN Mutation Rate")) +
      theme_classic() +
      theme(panel.background = element_blank(), 
            strip.background = element_rect(colour = "white", fill = "grey92")) +
      scale_color_manual(values=c( "cadetblue3", "navy")) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.ave.data$aveAge)))
    
 ### Inversion Length ####
 #######################################
    
    inv.Length <- ggplot(data = df.all.data[df.all.data$enVar == 0 & df.all.data$muInv == 0.001 &
                                               df.all.data$sigmaK == 0.75 & df.all.data$alpha == 0.002,], 
                          aes(x = sim_gen, y = aveLength, group = interaction(muBase, rep))) + 
      geom_line(aes(color = muBase), size = 0.75, alpha = 0.3) + 
      geom_line(data = df.ave.data[df.ave.data$enVar == 0 & df.ave.data$muInv == 0.001 &
                                     df.ave.data$sigmaK == 0.75 & df.ave.data$alpha == 0.002,], 
                aes(x = sim_gen, y = aveLength, group = interaction(muBase), color = muBase), size = 1.2) + 
      facet_wrap(~ mig1 , labeller = labeller(mig1 = mig.labels),
                 ncol = 6, nrow = 1) + 
      labs(title = "Average Inversion Length Through Time",
           y = "Average Inversion Length",
           x = "Generation") +
      guides(color = guide_legend(title = "QTN Mutation Rate")) +
          #   linetype = guide_legend(title = "Sigma Mu (Effect Size)")) +
      theme_classic() +
      theme(panel.background = element_blank(), 
            strip.background = element_rect(colour = "white", fill = "grey92")) +
      scale_color_manual(values=c( "cadetblue3", "navy")) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.ave.data$aveLength)))
    
    ## plot inversion length -- split by inversion FST just average
    ggplot(data = df.aveFSTsplit[df.aveFSTsplit$enVar == 0 & df.aveFSTsplit$muInv == 0.001 & 
                                   df.aveFSTsplit$sigmaK == 0.75 & df.aveFSTsplit$alpha == 0.002,], 
           aes(x = sim_gen, y = inv_lengthT10)) + 
      geom_line(aes(color = mig1), size = 0.75, alpha = 0.95) + 
      geom_line(data = df.aveFSTsplit[df.aveFSTsplit$enVar == 0 & df.aveFSTsplit$muInv == 0.001 & 
                                        df.aveFSTsplit$sigmaK == 0.75 & df.aveFSTsplit$alpha == 0.002,], 
                aes(x = sim_gen, y = inv_lengthB90, color = mig1), 
                size = 0.75, alpha = 0.95, linetype = 3) + 
      facet_wrap(~muBase + mig1, labeller = labeller(mig1 = mig.labels), ncol = 6, nrow = 2) + 
      labs(title = "Average Inversions Length Through Time",
           y = "Average Inversion Length",
           x = "Generation") +
      guides(color = guide_legend(title = "Migration Rate")) +
      scale_color_manual(values=c( "darkgrey","cadetblue1","cadetblue3", 
                                   "cornflowerblue", "navy", "black")) +
      theme_classic() +
      theme(panel.background = element_blank(), 
            strip.background = element_rect(colour = "white", fill = "grey92")) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.FSTsplit$inv_lengthT10)))
    
    
    ggplot(data = df.aveFSTsplit[df.aveFSTsplit$enVar == 0 & df.aveFSTsplit$muInv == 0.001 & 
                                   df.aveFSTsplit$sigmaK == 0.75 & df.aveFSTsplit$alpha == 0.002,], 
           aes(x = sim_gen, y = num_qtnsT10)) + 
      geom_line(aes(color = mig1), size = 0.75, alpha = 0.95) + 
      geom_line(data = df.aveFSTsplit[df.aveFSTsplit$enVar == 0 & df.aveFSTsplit$muInv == 0.001 & 
                                        df.aveFSTsplit$sigmaK == 0.75 & df.aveFSTsplit$alpha == 0.002,], 
                aes(x = sim_gen, y = numqtnsB90, color = mig1), 
                size = 0.75, alpha = 0.95, linetype = 3) + 
      facet_wrap(~muBase + mig1, labeller = labeller(mig1 = mig.labels), ncol = 6, nrow = 2) + 
      labs(title = "Average Number of Inversion QTNs Through Time",
           y = "Average Number of Inversion QTNs",
           x = "Generation") +
      guides(color = guide_legend(title = "Migration Rate")) +
      scale_color_manual(values=c( "darkgrey","cadetblue1","cadetblue3", 
                                   "cornflowerblue", "navy", "black")) +
      theme_classic() +
      theme(panel.background = element_blank(), 
            strip.background = element_rect(colour = "white", fill = "grey92")) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, 300))
   
### Inversion QTNs ####
#######################################   
    
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
    