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
  folderIn <- "results/Inversion/20210719_fullSummaryData/noInvControls/"#args[1] # 
  folderOut <- "figures/20210719_fullSummaryData/noInvControls/" #args[2] #
  seed <-"3385445"  #   args[3] # 
  
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

  
qtns <- df.muts[df.muts$type =="m2" & df.muts$FST != "NAN",]  
head(qtns[order(-as.numeric(as.character(qtns$FST))),], n=10)

manh.plot <- ggplot(qtns, aes(x = position, y = as.numeric(as.character(FST)))) +
  geom_point() +
  theme_classic() +
  labs(y = "FST") +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92"),
        text = element_text(size = 11)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.7))

df.muts$FST <- as.numeric(as.character(df.muts$FST))
freq <- df.muts$freq
altFreq <- 1-freq
df.calc <- cbind(freq, altFreq)
df.muts$MAF <- apply(df.calc, 1, min)
df.muts.MAF <- subset(df.muts, subset = MAF >= 0.01)
options(scipen = 999)
chrom_num <- 21
chrom_len <-  100000

df.muts.MAF$chrom <- NA
for(i in 1:chrom_num){
  if(i != chrom_num){
    chrom_end <- as.numeric(paste0(i, "00000"))
    chrom_start <- chrom_end - chrom_len
  } else {
    chrom_end <- chrom_num*chrom_len 
    chrom_start <- chrom_end - chrom_len
  }

  ## SELECTION
  for(j in 1:nrow(df.muts.MAF)){
    if(i != chrom_num & df.muts.MAF$position[j] <= chrom_end & df.muts.MAF$position[j] > chrom_start){
      df.muts.MAF$chrom[j] <- i
    } else if(i == chrom_num & df.muts.MAF$position[j] > chrom_start){
      df.muts.MAF$chrom[j] <- i
    }
  }
}  
head(df.muts.MAF)
df.muts.MAF$chrom <- as.factor(df.muts.MAF$chrom)
manh.plot <- ggplot(df.muts.MAF, aes(x = position, y = FST, 
                                        group =  chrom)) +
  geom_point(data = df.muts.MAF, aes(color = chrom, shape = chrom)) + 
  scale_shape_manual(name = "Chromosome", values = c(rep(21, 20), 15)) +
  scale_color_manual(name = "Chromosome", values = c(rep(c("navy", "lightblue"), 10), "darkgrey")) + 
  labs(title = expression(bold("No Inversion Control"))) + 
  theme_classic() +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92"),
        text = element_text(size = 11)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))
png(paste0(folderOut, seed, "_manh.png"), width = 7, height = 4, units = 'in', res = 300)
#ggarrange(manh.plot, manh.plot,  manh.plot, ncol = 1, nrow =3)
manh.plot
dev.off()
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
  #par(mfrow= c(1,2), mar = c(2,2,3,2), oma = c(1.5,1.5,3,0))
  df.invQTNs <- subset(df.muts.MAF, subset = df.muts.MAF$inOut == "in")
  df.invQTNs$FST <- as.numeric(as.character(df.invQTNs$FST))
  #hist(df.invQTNs$FST, main = "Selection",
       #xlab = "FST") 
  
  df.invQTNs.NS <- subset(df.muts.NS.MAF, subset = df.muts.NS.MAF$inOut == "in")
  df.invQTNs.NS$FST <- as.numeric(as.character(df.invQTNs.NS$FST))
  #hist(df.invQTNs.NS$FST, main = "No Selection",
       #xlab = "FST") 
 # mtext(expression(bold("Inversion Window QTN FST values")), side = 3, outer = TRUE, cex = 1.5)
 # mtext(expression(bold("FST")), side = 1, outer = TRUE)
 # mtext(expression(bold("Frequency")), side = 2, outer = TRUE)
  
  # create null distributions for outlier criteria
  null <- df.muts.NS.MAF$FST[!is.nan(df.muts.NS.MAF$FST) & df.muts.NS.MAF$type == "m2"] # criteria 1: compare to a null distribution of all QTNs in no-selec sim
  #null <- df.muts.NS.MAF$FST[df.muts.NS.MAF$inOut == "in"]
  null_neut <- df.muts.MAF$FST[df.muts.MAF$type == "m1"] # criteria 2: compare to a null distribution of neutral QTNs in selec sim
  
  # subset for just qtns
  df.qtnMuts.MAF <- df.muts.MAF[df.muts.MAF$type == "m2" | df.muts.MAF$type == "m1",]
  df.qtnMuts.NS.MAF <- df.muts.NS.MAF[df.muts.NS.MAF$type == "m2" | df.muts.NS.MAF$type == "m1",]
  
  # criteria calculation loop 
  df.qtnMuts.MAF$crit1_p.value <- NULL
  df.qtnMuts.MAF$crit2_p.value <- NULL
  df.qtnMuts.MAF$crit3_Va <- NULL
  for(i in 1:nrow(df.qtnMuts.MAF)){
    if(!is.nan(df.qtnMuts.MAF$FST[i])){
      obs <- df.qtnMuts.MAF$FST[i]
      df.qtnMuts.MAF$crit1_p.value[i] <- 1-rank(c(null, obs))[length(null)+1]/(length(null)+1)
      df.qtnMuts.MAF$crit2_p.value[i] <- 1-rank(c(null_neut, obs))[length(null_neut)+1]/(length(null_neut)+1)
      df.qtnMuts.MAF$crit3_Va[i] <- (df.qtnMuts.MAF$selCoef[i]^2)*df.qtnMuts.MAF$freq[i]*(1-df.qtnMuts.MAF$freq[i])
    } else {
      df.qtnMuts.MAF$crit1_p.value[i] <- NA
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
                  geom_histogram(position = "identity", alpha = 0.8, bins = 30) + 
                  scale_fill_manual(values = c("red", "blue")) +
                  labs(title =  "Criteria 1 - compare to null of all \nQTNs no-selection simulation",
                       x = "empirical p-value") +
                  theme(legend.position = "none")
    
  
  # plot for criteria 2
  crit2.plot <- ggplot(df.qtnMuts.MAF[df.qtnMuts.MAF$type == "m2",], aes(x = crit2_p.value, fill = inOut)) +
                  geom_histogram(position = "identity", alpha = 0.8, bins = 30) + 
                  scale_fill_manual(values = c("red", "blue")) +
                  labs(title =  "Criteria 2 - compare to null of neutral \nQTNs selection simulation",
                       x = "empirical p-value") +
                  theme(legend.position = "none") 
  
  # plot Va percent
  crit3.plot.leg <- ggplot(df.qtnMuts.MAF[df.qtnMuts.MAF$type == "m2",], aes(x = crit3_Va_perc, fill = inOut)) +
                       geom_histogram(position = "identity", alpha = 0.8, bins = 30) + 
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
                    geom_histogram(position = "identity", alpha = 0.8, bins = 30) + 
                    scale_fill_manual(values = c("red", "blue")) +
                    labs(title = "Criteria 3 - compare % Va within \nselection simulation subset for > 0.01",
                         x = "percent of Va explained") +
                    theme(legend.position = "none")
    
    
   crit.leg <- g_legend(crit3.plot.leg)
  
    pdf(paste0(folderOut, seed, "_adaptInvCriteria.pdf"), height = 5, width = 15)
      print(ggarrange(crit1.plot, crit2.plot, crit3.subplot, crit3.plot.leg, ncol = 4, widths = c(2.3,2.3,2.3,0.8)))
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
  center.bases <- NULL
  first.bases <- NULL
  final.bases <- NULL
  nonadapt.inv <- NULL
  center.bases.NA <- NULL
  first.bases.NA <- NULL
  final.bases.NA <- NULL
  for(i in 1:length(inv.IDs)){
    focalWindow <- invWindBases.MAF[invWindBases.MAF$invIDMAF == inv.IDs[i],]
    wind.length <- nrow(focalWindow)
    first.wind.base <- focalWindow$invWindBasesMAF[1]
    final.wind.base <- focalWindow$invWindBasesMAF[wind.length]
    center.base <-  (first.wind.base + final.wind.base)/2
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
      center.bases <- c(center.bases, center.base)
      first.bases <- c(first.bases, first.wind.base)
      final.bases <- c(final.bases, final.wind.base)
    } else {
      nonadapt.inv <- c(nonadapt.inv, inv.IDs[i])
      center.bases.NA <- c(center.bases.NA, center.base)
      first.bases.NA <- c(first.bases.NA, first.wind.base)
      final.bases.NA <- c(final.bases.NA, final.wind.base)
    }
  }
  df.crit.output <- as.data.frame(crit.output)
  df.crit.output$crit1_thres <- non.inv.crit1
  df.crit.output$crit2_thres <- non.inv.crit2
  df.crit.output$crit3_thres <- non.inv.crit3 
  
  df.crit.output
  inQTNs <- df.qtnMuts.MAF[df.qtnMuts.MAF$inOut == "Inside Inversion",]
  Va_perc_In <- sum(inQTNs$crit3_Va_perc)
  
  total.perc <- sum(Va_perc_In,Va_perc_Out) 
  total.perc 
  #sanity check
  
  inv.IDs.NS <- unique(invWindBasesNS.MAF$invIDMAF)
  center.bases.NS <- NULL
  first.bases.NS <- NULL
  final.bases.NS <- NULL
  for(i in 1:length(inv.IDs.NS)){
    focalWindow <- invWindBasesNS.MAF[invWindBasesNS.MAF$invIDMAF == inv.IDs.NS[i],]
    wind.length <- nrow(focalWindow)
    first.bases.NS <- c(first.bases.NS,focalWindow$invWindBasesMAF[1])
    final.bases.NS <- c(final.bases.NS,focalWindow$invWindBasesMAF[wind.length])
    center.bases.NS <-  c(center.bases.NS, (focalWindow$invWindBasesMAF[1] + focalWindow$invWindBasesMAF[wind.length])/2)
  }
 df.inv.NS <- data.frame(invWindID = inv.IDs.NS, first_base = first.bases.NS, final_base = final.bases.NS,
                            center_base = center.bases.NS)

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

  ## Nonadaptive Inversions
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

  # bind together datasets for boxplots 
  df.adaptSplitbox <- rbind(as.data.frame(adapt.inv.data.nosum), as.data.frame(non.adapt.inv.data.nosum))
  
  # Standard deviation
  if(length(adapt.inv) != 0){
    
  sd.Adapt <- aggregate(cbind(inv_age, mean_qtnSelCoef, num_qtns, inv_length, num_qtns_Lscaled)~sim_gen, 
                        data = adapt.inv.data.nosum, FUN = sd)
  colnames(sd.Adapt)[2:6] <- c("sd_inv_ageAdapt", "sd_qtnSelCoefAdapt", "sd_num_qtnsAdapt", 
                               "sd_inv_lengthAdapt", "sd_num_qtns_LscaledAdapt")
  } else {
    sd.Adapt <- as.data.frame(matrix(NA, ncol = 5, nrow = length(unique(non.adapt.inv.data$sim_gen))))
    sd.Adapt <- cbind(unique(non.adapt.inv.data$sim_gen), sd.Adapt)
    colnames(sd.Adapt)[1:6] <- c("sim_gen", "sd_inv_ageAdapt", "sd_qtnSelCoefAdapt", "sd_num_qtnsAdapt", 
                                   "sd_inv_lengthAdapt", "sd_num_qtns_LscaledAdapt")
  }
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
  
} ## CLOSE IF STATEMENT FOR IF THERE ARE NO INVERSIONS IN THE FINAL GENERATION
## end subset inversions
######################################################################################################
} ## CLOSE IF STATEMENT FOR IF THIS IS A NO INVERSION SIM

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

  
  pdf(paste0(folderOut, seed, "_LA.pdf"), height = 5, width = 7)
    LA.plot + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  dev.off()

  ## LA final for OUTPUT
  LA_final <- df.popDyn$localAdaptSA[nrow(df.popDyn)]

## Phenotypes ##
  df.pheno <- pivot_longer(df.popDyn[, c(1,10,14)], cols = c(meanPhenoP1, meanPhenoP2),
                            names_to = "pop", values_to = "meanPheno")
  pheno.plot.leg <- ggplot(data = df.pheno, 
                        aes(x = sim_gen, y = meanPheno, group = pop)) + 
    geom_line(aes(color = pop, linetype = pop), size = 0.75) + 
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
                       labels = c("Pop 1", "Pop 2")) +
    scale_linetype_manual(name = "Population",
                          values = c("solid", "dotted"),
                          labels = c("Pop 1", "Pop 2"))

  min_pheno <- min(c(df.popDyn$meanPhenoP2 - df.popDyn$sdPhenoP2 - 0.1,df.popDyn.NS$meanPhenoP2 - df.popDyn.NS$sdPhenoP2 - 0.1))
  max_pheno <- max(c(df.popDyn$meanPhenoP1 + df.popDyn$sdPhenoP1 + 0.1, df.popDyn.NS$meanPhenoP2 - df.popDyn.NS$sdPhenoP2 - 0.1))
  
  pheno.plot <- ggplot(data = df.popDyn, 
                         aes(x = sim_gen, y = meanPhenoP1)) + 
    geom_line(color = "cadetblue3", linetype = "solid", size = 0.75) + 
    geom_line(aes(y = meanPhenoP2, x = sim_gen), color = "navy", linetype = "dotted") + 
    geom_ribbon(aes(ymin= meanPhenoP1 - sdPhenoP1, ymax= meanPhenoP1 + sdPhenoP1), fill = "cadetblue3", alpha=0.2) +
    geom_ribbon(aes(ymin= meanPhenoP2 - sdPhenoP2, ymax= meanPhenoP2 + sdPhenoP2), fill = "navy", alpha=0.2) +
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = expression(bold("Selection")), y = "Phenotype", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(min_pheno, max_pheno)) + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 

  pheno.plot.NS <- ggplot(data = df.popDyn.NS, 
                            aes(x = sim_gen, y = meanPhenoP1)) + 
    geom_line(size = 0.75, color = "cadetblue3", linetype = "solid") + 
    geom_ribbon(aes(ymin= meanPhenoP1 - sdPhenoP1, ymax= meanPhenoP1 + sdPhenoP1), fill = "cadetblue3", alpha=0.2) +
    geom_line(aes(y = meanPhenoP2, x = sim_gen), color = "navy", linetype = "dotted") + 
    geom_ribbon(aes(ymin= meanPhenoP2 - sdPhenoP2, ymax= meanPhenoP2 + sdPhenoP2), fill = "navy", alpha=0.2) +
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = expression(bold("No Selection")), y = " ", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(min_pheno, max_pheno)) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  

  pheno.leg <- g_legend(pheno.plot.leg)

  pdf(paste0(folderOut, seed, "_pheno.pdf"), height = 5, width = 7)
    ggarrange(pheno.plot, pheno.plot.NS, pheno.leg, ncol = 3, widths = c(2.3,2.3,0.8))
  dev.off()

## Fitnesses ##
  
  fitP1.plot <- ggplot(data = df.popDyn, 
                       aes(x = sim_gen, y = meanFitP1)) + 
    geom_line(color = "cadetblue3", size = 0.75, linetype = "solid") + 
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
    geom_line(color = "navy", size = 0.75, linetype = "dotted") + 
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

## IF NO INVERSION SIM DO NOT EVALUATE THE REST OF THE CODE OTHER THAN OUTPUTS
if(df.params$muInv != 0){
  if(nrow(df.invDataFinalGen) > 0 & nrow(df.invDataFinalGen.NS) > 0){
## Average Inversion Age ##

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
  
  ## change the second (no selection) line to dotted line 
  inv.age.plot <- ggplot(data = df.invage, 
                         aes(x = sim_gen, y = inv_age, group = Adaptsplit)) + 
    geom_line(data = df.invage, aes(color = Adaptsplit, linetype = Adaptsplit), alpha = 0.9) + 
    geom_ribbon(data = df.AdaptSplit, aes(x = sim_gen, ymin=  inv_ageAdaptLower, 
                                           ymax= inv_ageAdapt + sd_inv_ageAdapt), 
               fill = inferno(4)[3], alpha=0.2, inherit.aes = FALSE) +
    geom_ribbon(data = df.AdaptSplit, aes(x = sim_gen, ymin = inv_ageNonAdaptLower,
                                           ymax = inv_ageNonAdapt + sd_inv_ageNonAdapt),
                fill = inferno(4)[2], alpha=0.2, inherit.aes = FALSE) +
    geom_ribbon(data = df.inv.data.NS, aes(x = sim_gen, ymin =  inv_ageNSLower,
                                          ymax = inv_age_NS + sd_inv_age_NS),
               fill = inferno(4)[1], alpha=0.2, inherit.aes = FALSE) +
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = " ", y = "Average Inversion Age", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_color_manual(name = "", labels = c("Adaptive", "Nonadaptive", "No Selection"), 
                       values=inferno(4)[3:1]) +
    scale_linetype_manual(name = "", labels = c("Adaptive", "Nonadaptive","No Selection"), 
                          values = c("solid", "solid", "dotted")) + 
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  
  pdf(paste0(folderOut, seed, "_invAge.pdf"), height = 5, width = 7)
   print(inv.age.plot)
  dev.off()

  
  # Subset for final generation
  final.inv <- df.adaptSplitboxNS[df.adaptSplitboxNS$sim_gen == 50000,]
  final.inv$adaptInv <- factor(final.inv$adaptInv, levels = c("Adaptive", "Nonadaptive", "No selection"))
  
  if(length(adapt.inv) != 0){
    color_scale <- inferno(4)[3:1]
  } else {
    color_scale <- inferno(4)[2:1]
  }
  
  pdf(paste0(folderOut, seed, "_invAgebox.pdf"), height = 5, width = 5)
    print(ggplot(data = final.inv, aes(x = adaptInv, y= inv_age, fill = adaptInv)) +
      geom_boxplot() + 
      scale_fill_manual(values = color_scale) + 
      theme_classic() +
      theme(legend.position = "none") + 
      theme(panel.background = element_blank(), 
            strip.background = element_rect(colour = "white", fill = "grey92")) +
      labs(title = " ", y = "Average Inversion Age", x = "Inversion Status") )
  dev.off()
  
  ## OUTPUT: dataframe of values split by adaptive and nonadaptive inversions
  df.invChar <- final.inv[,c("inv_id", "inv_age", "inv_length", "num_qtns_Lscaled", "adaptInv")]
  
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

  inv.length.plot <- ggplot(data = df.invlength, 
                            aes(x = sim_gen, y = inv_length, group = Adaptsplit)) + 
    geom_line(aes(color = Adaptsplit, linetype = Adaptsplit), alpha = 0.9) + 
    geom_ribbon(data = df.AdaptSplit, aes(x = sim_gen, ymin=  inv_lengthAdaptLower, 
                                            ymax= inv_lengthAdapt + sd_inv_lengthAdapt), 
                  fill = inferno(4)[3], alpha=0.2, inherit.aes = FALSE) +
    geom_ribbon(data = df.AdaptSplit, aes(x = sim_gen, ymin = inv_lengthNonAdaptLower,
                                          ymax = inv_lengthNonAdapt + sd_inv_lengthNonAdapt),
                fill = inferno(4)[2], alpha=0.2, inherit.aes = FALSE) +
    geom_ribbon(data = df.inv.data.NS, aes(x = sim_gen, ymin =  inv_lengthNSLower,
                                           ymax = inv_length_NS + sd_inv_length_NS),
                fill = inferno(4)[1], alpha=0.2, inherit.aes = FALSE) +
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = " ", y = "Average Inversion Length", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_color_manual(name = "", labels = c( "Adaptive", "Nonadaptive","No Selection"),
                       values=inferno(4)[3:1]) +
    scale_linetype_manual(name = "", labels = c("Adaptive", "Nonadaptive","No Selection"), 
                          values = c("solid", "solid", "dotted")) + 
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 50000)) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  
  pdf(paste0(folderOut, seed, "_invLength.pdf"), height = 5, width = 7)
    print(inv.length.plot)
  dev.off()

  pdf(paste0(folderOut, seed, "_invLengthbox.pdf"), height = 5, width = 5)
  print(ggplot(data = final.inv, aes(x = adaptInv, y= inv_length, fill = adaptInv)) +
    geom_boxplot() + 
    scale_fill_manual(values = color_scale) + 
    theme_classic() +
    theme(legend.position = "none") + 
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    labs(title = " ", y = "Average Inversion Length", x = "Inversion Status") )
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
  
  inv.qtns.Lscaled.plot <- ggplot(data = df.invQTNsLscaled, 
                                  aes(x = sim_gen, y = inv_numQTNs, group = Adaptsplit)) + 
    geom_line(aes(color = Adaptsplit, linetype = Adaptsplit), alpha = 0.9) + 
    geom_ribbon(data = df.AdaptSplit, aes(x = sim_gen, ymin= inv_numQTNsLscaledAdaptLower, 
                                          ymax= num_qtns_LscaledAdapt + sd_num_qtns_LscaledAdapt), 
                fill = inferno(4)[3], alpha=0.2, inherit.aes = FALSE) +
    geom_ribbon(data = df.AdaptSplit, aes(x = sim_gen, ymin = inv_numQTNsLscaledNonAdaptLower,
                                          ymax = num_qtns_LscaledNonAdapt + sd_num_qtns_LscaledNonAdapt),
                fill = inferno(4)[2], alpha=0.2, inherit.aes = FALSE) +
    geom_ribbon(data = df.inv.data.NS, aes(x = sim_gen, ymin =  inv_numQTNsLscaledNSLower,
                                           ymax = num_qtns_Lscaled_NS + sd_num_qtns_Lscaled_NS),
                fill = inferno(4)[1], alpha=0.2, inherit.aes = FALSE) +
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = "",
         y = "Average number of inversion QTNs \nscaled by inversion length",
         x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_color_manual(name = "", labels = c( "Adaptive", "Nonadaptive","No Selection"),
                       values=inferno(4)[3:1]) +
    scale_linetype_manual(name = "", labels = c("Adaptive", "Nonadaptive","No Selection"), values = c("solid", "solid", "dotted")) + 
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.invQTNsLscaled$inv_numQTNs))) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  
  pdf(paste0(folderOut, seed, "_invQTNsLscaled.pdf"), height = 5, width = 7)
    print(inv.qtns.Lscaled.plot)
  dev.off()
  
  pdf(paste0(folderOut, seed, "_invQTNsLscaledbox.pdf"), height = 5, width = 5)
    print(ggplot(data = final.inv, aes(x = adaptInv, y= num_qtns_Lscaled, fill = adaptInv)) +
      geom_boxplot() + 
      scale_fill_manual(values = color_scale) + 
      theme_classic() +
      theme(legend.position = "none") + 
      theme(panel.background = element_blank(), 
            strip.background = element_rect(colour = "white", fill = "grey92")) +
      labs(title = " ", y = "Average number of inversion QTNs \nscaled by inversion length", x = "Inversion Status"))
  dev.off()
### end initial full data plotting
######################################################################################################
  
  
  
 
######################################################################################################  
#### Add chromosome number ####
  
  options(scipen = 999)
  chrom_num <- 21
  chrom_len <-  100000
  
  df.muts.MAF$chrom <- NA
  df.muts.NS.MAF$chrom <- NA
  for(i in 1:chrom_num){
    if(i != chrom_num){
      chrom_end <- as.numeric(paste0(i, "00000"))
      chrom_start <- chrom_end - chrom_len
    } else {
      chrom_end <- chrom_num*chrom_len 
      chrom_start <- chrom_end - chrom_len
    }
    
    ## SELECTION
    for(j in 1:nrow(df.muts.MAF)){
      if(i != chrom_num & df.muts.MAF$position[j] <= chrom_end & df.muts.MAF$position[j] > chrom_start){
        df.muts.MAF$chrom[j] <- i
      } else if(i == chrom_num & df.muts.MAF$position[j] > chrom_start){
        df.muts.MAF$chrom[j] <- i
      }
    }
    
    ## NO SELECTION
    for(j in 1:nrow(df.muts.NS.MAF)){
      if(i != chrom_num & df.muts.NS.MAF$position[j] <= chrom_end & df.muts.NS.MAF$position[j] > chrom_start){
        df.muts.NS.MAF$chrom[j] <- i
      } else if(i == chrom_num & df.muts.NS.MAF$position[j] > chrom_start){
        df.muts.NS.MAF$chrom[j] <- i
      }
    }
  }
  
 
#### end add chromosome number
######################################################################################################  

  
  
######################################################################################################   
#### plot inversion origin dynamics ####
  # Subset for final generation
  df.invFinalGen <- subset(df.invTime, subset = sim_gen == 50000)
  df.InvDataOrigin <- left_join(df.invFinalGen, df.invData, by = "inv_id")
  df.invFinalGen.NS <- subset(df.invTime.NS, subset = sim_gen == 50000)
  df.InvDataOrigin.NS <- left_join(df.invFinalGen.NS, df.invData.NS, by = "inv_id")
  
  ## SELECTION 
  # subset for MAF > 0.01
  freq <- df.InvDataOrigin$freq
  altFreq <- 1-freq
  df.calc <- cbind(freq, altFreq)
  df.InvDataOrigin$MAF <- apply(df.calc, 1, min)
  df.InvDataOriginMAF <- subset(df.InvDataOrigin, subset = MAF >= 0.01)
  df.InvDataOriginMAF$qtnSelCoefsum <- df.InvDataOriginMAF$mean_qtnSelCoef*df.InvDataOriginMAF$num_qtns
  df.InvDataOriginMAF$evHistory <- NULL
  for(i in 1:nrow(df.InvDataOriginMAF)){
    
    ## identify population that the inversions originated in
    if(df.InvDataOriginMAF$freq_p1[i] > df.InvDataOriginMAF$freq_p2[i]){
      df.InvDataOriginMAF$pop[i] <- "Pop 1"
    } else {
      df.InvDataOriginMAF$pop[i] <- "Pop 2"
    }
    
    # idetnify the evolutionary history of the inversion by documenting whether it was in one of
    # our four categories of capturing QTNs initially or being neutral and whether or not it gained
    # QTNs over time
    if(df.InvDataOriginMAF$invOrigNumQTNs[i] > 0){
      if(df.InvDataOriginMAF$num_qtns[i] - df.InvDataOriginMAF$invOrigNumQTNs[i] == 0){
        df.InvDataOriginMAF$evHistory[i] <- "capture_nogain"
      } else {
        df.InvDataOriginMAF$evHistory[i] <- "capture_gain"
      }
    } else {
      if(df.InvDataOriginMAF$num_qtns[i] - df.InvDataOriginMAF$invOrigNumQTNs[i] == 0){
        df.InvDataOriginMAF$evHistory[i] <- "neutral_nogain"
      } else {
        df.InvDataOriginMAF$evHistory[i] <- "neutral_gain"
      }
    }
    
    # get information about the origin FST values of the inversion to compare to the initial
    df.temp <- df.invTime[df.invTime$inv_id == df.InvDataOriginMAF$inv_id[i],]
    df.first <- df.temp[order(df.temp$sim_gen),][1,]
    df.InvDataOriginMAF$originFST[i] <- df.first$inv_FST
  }
  
  ## NO SELECTION
  freq <- df.InvDataOrigin.NS$freq
  altFreq <- 1-freq
  df.calc <- cbind(freq, altFreq)
  df.InvDataOrigin.NS$MAF <- apply(df.calc, 1, min)
  df.InvDataOriginMAF.NS <- subset(df.InvDataOrigin.NS, subset = MAF >= 0.01)
  df.InvDataOriginMAF.NS$qtnSelCoefsum <- df.InvDataOriginMAF.NS$mean_qtnSelCoef*df.InvDataOriginMAF.NS$num_qtns
  df.InvDataOriginMAF.NS$evHistory <- NULL
  for(i in 1:nrow(df.InvDataOriginMAF.NS)){
    
    ## identify population that the inversions originated in
    if(df.InvDataOriginMAF.NS$freq_p1[i] > df.InvDataOriginMAF.NS$freq_p2[i]){
      df.InvDataOriginMAF.NS$pop[i] <- "Pop 1"
    } else {
      df.InvDataOriginMAF.NS$pop[i] <- "Pop 2"
    }
    
    # idetnify the evolutionary history of the inversion by documenting whether it was in one of
    # our four categories of capturing QTNs initially or being neutral and whether or not it gained
    # QTNs over time
    if(df.InvDataOriginMAF.NS$invOrigNumQTNs[i] > 0){
      if(df.InvDataOriginMAF.NS$num_qtns[i] - df.InvDataOriginMAF.NS$invOrigNumQTNs[i] == 0){
        df.InvDataOriginMAF.NS$evHistory[i] <- "capture_nogain"
      } else {
        df.InvDataOriginMAF.NS$evHistory[i] <- "capture_gain"
      }
    } else {
      if(df.InvDataOriginMAF.NS$num_qtns[i] - df.InvDataOriginMAF.NS$invOrigNumQTNs[i] == 0){
        df.InvDataOriginMAF.NS$evHistory[i] <- "neutral_nogain"
      } else {
        df.InvDataOriginMAF.NS$evHistory[i] <- "neutral_gain"
      }
    }
    
    # get information about the origin FST values of the inversion to compare to the initial
    df.temp <- df.invTime.NS[df.invTime.NS$inv_id == df.InvDataOriginMAF.NS$inv_id[i],]
    df.first <- df.temp[order(df.temp$sim_gen),][1,]
    df.InvDataOriginMAF.NS$originFST[i] <- df.first$inv_FST
  }
  
  ## OUTPUT: Calculate propotions 
  ## adaptive inversions
  capture_gain_p <- sum((df.InvDataOriginMAF$inv_id %in% adapt.inv) & df.InvDataOriginMAF$evHistory == "capture_gain")/
                       sum((df.InvDataOriginMAF$inv_id %in% adapt.inv))
  capture_no_gain_p <- sum((df.InvDataOriginMAF$inv_id %in% adapt.inv) & df.InvDataOriginMAF$evHistory == "capture_nogain")/
                           sum((df.InvDataOriginMAF$inv_id %in% adapt.inv))
  neutral_gain_p <- sum((df.InvDataOriginMAF$inv_id %in% adapt.inv) & df.InvDataOriginMAF$evHistory == "neutral_gain")/
                        sum((df.InvDataOriginMAF$inv_id %in% adapt.inv))
  neutral_no_gain_p <- sum((df.InvDataOriginMAF$inv_id %in% adapt.inv) & df.InvDataOriginMAF$evHistory == "neutral_nogain")/
                            sum((df.InvDataOriginMAF$inv_id %in% adapt.inv))
  
  ## nonadaptive inversions
  capture_gain_p.NA <- sum(!(df.InvDataOriginMAF$inv_id %in% adapt.inv) & df.InvDataOriginMAF$evHistory == "capture_gain")/
                           sum(!(df.InvDataOriginMAF$inv_id %in% adapt.inv))
  capture_no_gain_p.NA <- sum(!(df.InvDataOriginMAF$inv_id %in% adapt.inv) & df.InvDataOriginMAF$evHistory == "capture_nogain")/
                              sum(!(df.InvDataOriginMAF$inv_id %in% adapt.inv))
  neutral_gain_p.NA <- sum(!(df.InvDataOriginMAF$inv_id %in% adapt.inv) & df.InvDataOriginMAF$evHistory == "neutral_gain")/
                           sum(!(df.InvDataOriginMAF$inv_id %in% adapt.inv))
  neutral_no_gain_p.NA <- sum(!(df.InvDataOriginMAF$inv_id %in% adapt.inv) & df.InvDataOriginMAF$evHistory == "neutral_nogain")/
                              sum(!(df.InvDataOriginMAF$inv_id %in% adapt.inv))
  
  ## no selection inversions
  capture_gain_p.NS <- sum(df.InvDataOriginMAF.NS$evHistory == "capture_gain")/nrow(df.InvDataOriginMAF.NS)
  capture_no_gain_p.NS <- sum(df.InvDataOriginMAF.NS$evHistory == "capture_nogain")/nrow(df.InvDataOriginMAF.NS)
  neutral_gain_p.NS <- sum(df.InvDataOriginMAF.NS$evHistory == "neutral_gain")/nrow(df.InvDataOriginMAF.NS)
  neutral_no_gain_p.NS <- sum(df.InvDataOriginMAF.NS$evHistory == "neutral_nogain")/nrow(df.InvDataOriginMAF.NS)
  mat <- matrix(c(capture_gain_p, capture_no_gain_p, neutral_gain_p, neutral_no_gain_p, 
                  capture_gain_p.NA, capture_no_gain_p.NA, neutral_gain_p.NA, neutral_no_gain_p.NA, 
                  capture_gain_p.NS, capture_no_gain_p.NS, neutral_gain_p.NS, neutral_no_gain_p.NS), ncol = 3, nrow = 4)
  
  colnames(mat) <- c("Adaptive", "Nonadaptive", "No Selection")
  rownames(mat) <- c("Capture Gain", "Capture No Gain", "Neutral Gain", "Neutral No Gain")
  
  pdf(paste0(folderOut, seed, "_invOriginBarplot.pdf"), height = 6, width = 5)
    par(mfrow=c(1,1), mar = c(2.5,4,2.5,2.5))
    barplot(mat, col = viridis(5)[4:1], ylim = c(0, 1.2),
          args.legend = list(bty = "n", x = "top", ncol = 4), ylab = "Proportion")
    legend("top", legend = rownames(mat), fill = viridis(5)[4:1], ncol = 2, cex = 0.85)
  dev.off()
  
  ## OUTPUT: for each of the three groups, adaptive, nonadaptive and no selection inversions, output
  # Number of inversions, starting average QTNs, average FST & absolute value of the QTN effects, and then average 
  # ending average QTNs, average FST & absolute value of the QTN effects
  num_inv <- sum(df.InvDataOriginMAF$inv_id %in% adapt.inv)
  num_inv_NA <- sum(!df.InvDataOriginMAF$inv_id %in% adapt.inv)
  num_inv_NS <- nrow(df.InvDataOriginMAF.NS)
  ave_start_QTNs <- mean(df.InvDataOriginMAF$invOrigNumQTNs[df.InvDataOriginMAF$inv_id %in% adapt.inv])
  ave_start_QTNs_NA <- mean(df.InvDataOriginMAF$invOrigNumQTNs[!(df.InvDataOriginMAF$inv_id %in% adapt.inv)])
  ave_start_QTNs_NS <- mean(df.InvDataOriginMAF.NS$invOrigNumQTNs)
  ave_end_QTNs <- mean(df.InvDataOriginMAF$num_qtns[df.InvDataOriginMAF$inv_id %in% adapt.inv])
  ave_end_QTNs_NA <- mean(df.InvDataOriginMAF$num_qtns[!(df.InvDataOriginMAF$inv_id %in% adapt.inv)])
  ave_end_QTNs_NS <- mean(df.InvDataOriginMAF.NS$num_qtns)
  ave_start_FST <- mean(df.InvDataOriginMAF$originFST[df.InvDataOriginMAF$inv_id %in% adapt.inv])
  ave_start_FST_NA <- mean(df.InvDataOriginMAF$originFST[!(df.InvDataOriginMAF$inv_id %in% adapt.inv)])
  ave_start_FST_NS <- mean(df.InvDataOriginMAF.NS$originFST)
  ave_end_FST <- mean(df.InvDataOriginMAF$inv_FST[df.InvDataOriginMAF$inv_id %in% adapt.inv])
  ave_end_FST_NA <- mean(df.InvDataOriginMAF$inv_FST[!(df.InvDataOriginMAF$inv_id %in% adapt.inv)])
  ave_end_FST_NS <- mean(df.InvDataOriginMAF.NS$inv_FST)
  ave_abV_start_qtnSelCoef <- mean(abs(df.InvDataOriginMAF$invOrigQTNSelCoef[df.InvDataOriginMAF$inv_id %in% adapt.inv]))
  ave_abV_start_qtnSelCoef_NA <- mean(abs(df.InvDataOriginMAF$invOrigQTNSelCoef[!(df.InvDataOriginMAF$inv_id %in% adapt.inv)]))
  ave_abV_start_qtnSelCoef_NS <- mean(abs(df.InvDataOriginMAF.NS$invOrigQTNSelCoef))
  ave_abV_end_qtnSelCoef <- mean(abs(df.InvDataOriginMAF$qtnSelCoefsum[df.InvDataOriginMAF$inv_id %in% adapt.inv]))
  ave_abV_end_qtnSelCoef_NA <- mean(abs(df.InvDataOriginMAF$qtnSelCoefsum[!(df.InvDataOriginMAF$inv_id %in% adapt.inv)]))
  ave_abV_end_qtnSelCoef_NS <- mean(abs(df.InvDataOriginMAF.NS$qtnSelCoefsum))
  
  
  ## Subset dataframe to get how the MAF filtered inversions change through time
  ## Then subset for each population for plotting separate color scales
  ## SELECTION
  inv.IDs <- as.vector(df.InvDataOriginMAF$inv_id)
  df.invFinalAllData <- df.invTime[df.invTime$inv_id %in% inv.IDs, ]
  df.invFinalAllData$qtnSelCoefsum <- df.invFinalAllData$mean_qtnSelCoef*df.invFinalAllData$num_qtns
  df.invFinalAllDataPop <- left_join(df.invFinalAllData,
                                     df.InvDataOriginMAF[c(2,13:19,21)], 
                                     by = "inv_id")
  df.invFinalsubset <- df.invFinalAllDataPop %>% filter(sim_gen %in% seq(0, 50000, by = 1000))
  df.invFinalsubset$inv_id <- as.factor(df.invFinalsubset$inv_id)
  df.invFinalsubset$adaptInv <- "Nonadaptive"
  df.invFinalsubset$adaptInv[df.invFinalsubset$inv_id %in% adapt.inv] <- "Adaptive"
  df.invFinalsubset$adaptInv <- as.factor(df.invFinalsubset$adaptInv)
  df.pop1 <- df.invFinalsubset[df.invFinalsubset$pop == "Pop 1",]
  df.pop2 <- df.invFinalsubset[df.invFinalsubset$pop == "Pop 2",]
  df.pop1$pop <- as.factor(df.pop1$pop)
  df.pop2$pop <- as.factor(df.pop2$pop)
  
  ## NO SELECTION
  inv.IDs.NS <- as.vector(df.InvDataOriginMAF.NS$inv_id)
  df.invFinalAllData.NS <- df.invTime.NS[df.invTime.NS$inv_id %in% inv.IDs.NS, ]
  df.invFinalAllData.NS$qtnSelCoefsum <- df.invFinalAllData.NS$mean_qtnSelCoef*df.invFinalAllData.NS$num_qtns
  df.invFinalAllDataPop.NS <- left_join(df.invFinalAllData.NS,
                                     df.InvDataOriginMAF.NS[c(2,13:19,21)], 
                                     by = "inv_id")
  df.invFinalsubset.NS <- df.invFinalAllDataPop.NS %>% filter(sim_gen %in% seq(0, 50000, by = 1000)) 
  df.invFinalsubset.NS$inv_id <- as.factor(df.invFinalsubset.NS$inv_id)
  df.pop1.NS <- df.invFinalsubset.NS[df.invFinalsubset.NS$pop == "Pop 1",]
  df.pop2.NS <- df.invFinalsubset.NS[df.invFinalsubset.NS$pop == "Pop 2",]
  df.pop1.NS$pop <- as.factor(df.pop1.NS$pop)
  df.pop2.NS$pop <- as.factor(df.pop2.NS$pop)


  ## plot using pop as separate dataframe
  colorCountBlue <- length(unique(df.pop1$inv_id[df.pop1$adaptInv == "Adaptive"]))
  getPaletteBlue <-  colorRampPalette(brewer.pal(9, "Blues"))
  colorCountRed <- length(unique(df.pop2$inv_id[df.pop2$adaptInv == "Adaptive"]))
  getPaletteRed <-  colorRampPalette(brewer.pal(9, "Reds"))
  
  max_value_origin <- max(c(abs(df.invFinalsubset$qtnSelCoefsum), abs(df.invFinalsubset.NS$qtnSelCoefsum)))
 
  plot.inv.orig.col <- ggplot(df.pop1[df.pop1$adaptInv == "Adaptive",], aes(x = sim_gen, y = qtnSelCoefsum)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_point(aes(color = inv_id, size = inv_FST, group = interaction(inv_id, inv_FST)), alpha = 0.8) + 
    geom_line(aes(color = inv_id, group =  inv_id), alpha = 0.8) + 
    scale_color_manual(values = getPaletteBlue(colorCountBlue)[colorCountBlue:1]) + 
    new_scale_color() +
    geom_point(data = df.pop2[df.pop2$adaptInv == "Adaptive",], aes(x = sim_gen, y = qtnSelCoefsum, color = inv_id, 
                                   size = inv_FST), alpha = 0.8, inherit.aes = FALSE) + 
    geom_line(data = df.pop2[df.pop2$adaptInv == "Adaptive",], aes(x = sim_gen, y = qtnSelCoefsum, color = inv_id, 
                                  group = inv_id), alpha = 0.8, inherit.aes = FALSE) + 
    scale_color_manual(values = getPaletteRed(colorCountRed)[colorCountRed:1]) + 
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92"),
          text = element_text(size = 15)) +
    theme(axis.line.x = element_line(color="black", size = 0.2),
          axis.line.y = element_line(color="black", size = 0.2)) + 
    labs(title = expression(bold("Selection - Adaptive Inversions")),
         y = " \n ",
         x = " ") +
    ylim(-max_value_origin, max_value_origin) +
    xlim(0, 50000) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
    theme(legend.position = "none")
  
 
  colorCountBlue.NA <- length(unique(df.pop1$inv_id[df.pop1$adaptInv == "Nonadaptive"]))+4
  colorCountRed.NA <- length(unique(df.pop2$inv_id[df.pop2$adaptInv == "Nonadaptive"]))+4
  
  plot.inv.orig.col.NA <- ggplot(df.pop1[df.pop1$adaptInv == "Nonadaptive",], aes(x = sim_gen, y = qtnSelCoefsum)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_point(aes(color = inv_id, size = inv_FST), alpha = 0.8) + 
    geom_line(aes(color = inv_id, group =  inv_id), alpha = 0.8) + 
    scale_color_manual(values = getPaletteBlue(colorCountBlue.NA)[colorCountBlue.NA:1]) + 
    scale_alpha_manual(values = c(0.95, 0.55)) +
    new_scale_color() +
    geom_point(data = df.pop2[df.pop2$adaptInv == "Nonadaptive",], aes(x = sim_gen, y = qtnSelCoefsum, color = inv_id, 
                                   size = inv_FST), alpha = 0.8, inherit.aes = FALSE) + 
    geom_line(data = df.pop2[df.pop2$adaptInv == "Nonadaptive",], aes(x = sim_gen, y = qtnSelCoefsum, color = inv_id, 
                                  group = inv_id), alpha = 0.8, inherit.aes = FALSE) + 
    scale_color_manual(values = getPaletteRed(colorCountRed.NA)[colorCountRed.NA:1]) + 
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92"),
          text = element_text(size = 15)) +
    theme(axis.line.x = element_line(color="black", size = 0.2),
          axis.line.y = element_line(color="black", size = 0.2)) + 
    labs(title = expression(bold("Selection - Nonadaptive Inversions")),
         y = "sum of each Inversion QTNs \neffects on phenotype",
         x = " ") +
    ylim(-max_value_origin, max_value_origin) +
    xlim(0, 50000) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
    theme(legend.position = "none")
  
  colorCountBlue.NS <- length(unique(df.pop1.NS$inv_id))
  colorCountRed.NS <- length(unique(df.pop2.NS$inv_id))
  
  plot.inv.orig.col.NS <- ggplot(df.pop1.NS, aes(x = sim_gen, y = qtnSelCoefsum)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_point(aes(color = inv_id, size = inv_FST), alpha = 0.8) + 
    geom_line(aes(color = inv_id, group = inv_id), alpha = 0.8) + 
    scale_color_manual(values = getPaletteBlue(colorCountBlue.NS)[colorCountBlue.NS:1]) + 
    new_scale_color() +
    geom_point(data = df.pop2.NS, aes(x = sim_gen, y = qtnSelCoefsum, color = inv_id, size = inv_FST), 
               alpha = 0.8, inherit.aes = FALSE) + 
    geom_line(data = df.pop2.NS, aes(x = sim_gen, y = qtnSelCoefsum, color = inv_id, group = inv_id), 
              alpha = 0.8, inherit.aes = FALSE) + 
    scale_color_manual(values = getPaletteRed(colorCountRed.NS)[colorCountRed.NS:1]) + 
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92"),
          text = element_text(size = 15)) +
    theme(axis.line.x = element_line(color="black", size = 0.2),
          axis.line.y = element_line(color="black", size = 0.2)) + 
    labs(title = expression(bold("No Selection")),
         y = " \n ",
         x = "Generation") +
    ylim(-max_value_origin, max_value_origin) +
    xlim(0, 50000) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
    theme(legend.position = "none")
  
  plot.inv.orig <- ggplot(df.invFinalsubset, aes(x = sim_gen, y = qtnSelCoefsum)) +
    geom_point(aes(color = pop, size = inv_FST)) +
    geom_line(aes(color = pop, group = inv_id)) +
    scale_color_manual(values=c("navy", "red")) + 
    theme_classic() +
    theme(panel.background = element_blank(),
          strip.background = element_rect(colour = "white", fill = "grey92"),
          text = element_text(size = 15)) +
    labs(title = expression(bold("Selection")),
         y = " \n ",
         x = "Generation") +
    ylim(-max(df.invFinalsubset$qtnSelCoefsum), max(df.invFinalsubset$qtnSelCoefsum)) +
    xlim(0, 50000) +
    guides(color = guide_legend(title = "Pop with Highest\nFrequency of \nInversion"),
           names = c("Pop 1", "Pop 2")) +
    guides(size = guide_legend(title = "Inversion FST")) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
 
  legInvOrig <- g_legend(plot.inv.orig)
  
  ## blank graph
  blank <- ggplot() + theme_void()
  
  pdf(paste0(folderOut, seed, "_invOrigin.pdf"), height = 12, width = 9)
    print(ggarrange(plot.inv.orig.col, blank, plot.inv.orig.col.NA, legInvOrig, plot.inv.orig.col.NS, blank, nrow = 3, ncol = 2, widths = c(2.3,0.8,2.3,0.8,2.3,0.8)))
  dev.off()
  
#### end plot origin dynamics
######################################################################################################  

######################################################################################################  
#### process VCF ####
## Selection Simulation
  
  # VCF file
  vcf.MAF <- read.vcfR(paste0(folderIn,seed, "_InversionVCF_MAF01.recode.vcf"))
  vcf.MAF.NS <- read.vcfR(paste0(folderIn,seed, "noSel_InversionVCF_MAF01.recode.vcf"))
  
  head(vcf.MAF)
  head(vcf.MAF@fix, 50)
  dim(vcf.MAF@fix)
  
  df.muts.MAF$FST <- as.numeric(as.character(df.muts.MAF$FST))
  df.muts.NS.MAF$FST <- as.numeric(as.character(df.muts.NS.MAF$FST))
  colnames(df.muts.MAF)[2] <- "LocusName"
  colnames(df.muts.NS.MAF)[2] <- "LocusName"
  
  # example of how to find a specific mutation in the vcf file
  df.muts.MAF[2,]
  vcf.MAF@fix[grep(df.muts.MAF$LocusName[1], vcf.MAF@fix[,"INFO"]),]
  
  geno <- vcf.MAF@gt[,-1] # this gets individual ids and genotypes
  geno.NS <- vcf.MAF.NS@gt[,-1]
  position <- getPOS(vcf.MAF) # this gets position of mutations
  position.NS <- getPOS(vcf.MAF.NS)
  chromosome <- getCHROM(vcf.MAF) # Identify chromosome ID for each mutation
  chromosome.NS <- getCHROM(vcf.MAF.NS)
  
  if (sum(duplicated(position)) != 0){
    print("Selection simulation needs to be checked for duplicated locus positions")
  }
  if (sum(duplicated(position.NS)) != 0){
    print("No selection simulation needs to be checked for duplicated locus positions")
  }
  
  # convert the genome data to only have three options 0 (homoz), 1 (heteroz), or 2 (homoz alt)
  G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
  G.NS <- matrix(NA, nrow = nrow(geno.NS), ncol = ncol(geno.NS))
  G[geno %in% c("0/0", "0|0")] <- 0
  G.NS[geno.NS %in% c("0/0", "0|0")] <- 0
  G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
  G.NS[geno.NS  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
  G[geno %in% c("1/1", "1|1")] <- 2
  G.NS[geno.NS %in% c("1/1", "1|1")] <- 2
  
  # calculate allele frequencies
  # sum across rows to get how many copies of the allele then divide by 2 times the number of columns (ind x diploid)
  a_freq <- rowSums(G)/(2*ncol(G))
  a_freq.NS <- rowSums(G.NS)/(2*ncol(G.NS))
  #hist(a_freq) 
  #hist(a_freq.NS) 
  
  # get individual name information from vcf file (e.g., i999) 
  vcf_ind <- data.frame(vcf_ind=colnames(vcf.MAF@gt)[-1])
  vcf_ind.NS <- data.frame(vcf_ind=colnames(vcf.MAF.NS@gt)[-1])
  
  # store mutation metadata which includes mut ID (MID), selection coef (S)
  # dominance (DOM), population origin (PO), ? (GO), mutation type (MT), ? (AC), read depth (DP)
  meta <- vcf.MAF@fix[,"INFO"]
  meta.NS <- vcf.MAF.NS@fix[,"INFO"]
  head(meta)
  length(meta)
  head(meta.NS)
  length(meta.NS)
  
  # identify mutation ids and make sure its the number of muts you expect
  LocusNameOGorder <- regmatches(meta, regexpr("[0-9]+[0-9]", meta))
  LocusNameOGorder.NS <- regmatches(meta.NS, regexpr("[0-9]+[0-9]", meta.NS))
  length(LocusNameOGorder)
  length(LocusNameOGorder.NS)
  
  ## make a data frame where we have the locus and position in the same dataframe (useful later)
  df.OGlocusInfo <- data.frame(LocusName = LocusNameOGorder, position = position, chromosome = chromosome)
  df.OGlocusInfo.NS <- data.frame(LocusName = LocusNameOGorder.NS, position = position.NS, chromosome = chromosome.NS)
  df.ord <- df.OGlocusInfo[order(as.numeric(df.OGlocusInfo$position)),]
  df.ord.NS <- df.OGlocusInfo.NS[order(as.numeric(df.OGlocusInfo.NS$position)),]
  
  # sort mutations and make sure vcf matches slim output
  # slim output is one base off so add a base to match vcf position
  df.muts.MAF$position_vcf <- df.muts.MAF$position + 1
  df.muts.NS.MAF$position_vcf <- df.muts.NS.MAF$position + 1
  df.mutsMAFord <- df.muts.MAF[order(df.muts.MAF$position_vcf),]
  df.mutsMAFord.NS <- df.muts.NS.MAF[order(df.muts.NS.MAF$position_vcf),]
  head(df.mutsMAFord)
  vcf_pos <- as.numeric(vcf.MAF@fix[,"POS"])
  vcf_pos.NS <- as.numeric(vcf.MAF.NS@fix[,"POS"])
  head(vcf_pos)
 
  # 
  vcf_muts <- data.frame(vcf_muts=vcf.MAF@fix[,8])
  vcf_muts.NS <- data.frame(vcf_muts.NS=vcf.MAF.NS@fix[,8])
  colnames(G) <- vcf_ind$vcf_ind # adds individual ids as column names
  colnames(G.NS) <- vcf_ind.NS$vcf_ind # adds individual ids as column names
  rownames(G) <- regmatches(meta, regexpr("[0-9]+[0-9]", meta)) #ADD MUTATION NAMES TO G
  rownames(G.NS) <- regmatches(meta.NS, regexpr("[0-9]+[0-9]", meta.NS)) #ADD MUTATION NAMES TO G
  
  # head(G[,1:5])
  # dim(G)
  # head(G.NS[,1:5])
  # dim(G.NS)
  
  dim(vcf_ind)
  head(vcf_ind)
  head(df.indPheno)
  df.indPheno$vcf_ind <- paste0("i",0:1999) # hard coding
  df.indPheno.NS$vcf_ind <- paste0("i",0:1999) # hard coding
  
  # The individual IDs in Slim do not match the IDs in the VCF file. 
  # I will assume they are in the same order
  tail(df.indPheno)
  tail(df.indPheno.NS)
  
  # Add vcf individual IDs with slim individual IDs and pop numbers 
  indPhen_df_vcf <- merge(vcf_ind, df.indPheno, by="vcf_ind")
  indPhen_df_vcf.NS <- merge(vcf_ind.NS, df.indPheno.NS, by="vcf_ind")
  
  dim(df.indPheno)
  dim(indPhen_df_vcf)
  
  dim(df.indPheno.NS)
  dim(indPhen_df_vcf.NS)
  
  # reorder the dataframe so that it is by subpop first then ind ID
  indPhen_df_vcf <- indPhen_df_vcf[order(indPhen_df_vcf$subpop, indPhen_df_vcf$id),]
  indPhen_df_vcf.NS <- indPhen_df_vcf.NS[order(indPhen_df_vcf.NS$subpop, indPhen_df_vcf.NS$id),]
  
  head(indPhen_df_vcf)
  tail(indPhen_df_vcf)
  head(indPhen_df_vcf.NS)
  tail(indPhen_df_vcf.NS)
  
  # split individuals by subpop
  pop1_ids <- which(indPhen_df_vcf$subpop==1)
  pop2_ids <- which(indPhen_df_vcf$subpop==2)
  pop1_ids.NS <- which(indPhen_df_vcf.NS$subpop==1)
  pop2_ids.NS <- which(indPhen_df_vcf.NS$subpop==2)
  
  # split G by subpop
  G_pop1 <- G[, pop1_ids]
  G_pop2 <- G[, pop2_ids]
  dim(G_pop1)
  # head(G_pop1[,1:5])
  # head(G_pop2[,1:5])
  G_pop1.NS <- G.NS[, pop1_ids.NS]
  G_pop2.NS <- G.NS[, pop2_ids.NS]
  dim(G_pop1.NS)
  # head(G_pop1.NS[,1:5])
  # head(G_pop2.NS[,1:5])
  
  ## clustering ##
  # first transform matrix so that columns are mutations
  fordistp1 <- as.data.frame(t(G_pop1))
  fordistp2 <- as.data.frame(t(G_pop2))
  fordistp1.NS <- as.data.frame(t(G_pop1.NS))
  fordistp2.NS <- as.data.frame(t(G_pop2.NS))
  # this calculates the distance matrix for individuals based on similar
  # mutation combinations and euclidean is the square root of sum of squares 
  # of mutational dissimilarities
  dist_matp1 <- dist(fordistp1, method="euclidean")
  dist_matp2 <- dist(fordistp2, method="euclidean")
  dist_matp1.NS <- dist(fordistp1.NS, method="euclidean")
  dist_matp2.NS <- dist(fordistp2.NS, method="euclidean")
  # This then clusters individuals into distinct groups based on genetic 
  # distances calculated previously 
  pop1_clust <- hclust(dist_matp1, method = "ward.D")
  pop2_clust <- hclust(dist_matp2, method = "ward.D")
  pop1_clust.NS <- hclust(dist_matp1.NS, method = "ward.D")
  pop2_clust.NS <- hclust(dist_matp2.NS, method = "ward.D")
  str(pop1_clust)
  str(pop2_clust)
  str(pop1_clust.NS)
  str(pop2_clust.NS)
  pop1_order <- pop1_clust$order
  pop2_order <- pop2_clust$order
  pop1_order.NS <- pop1_clust.NS$order
  pop2_order.nS <- pop2_clust.NS$order

  # Why are we doing this for all allele affect sizes? shouldn't we just be looking at m2 mutations?
  # create variables that identify which mutations are which
  whichinversionmuts <- grep("MT=3", vcf.MAF@fix[,"INFO"]) #inversions
  whichqtnmuts <- grep("MT=2", vcf.MAF@fix[,"INFO"]) #qtns
  whichneutmuts <- grep("MT=1", vcf.MAF@fix[,"INFO"]) #neut
  whichinversionmuts.NS <- grep("MT=3", vcf.MAF.NS@fix[,"INFO"]) #inversions
  whichqtnmuts.NS <- grep("MT=2", vcf.MAF.NS@fix[,"INFO"]) #qtns
  whichneutmuts.NS <- grep("MT=1", vcf.MAF.NS@fix[,"INFO"]) #neut
  
  #vcf.MAF@fix[whichinversionmuts,"INFO"] # if you want info for specific mutation types
  # store info for mutation sin another variable and split data into columns
  info <- str_split(vcf.MAF@fix[,"INFO"], pattern =";", simplify=TRUE)
  info.NS <- str_split(vcf.MAF.NS@fix[,"INFO"], pattern =";", simplify=TRUE)
  head(info)
  dim(info)
  head(info.NS)
  dim(info.NS)
  
  # find allele effect sizes
  a <- as.numeric(substring(info[,2], first=3)) 
  a.NS <- as.numeric(substring(info.NS[,2], first=3)) 
  # head(a)
  # hist(a, breaks=seq(-0.01, 0.01, length.out=101))
  # head(a.NS)
  # hist(a.NS, breaks=seq(-0.01, 0.01, length.out=101))
  summary(a)
  length(a)
  dim(G)
  summary(a.NS)
  length(a.NS)
  dim(G.NS)

  #G * a gives the overall effect size of the mutations on the phenotype
  G1_alpha <- G_pop1*a 
  #head(G1_alpha[,1:10])
  # hist(G1_alpha, breaks=seq(-0.02, 0.02, length.out=101))
  
  G2_alpha <- G_pop2*a 
  #head(G2_alpha[,1:10])
  #hist(G2_alpha, breaks=seq(-0.02, 0.02, length.out=101))
  
  G1_alpha.NS <- G_pop1.NS*a.NS
  #head(G1_alpha.NS[,1:10])
  #hist(G1_alpha.NS, breaks=seq(-0.02, 0.02, length.out=101))
  
  G2_alpha.NS <- G_pop2.NS*a.NS
  #head(G2_alpha.NS[,1:10])
  #hist(G2_alpha.NS, breaks=seq(-0.02, 0.02, length.out=101))
  
  
  # this gives the distribution of phenotypes in each population 
  # population 1 is evolving to an optimum of 1
  # population 2 is evolving to an optimum of -1
  #hist(colSums(G1_alpha))
  #hist(colSums(G2_alpha))
  #hist(colSums(G1_alpha.NS))
  #hist(colSums(G2_alpha.NS))
  
  # Sanity check - mutations in rows
  #head(G[1:100,1:10])
  #t(G1_alpha[1:100,1:2])
  
  # Sanity check - mutations in rows
  #head(G.NS[1:100,1:10])
  #t(G1_alpha.NS[1:100,1:2])
  
  # get position of all mutations and plot a hist of mutation locations
  # the regions with higher frequencies of mutations are potentially inverted regions
  # hist(vcf_pos, breaks=seq(0,2600000, length.out=100))
  # hist(df.mutsMAFord$position, breaks=seq(0,2600000, length.out=100))
  # hist(vcf_pos.NS, breaks=seq(0,2600000, length.out=100))
  # hist(df.mutsMAFord.NS$position, breaks=seq(0,2600000, length.out=100))
  # 

  # df.mutsMAFord$is_vcf <- NA
  # df.mutsMAFord.NS$is_vcf <- NA
  # # Check to make sure FST values match between files
  # # this is slow, but correct
  # G_FST <- rep(NA, nrow(G)) 
  # count <- 0
  # for (i in 1:nrow(df.mutsMAFord)){
  #   x <- grep(df.mutsMAFord$LocusName[i], vcf.MAF@fix[,"INFO"])
  #   if (length(x)==1){
  #     G_FST[x] <- df.mutsMAFord$FST[i]
  #   }
  #   count <- count+1
  #   print(count)
  # }
  # 
  # hist(df.mutsMAFord$FST, breaks = 50)
  # head(G_FST)
  # hist(G_FST, breaks = 50)
  # 
  # sum(is.na(G_FST))
  # mutations that don't match up
  
  # G_FST.NS <- rep(NA, nrow(G.NS)) 
  # count.NS <- 0
  # for (i in 1:nrow(df.mutsMAFord.NS)){
  #   x <- grep(df.mutsMAFord.NS$LocusName[i], vcf.MAF.NS@fix[,"INFO"])
  #   if (length(x)==1){
  #     G_FST.NS[x] <- df.mutsMAFord.NS$FST[i]
  #   }
  #   count.NS <- count.NS+1
  #   print(count.NS)
  # }
  # 
  # hist(df.mutsMAFord.NS$FST, breaks = 50)
  # head(G_FST.NS)
  # hist(G_FST.NS, breaks = 50)
  # 
  # sum(is.na(G_FST.NS))
  # # mutations that don't match up
  # 
  # length(a)
  # length(G_FST)
  # 
  # length(a.NS)
  # length(G_FST.NS)
  # 
  # head(df.mutsMAFord)
  # head(df.mutsMAFord.NS)
  #df.mutsMAFord<- df.mutsMAF[order(df.mutsMAF$position),]
  #G.ord <- G[order(G)]
  
#### end process VCF
######################################################################################################    

  
  
  
######################################################################################################      
#### plot heatmaps

  #hist(a)
  a2 <- a
  # make an arbitary cutoff to visualize loci effect on phenotype
  a2[a>0.001] <- 1
  a2[a<0.001] <- -1
  
  G1_alpha <- G_pop1*a2*df.mutsMAFord$FST # make sure G and a line up
  G2_alpha <- G_pop2*a2*df.mutsMAFord$FST # make sure G and a line up
  
  #hist(G_pop1*a2)

  png(paste0(folderOut, seed, "_heatmapPop1alphaFST.png"), type = "cairo", width = 7, height = 7, units = 'in', res = 300)
    heatmap(t(G1_alpha[, pop1_order]),   
            main="Pop1 G*a+-*FST",cexCol = 0.3,
            Colv = NA, useRaster=TRUE,
            scale="none",
            col=two.colors(100, start = "blue", end="red", middle="white"))
  dev.off()
  
  png(paste0(folderOut, seed, "_heatmapPop2alphaFST.png"), type = "cairo", width = 7, height = 7, units = 'in', res = 300)
    heatmap(t(G2_alpha[, pop2_order]),   
            main="Pop2 G*a+-*FST",cexCol = 0.3,
            Colv = NA, useRaster=TRUE,
            scale="none",
            col=two.colors(100, start = "blue", end="red", middle="white") )
  dev.off()
  
  G_ref1 <- G_pop1
  G_ref2 <- G_pop2
  
  af_pop1 <- rowSums(G_pop1)/(2*ncol(G_pop1))
  af_pop2 <- rowSums(G_pop2)/(2*ncol(G_pop2))
  #ist(af_pop1)
  #hist(af_pop2)
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

  png(paste0(folderOut, seed, "_heatmapPop1geno.png"), type = "cairo", width = 7, height = 7, units = 'in', res = 300)
    heatmap(t(G_ref1[,pop1_order]), Rowv = NA,  main="Pop1 genotypes",cexCol = 0.3,
            Colv = NA, useRaster=TRUE,
            scale="none")
  dev.off()
  
  png(paste0(folderOut, seed, "_heatmapPop2geno.png"), type = "cairo", width = 7, height = 7, units = 'in', res = 300)
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
  
if(length(adapt.inv) != 0){
# list the G matrix with the position of each mutation and the chromosome 
# it is on. 
  training <- list(G = G, 
                  position = df.ord$position,
                  chromosome = df.mutsMAFord$chrom)

  training.NS <- list(G = G.NS, 
                   position = df.ord.NS$position,
                   chromosome = df.mutsMAFord.NS$chrom)
  # confirm that these are all in the proper order 
  
# puts it in the raw format and stores likelihood genotype probability
  G_coded <- add_code256(big_copy(t(training$G),
                                  type = "raw"), 
                        code=bigsnpr:::CODE_012)
  G_coded.NS <- add_code256(big_copy(t(training.NS$G),
                                  type = "raw"), 
                         code=bigsnpr:::CODE_012)

# this is doing SNP pruning - removing correlated SNPs
  newpc <- snp_autoSVD(G=G_coded, infos.chr = as.integer(training$chromosome),
                       infos.pos = training$position, roll.size = 0)
 
  newpc.NS <- try(snp_autoSVD(G=G_coded.NS, infos.chr = as.integer(training.NS$chromosome),
                       infos.pos = training.NS$position, roll.size = 0))
  if("try-error" %in% class(newpc.NS)){
    pcadapt.fail <- TRUE
  } else {
  # take snps with highest MAF and correlate snps around it
  # Snps with R^2 > 0.2 are removed
  # the subset is the indexes of the remaining SNPs
  str(newpc) #2760
  str(newpc.NS) #4056

# These are the indexes of the quasi-independent 
# set of loci that are kept after pruning for LD
  which_pruned <-  attr(newpc, which = "subset")
  which_pruned.NS <-  attr(newpc.NS, which = "subset")
  length(which_pruned)
  length(which_pruned.NS)

  training$G_coded <- G_coded
  training$G_pruned <- training$G[which_pruned,]
  training$which_pruned <- which_pruned
  training.NS$G_coded <- G_coded.NS
  training.NS$G_pruned <- training.NS$G[which_pruned.NS,]
  training.NS$which_pruned <- which_pruned.NS

  df.muts.MAF$quasi_indep <- FALSE
  df.muts.MAF$quasi_indep[training$which_pruned] <- TRUE
  df.muts.NS.MAF$quasi_indep <- FALSE
  df.muts.NS.MAF$quasi_indep[training.NS$which_pruned] <- TRUE

#### end bigsnpr
######################################################################################################
  
######################################################################################################
#### PCADAPT

  # check structure of data
  #training$G[1:6, 1:6]
 # training.NS$G[1:6, 1:6]
  head(df.mutsMAFord)
  head(df.mutsMAFord[order(df.mutsMAFord$position_vcf), 1:ncol(df.mutsMAFord)])
  head(df.mutsMAFord.NS[order(df.mutsMAFord.NS$position_vcf), 1:ncol(df.mutsMAFord.NS)])
  
  # calculate pca loadings from genotype matrix 
  # first you need to convert to a lfmm file
  gename <- paste0(seed, "_genotypes.lfmm")
  gename.NS <- paste0(seed, "noSel_genotypes.lfmm")
  write.lfmm(t(training$G), paste0(folderIn, gename))
  write.lfmm(t(training.NS$G), paste0(folderIn, gename.NS))
  pcafile <- read.pcadapt(paste0(folderIn, gename), type="lfmm")
  pcafile.NS <- read.pcadapt(paste0(folderIn,gename.NS), type="lfmm")
  pca_all <- pcadapt(pcafile,K=2)
  pca_all.NS <- pcadapt(pcafile.NS,K=2)
  head(pca_all$loadings)
  head(pca_all.NS$loadings)
  str(pca_all)
  
  pdf(paste0(folderOut, seed, "_pcaScores.pdf"), width = 5, height = 5)
    plot(pca_all$scores[,1], pca_all$scores[,2], col = c(rep("red", 1000), rep("blue", 1000)))
  dev.off()
  
  pdf(paste0(folderOut, seed, "noSel_pcaScores.pdf"), width = 5, height = 5)
    plot(pca_all.NS$scores[,1], pca_all.NS$scores[,2], col = c(rep("red", 1000), rep("blue", 1000)))
  dev.off()
  
  head(df.muts.MAF)
  head(pca_all)
  
  # add pca loadings to muts dataframe 
  df.mutsMAFord$pca_ALL_PC1_loadings <- pca_all$loadings[,1]
  df.mutsMAFord$pca_ALL_PC2_loadings <- pca_all$loadings[,2]
  df.mutsMAFord.NS$pca_ALL_PC1_loadings <- pca_all.NS$loadings[,1]
  df.mutsMAFord.NS$pca_ALL_PC2_loadings <- pca_all.NS$loadings[,2]
  head(df.mutsMAFord)
  
  pdf(paste0(folderOut, seed, "_pcaLoadingsPos.pdf"), width = 5, height = 5)
    plot(df.mutsMAFord$position_vcf, df.mutsMAFord$pca_ALL_PC1_loadings)
  dev.off()
  
  pdf(paste0(folderOut, seed, "noSel_pcaLoadingsPos.pdf"), width = 5, height = 5)
    plot(df.mutsMAFord.NS$position_vcf, df.mutsMAFord.NS$pca_ALL_PC1_loadings)
  dev.off()
  
### PCA loadings for pruned data ####
  gename2 <- paste0(folderIn, seed, "_genotypes_pruned.lfmm")
  gename2.NS <- paste0(folderIn, seed, "noSel_genotypes_pruned.lfmm")
  write.lfmm(t(training$G_pruned), gename2)
  write.lfmm(t(training.NS$G_pruned), gename2.NS)
  pcafile2 <- read.pcadapt(gename2, type="lfmm")
  pcafile2.NS <- read.pcadapt(gename2.NS, type="lfmm")
  pca_pruned <- pcadapt(pcafile2,K=2)
  pca_pruned.NS <- pcadapt(pcafile2.NS,K=2)
  
  pdf(paste0(folderOut, seed, "_pcaScoresPruned.pdf"), width =5, height = 5)
    plot(pca_pruned$scores[,1], pca_pruned$scores[,2], col = c(rep("red", 1000), rep("blue", 1000)))
  dev.off()
  
  pdf(paste0(folderOut, seed, "noSel_pcaScoresPruned.pdf"), width =5, height = 5)
    plot(pca_pruned.NS$scores[,1], pca_pruned.NS$scores[,2], col = c(rep("red", 1000), rep("blue", 1000)))
  dev.off()
  
  
  # add pca loadings to muts dataframe
  df.mutsMAFord$pca_PRUNED_PC1_loadings <- NA 
  df.mutsMAFord$pca_PRUNED_PC2_loadings <- NA
  df.mutsMAFord$pca_PRUNED_PC1_loadings[training$which_pruned] <- pca_pruned$loadings[,1]
  df.mutsMAFord$pca_PRUNED_PC2_loadings[training$which_pruned] <- pca_pruned$loadings[,2] 
  df.mutsMAFord.NS$pca_PRUNED_PC1_loadings <- NA 
  df.mutsMAFord.NS$pca_PRUNED_PC2_loadings <- NA
  df.mutsMAFord.NS$pca_PRUNED_PC1_loadings[training.NS$which_pruned] <- pca_pruned.NS$loadings[,1]
  df.mutsMAFord.NS$pca_PRUNED_PC2_loadings[training.NS$which_pruned] <- pca_pruned.NS$loadings[,2] 
    
  ### outlier detection ###
  ## all data
  df.mutsMAFord$pcadapt_4.3.3_ALL_chisq <- as.numeric(pca_all$chi2.stat)
  df.mutsMAFord$pcadapt_4.3.3_ALL_log10p <- -log10(pca_all$pvalues)
  
  df.mutsMAFord.NS$pcadapt_4.3.3_ALL_chisq <- as.numeric(pca_all.NS$chi2.stat)
  df.mutsMAFord.NS$pcadapt_4.3.3_ALL_log10p <- -log10(pca_all.NS$pvalues)
 
  ### outlier detection ### 
  ## pruned data
  test <- try(snp_gc(snp_pcadapt(training$G_coded, U.row = newpc$u[,1])))
  if("try-error" %in% class(test)){
    df.mutsMAFord$pcadapt_4.3.3_PRUNED_log10p <- 0
    df.mutsMAFord$pcadapt_4.3.3_PRUNED_pvalue <- NA
    df.mutsMAFord$qvalue <- NA
    # call outliers based on qvalues
    df.mutsMAFord$pcadapt_outlier <- FALSE
    df.mutsMAFord$pcadapt_outlier <- as.factor(df.mutsMAFord$pcadapt_outlier)
  } else {
    df.mutsMAFord$pcadapt_4.3.3_PRUNED_log10p <- -predict(test,log10=T)
    df.mutsMAFord$pcadapt_4.3.3_PRUNED_pvalue <- 10^(-df.mutsMAFord$pcadapt_4.3.3_PRUNED_log10p)
    df.mutsMAFord$qvalue <- qvalue(df.mutsMAFord$pcadapt_4.3.3_PRUNED_pvalue)$qvalues  
    df.mutsMAFord$pcadapt_outlier <- ifelse(df.mutsMAFord$qvalue > 0.01, FALSE, TRUE)
    df.mutsMAFord$pcadapt_outlier <- as.factor(df.mutsMAFord$pcadapt_outlier)
  }
  
  test.NS <- snp_gc(snp_pcadapt(training.NS$G_coded, U.row = newpc.NS$u[,1]))
  if("try-error" %in% class(test.NS)){
    df.mutsMAFord.NS$pcadapt_4.3.3_PRUNED_log10p <- NA
    df.mutsMAFord.NS$pcadapt_4.3.3_PRUNED_pvalue <- NA
    df.mutsMAFord.NS$qvalue <- NA
    df.mutsMAFord.NS$pcadapt_outlier <- FALSE
    df.mutsMAFord.NS$pcadapt_outlier <- as.factor(df.mutsMAFord.NS$pcadapt_outlier)
  } else {
    df.mutsMAFord.NS$pcadapt_4.3.3_PRUNED_log10p <- -predict(test.NS,log10=T)
    df.mutsMAFord.NS$pcadapt_4.3.3_PRUNED_pvalue <- 10^(-df.mutsMAFord.NS$pcadapt_4.3.3_PRUNED_log10p)
    df.mutsMAFord.NS$qvalue <- qvalue(df.mutsMAFord.NS$pcadapt_4.3.3_PRUNED_pvalue)$qvalues
    df.mutsMAFord.NS$pcadapt_outlier <- ifelse(df.mutsMAFord.NS$qvalue > 0.01, FALSE, TRUE)
    df.mutsMAFord.NS$pcadapt_outlier <- as.factor(df.mutsMAFord.NS$pcadapt_outlier)
  }
  
  plot(newpc$u[,1], newpc$u[,2], col = c(rep("red", 1000), rep("blue", 1000)))
  plot(newpc.NS$u[,1], newpc.NS$u[,2], col = c(rep("red", 1000), rep("blue", 1000)))
  
  dim(df.mutsMAFord)
  dim(df.mutsMAFord.NS)
  pcadapt.fail <- FALSE
  }
} 
  #### end PCADAPT
######################################################################################################



######################################################################################################
#### OutFLANK

if(length(adapt.inv) != 0){
  if("try-error" %in% class(newpc.NS)){
    outflank.fail <- TRUE
  } else {
   
  FstDataFrame <- MakeDiploidFSTMat(t(training$G),rownames(training$G),
                                    c(rep("Pop1", 1000), rep("Pop2", 1000)))
  colnames(FstDataFrame)[3] <- "FST_outflank"
  
  
  # ask about qthreshold and trim fraction
  out_pruned <- OutFLANK(FstDataFrame[training$which_pruned,], NumberOfSamples=2, 
                         LeftTrimFraction=0.05, RightTrimFraction=0.05,
                          Hmin=0.1, qthreshold=0.01)
  outflank.test <- try(out_pruned == 0)
  if("try-error" %in% class(outflank.test)){
  
  head(df.mutsMAFord)
  colnames(df.mutsMAFord)[10] <- "FST_slim"
  df.FST.temp <- merge(df.mutsMAFord, FstDataFrame, by = "LocusName")
  df.FST <- df.FST.temp[order(df.FST.temp$position_vcf),]
  dim(df.FST)
  
  df.out <- pOutlierFinderChiSqNoCorr(df.FST, 
                                  Fstbar = out_pruned$FSTNoCorrbar, 
                                  dfInferred = out_pruned$dfInferred, Hmin=0.1)
  
  df.out$OutFLANK_0.2_PRUNED_log10p <- -log10(df.out$pvaluesRightTail)
  df.out$OutFLANK_0.2_PRUNED_log10p_add <- -log10(df.out$pvaluesRightTail + 1/10000000000000000000)
  
  pdf(paste0(folderOut, seed, "_outflankFstHist.pdf"), width =5, height = 5)
    OutFLANKResultsPlotter(out_pruned, Zoom = T)
  dev.off()
  
  
  FstDataFrame.NS <- MakeDiploidFSTMat(t(training.NS$G),rownames(training.NS$G),
                                       c(rep("Pop1", 1000), rep("Pop2", 1000)))
  colnames(FstDataFrame.NS)[3] <- "FST_outflank"
  
  out_pruned.NS <- OutFLANK(FstDataFrame.NS[training.NS$which_pruned,], NumberOfSamples=2, 
                            LeftTrimFraction=0.05, RightTrimFraction=0.05,
                            Hmin=0.1, qthreshold=0.01)     
  str(out_pruned.NS)
  head(df.mutsMAFord.NS)
  colnames(df.mutsMAFord.NS)[10] <- "FST_slim"
  df.FST.temp.NS <- merge(df.mutsMAFord.NS, FstDataFrame.NS, by = "LocusName")
  df.FST.NS <- df.FST.temp.NS[order(df.FST.temp.NS$position_vcf),]
  dim(df.FST.NS)
  df.out.NS <- pOutlierFinderChiSqNoCorr(df.FST.NS, 
                                         Fstbar = out_pruned.NS$FSTNoCorrbar, 
                                         dfInferred = out_pruned.NS$dfInferred, Hmin=0.1)
  df.out.NS$OutFLANK_0.2_PRUNED_log10p <- -log10(df.out.NS$pvaluesRightTail)
  df.out.NS$OutFLANK_0.2_PRUNED_log10p_add <- -log10(df.out.NS$pvaluesRightTail + 1/10000000000000000000)
  
  
  pdf(paste0(folderOut, seed, "noSel_outflankFstHist.pdf"), width =5, height = 5)
    OutFLANKResultsPlotter(out_pruned.NS, Zoom = T)
  dev.off()
  
  outflank.fail <- FALSE
   } else {
    outflank.fail <- TRUE # if statement for outflank error
   }
  }# if statement for error in pcadapt
}
#### end OutFLANK
######################################################################################################



######################################################################################################
#### start Outlier Plotting

#########################
#### Manhattan plots ####
  df.adaptInv <- as.data.frame(cbind(ID = adapt.inv, center_bases = center.bases, first_bases = first.bases, 
                                     final_bases = final.bases))
  df.nonadaptInv <- as.data.frame(cbind(ID = nonadapt.inv, center_bases = center.bases.NA, first_bases = first.bases.NA, 
                                        final_bases = final.bases.NA))
  # SELECTION #
  df.neutQTNmuts <- df.muts.MAF[df.muts.MAF$inOut != "inv",]
  df.neutQTNmuts$inOut <- factor(df.neutQTNmuts$inOut)
  df.neutQTNmuts$chrom <- factor(df.neutQTNmuts$chrom)
  
  df.neutQTNmuts.NS <- df.muts.NS.MAF[df.muts.NS.MAF$inOut != "inv",]
  df.neutQTNmuts.NS$inOut <- factor(df.neutQTNmuts.NS$inOut)
  df.neutQTNmuts.NS$chrom <- factor(df.neutQTNmuts.NS$chrom)
  
  manh.plot <- ggplot(df.neutQTNmuts, aes(x = position, y = FST, 
                                          group = interaction(inOut, chrom))) 
   if(length(adapt.inv) != 0){
     manh.plot <- manh.plot + geom_rect(data=df.adaptInv, mapping=aes(xmin=first_bases, xmax=final_bases, ymin=0,
                                             ymax=max(c(df.neutQTNmuts$FST, df.neutQTNmuts.NS$FST))*1.1), 
                                        fill = "tan1", color="black", alpha=0.5, inherit.aes = FALSE) 
   }
    manh.plot <- manh.plot +
    geom_point(data = df.neutQTNmuts, aes(color = chrom, shape = inOut)) + 
    scale_shape_manual(name = "QTN location", 
                       labels = c("Inside Inversion", "Neutral", "Outside Inversion"), 
                       values=c(19, 15, 1)) +
    scale_color_manual(name = "Chromosome", values = c(rep(c("navy", "lightblue"), 10), "darkgrey")) + 
    labs(title = expression(bold("Selection - Adaptive Inversions"))) + 
    theme(legend.position = "none") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92"),
          text = element_text(size = 11)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(c(df.neutQTNmuts$FST, df.neutQTNmuts.NS$FST)))*1.1) 
  

  manh.plot.NA <- ggplot(df.neutQTNmuts, aes(x = position, y = FST, 
                                          group = interaction(inOut, chrom))) 
    if(length(nonadapt.inv) != 0){
      manh.plot.NA <- manh.plot.NA + geom_rect(data=df.nonadaptInv, mapping=aes(xmin=first_bases, xmax=final_bases, ymin=0,
                                            ymax=max(c(df.neutQTNmuts$FST, df.neutQTNmuts.NS$FST))*1.1), 
                                            fill = "tan1", color="black", alpha=0.5, inherit.aes = FALSE)
    }
    manh.plot.NA <- manh.plot.NA + geom_point(data = df.neutQTNmuts, aes(color = chrom, shape = inOut)) + 
    scale_shape_manual(name = "QTN location", 
                       labels = c("Inside Inversion", "Neutral", "Outside Inversion"), 
                       values=c(19, 15, 1)) +
    scale_color_manual(name = "Chromosome", values = c(rep(c("navy", "lightblue"), 10), "darkgrey")) + 
    labs(title = expression(bold("Selection - Nonadaptive Inversions"))) + 
    theme(legend.position = "none") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92"),
          text = element_text(size = 11)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(c(df.neutQTNmuts$FST, df.neutQTNmuts.NS$FST))*1.1))
  
  # NO SELECTION #
   manh.plot.NS <- ggplot(df.neutQTNmuts.NS, aes(x = position, y = FST, 
                                                group = interaction(inOut, chrom)))
    if(nrow(df.inv.NS) != 0){
      manh.plot.NS <- manh.plot.NS + geom_rect(data=df.inv.NS, mapping=aes(xmin=first_base, xmax=final_base, ymin=0,
                                               ymax=max(c(df.neutQTNmuts$FST, df.neutQTNmuts.NS$FST))*1.1), 
                                      fill = "tan1", color="black", alpha=0.5, inherit.aes = FALSE)
    }
    manh.plot.NS <- manh.plot.NS + geom_point(aes(color = chrom, shape = inOut)) + 
    scale_color_manual(name = "Chromosome",
                       values = c(rep(c("navy", "lightblue"), 10), "darkgrey")) +
    scale_shape_manual(name = "QTN location", labels = c("Inside Inversion", "Neutral", "Outside Inversion"), 
                       values=c(19, 15, 1)) +
    labs(title = expression(bold("No Selection"))) + 
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92"),
          text = element_text(size = 11)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(c(df.neutQTNmuts$FST, df.neutQTNmuts.NS$FST))*1.1))
  
  ## No legend
  manh.plot.NS.noleg <- manh.plot.NS + theme(legend.position = "none")
  manh.plot.noleg <-  manh.plot+ theme(legend.position = "none")
  manh.plot.NA.noleg <- manh.plot.NA  + theme(legend.position = "none")
  legManh <- g_legend(manh.plot)
  
  ## Print plot
  png(paste0(folderOut, seed, "_manh.png"), type = "cairo", width = 7, height = 7, units = 'in', res = 300)
    print(ggarrange(manh.plot.noleg, blank, manh.plot.NA.noleg, legManh, manh.plot.NS.noleg, blank, nrow = 3, ncol = 2, widths = c(2.3,0.8,2.3,0.8,2.3,0.8)))
  dev.off()

if(length(adapt.inv) != 0 ){
  if(pcadapt.fail == FALSE & outflank.fail == FALSE){
   
    
#######################
#### Outlier Plots ####
## Adaptive Inversions
  ## Get max value for scales on outlier plots
  max_value_log10 <- max(c(df.out$OutFLANK_0.2_PRUNED_log10p_add[!is.na(df.out$OutFLANK_0.2_PRUNED_log10p_add)], 
                           df.out$pcadapt_4.3.3_PRUNED_log10p[!is.na(df.out$pcadapt_4.3.3_PRUNED_log10p)]))*1.1
  
  ## Outflank
  outflank.log10p <- ggplot(df.out, aes(x = position_vcf, y = OutFLANK_0.2_PRUNED_log10p_add))  
    if(length(adapt.inv) != 0){
      outflank.log10p <- outflank.log10p + geom_rect(data=df.adaptInv, mapping=aes(xmin=first_bases, xmax=final_bases, ymin=0,
                                            ymax=max_value_log10), fill = "tan1", 
                                            color="black", alpha=0.5, inherit.aes = FALSE) 
    }
  outflank.log10p <- outflank.log10p + geom_point(data = df.out[!is.na(df.out$OutlierFlag),], aes(color = OutlierFlag, 
                                                                     shape = OutlierFlag))+
    scale_color_manual(values = c("black", "red")) +
    scale_shape_manual(values = c(19, 1)) +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92"),
          text = element_text(size = 11)) +
    labs(title = "OutFLANK",
         y = "-log10(p-values)",
         x = " ") + 
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,max_value_log10))

  ## PCAdapt
  pcadapt.log10p <- ggplot(df.mutsMAFord, aes(x = position_vcf, y = pcadapt_4.3.3_PRUNED_log10p))
    if(length(adapt.inv) != 0){
      pcadapt.log10p <- pcadapt.log10p + geom_rect(data=df.adaptInv, mapping=aes(xmin=first_bases, xmax=final_bases, ymin=0,
                                              ymax=max_value_log10), fill = "tan1", color="black", 
                                              alpha=0.5, inherit.aes = FALSE)
    } 
  
  pcadapt.log10p <- pcadapt.log10p + geom_point(data = df.mutsMAFord, aes(color = pcadapt_outlier, 
                                                                          shape = pcadapt_outlier)) +
    scale_color_manual(name = "Outlier", 
                       values = c("black",  "red")) +
    scale_shape_manual(name = "Outlier", 
                       values = c(19, 1)) +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "dimgrey"),
          text = element_text(size = 11)) +
    labs(title = "PCAdapt",
         y = "-log10(p-value)",
         x = "Genome Position") + 
    theme(legend.position = "none") + 
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value_log10))  
  
  ## Manhattan plot
  manh.plot.slim <- manh.plot + labs(title = "Adaptive Inversions - Simulation Output",
                                     x = " ") +
    theme(legend.position = "none") 
  
  outflank.log10pNoleg <- outflank.log10p + theme(legend.position = "none")
  outlierleg <- get_legend(outflank.log10p)
  
  ## Print plots
  png(paste0(folderOut, seed, "_outlierTests.png"), type = "cairo", width = 7, height = 7, units = 'in', res = 300)
    print(ggarrange(manh.plot.slim, blank, outflank.log10pNoleg, outlierleg, pcadapt.log10p, blank, nrow = 3, ncol = 2, widths = c(2.3,0.8,2.3,0.8,2.3,0.8)))
  dev.off()
  
  ## Nonadaptive Inversions
  
  ## Outflank
  outflank.log10p.NA <- ggplot(df.out[!is.na(df.out$OutlierFlag),], aes(x = position_vcf, y = OutFLANK_0.2_PRUNED_log10p_add))
    if(length(nonadapt.inv) != 0){
      outflank.log10p.NA <- outflank.log10p.NA + geom_rect(data=df.nonadaptInv, mapping=aes(xmin=first_bases, xmax=final_bases, ymin=0,
                                            ymax=max_value_log10), fill = "tan1", 
              color="black", alpha=0.5, inherit.aes = FALSE)
    }
    outflank.log10p.NA <- outflank.log10p.NA + geom_point(data = df.out[!is.na(df.out$OutlierFlag),], aes(color = OutlierFlag, shape = OutlierFlag)) +
    scale_color_manual(values = c("black", "red")) +
    scale_shape_manual(values = c(19, 1)) +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92"),
          text = element_text(size = 11)) +
    labs(title = "OutFLANK",
         y = "-log10(p-values)",
         x = " ") + 
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,max_value_log10))
  
  ## PCAdapt
  pcadapt.log10p.NA <- ggplot(df.mutsMAFord, aes(x = position_vcf, y = pcadapt_4.3.3_PRUNED_log10p))
    if(length(nonadapt.inv) != 0){
      pcadapt.log10p.NA <- pcadapt.log10p.NA + geom_rect(data=df.nonadaptInv, mapping=aes(xmin=first_bases, xmax=final_bases, ymin=0,
                                            ymax=max_value_log10), fill = "tan1", color="black", alpha=0.5, inherit.aes = FALSE)
    }
    pcadapt.log10p.NA <- pcadapt.log10p.NA + geom_point(data = df.mutsMAFord, aes(color = pcadapt_outlier, shape = pcadapt_outlier)) +
    scale_color_manual(name = "Outlier", 
                       values = c("black",  "red")) +
    scale_shape_manual(name = "Outlier", 
                       values = c(19, 1)) +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "dimgrey"),
          text = element_text(size = 11)) +
    labs(title = "PCAdapt",
         y = "-log10(p-value)",
         x = "Genome Position") + 
    theme(legend.position = "none") + 
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value_log10)) 
   
  ## Manhattan plot
  manh.plot.slim.NA <- manh.plot.NA + labs(title = "Nonadaptive Inversions - Simulation Output",
                                           x = " ") +
    theme(legend.position = "none") 
  
  outflank.log10pNoleg.NA <- outflank.log10p.NA + theme(legend.position = "none")
  outlierlegNA <- get_legend(outflank.log10p.NA)
  
  ## Print plots
  png(paste0(folderOut, seed, "_outlierTestsNA.png"), type = "cairo", width = 7, height = 7, units = 'in', res = 300)
   print(ggarrange(manh.plot.slim.NA, blank, outflank.log10pNoleg.NA, outlierlegNA, pcadapt.log10p.NA, blank, nrow = 3, ncol = 2, widths = c(2.3,0.8,2.3,0.8,2.3,0.8)))
  dev.off()
  
  
  ## NO SELECTION ##
 max_value_log10.NS <- max(c(df.out.NS$OutFLANK_0.2_PRUNED_log10p_add[!is.na(df.out.NS$OutFLANK_0.2_PRUNED_log10p_add)], 
                             df.mutsMAFord.NS$pcadapt_4.3.3_PRUNED_log10p[!is.na(df.out.NS$pcadapt_4.3.3_PRUNED_log10p)]))*1.1
 
 pcadapt.log10p.NS <- ggplot(df.mutsMAFord.NS, aes(x = position_vcf, y = pcadapt_4.3.3_PRUNED_log10p))
    if(nrow(df.inv.NS) != 0){
      pcadapt.log10p.NS <- pcadapt.log10p.NS + geom_rect(data=df.inv.NS, mapping=aes(xmin=first_base, xmax=final_base, ymin=0,
                                               ymax=max_value_log10.NS), fill = "tan1", 
              color="black", alpha=0.5, inherit.aes = FALSE)
    }
    pcadapt.log10p.NS <- pcadapt.log10p.NS + geom_point(data = df.mutsMAFord.NS, aes(color = pcadapt_outlier, shape = pcadapt_outlier)) +
    scale_color_manual(name = "Outlier", 
                       values = c("black",  "red")) +
    scale_shape_manual(name = "Outlier", 
                       values = c(19, 1)) +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "dimgrey"),
          text = element_text(size = 11)) +
    labs(title = "PCAdapt",
         y = "-log10(p-value)",
         x = "Genome Position") + 
    theme(legend.position = "none") + 
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_value_log10.NS)) 
  
  
  ## Outflank
  outflank.log10p.NS <- ggplot(df.out.NS[!is.na(df.out.NS$OutlierFlag),], aes(x = position_vcf, y = OutFLANK_0.2_PRUNED_log10p_add))
    if(nrow(df.inv.NS) != 0){
      outflank.log10p.NS <- outflank.log10p.NS + geom_rect(data=df.inv.NS, mapping=aes(xmin=first_base, xmax=final_base, ymin=0,
                                                                                       ymax=max_value_log10.NS), fill = "tan1", 
                                                            color="black", alpha=0.5, inherit.aes = FALSE)
    }
    outflank.log10p.NS <- outflank.log10p.NS + geom_point(data = df.out.NS[!is.na(df.out.NS$OutlierFlag),], aes(color = OutlierFlag, shape = OutlierFlag))+
    scale_color_manual(values = c("black", "red")) +
    scale_shape_manual(values = c(19, 1)) +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92"),
          text = element_text(size = 11)) +
    labs(title = "OutFLANK",
         y = "-log10(p-values)",
         x = " ") + 
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,max_value_log10.NS))
  
  manh.plot.NS.noleg <- manh.plot.NS + labs(title = "No Selection - Simulation Output",
                                           x = " ") +
    theme(legend.position = "none") 

  outflank.log10pNoleg.NS <- outflank.log10p.NS + theme(legend.position = "none")
  outlierleg.NS <- g_legend(outflank.log10p.NS)
  
  png(paste0(folderOut, seed, "noSel_outlierTests.png"), type = "cairo", width = 7, height = 7, units = 'in', res = 300)
    print(ggarrange(manh.plot.NS.noleg, blank, outflank.log10pNoleg.NS, outlierleg.NS, pcadapt.log10p.NS, blank, nrow = 3, ncol = 2, widths = c(2.3,0.8,2.3,0.8,2.3,0.8)))
  dev.off()

#### end outlier plotting
######################################################################################################


######################################################################################################    
## Outlier identification from genome scans

  
  colnames(df.qtnMuts.MAF)[2] <- "LocusName"
  df.outlierComp <- left_join(df.qtnMuts.MAF[c(2,12:ncol(df.qtnMuts.MAF))], df.out, by = "LocusName")
  dim(df.outlierComp)
  
  
  
  colnames(df.qtnMuts.NS.MAF)[2] <- "LocusName"
  df.outlierComp.NS <- left_join(df.qtnMuts.NS.MAF[c(2,12:ncol(df.qtnMuts.NS.MAF))], df.out.NS, by = "LocusName")
  dim(df.outlierComp.NS)
  } else {
  df.outlierComp <- NULL
  df.outlierComp.NS <- NULL
  } 
  
  } else {
  df.outlierComp <- NULL
  df.outlierComp.NS <- NULL
  }
  
  
if(length(adapt.inv) != 0){
  if(pcadapt.fail == FALSE & outflank.fail == FALSE){
  
  false_neg_pcadapt <- 0
  true_pos_pcadapt <- 0
  true_pos_outflank <- 0
  false_neg_outflank <- 0
  for(i in 1:nrow(df.adaptInv)){
  
    # subset dataframe for the qtns in the focal inversion
    focal.inv <- df.outlierComp[df.outlierComp$position >= df.adaptInv$first_bases[i] & 
                     df.outlierComp$position <= df.adaptInv$final_bases[i], ]
    # count number of qtns
    var <- nrow(focal.inv)
    
    # get top 10% of values by each test
    outflank.top <- focal.inv[order(focal.inv$OutFLANK_0.2_PRUNED_log10p_add, decreasing = TRUE),][1:round(var*0.1),]
    pcadapt.top <- focal.inv[order(focal.inv$pcadapt_4.3.3_PRUNED_log10p, decreasing = TRUE),][1:round(var*0.1),]
    # sum the trues
    outflank.T <- sum(as.logical(outflank.top$OutlierFlag), na.rm = TRUE)
    pcadapt.T <- sum(as.logical(pcadapt.top$pcadapt_outlier), na.rm = TRUE)
    # if they are all trues: true positive else false negative
    ifelse(outflank.T == nrow(outflank.top), true_pos_outflank <- true_pos_outflank + 1, false_neg_outflank <- false_neg_outflank + 1)
    ifelse(pcadapt.T == nrow(pcadapt.top), true_pos_pcadapt <- true_pos_pcadapt + 1, false_neg_pcadapt <- false_neg_pcadapt + 1)
    
  }
  
} else {
  false_neg_pcadapt <- NA
  true_pos_pcadapt <- NA
  true_pos_outflank <- NA
  false_neg_outflank <- NA
  }
} else {
  false_neg_pcadapt <- NA
  true_pos_pcadapt <- NA
  true_pos_outflank <- NA
  false_neg_outflank <- NA
}


if(length(adapt.inv) != 0){
  if(pcadapt.fail == FALSE & outflank.fail == FALSE){
  true_neg_pcadapt <- 0
  false_pos_pcadapt <- 0
  true_neg_outflank <- 0
  false_pos_outflank <- 0
  for(i in 1:nrow(df.nonadaptInv)){
    focal.inv <- df.outlierComp[df.outlierComp$position >= df.nonadaptInv$first_bases[i] & 
                                df.outlierComp$position <= df.nonadaptInv$final_bases[i], ]
    var <- nrow(focal.inv)
  
    # get top 10% of values by each test
    outflank.top <- focal.inv[order(focal.inv$OutFLANK_0.2_PRUNED_log10p_add, decreasing = TRUE),][1:round(var*0.1),]
    pcadapt.top <- focal.inv[order(focal.inv$pcadapt_4.3.3_PRUNED_log10p, decreasing = TRUE),][1:round(var*0.1),]
  
    # sum the trues
    outflank.T <- sum(as.logical(outflank.top$OutlierFlag), na.rm = TRUE)
    pcadapt.T <- sum(as.logical(pcadapt.top$pcadapt_outlier), na.rm = TRUE)
  
    #true negative and false positive for each test
    ifelse(outflank.T != nrow(outflank.top), true_neg_outflank <- true_neg_outflank + 1, false_pos_outflank <- false_pos_outflank + 1)
    ifelse(pcadapt.T != nrow(pcadapt.top), true_neg_pcadapt <- true_neg_pcadapt + 1, false_pos_pcadapt <- false_pos_pcadapt + 1)
  
  }
} else {
    true_neg_pcadapt <- NA
    false_pos_pcadapt <- NA
    true_neg_outflank <- NA
    false_pos_outflank <- NA
  }
} else {
    true_neg_pcadapt <- NA
    false_pos_pcadapt <- NA
    true_neg_outflank <- NA
    false_pos_outflank <- NA
}


if(length(adapt.inv) != 0){
  if(pcadapt.fail == FALSE & outflank.fail == FALSE){

  true_neg_outflank_NS <- 0
  true_neg_pcadapt_NS <- 0
  false_pos_outflank_NS <- 0
  false_pos_pcadapt_NS <- 0
  for(i in 1:nrow(df.inv.NS)){
    focal.inv.NS <- df.outlierComp.NS[df.outlierComp.NS$position >= df.inv.NS$first_base[i] & 
                                  df.outlierComp.NS$position <= df.inv.NS$final_base[i], ]
    var <- nrow(focal.inv.NS)
  
    # get top 10% of values by each test
    outflank.top <- focal.inv.NS[order(focal.inv.NS$OutFLANK_0.2_PRUNED_log10p_add, decreasing = TRUE),][1:round(var*0.1),]
    pcadapt.top <- focal.inv.NS[order(focal.inv.NS$pcadapt_4.3.3_PRUNED_log10p, decreasing = TRUE),][1:round(var*0.1),]
  
    # sum the trues
    pcadapt.T <- sum(as.logical(pcadapt.top$pcadapt_outlier), na.rm = TRUE)
    outflank.T <- sum(as.logical(outflank.top$OutlierFlag), na.rm = TRUE)
  
    #true negative and false positive for each test
    ifelse(outflank.T != nrow(outflank.top), true_neg_outflank_NS <- true_neg_outflank_NS + 1, false_pos_outflank_NS <- false_pos_outflank_NS + 1)
    ifelse(pcadapt.T != nrow(pcadapt.top), true_neg_pcadapt_NS <- true_neg_pcadapt_NS + 1, false_pos_pcadapt_NS <- false_pos_pcadapt_NS + 1)
  
    }
  } else {
    true_neg_outflank_NS <- NA
    true_neg_pcadapt_NS <- NA
    false_pos_outflank_NS <- NA
    false_pos_pcadapt_NS <- NA
  }
} else {
  true_neg_outflank_NS <- NA
  true_neg_pcadapt_NS <- NA
  false_pos_outflank_NS <- NA
  false_pos_pcadapt_NS <- NA
}

if(length(adapt.inv) != 0 ){
  if("try-error" %in% class(test) | pcadapt.fail == TRUE | outflank.fail == TRUE){
    false_neg_pcadapt <- NA
    true_pos_pcadapt <- NA
    true_neg_pcadapt <- NA
    false_pos_pcadapt <- NA
  } 
} else {
  false_neg_pcadapt <- NA
  true_pos_pcadapt <- NA
  true_neg_pcadapt <- NA
  false_pos_pcadapt <- NA
}
  
if(length(adapt.inv) != 0){  
  if("try-error" %in% class(test.NS) | pcadapt.fail == TRUE | outflank.fail == TRUE){
    true_neg_pcadapt_NS <- NA
    false_pos_pcadapt_NS <- NA
  } 
} else {
  true_neg_pcadapt_NS <- NA
  false_pos_pcadapt_NS <- NA
}

#### end Outlier identification
######################################################################################################

######################################################################################################    
## Data to output to table

## Q1&2: VA in inversions; final LA
## Q3: # QTNs start vs. end for adaptive inversions;  in FST start vs end
## Q4: average inversion characterists (age, length, num qtns/length) between adaptive, non-adaptive and no selection
## Q5: outlier data on false positives, true positives, etc.


  output.data <- c(seed, Va_perc_In, LA_final, num_inv, num_inv_NA, num_inv_NS,
                 capture_gain_p, capture_no_gain_p, neutral_gain_p, neutral_no_gain_p,
                 capture_gain_p.NA,  capture_no_gain_p.NA, neutral_gain_p.NA, neutral_no_gain_p.NA,
                 capture_gain_p.NS,  capture_no_gain_p.NS, neutral_gain_p.NS, neutral_no_gain_p.NS,
                 ave_start_QTNs, ave_start_QTNs_NA, ave_start_QTNs_NS,
                 ave_end_QTNs, ave_end_QTNs_NA, ave_end_QTNs_NS, 
                 ave_start_FST, ave_start_FST_NA, ave_start_FST_NS, 
                 ave_end_FST, ave_end_FST_NA, ave_end_FST_NS,  
                 ave_abV_start_qtnSelCoef, ave_abV_start_qtnSelCoef_NA, ave_abV_start_qtnSelCoef_NS,
                 ave_abV_end_qtnSelCoef, ave_abV_end_qtnSelCoef_NA, ave_abV_end_qtnSelCoef_NS, 
                 true_pos_pcadapt, true_pos_outflank, false_neg_pcadapt, false_neg_outflank, 
                 true_neg_pcadapt, true_neg_outflank, false_pos_pcadapt, false_pos_outflank,
                 true_neg_pcadapt_NS,  false_pos_pcadapt_NS, 
                 true_neg_outflank_NS, false_pos_outflank_NS)



  write.table(cbind(seed, df.invChar), paste0(folderOut, "outputInvChar_finalGen.txt"), col.names = F, row.names = F, append = T)
  write.table(cbind(seed, df.AdaptSplit), paste0(folderOut, "outputInvChar_allData.txt"), col.names = F, row.names = F, append = T)
  write.table(cbind(seed, df.crit.output), paste0(folderOut, "outputAdaptInvCrit.txt"), col.names = F, row.names = F, append = TRUE)

#### end Data output
######################################################################################################
## END IF STATEMENT FOR IF THIS IS A NO INVERSION SIMULATION
} else {
    output.data <- as.data.frame(c(seed, NA, LA_final, rep(NA, 45)))
  } 
} else {
  output.data <- as.data.frame(c(seed, NA, LA_final, rep(NA, 45)))
}
  write.table(t(output.data),  paste0(folderOut, "outputSumData.txt"), col.names = F, row.names = F, append = T)
  
  