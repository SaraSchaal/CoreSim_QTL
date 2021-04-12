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
#### Subset Inversions ####

### SELECTION SIMS ###
## Top 10 percent of FST values
top10.data <- df.invAllData %>%
  group_by(sim_gen) %>%
  filter(inv_FST>=quantile(inv_FST, 0.9)) %>%
  summarise_at(c("inv_age", "mean_qtnSelCoef", "num_qtns", "inv_length", "num_qtns_Lscaled"), 
               mean, .groups = "keep") %>%
  rename(inv_ageT10 = inv_age, mean_qtnSelCoefT10 = mean_qtnSelCoef, 
         num_qtnsT10 = num_qtns, inv_lengthT10 = inv_length, 
         num_qtns_LscaledT10 = num_qtns_Lscaled)

top10.data.nosum <- df.invAllData %>%
  group_by(sim_gen) %>%
  filter(inv_FST>=quantile(inv_FST, 0.9))

## Bottom 90 percent of FST values
bottom90.data <- df.invAllData %>%
  group_by(sim_gen) %>%
  filter(inv_FST<quantile(inv_FST, 0.9)) %>%
  summarise_at(c("inv_age", "mean_qtnSelCoef", "num_qtns", "inv_length", "num_qtns_Lscaled"), 
               mean, .groups = "keep") %>%
  rename(inv_ageB90 = inv_age, mean_qtnSelCoefB90 = mean_qtnSelCoef, 
         numqtnsB90 = num_qtns, inv_lengthB90 = inv_length,
         num_qtns_LscaledB90 = num_qtns_Lscaled) 

bot90.data.nosum <- df.invAllData %>%
  group_by(sim_gen) %>%
  filter(inv_FST<quantile(inv_FST, 0.9))


sd.Top10 <- aggregate(cbind(num_qtns, mean_qtnSelCoef, num_qtns, inv_length, num_qtns_Lscaled)~sim_gen, 
                      data = top10.data.nosum, FUN = sd)
colnames(sd.Top10)[2:6] <- c("sd_inv_ageT10", "sd_qtnSelCoefT10", "sd_num_qtnsT10", 
                             "sd_inv_lengthT10", "sd_num_qtns_LscaledT10")
sd.Bottom90 <- aggregate(cbind(num_qtns, mean_qtnSelCoef, num_qtns, inv_length, 
                               num_qtns_Lscaled)~sim_gen, data = bot90.data.nosum, FUN = sd)
colnames(sd.Bottom90)[2:6] <- c("sd_inv_ageB90", "sd_qtnSelCoefB90", "sd_num_qtnsB90", 
                                "sd_inv_lengthB90", "sd_num_qtns_LscaledB90")              

## Join dataframes with parameters
df.FSTsplitTb <- full_join(top10.data, bottom90.data, by = "sim_gen")

## convert to data frame and factor parameter columns
df.FSTsplitTemp <- as.data.frame(df.FSTsplitTb)  
df.FSTsplitTemp2 <- left_join(df.FSTsplitTemp, sd.Top10, by = "sim_gen")
df.FSTsplit <- left_join(df.FSTsplitTemp2, sd.Bottom90, by = "sim_gen")

### NO SELECTION SIMS ###
## Top 10 percent of FST values
top10.data.NS <- df.invAllData.NS %>%
  group_by(sim_gen) %>%
  filter(inv_FST>=quantile(inv_FST, 0.9)) %>%
  summarise_at(c("inv_age", "mean_qtnSelCoef", "num_qtns", "inv_length", "num_qtns_Lscaled"), 
               mean, .groups = "keep") %>%
  rename(inv_ageT10 = inv_age, mean_qtnSelCoefT10 = mean_qtnSelCoef, 
         num_qtnsT10 = num_qtns, inv_lengthT10 = inv_length, 
         num_qtns_LscaledT10 = num_qtns_Lscaled)

top10.data.nosum.NS <- df.invAllData.NS %>%
  group_by(sim_gen) %>%
  filter(inv_FST>=quantile(inv_FST, 0.9)) 

## Bottom 90 percent of FST values
bottom90.data.NS <- df.invAllData.NS %>%
  group_by(sim_gen) %>%
  filter(inv_FST<quantile(inv_FST, 0.9)) %>%
  summarise_at(c("inv_age", "mean_qtnSelCoef", "num_qtns", "inv_length", "num_qtns_Lscaled"), 
               mean, .groups = "keep") %>%
  rename(inv_ageB90 = inv_age, mean_qtnSelCoefB90 = mean_qtnSelCoef, 
         numqtnsB90 = num_qtns, inv_lengthB90 = inv_length, 
         num_qtns_LscaledB90 = num_qtns_Lscaled) 

bot90.data.nosum.NS <- df.invAllData.NS %>%
  group_by(sim_gen) %>%
  filter(inv_FST<quantile(inv_FST, 0.9)) 

sd.Top10.NS <- aggregate(cbind(inv_age, mean_qtnSelCoef, num_qtns, inv_length, 
                               num_qtns_Lscaled)~sim_gen, data = top10.data.nosum.NS, FUN = sd)
colnames(sd.Top10.NS)[2:6] <- c("sd_inv_ageT10", "sd_qtnSelCoefT10", "sd_num_qtnsT10", 
                                "sd_inv_lengthT10", "sd_num_qtns_LscaledT10")
sd.Bottom90.NS <- aggregate(cbind(inv_age, mean_qtnSelCoef, num_qtns, inv_length,
                                  num_qtns_Lscaled)~sim_gen, data = bot90.data.nosum.NS, FUN = sd)
colnames(sd.Bottom90.NS)[2:6] <- c("sd_inv_ageB90", "sd_qtnSelCoefB90", "sd_num_qtnsB90",
                                   "sd_inv_lengthB90", "sd_num_qtns_LscaledB90")

## Join dataframes with parameters
df.FSTsplit.NSTb <- full_join(top10.data.NS, bottom90.data.NS, by = "sim_gen")

## convert to data frame and factor parameter columns
df.FSTsplit.NStemp <- as.data.frame(df.FSTsplit.NSTb)  
df.FSTsplit.NStemp2 <- left_join(df.FSTsplit.NStemp, sd.Top10.NS, by = "sim_gen")
df.FSTsplit.NS <- left_join(df.FSTsplit.NStemp2, sd.Bottom90.NS, by = "sim_gen")

## end subset inversions
######################################################################################################




######################################################################################################
#### Inital full data plotting #####

# This function saves the legend into an object to use for plotting as an element in ggarrange
g_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

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
 # geom_line(data = df.popDyn, aes(y = meanPhenoP1 + sdPhenoP1, x = sim_gen), color = "cadetblue3", linetype = "dotted", alpha = 0.75) + 
 # geom_line(data = df.popDyn, aes(y = meanPhenoP1 - sdPhenoP1, x = sim_gen), color = "cadetblue3", linetype = "dotted", alpha = 0.75) +
  #geom_line(aes(y = meanPhenoP2, x = sim_gen), color = "navy") + 
 # geom_line(aes(y = meanPhenoP2 + sdPhenoP2, x = sim_gen), color = "navy", linetype = "dotted", alpha = 0.5) + 
 # geom_line(aes(y = meanPhenoP2 - sdPhenoP2, x = sim_gen), color = "navy", linetype = "dotted", alpha = 0.5) +
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
    geom_line(aes(y = meanPhenoP1 + sdPhenoP1, x = sim_gen), color = "cadetblue3", linetype = "dotted", alpha = 0.75) + 
    geom_line(aes(y = meanPhenoP1 - sdPhenoP1, x = sim_gen), color = "cadetblue3", linetype = "dotted", alpha = 0.75) +
    geom_line(aes(y = meanPhenoP2, x = sim_gen), color = "navy") + 
    geom_line(aes(y = meanPhenoP2 + sdPhenoP2, x = sim_gen), color = "navy", linetype = "dotted", alpha = 0.5) + 
    geom_line(aes(y = meanPhenoP2 - sdPhenoP2, x = sim_gen), color = "navy", linetype = "dotted", alpha = 0.5) +
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
    geom_line(aes(y = meanPhenoP1 + sdPhenoP1, x = sim_gen), color = "cadetblue3", linetype = "dotted", alpha = 0.75) + 
    geom_line(aes(y = meanPhenoP1 - sdPhenoP1, x = sim_gen), color = "cadetblue3", linetype = "dotted", alpha = 0.75) + 
    geom_line(aes(y = meanPhenoP2, x = sim_gen), color = "navy") + 
    geom_line(aes(y = meanPhenoP2 + sdPhenoP2, x = sim_gen), color = "navy", linetype = "dotted", alpha = 0.5) + 
    geom_line(aes(y = meanPhenoP2 - sdPhenoP2, x = sim_gen), color = "navy", linetype = "dotted", alpha = 0.5) +
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
    geom_line(aes(y = meanFitP1 + sdFitP1, x = sim_gen), color = "lightgrey") + 
    geom_line(aes(y = meanFitP1 - sdFitP1, x = sim_gen), color = "lightgrey") + 
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = "Population 1", y = "Fitness", x = "Generation") +
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
    geom_line(aes(y = meanFitP2 + sdFitP2, x = sim_gen), color = "lightgrey") + 
    geom_line(aes(y = meanFitP2 - sdFitP2, x = sim_gen), color = "lightgrey") + 
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = "Population 2", y = " ", x = "Generation") +
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
  # SELECTION #
  df.invage <- pivot_longer(df.FSTsplit[, c(1,2,7)], cols = c(inv_ageT10, inv_ageB90),
                            names_to = "FSTsplit", values_to = "inv_age")
  inv.age.plot <- ggplot(data = df.invage, 
                         aes(x = sim_gen, y = inv_age, group = FSTsplit)) + 
    geom_line(aes(color = FSTsplit), size = 0.75) + 
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = " ", y = "Average Inversion Age", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_color_manual(labels = c("Bot 90%", "Top 10%"), 
                       values=c("darkseagreen4", "darkseagreen3", "darkgrey", "darkgrey")) +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.FSTsplit$inv_ageT10))) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  
  # NO SELECTION #
  df.invage.NS <- pivot_longer(df.FSTsplit.NS[, c(1,2,7)], cols = c(inv_ageT10, inv_ageB90),
                               names_to = "FSTsplit", values_to = "inv_age")
  
  inv.age.plot.NS <- ggplot(data = df.invage.NS, 
                            aes(x = sim_gen, y = inv_age, group = FSTsplit)) + 
    geom_line(aes(color = FSTsplit), size = 0.75) + 
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = " ", y = "", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    
    scale_color_manual(labels = c("Top 10%", "Bot 90%"), 
                       values=c("darkseagreen3", "darkseagreen4")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.FSTsplit$inv_ageT10))) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  
  leg <- g_legend(inv.age.plot.NS)
  
  inv.age.plot.NS.noleg <- inv.age.plot.NS + theme(legend.position = "none")
  
  pdf(paste0(folderOut, seed, "_invAge.pdf"), height = 5, width = 7)
  ggarrange(inv.age.plot, inv.age.plot.NS.noleg, leg, 
            labels = c("Selection", "No Selection"), ncol = 3, widths = c(2.3,2.3,0.8))
  dev.off()

## Average Inversion Length ##
  # SELECTION #  
  df.invlength <- pivot_longer(df.FSTsplit[, c(1,5,10)], cols = c(inv_lengthT10, inv_lengthB90),
                               names_to = "FSTsplit", values_to = "inv_length")
  inv.length.plot <- ggplot(data = df.invlength, 
                            aes(x = sim_gen, y = inv_length, group = FSTsplit)) + 
    geom_line(aes(color = FSTsplit), size = 0.75, alpha = 0.9) + 
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = " ", y = "Average Inversion Length", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_color_manual(labels = c( "Top 10%", "Bot 90%"),
                       values=c("lightsalmon3", "lightsalmon1")) +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 50000)) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  
  
  # NO SELECTION #
  df.invlength.NS <- pivot_longer(df.FSTsplit.NS[, c(1,5,10)], cols = c(inv_lengthT10, inv_lengthB90),
                                  names_to = "FSTsplit", values_to = "inv_length")
  
  inv.length.plot.NS <- ggplot(data = df.invlength.NS, 
                               aes(x = sim_gen, y = inv_length, group = FSTsplit)) + 
    geom_line(aes(color = FSTsplit), size = 0.75, alpha= 0.9) + 
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = " ", y = "", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_color_manual(labels = c("Top 10%", "Bot 90%"),
                       values=c("lightsalmon1", "lightsalmon3")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 50000)) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  
  legLeng <- g_legend(inv.length.plot.NS)
  
  inv.length.plot.NS.noLeg <- inv.length.plot.NS + theme(legend.position = "none")
  
  pdf(paste0(folderOut, seed, "_invLength.pdf"), height = 5, width = 7)
  ggarrange(inv.length.plot, inv.length.plot.NS.noLeg, legLeng, labels = c("Selection", "No Selection"),
            ncol = 3, widths = c(2.3,2.3,0.8))
  dev.off()

## Average number of QTNs in Inverison ##
  # SELECTION #
  df.invQTNs <- pivot_longer(df.FSTsplit[, c(1,4,9)], cols = c(num_qtnsT10, numqtnsB90),
                             names_to = "FSTsplit", values_to = "inv_qtnNum")
  
  inv.qtns.plot <- ggplot(data = df.invQTNs, 
                          aes(x = sim_gen, y = inv_qtnNum, group = FSTsplit)) + 
    geom_line(aes(color = FSTsplit), size = 0.75) + 
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = "",
         y = "Average Number of inversion QTNs",
         x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_color_manual(labels = c( "Top 10%", "Bot 90%"), 
                       values=c( "thistle", "plum4")) +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.invQTNs$inv_qtnNum))) + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  
  
  # NO SELECTION #    
  df.invQTNs.NS <- pivot_longer(df.FSTsplit.NS[, c(1,4,9)], cols = c(num_qtnsT10, numqtnsB90),
                                names_to = "FSTsplit", values_to = "inv_qtnNum")
  
  inv.qtns.plot.NS <- ggplot(data = df.invQTNs.NS, 
                             aes(x = sim_gen, y = inv_qtnNum, group = FSTsplit)) + 
    geom_line(aes(color = FSTsplit), size = 0.75) + 
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = " ", y = "", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_color_manual(labels = c("Top 10%", "Bot 90%"), 
                       values=c("thistle", "plum4")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.invQTNs$inv_qtnNum))) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  
  legQTNs <- g_legend(inv.qtns.plot.NS)
  
  inv.qtns.plot.NS.noleg <- inv.qtns.plot.NS + theme(legend.position = "none")
  
  pdf(paste0(folderOut, seed, "_invQTNs.pdf"), height = 5, width = 7)
  ggarrange(inv.qtns.plot, inv.qtns.plot.NS.noleg, legQTNs, labels = c("Selection", "No Selection"),
            ncol = 3, widths = c(2.3,2.3,0.8))
  dev.off()
  
## Inversion QTNs Length Scaled ##
  # Selection
  df.invQTNsLscaled <- pivot_longer(df.FSTsplit[, c(1,6,11)], cols = c(num_qtns_LscaledT10, num_qtns_LscaledB90),
                                    names_to = "FSTsplit", values_to = "inv_qtnNum")
  
  inv.qtns.Lscaled.plot <- ggplot(data = df.invQTNsLscaled, 
                                  aes(x = sim_gen, y = inv_qtnNum, group = FSTsplit)) + 
    geom_line(aes(color = FSTsplit), size = 0.75) + 
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = "",
         y = "Average number of inversion QTNs \nscaled by inversion length",
         x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_color_manual(labels = c("Top 10%", "Bot 90%"), 
                       values=c(   "plum4", "thistle")) +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.invQTNsLscaled$inv_qtnNum))) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  
  # No Selection    
  df.invQTNs.Lscaled.NS <- pivot_longer(df.FSTsplit.NS[, c(1,6,11)], cols = c(num_qtns_LscaledT10, num_qtns_LscaledB90),
                                        names_to = "FSTsplit", values_to = "inv_qtnNum")
  
  inv.qtns.Lscaled.plot.NS <- ggplot(data = df.invQTNs.Lscaled.NS, 
                                     aes(x = sim_gen, y = inv_qtnNum, group = FSTsplit)) + 
    geom_line(aes(color = FSTsplit), size = 0.75) + 
    geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
    labs(title = " ", y = "", x = "Generation") +
    theme_classic() +
    theme(panel.background = element_blank(), 
          strip.background = element_rect(colour = "white", fill = "grey92")) +
    scale_color_manual(labels = c("Top 10%", "Bot 90%"), 
                       values=c("plum4", "thistle")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.invQTNsLscaled$inv_qtnNum))) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  
  legQTNsLscaled <- g_legend(inv.qtns.Lscaled.plot.NS)
  
  inv.qtns.Lscaled.plot.NS.noleg <- inv.qtns.Lscaled.plot.NS + theme(legend.position = "none")
  
  pdf(paste0(folderOut, seed, "_invQTNsLscaled.pdf"), height = 5, width = 7)
  ggarrange(inv.qtns.Lscaled.plot, inv.qtns.Lscaled.plot.NS.noleg, legQTNsLscaled, labels = c("Selection", "No Selection"),
            ncol = 3, widths = c(2.3,2.3,0.8))
  dev.off()
  
### end initial full data plotting
######################################################################################################
  
  
  
######################################################################################################
#### Calculate Inversion Windows ####
  
  # merge inversion data files together and remove the extra sim-gen column from metadata file
  df.invAllData <- left_join(df.invTime, df.invData[,c(1,3:ncol(df.invData))], by = "inv_id")
  df.invAllDataNS <- left_join(df.invTime.NS, df.invData.NS[,c(1,3:ncol(df.invData))], by = "inv_id")
  
  # filter the inversion data to get the final generation
  df.invDataFinalGen <- filter(df.invAllData, sim_gen == 50000)
  df.invDataFinalGenNS <- filter(df.invAllDataNS, sim_gen == 50000)
  
  # identify all the positions in the genome that are within the inversion windows
  invWind <- function(x){
    invWindowBases <- NULL
    numSegInv <- 0
    for(i in 1:nrow(x)) {
      invWindowBases <- c(invWindowBases, seq(from = x$inv_pos[i], 
                                              to = x$inv_end[i], by = 1))
      numSegInv <- numSegInv + 1
    }
    
    # identify all the positions in the genome that are within the inversion windows
    # filtered for MAF of 0.01
    invWindowBasesMAF <- NULL
    numSegInvMAF <- 0
    for(i in 1:nrow(x)) {
      Inv_wtype_freq <-  1 - x$freq[i]					
      MAF <-  min(x$freq[i], Inv_wtype_freq)				
      if(MAF > 0.01){
        invWindowBasesMAF <- c(invWindowBasesMAF, seq(from = x$inv_pos[i], 
                                                      to = x$inv_end[i], by = 1))
        numSegInvMAF <- numSegInvMAF + 1
      }
    }
    #segInv <- paste("number of segregating inversions =",  numSegInv, sep = " ")
    #segInvMAF <- paste("number of segregating inversions filtered for MAF 0.01 =",  numSegInvMAF, sep = " ")
    obj_list <- list("invWindBases" = invWindowBases, "invWindBasesMAF" = invWindowBasesMAF, 
                     "numSeqInv" = numSegInv, "numSeqInvMAF" = numSegInvMAF)
    return(obj_list)
  }
  
  invWindBases <- invWind(df.invDataFinalGen)[[1]]
  invWindBasesNS <- invWind(df.invDataFinalGenNS)[[1]]  
  
### end inversion window calculation
######################################################################################################  

  

######################################################################################################  
#### Identify which QTNs fall within inversion windows ####
  
  qtnMuts <- filter(df.muts, type == "m2")
  qtnMuts.NS <- filter(df.muts.NS, type == "m2")
  
  # identify which qtns overlap with inv window locations
  invWinQTNrows <- which(qtnMuts$position %in% invWindBases)
  invWinQTNrows.NS <- which(qtnMuts.NS$position %in% invWindBasesNS)
  
  # use this to identify how many positions have multiple qtns 
  mult.qtns <- length(which(qtnMuts$position %in% invWindBases)) - length(intersect(invWindBases, qtnMuts$position))
  mult.qtns.NS <- length(which(qtnMuts.NS$position %in% invWindBasesNS)) - length(intersect(invWindBasesNS, qtnMuts.NS$position))
  
  ## Selection
  for(i in 1:nrow(qtnMuts)){
    if(i %in% invWinQTNrows){
      qtnMuts$inOut[i] <- "in"
    } else {
      qtnMuts$inOut[i] <- "out"
    } 
  }
  
  neutMuts <- filter(df.muts, type == "m1")
  neutMuts$inOut <- "neut"
  
  df.neutQtnMutstemp <- rbind(qtnMuts, neutMuts)
  df.neutQtnMutstemp$FST <- as.numeric(as.character(df.neutQtnMutstemp$FST))
  df.neutQtnMuts <- df.neutQtnMutstemp[!is.nan(df.neutQtnMutstemp$FST),]
  
  ## No selection
  for(i in 1:nrow(qtnMuts.NS)){
    if(i %in% invWinQTNrows.NS){
      qtnMuts.NS$inOut[i] <- "in"
    } else {
      qtnMuts.NS$inOut[i] <- "out"
    } 
  }
  
  neutMuts.NS <- filter(df.muts.NS, type == "m1")
  neutMuts.NS$inOut <- "neut"
  
  
  df.neutQtnMuts.NStemp <- rbind(qtnMuts.NS, neutMuts.NS)
  df.neutQtnMuts.NStemp$FST <- as.numeric(as.character(df.neutQtnMuts.NStemp$FST))
  df.neutQtnMuts.NS <- df.neutQtnMuts.NStemp[!is.nan(df.neutQtnMuts.NStemp$FST),]

### end identifying inversion window QTNs
######################################################################################################  
  
######################################################################################################  
#### Manhattan plot ####
  # SELECTION #
  manh.plot <- ggplot(df.neutQtnMuts, aes(x = position, y = FST, group = inOut)) + 
    geom_point(aes(color = inOut)) + 
    scale_color_manual(values=c( "red", "goldenrod", "navy")) +
    xlim(0, 2100000) + 
    ylim(0, max(df.neutQtnMuts$FST)) +
    labs(title = "Selection") + 
    theme(legend.position = "none")
  
  # NO SELECTION #
  manh.plot.NS <- ggplot(df.neutQtnMuts.NS, aes(x = position, y = FST, group = inOut)) + 
    geom_point(aes(color = inOut)) + 
    scale_color_manual(values=c( "red", "goldenrod", "navy")) +
    xlim(0, 2100000) + 
    ylim(0, max(df.neutQtnMuts$FST)) +
    labs(title = "No Selection") 
  
  # No legend
  manh.plot.NS.noleg <- manh.plot.NS + theme(legend.position = "none")
  
  legManh <- g_legend(manh.plot.NS)
  
  ## TO DO: ADD CHROMOSOME COLORS
  png(paste0(folderOut, seed, "_manh.png"), width = 700, height = 700, units = "px")
  ggarrange(manh.plot, manh.plot.NS.noleg, legManh, ncol = 3, widths = c(2.3,2.3,0.8))
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
    geom_point(aes(color = pop, size = inv_FST), alpha = 0.8) + 
    geom_line(aes(color = pop, group = inv_id), alpha = 0.8) + 
    scale_color_manual(values=c("navy", "red")) + 
    scale_size(range = c(0.5, 4), breaks = c(0.00001, 0.05, 0.15, 0.2)) + 
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
  
  # example of how to find a specific mutation in the vcf file
  df.mutsMAF[2,]
  vcf.MAF@fix[grep(df.mutsMAF$mutID[1], vcf.MAF@fix[,"INFO"]),]
  
  geno <- vcf.MAF@gt[,-1] # this gets individual ids and genotypes
  position <- getPOS(vcf.MAF) # this gets position of mutations
  
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
  length(regmatches(meta, regexpr("[0-9]+[0-9]", meta)))
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
  vcf_pos <- as.numeric(vcf.MAF@fix[,"POS"])
  hist(vcf_pos, breaks=seq(0,2100000, length.out=100))
  hist(df.mutsMAF$position, breaks=seq(0,2100000, length.out=100))
  
  # sort mutations and make sure vcf matches slim output
  # slim output is one base off so add a base to match vcf position
  head(sort(df.mutsMAF$position))
  head(vcf_pos)
  df.mutsMAF$position_vcf <- df.mutsMAF$position + 1
  df.mutsMAF$is_vcf <- NA
  
  # Check to make sure FST values match between files
  # this is slow, but correct
  G_FST <- rep(NA, nrow(G)) 
  count <- 0
  for (i in 1:nrow(df.mutsMAF)){
    x <- grep(df.mutsMAF$mutID[i], vcf.MAF@fix[,"INFO"])
    if (length(x)==1){
      G_FST[x] <- df.mutsMAF$FST[i]
    }
    count <- count+1
    print(count)
  }
  
  hist(df.mutsMAF$FST, breaks = 50)
  head(G_FST)
  hist(G_FST, breaks = 50)
  
  sum(is.na(G_FST))
  # mutations that don't match up
  
  length(a)
  length(G_FST)
  
  head(df.mutsMAF)
  df.mutsMAFord<- df.mutsMAF[order(df.mutsMAF$position),]
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
  
  pdf(paste0("figures/", seed, "_heatmapPop1alphaFST.pdf"), height = 5, width = 7)
  heatmap(t(G1_alpha[, pop1_order]),   
          main="Pop1 G*a+-*FST",cexCol = 0.3,
          Colv = NA, useRaster=TRUE,
          scale="none",
          # Rowv = NA, 
          col=two.colors(100, start = "blue", end="red", middle="white"))
  # ADDING BREAKS SCREWS UP EVERYTHING
  #breaks=seq(-0.005, 0.005, length.out = 101))
  dev.off()
  
  pdf(paste0("figures/", seed, "_heatmapPop2alphaFST.pdf"), height = 5, width = 7)
  
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

  pdf(paste0("figures/", seed, "_heatmapPop1geno.pdf"), height = 5, width = 7)
  
  heatmap(t(G_ref1[,pop1_order]), Rowv = NA,  main="Pop1 genotypes",cexCol = 0.3,
          Colv = NA, useRaster=TRUE,
          scale="none")
  dev.off()
  
  pdf(paste0("figures/", seed, "_heatmapPop1geno.pdf"), height = 5, width = 7)
  heatmap(t(G_ref2[,pop2_order]), Rowv = NA,  main="Pop2 genotypes",cexCol = 0.3,
          Colv = NA, useRaster=TRUE,
          scale="none")
  dev.off()
#### end plot heatmanps
######################################################################################################    

######################################################################################################
#### BigSnpr ####
## This chunk of code finds a quasi-independent set of SNPs
# Identify chromosome ID for each mutation
chromosome <- getCHROM(vcf.MAF)

# list the G matrix with the position of each mutation and the chromosome 
# it is on. This is used for
training <- list(G = G, 
                 position = position,
                 chromosome = chromosome)

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

# These are the indexes of the quasi-independent 
# set of loci that are kept after pruning for LD
which_pruned = attr(newpc, which = "subset")

training$G_coded <- G_coded
training$G_pruned <- training$G[which_pruned,]
training$which_pruned <- which_pruned

df.mutsMAF$quasi_indep <- FALSE
df.mutsMAF$quasi_indep[training$which_pruned] <- TRUE

#### end bigsnpr
######################################################################################################
  
######################################################################################################
#### PCADAPT
gename <- paste0(seed, "_genotypes.lfmm")
write.lfmm(t(training$G), gename)
pcafile <- read.pcadapt(gename, type="lfmm")
pca_all <- pcadapt(pcafile,K=3)
head(pca_all$loadings)
df.mutsMAF$pca_ALL_PC1_loadings <- pca_all$loadings[,1]
df.mutsMAF$pca_ALL_PC2_loadings <- pca_all$loadings[,2]
df.mutsMAF$pca_ALL_PC3_loadings <- pca_all$loadings[,3]
head(df.mutsMAF)


### PCA loadings if pruned data is used ####
gename2 <- paste0(seed, "_genotypes_pruned.lfmm")
write.lfmm(t(training$G_pruned), gename2)
pcafile2 <- read.pcadapt(gename2, type="lfmm")
pca_pruned <- pcadapt(pcafile2,K=3)

df.mutsMAF$pca_PRUNED_PC1_loadings <- NA 
df.mutsMAF$pca_PRUNED_PC2_loadings <- NA
df.mutsMAF$pca_PRUNED_PC1_loadings[training$which_pruned] <- pca_pruned$loadings[,1]
df.mutsMAF$pca_PRUNED_PC2_loadings[training$which_pruned] <- pca_pruned$loadings[,2]  

cor.test(df.mutsMAF$pca_ALL_PC1_loadings, df.mutsMAF$pca_PRUNED_PC1_loadings)
cor.test(df.mutsMAF$pca_ALL_PC2_loadings, df.mutsMAF$pca_PRUNED_PC2_loadings)

df.mutsMAF$pcadapt_4.3.3_ALL_chisq <- as.numeric(pca_all$chi2.stat)
df.mutsMAF$pcadapt_4.3.3_ALL_log10p <- -log10(pca_all$pvalues)

plot(df.mutsMAF$position, df.mutsMAF$pcadapt_4.3.3_ALL_chisq)
plot(df.mutsMAF$position, df.mutsMAF$pcadapt_4.3.3_ALL_log10p)

test <- snp_gc(snp_pcadapt(training$G_coded, U.row = newpc$u[,1]))
df.mutsMAF$pcadapt_4.3.3_PRUNED_log10p <- -predict(test,log10=T)

plot(df.mutsMAF$position, df.mutsMAF$pcadapt_4.3.3_PRUNED_log10p )
cor.test(df.mutsMAF$pcadapt_4.3.3_ALL_log10p, df.mutsMAF$pcadapt_4.3.3_PRUNED_log10p, method = "spearman")
plot(df.mutsMAF$position,df.mutsMAF$pcadapt_4.3.3_PRUNED_log10p)

#### end PCADAPT
######################################################################################################



######################################################################################################
#### OutFLANK

FstDataFrame <- MakeDiploidFSTMat(t(training$G),final_df$vcf_ord,
                                  ind$group[ind$infinal])
out_ini <- OutFLANK(FstDataFrame, NumberOfSamples=k) 
str(out_ini)

out_pruned <- OutFLANK(FstDataFrame[training$which_pruned,], NumberOfSamples=k)     
str(out_pruned)

P1 <- pOutlierFinderChiSqNoCorr(FstDataFrame, 
                                Fstbar = out_pruned$FSTNoCorrbar, 
                                dfInferred = out_pruned$dfInferred, Hmin=0.1)
P1 <- P1[order(P1$LocusName),]
identical(P1$LocusName, df.mutsMAF$vcf_ord) # should be TRUE
identical(out_ini$results$LocusName, df.mutsMAF$vcf_ord) # should be TRUE
forfinal <- data.frame(vcf_ord = out_ini$results$LocusName,
                       OutFLANK_0.2_FST = out_ini$results$FST,
                       OutFLANK_0.2_He = out_ini$results$He,
                       OutFLANK_0.2_ALL_log10p = -log10(out_ini$results$pvaluesRightTail),
                       OutFLANK_0.2_PRUNED_log10p = -log10(P1$pvaluesRightTail))
df.mutsMAF_temp <- merge(df.mutsMAF, forfinal, by="vcf_ord", all.x=TRUE)
dim(df.mutsMAF)
dim(df.mutsMAF_temp)
head(df.mutsMAF_temp)

#### end OutFLANK
######################################################################################################





######################################################################################################    
## COPY AND PASTE WHERE NEEDED
pdf(paste0("figures/", seed, "_heatmapPop1geno.pdf"), height = 5, width = 7)

dev.off()

png(paste0("figures/", seed, "XXXX.png"), width = 480, height = 480, units = "px")

dev.off()