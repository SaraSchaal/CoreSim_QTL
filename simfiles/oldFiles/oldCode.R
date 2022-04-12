#########################################################################################################
### start outlier detection
# df.outlierComp$testNum <- NA
# df.outlierComp$whichTests <- NA
# for(i in 1:nrow(df.outlierComp)){
#   if(!is.na(df.outlierComp$OutlierFlag[i]) & !is.na(df.outlierComp$pcadapt_outlier[i])){
#     if(df.outlierComp$crit1_p.value[i] != 0 | df.outlierComp$crit2_p.value[i] != 0 & df.outlierComp$OutlierFlag[i] == FALSE & df.outlierComp$pcadapt_outlier[i] == FALSE){
#       df.outlierComp$testNum[i] <- 0
#       df.outlierComp$whichTests[i] <- "none"
#     } else if(df.outlierComp$crit1_p.value[i] == 0 & df.outlierComp$crit2_p.value[i] == 0 & df.outlierComp$OutlierFlag[i] == TRUE & df.outlierComp$pcadapt_outlier[i] == TRUE){
#       df.outlierComp$testNum[i] <- 3
#       df.outlierComp$whichTests[i] <- "crit_OutFLANK_PCAdapt"
#     } else if(df.outlierComp$crit1_p.value[i] == 0 & df.outlierComp$crit2_p.value[i] == 0 & df.outlierComp$OutlierFlag[i] == TRUE & df.outlierComp$pcadapt_outlier[i] == FALSE){
#       df.outlierComp$testNum[i] <- 2
#       df.outlierComp$whichTests[i] <- "crit_OutFLANK"
#     } else if(df.outlierComp$crit1_p.value[i] == 0 & df.outlierComp$crit2_p.value[i] == 0 & df.outlierComp$OutlierFlag[i] == FALSE & df.outlierComp$pcadapt_outlier[i] == TRUE){
#       df.outlierComp$testNum[i] <- 2
#       df.outlierComp$whichTests[i] <- "crit_PCAdapt"
#     } else if(df.outlierComp$crit1_p.value[i] != 0 | df.outlierComp$crit2_p.value[i] != 0 & df.outlierComp$OutlierFlag[i] == TRUE & df.outlierComp$pcadapt_outlier[i] == FALSE){
#       df.outlierComp$testNum[i] <- 1
#       df.outlierComp$whichTests[i] <- "OutFLANK"
#     } else if(df.outlierComp$crit1_p.value[i] != 0 | df.outlierComp$crit2_p.value[i] != 0 & df.outlierComp$OutlierFlag[i] == FALSE & df.outlierComp$pcadapt_outlier[i] == TRUE){ 
#       df.outlierComp$testNum[i] <- 1
#       df.outlierComp$whichTests[i] <- "PCAdapt"
#     }
#   }
# }
# 
# df.outlierComp$whichTests <- as.factor(df.outlierComp$whichTests)
# df.outlierComp$whichTests <- relevel(df.outlierComp$whichTests, "none")
# levels(df.outlierComp$inOut)
# overall.outliers <- ggplot(df.outlierComp[!is.na(df.outlierComp$whichTests),], aes(x = position_vcf, y = FST_slim)) + 
#   geom_point(aes(color = whichTests), alpha = 0.8) +
#   facet_wrap(~inOut, nrow = 4, ncol = 1) + 
#   scale_color_manual(values = c("black", "blue", "red", "goldenrod")) +
#   theme_classic() +
#   theme(panel.background = element_blank(), 
#         strip.background = element_rect(colour = "white", fill = "grey92"),
#         text = element_text(size = 11)) +
#   labs(title = "Outlier Comparison",
#        y = "FST",
#        x = "Genome Position") + 
#   ylim(c(0, 0.3))
### end outlier detection
#########################################################################################################





#########################################################################################################
### start outflank fst plot
# 
# outflank.fst <- ggplot(df.out, aes(x = position_vcf, y = FST_outflank)) +
#   geom_rect(data=df.adaptInv, mapping=aes(xmin=first_bases, xmax=final_bases, ymin=0,
#                                           ymax=max(c(df.out$FST_outflank, df.neutQTNmuts$FST)) + 0.05), fill = "tan1", 
#             color="black", alpha=0.5, inherit.aes = FALSE) +
#   geom_point(aes(color = OutlierFlag, shape = OutlierFlag))+
#   scale_color_manual(values = c("black", "red")) +
#   scale_shape_manual(values = c(19, 1)) +
#   theme_classic() +
#   theme(panel.background = element_blank(), 
#         strip.background = element_rect(colour = "white", fill = "grey92"),
#         text = element_text(size = 11)) +
#   labs(title ="OutFLANK",
#        y = "FST",
#        x = "Genome Position") + 
#   theme(legend.position = "none") + 
#   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
#   scale_y_continuous(expand = c(0, 0), limits = c(0,max(c(df.out$FST_outflank, df.neutQTNmuts$FST)) + 0.05))

### end outflank fst plot
#########################################################################################################


#########################################################################################################
### start arrows on plot

# annotate("segment", x = center.bases, y = rep(max(c(df.out$FST_outflank, df.neutQTNmuts$FST)), length(center.bases)), 
#xend = center.bases, yend = rep(max(df.neutQTNmuts$FST + 0.01), length(center.bases)),
#arrow = arrow(length = unit(0.5, "cm")))

### end arrows on plot
#########################################################################################################


#########################################################################################################
### Plots on inversions split by FST 
# ## Average Inversion Length ##
# # SELECTION #  
# df.invlength <- pivot_longer(df.FSTsplit[, c(1,5,10)], cols = c(inv_lengthAdapt, inv_lengthNonAdapt),
#                              names_to = "FSTsplit", values_to = "inv_length")
# 
# inv.length.plot <- ggplot(data = df.invlength, 
#                           aes(x = sim_gen, y = inv_length, group = FSTsplit)) + 
#   geom_line(aes(color = FSTsplit), size = 0.75, alpha = 0.9) + 
#   geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
#   labs(title = " ", y = "Average Inversion Length", x = "Generation") +
#   theme_classic() +
#   theme(panel.background = element_blank(), 
#         strip.background = element_rect(colour = "white", fill = "grey92")) +
#   scale_color_manual(labels = c( "NonAdaptive", "Adaptive"),
#                      values=c("lightsalmon3", "lightsalmon1")) +
#   theme(legend.position = "none") +
#   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 50000)) +
#   theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
# 
# 
# # NO SELECTION #
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
# pdf(paste0(folderOut, seed, "_invLength.pdf"), height = 5, width = 7)
# ggarrange(inv.length.plot, inv.length.plot.NS.noLeg, legLeng, labels = c("Selection", "No Selection"),
#           ncol = 3, widths = c(2.3,2.3,0.8))
# dev.off()
# 
# 
# ggplot(data = final.inv, aes(x = adaptInv, y= inv_length, fill = adaptInv)) +
#   geom_boxplot() + 
#   scale_fill_manual(values = c("paleturquoise2","skyblue", "skyblue4")) + 
#   theme_classic() +
#   theme(panel.background = element_blank(), 
#         strip.background = element_rect(colour = "white", fill = "grey92")) +
#   labs(title = " ", y = "Average Inversion Length", x = "Inversion Status") 


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

# ggarrange(inv.length.plot, inv.length.plot.NS.noLeg, legLeng, labels = c("Selection", "No Selection"),
#           ncol = 3, widths = c(2.3,2.3,0.8))


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

# ggarrange(inv.qtns.Lscaled.plot, inv.qtns.Lscaled.plot.NS.noleg, legQTNsLscaled, labels = c("Selection", "No Selection"),
#           ncol = 3, widths = c(2.3,2.3,0.8))

### End Plots on inversions split by FST
#########################################################################################################



#########################################################################################################
### Plots on inversion outlier
# crit 1
# ggplot(df.qtnMuts.MAF[df.qtnMuts.MAF$type == "m2",], aes(x = crit1_p.value, fill = inOut)) +
#   geom_histogram(color = "black") +
#   facet_grid(inOut~.) + 
#   theme(legend.position = "n") +
#   scale_fill_manual(values = c("red", "blue", "gold")) + 
#   labs(title = "Criteria 1 - compare to null of inversion QTNs \nno-selection simulation",
#        x = "empirical p-value")
# 
# # crit 2
# ggplot(df.qtnMuts.MAF[df.qtnMuts.MAF$type == "m2",], aes(x = crit2_p.value, fill = inOut)) +
#   geom_histogram(color = "black") +
#   facet_grid(inOut~.) + 
#   theme(legend.position = "n") +
#   scale_fill_manual(values = c("red", "blue")) + 
#   labs(title = "Criteria 2 - compare to null of neutral QTNs \nselection simulation",
#        x = "empirical p-value")
# 
# # crit 3
# ggplot(df.qtnMuts.MAF[df.qtnMuts.MAF$type == "m2",], aes(x = crit3_Va_perc, fill = inOut)) +
#   geom_histogram( colour = "black") +
#   facet_grid(inOut~.) +
#   theme(legend.position = "n") +
#   scale_fill_manual(values = c("red", "blue")) +
#   labs(title = "Criteria 3 - compare additive genetic variance percent within \nselection simulation",
#        x = "percent of Va explained")
# # crit 3
# ggplot(df.qtnMuts.MAF[df.qtnMuts.MAF$type == "m2" & df.qtnMuts.MAF$crit3_Va_perc >= 0.01, ],
#        aes(x = crit3_Va_perc, fill = inOut)) +
#   geom_histogram( colour = "black") +
#   facet_grid(inOut~.) +
#   theme(legend.position = "n") +
#   scale_fill_manual(values = c("red", "blue")) +
#   labs(title = "Criteria 3 - compare additive genetic variance percent within \nselection simulation subset for > 0.01",
#        x = "percent of Va explained")
# 


### End Plots
#########################################################################################################


#########################################################################################################
### Notes on adaptive inversions


## First we need to identify adaptive QTNs
# 1st QTN criteria -- addresses genetic drift causing false positive inversion outliers 
# (comparing selection sim to no selec)
# a) get the qtns from the no selection simulation inside inversions and record the FST values as our
# null distribution
# b) loop through all qtns in the final generation of selection simulation to get empirical p


# #) 2nd criteria -- addresses genetic drift in the same simulation (within selection sim)
# a) get all the neutral loci on the final chromosome and get the distribution of their FST values (calculate this outside the loop that I am stepping through for QTNs)
# b) for the qtns in the focal qtn (selection sim) 
# c) record the empirical p value :  1-rank(c(null, obs))[length(null)+1]/(length(null)+1)

# 3rd criteria -- added genetic variation for each mutation and the proportion/percent of added genetic variation
# #) sanity check: this should be zero for neutral mutations

# summarize these criteria for each inversion
# 1st inv criteria 
#) loop through the inversions themselves 
# average -log10p in each inversion compared to the average -log10p for everything that is non-inverted (blue manh)
# count the number of FST outliers with p = 0 / inversion window size is > the number of FST outliers in blue / total
# blue region (do not include the neutral chromosome in that calculation)

# 2nd criteria
# average -log10p in each inversion compared to the average -log10p for qtns on neutral chromosome
# count the number of FST outliers with p = 0 / inversion window size is > the number of FST outliers in neut / total
# neut region 



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


#############################################################################################################
### adaptive inversion split for No selection sims 
#
# adapt.inv.data.NS <- df.invAllData.NS %>%
#    group_by(sim_gen) %>%
#    filter(inv_id %in% adapt.inv.NS) %>%
#    summarise_at(c("inv_age", "mean_qtnSelCoef", "num_qtns", "inv_length", "num_qtns_Lscaled"), 
#                 mean, .groups = "keep") %>%
#    rename(inv_ageAdapt = inv_age, mean_qtnSelCoefAdapt = mean_qtnSelCoef, 
#           num_qtnsAdapt = num_qtns, inv_lengthAdapt = inv_length, 
#           num_qtns_LscaledAdapt = num_qtns_Lscaled)
#  
#  adap.inv.data.nosum.NS <- df.invAllData.NS %>%
#    group_by(sim_gen) %>%
#    filter(inv_id %in% adapt.inv)
#  
#  ## Bottom 90 percent of FST values
#  non.adap.inv.data.NS <- df.invAllData.NS %>%
#    group_by(sim_gen) %>%
#    filter(!inv_id %in% adapt.inv) %>%
#    summarise_at(c("inv_age", "mean_qtnSelCoef", "num_qtns", "inv_length", "num_qtns_Lscaled"), 
#                 mean, .groups = "keep") %>%
#    rename(inv_ageB90 = inv_age, mean_qtnSelCoefB90 = mean_qtnSelCoef, 
#           numqtnsB90 = num_qtns, inv_lengthB90 = inv_length, 
#           num_qtns_LscaledB90 = num_qtns_Lscaled) 
#  
# non.adap.inv.data.nosum.NS <- df.invAllData.NS %>%
#   group_by(sim_gen) %>%
#   filter(!inv_id %in% adapt.inv) 
# 
# sd.Top10.NS <- aggregate(cbind(inv_age, mean_qtnSelCoef, num_qtns, inv_length, 
#                                num_qtns_Lscaled)~sim_gen, data = top10.data.nosum.NS, FUN = sd)
# colnames(sd.Top10.NS)[2:6] <- c("sd_inv_ageT10", "sd_qtnSelCoefT10", "sd_num_qtnsT10", 
#                                 "sd_inv_lengthT10", "sd_num_qtns_LscaledT10")
# sd.Bottom90.NS <- aggregate(cbind(inv_age, mean_qtnSelCoef, num_qtns, inv_length,
#                                   num_qtns_Lscaled)~sim_gen, data = bot90.data.nosum.NS, FUN = sd)
# colnames(sd.Bottom90.NS)[2:6] <- c("sd_inv_ageB90", "sd_qtnSelCoefB90", "sd_num_qtnsB90",
#                                    "sd_inv_lengthB90", "sd_num_qtns_LscaledB90")
# 
# ## Join dataframes with parameters
# df.FSTsplit.NSTb <- full_join(top10.data.NS, bottom90.data.NS, by = "sim_gen")
# 
# ## convert to data frame and factor parameter columns
# df.FSTsplit.NStemp <- as.data.frame(df.FSTsplit.NSTb)  
# df.FSTsplit.NStemp2 <- left_join(df.FSTsplit.NStemp, sd.Top10.NS, by = "sim_gen")
# df.FSTsplit.NS <- left_join(df.FSTsplit.NStemp2, sd.Bottom90.NS, by = "sim_gen")

### CLOSE CODE FOR AVERAGING REPS ###
#############################################################################################################


#############################################################################################################
## Average number of QTNs in Inverison ##
# SELECTION #
# df.invQTNs <- pivot_longer(df.FSTsplit[, c(1,4,9)], cols = c(num_qtnsT10, numqtnsB90),
#                            names_to = "FSTsplit", values_to = "inv_qtnNum")
# 
# inv.qtns.plot <- ggplot(data = df.invQTNs, 
#                         aes(x = sim_gen, y = inv_qtnNum, group = FSTsplit)) + 
#   geom_line(aes(color = FSTsplit), size = 0.75) + 
#   geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
#   labs(title = "",
#        y = "Average Number of inversion QTNs",
#        x = "Generation") +
#   theme_classic() +
#   theme(panel.background = element_blank(), 
#         strip.background = element_rect(colour = "white", fill = "grey92")) +
#   scale_color_manual(labels = c( "Top 10%", "Bot 90%"), 
#                      values=c( "thistle", "plum4")) +
#   theme(legend.position = "none") +
#   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
#   scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.invQTNs$inv_qtnNum))) + 
#   theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
# 
# 
# # NO SELECTION #    
# df.invQTNs.NS <- pivot_longer(df.FSTsplit.NS[, c(1,4,9)], cols = c(num_qtnsT10, numqtnsB90),
#                               names_to = "FSTsplit", values_to = "inv_qtnNum")
# 
# inv.qtns.plot.NS <- ggplot(data = df.invQTNs.NS, 
#                            aes(x = sim_gen, y = inv_qtnNum, group = FSTsplit)) + 
#   geom_line(aes(color = FSTsplit), size = 0.75) + 
#   geom_vline(xintercept = 10000, linetype = "dashed", color = "black") +
#   labs(title = " ", y = "", x = "Generation") +
#   theme_classic() +
#   theme(panel.background = element_blank(), 
#         strip.background = element_rect(colour = "white", fill = "grey92")) +
#   scale_color_manual(labels = c("Top 10%", "Bot 90%"), 
#                      values=c("thistle", "plum4")) +
#   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
#   scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.invQTNs$inv_qtnNum))) +
#   theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
# 
# legQTNs <- g_legend(inv.qtns.plot.NS)
# 
# inv.qtns.plot.NS.noleg <- inv.qtns.plot.NS + theme(legend.position = "none")
# 
# pdf(paste0(folderOut, seed, "_invQTNs.pdf"), height = 5, width = 7)
# ggarrange(inv.qtns.plot, inv.qtns.plot.NS.noleg, legQTNs, labels = c("Selection", "No Selection"),
#           ncol = 3, widths = c(2.3,2.3,0.8))
# dev.off()
### END num qtns in inversion 
#############################################################################################################

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
#### Run the following code chunk once to get full data files to do further analyses on ####
# df <- read.table("figures/20210514_noAdaptInv/outputAdaptInvCrit.txt")
# ## Inversion Through Time
# df.invTime <- NULL
# count <- 0
# for(i in 1:nrow(df.params)){
#   seed <- df.params$seed[i]
#   invTimeNewFile <- read.table(paste(folder, seed, "_outputInvTime.txt", sep=""), header = TRUE, stringsAsFactors = FALSE)
#   if(nrow(invTimeNewFile) > 0){
#     invTimeNewFile$seed <- seed
#     df.invTime <-  rbind(df.invTime, invTimeNewFile)
#   }
#   count <- count + 1
#   print(count)
# }
# write.table(df.invTime, "FullSet_invTime.txt", row.names = FALSE)
# 
# ## Inversion Summary Data
# df.invData <- NULL
# count <- 0
# for(i in 1:nrow(df.params)){
#   seed <- df.params$seed[i]
#   invData <- read.table(paste(folder, seed, "_outputInvSumInfo.txt", sep=""), header = TRUE,
#                         stringsAsFactors = FALSE)
#   if(nrow(invData) > 0){
#     invData$seed <- seed
#     df.invData <- rbind(df.invData, invData)
#   }
#   count <- count + 1
#   print(count)
# }
# write.table(df.invData, "FullSet_invData.txt", row.names = FALSE)
# 
# ## QTN Inversion Through Time
# df.invQTNTime <- NULL
# count <- 0
# for(i in 1:nrow(df.params)){
#   seed <- df.params$seed[i]
#   invQTNTimeNewFile <- read.table(paste(folder, seed, "_outputInvQtn.txt", sep=""), header = TRUE, stringsAsFactors = FALSE)
#   if(nrow(invQTNTimeNewFile) > 0){
#     invQTNTimeNewFile$seed <- seed
#     df.invQTNTime <-  rbind(df.invQTNTime, invQTNTimeNewFile)
#   }
#   count <- count + 1
#   print(count)
# }
# write.table(df.invQTNTime, "FullSet_invQTNTime.txt", row.names = FALSE)
# 
# ## QTN Inversion Summary Data
# df.invQTNData <- NULL
# count <- 0
# for(i in 1:nrow(df.params)){
#   seed <- df.params$seed[i]
#   invQTNDataNewFile <- read.table(paste(folder, seed, "_outputInvQtnSumInfo.txt", sep=""), header = TRUE, stringsAsFactors = FALSE)
#   if(nrow(invQTNDataNewFile) > 0){
#     invQTNDataNewFile$seed <- seed
#     df.invQTNData <-  rbind(df.invQTNData, invQTNDataNewFile)
#   }
#   count <- count + 1
#   print(count)
# }
# write.table(df.invQTNData, "FullSet_invQTNSumInfo.txt", row.names = FALSE)
# 
## Pop Dynamics
# df.popDyn <- NULL
# count <- 0
# for(i in 1:nrow(df.params)){
#   seed <- df.params$seed[i]
#   popDynNewFile <- read.table(paste(folder, seed, "_outputPopDynam.txt", sep=""), header = TRUE, stringsAsFactors = FALSE)
#   if(nrow(popDynNewFile) > 0){
#     popDynNewFile$seed <- seed
#     df.popDyn <-  rbind(df.popDyn, popDynNewFile)
#   }
#   count <- count + 1
#   print(count)
# }
# write.table(df.popDyn, "FullSet_popDyn.txt", row.names = FALSE)
# # 
# # ## SimStats
# df.simStats <- NULL
# count <- 0
# for(i in 1:nrow(df.params)){
#    seed <- df.params$seed[i]
#    simStatsNewFile <- read.table(paste(folder, seed, "_outputSimStats.txt", sep=""), stringsAsFactors = FALSE)
#    simStatsNewFile$seed <- seed
#    df.simStats <-  rbind(df.simStats, simStatsNewFile)
#    count <- count + 1
#    print(count)
#  }
#  df.simStats <- df.simStats[,2:ncol(df.simStats)]
#  colnames(df.simStats) <- c("mig1", "mig2", "pop1N", "pop2N", "mu_base", "mu_inv", "r", "alpha", "sigmaK", "burnin", "dom", "enVar", "Seed")
#  write.table(df.simStats, "FullSet_simStats.txt", row.names = FALSE)
# 
# ## Mutations File
# df.finalMuts <- NULL
# count <- 0
# no.Data <- NULL
# for(i in 1:nrow(df.params)){
#   seed <- df.params$seed[i]
#   if(file.exists(paste(folder, seed, "_outputMutations.txt", sep=""))){
#     finalMutsNewFile <- read.table(paste(folder, seed, "_outputMutations.txt", sep=""), stringsAsFactors = FALSE, header = TRUE)
#     finalMutsNewFile$seed <- seed
#     df.finalMuts <-  rbind(df.finalMuts, finalMutsNewFile)
#   } else {
#     no.Data <- c(no.Data, seed)
#   }
#   count <- count + 1
#   print(count)
# }
# colnames(df.finalMuts)[2] <- "mut_id" 
# write.table(df.finalMuts, "FullSet_finalMuts.txt", row.names = FALSE)


