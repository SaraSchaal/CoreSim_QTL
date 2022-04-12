DATADIR <- "/scratch/schaal.s/InversionSimulations/results/20220220/10gen"
folder <- "./results/Inversion/20220302_reviews/10gen/"
PARAMSDIR <- "/scratch/schaal.s/InversionSimulations/src/"
DATAOUT <- "/scratch/schaal.s/InversionSimulations/figures/"

df.params <- read.table(paste0(PARAMSDIR, "FullSet_dfparams.txt"))

df.popDyn <- NULL
count <- 0
for(i in 1:nrow(df.params)){
  seed <- df.params$seed[i]
  popDynNewFile <- read.table(paste(DATADIR, seed, "_outputPopDynam.txt", sep=""), header = TRUE, stringsAsFactors = FALSE)
  if(nrow(popDynNewFile) > 0){
    popDynNewFile$seed <- seed
    df.popDyn <-  rbind(df.popDyn, popDynNewFile)
  }
  count <- count + 1
  print(count)
}
write.table(df.popDyn, paste0(DATAOUT,"FullSet_popDyn.txt"), row.names = FALSE)

df.popDyn_stats <- left_join(df.popDyn, df.params, by ="seed")

popDyn_plotting_df <- aggregate(cbind(localAdaptSA, invVA, outVA)~sim_gen +muProp + muInv + sigmaK + alpha + rep + mig1 + enVar, FUN = mean, data = df.popDyn_stats)
write.table(popDyn_plotting_df, paste0(DATAOUT, "FullSet_popDynam.txt"))
# 
# for(i in 2:(ncol(popDyn_plotting_df)-3)){
#   popDyn_plotting_df[,i] <- as.factor(popDyn_plotting_df[,i])
# }
# popDyn_plotting_df$sigmaK <- factor(popDyn_plotting_df$sigmaK, c(3, 1.5, 0.75))
# popDyn_plotting_df$sigmaK <- recode_factor(popDyn_plotting_df$sigmaK, '3' = 'Weak Selection', '1.5' = 'Moderate Selection','0.75' = 'Strong Selection')
# popDyn_plotting_df$muInv <- recode_factor(popDyn_plotting_df$muInv, '0.001' = 'Inversions', '0' = 'No Inversions')
# 
# 
# 
# ggplot(data = popDyn_plotting_df, 
#        aes(x = sim_gen, y = localAdaptSA, group = interaction(mig1, rep, muInv, sigmaK))) + 
#   geom_line(aes(color = mig1, linetype = rep)) + 
#   scale_color_manual(values = viridis(6)) + 
#  # facet_grid(muInv~sigmaK) + 
#   facet_wrap(~muInv + sigmaK) +
#   labs(x = "Generation", y = "Local Adaptation") + 
#   guides(color = guide_legend(title = "Migration Rate"),
#          linetype = guide_legend(title = "Replicate")) +
#   theme_classic() +
#   theme(panel.background = element_blank(), 
#         strip.background = element_rect(colour = "white", fill = "grey92")) +
#   ylim(-0.05, 0.6)
# 
# ggplot(data = popDyn_plotting_df[popDyn_plotting_df$muProp ==0.01 & 
#                                    popDyn_plotting_df$alpha ==0.2 &
#                                    popDyn_plotting_df$enVar == 0.1,], 
#        aes(x = sim_gen, y = localAdaptSA, group = interaction(mig1, rep, muInv, sigmaK))) + 
#   geom_line(aes(color = mig1, linetype = rep)) +
#   scale_color_manual(values = viridis(6)) + 
#   #facet_grid(muInv~sigmaK) + 
#   facet_wrap(~muInv + sigmaK) + 
#   labs(x = "Generation", y = "Local Adaptation") + 
#   guides(color = guide_legend(title = "Migration Rate"),
#          linetype = guide_legend(title = "Replicate")) +
#   theme_classic() +
#   theme(panel.background = element_blank(), 
#         strip.background = element_rect(colour = "white", fill = "grey92")) +
#   ylim(-0.05, 1)
# 
#   



