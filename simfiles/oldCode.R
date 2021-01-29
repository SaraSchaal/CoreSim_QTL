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
