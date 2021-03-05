### Code for analyzing inversion characteristics ###
# Settings
options(scipen = 999)
# Libraries
library(dplyr)
library(ggplot2)

## Download Data
folder <- "./results/Origin/"
df.InvOrigin <- read.table(paste(folder, "3384725_outputOriginInvQtn.txt", sep = ""), header = TRUE)
df.invData <- read.table(paste(folder, "3384725_outputInvSumInfo.txt", sep = ""), header = TRUE)
df.invTime <- read.table(paste(folder, "3384725_outputInvTime.txt", sep = ""), header = TRUE)
#df.params <- read.table(paste(folder, "invSimParams.txt", sep = ""), header = TRUE)
#dim(df.params)

# Subset for final generation
df.invFinalGen <- subset(df.invTime, subset = sim_gen == 50000)


# get sum of effect size on the phenotype
df.InvQTNsum <- aggregate(qtnSelCoef~inv_id, FUN = sum, data = df.InvOrigin)
df.InvQTNnum <- aggregate(qtn_id~inv_id, FUN = length, data = df.InvOrigin)
df.InvOriginGen <- aggregate(sim_gen~inv_id, FUN = mean, data = df.InvOrigin)
df.InvQTNsumNum <- left_join(df.InvQTNsum, df.InvQTNnum, by = "inv_id")
colnames(df.InvQTNsumNum)[3] <- "origin_num_qtns"
df.InvQTNsumNumGen <- left_join(df.InvQTNsumNum, df.InvOriginGen, by = "inv_id")
colnames(df.InvQTNsumNumGen)[4] <- "origin_gen"

# merge inversion data with origin dynamics
df.InvDataOrigin <- left_join(df.invFinalGen, df.InvQTNsumNumGen, by = "inv_id")
#check
dim(df.InvOrigin[df.InvOrigin$inv_id == 40748419,])


# subset for MAF > 0.01
df.InvDataOriginMAF <- df.InvDataOrigin[df.InvDataOrigin$freq > 0.01,]
for(i in 1:nrow(df.InvDataOriginMAF)){
  if(df.InvDataOriginMAF$freq_p1[i] > df.InvDataOriginMAF$freq_p2[i]){
    df.InvDataOriginMAF$pop[i] <- "pop1"
  } else {
    df.InvDataOriginMAF$pop[i] <- "pop2"
    
  }
}



## Subset dataframe to get how the inversions change through time
inv.IDs <- as.vector(df.InvDataOriginMAF$inv_id)
df.invFinalAllData <- df.invTime[df.invTime$inv_id %in% inv.IDs, ]
df.invFinalAllData$qtnSelCoefsum <- df.invFinalAllData$mean_qtnSelCoef*df.invFinalAllData$num_qtns
df.invFinalAllDataPop <- left_join(df.invFinalAllData, df.InvDataOriginMAF[c(2,15)], by = "inv_id")

## Plotting
par(mfrow = c(1,2))
plot(origin_num_qtns~origin_gen, data = df.InvDataOriginMAF, xlab = "Inversion Origin Generation",
     ylab = "Origin number of QTNs in Inversion", col = "cornflowerblue", pch = 19, ylim = c(0, 400))
plot(num_qtns~origin_gen, data = df.InvDataOriginMAF, xlab = "Inversion Origin Generation",
     ylab = "Final number of QTNs in Inversion", col = "dodgerblue4", pch = 19, ylim = c(0, 400))

ggplot(df.InvDataOriginMAF, aes(x = origin_gen, y = qtnSelCoef, group = pop)) + 
  geom_point(aes(color = pop, size = inv_FST), alpha = 0.8) + 
  scale_color_manual(values=c("navy", "red")) + 
  theme_classic() +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92"),
        text = element_text(size = 15)) +
  labs(title = "Total Effect of Inversion QTNs on Phenotype",
       y = "sum of each Inversion QTNs effects on phenotype",
       x = "Inversion Origin Generation") +
  guides(color = guide_legend(title = "Pop with Highest\nFrequency of Inv")) +
  guides(size = guide_legend(title = "Inversion FST")) 


plot(qtnSelCoef~origin_gen, data = df.InvDataOriginMAF, xlab = "Inversion Origin Generation",
     ylab = "Sum of QTN effects on phenotype in Inversion", col = "cyan4", pch = 19)


ggplot(df.invFinalAllDataPop, aes(x = sim_gen, y = qtnSelCoefsum, group = pop)) + 
  geom_point(aes(color = pop, size = inv_FST), alpha = 0.8) + 
  geom_line(aes(color = pop), alpha = 0.8) + 
  scale_color_manual(values=c("navy", "red")) + 
  scale_size(range = c(0.5, 4), breaks = c(0.00001, 0.05, 0.15, 0.2)) + 
  theme_classic() +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92"),
        text = element_text(size = 15)) +
  labs(title = "Total Effect of Inversion QTNs on Phenotype",
       y = "sum of each Inversion QTNs effects on phenotype",
       x = "Generation") +
  guides(color = guide_legend(title = "Pop with Highest\nFrequency of Inv")) +
  guides(size = guide_legend(title = "Inversion FST")) 

ggplot(df.invFinalAllDataPop, aes(x = sim_gen, y = qtnSelCoefsum, group = pop)) + 
  geom_point(aes(color = pop, size = inv_FST), alpha = 0.8) + 
#  geom_line(aes(color = pop), alpha = 0.8) + 
  scale_color_manual(values=c("navy", "red")) + 
  scale_size(range = c(0.5, 4), breaks = c(0.00001, 0.05, 0.15, 0.2)) + 
  theme_classic() +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92"),
        text = element_text(size = 15)) +
  labs(title = "Total Effect of Inversion QTNs on Phenotype",
       y = "sum of each Inversion QTNs effects on phenotype",
       x = "Generation") +
  guides(color = guide_legend(title = "Pop with Highest\nFrequency of Inv")) +
  guides(size = guide_legend(title = "Inversion FST")) 


example <- df.invFinalAllDataPop[df.invFinalAllDataPop$sim_gen == 50000,][1,2]
df.example <- df.invFinalAllDataPop[df.invFinalAllDataPop$inv_id == example,]

ggplot(df.example, aes(x = sim_gen, y = qtnSelCoefsum)) + 
  geom_point(aes(color = pop, size = inv_FST), alpha = 0.8) + 
  geom_line(aes(color = pop), alpha = 0.8) + 
  scale_color_manual(values=c("navy", "red")) + 
  scale_size(range = c(0.5, 4), breaks = c(0.00001, 0.05, 0.15, 0.2)) + 
  theme_classic() +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92"),
        text = element_text(size = 15)) +
  labs(title = "Total Effect of Inversion QTNs on Phenotype",
       y = "sum of each Inversion QTNs effects on phenotype",
       x = "Generation") +
  guides(color = guide_legend(title = "Pop with Highest\nFrequency of Inv")) +
  guides(size = guide_legend(title = "Inversion FST")) 


### How often does overlap occur in this simulation
inv.data <- left_join(df.invTime, df.invData[c(1, 3:6)], by = "inv_id")

gen200 <- subset(inv.data, subset = inv.data$sim_gen == 200)

par(mfrow = c(1,1))
plot(x =  gen200$inv_pos, y = gen200$inv_id, data = gen200, col = "red", pch = 19) 
points(x = gen200$inv_end, y = gen200$inv_id, data = gen200, col = "blue", pch = 19) 


unique.gens <- unique(inv.data$sim_gen)
overlap <- 0 
df.overlap <- NULL 
for(i in 1:length(unique(inv.data$sim_gen))){
  genData <- subset(inv.data, subset = inv.data$sim_gen == unique.gens[i])
  for(j in 1:nrow(genData)){
    for(k in 1:nrow(genData)){
      if(j != k){
        if(genData$inv_pos[j] < genData$inv_pos[k] & genData$inv_end[j] > genData$inv_end[k]){
          overlap = overlap + 1
          df.overlap <- rbind(df.overlap, c(unique.gens[i], genData$inv_id[j], 
                                            genData$inv_id[k], genData$inv_pos[j], 
                                            genData$inv_end[j], genData$inv_pos[k], 
                                            genData$inv_end[k]))
        }
      }
    }
  }
}
colnames(df.overlap) <- c("sim_gen", "inv_id1", "inv_id2", "inv_id1_start", "inv_id1_end", "inv_id2_start", "inv_id2_end")
df.overlap <- as.data.frame(df.overlap)
head(df.overlap)
nrow(df.overlap)

length(unique(paste(df.overlap$inv_id1, df.overlap$inv_id2)))


df.overlap$IDpaste <- paste(df.overlap$inv_id1, df.overlap$inv_id2)
df.unique <- as.vector(unique(df.overlap$IDpaste))
df.overlap[df.overlap$IDpaste %in% df.unique, ]
