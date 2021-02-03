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




