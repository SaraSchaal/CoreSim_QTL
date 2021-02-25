## analyze inside and outside inversion FST
# load libraries
library(ggplot2)
library(ggpubr)

# set folder to work from
folder <- "results/Inversion/20210220_inOutInvFST/"
samp <- 3384725

# Read in data
## No inversion window
### mutation files
df.noSelMuts <- read.table(paste0(folder, samp, "noSel_outputMutations.txt"), header = TRUE, 
                           stringsAsFactors = TRUE, fill = TRUE)
df.SelMuts <- read.table(paste0(folder, samp, "_outputMutations.txt"), header = TRUE, 
                         stringsAsFactors = TRUE, fill = TRUE)

### inversion qtn files
df.noSelInvQTN <- read.table(paste0(folder, samp, "noSel_outputInvQtn.txt"), header = TRUE, 
                             stringsAsFactors = TRUE)
df.SelInvQTN <- read.table(paste0(folder, samp, "_outputInvQtn.txt"), header = TRUE, 
                           stringsAsFactors = TRUE)
## Inversion Window
### mutation files
df.noSeliWMuts <- read.table(paste0(folder, samp, "noSelinvWind_outputMutations.txt"), header = TRUE, 
                           stringsAsFactors = TRUE, fill = TRUE)
df.SeliWMuts <- read.table(paste0(folder, samp, "invWind_outputMutations.txt"), header = TRUE, 
                         stringsAsFactors = TRUE, fill = TRUE)

### inversion qtn files
df.noSeliWInvQTN <- read.table(paste0(folder, samp, "noSelinvWind_outputInvQtn.txt"), header = TRUE, 
                             stringsAsFactors = TRUE)
df.SeliWInvQTN <- read.table(paste0(folder, samp, "invWind_outputInvQtn.txt"), header = TRUE, 
                           stringsAsFactors = TRUE)

# Select only those in the final generation
df.noSelInvQTNfinalGen <- df.noSelInvQTN[df.noSelInvQTN$sim_gen == 50000,]
df.SelInvQTNfinalGen <- df.SelInvQTN[df.SelInvQTN$sim_gen == 50000,]
df.noSeliWInvQTNfinalGen <- df.noSeliWInvQTN[df.noSeliWInvQTN$sim_gen == 50000,]
df.SeliWInvQTNfinalGen <- df.SeliWInvQTN[df.SeliWInvQTN$sim_gen == 50000,]

# Select only QTN mutations
df.noSelMutsQTNs <- df.noSelMuts[df.noSelMuts$type == "m2",]
df.SelMutsQTNs <- df.SelMuts[df.SelMuts$type == "m2",]
df.noSeliWMutsQTNs <- df.noSeliWMuts[df.noSeliWMuts$type == "m2",]
df.SeliWMutsQTNs <- df.SeliWMuts[df.SeliWMuts$type == "m2",]

# Select only neut mutations
df.noSelMutsNeut <- df.noSelMuts[df.noSelMuts$type == "m1",]
df.SelMutsNeut <- df.SelMuts[df.SelMuts$type == "m1",]
df.noSeliWMutsNeut <- df.noSeliWMuts[df.noSeliWMuts$type == "m1",]
df.SeliWMutsNeut <- df.SeliWMuts[df.SeliWMuts$type == "m1",]
df.noSelMutsNeut$inOut <- "neut"
df.SelMutsNeut$inOut <- "neut"
df.noSeliWMutsNeut$inOut <- "neut"
df.SeliWMutsNeut$inOut <- "neut"
colnames(df.noSelMutsNeut)[10] <- "adapPhenoDiv"
colnames(df.SelMutsNeut)[10] <- "adapPhenoDiv"

# Identify inside vs. outside inversions
df.noSelMutsQTNs$inOut <- NULL
for(i in 1:nrow(df.noSelMutsQTNs)){
  if(length(intersect(df.noSelInvQTNfinalGen$qtn_id, df.noSelMutsQTNs$mutID[i])) > 0){
    df.noSelMutsQTNs$inOut[i] <- "in"
  } else {
    df.noSelMutsQTNs$inOut[i] <- "out"
  }
}

for(i in 1:nrow(df.SelMutsQTNs)){
  if(length(intersect(df.SelInvQTNfinalGen$qtn_id, df.SelMutsQTNs$mutID[i])) > 0){
    df.SelMutsQTNs$inOut[i] <- "in"
  } else {
    df.SelMutsQTNs$inOut[i] <- "out"
  }
}

for(i in 1:nrow(df.noSeliWMutsQTNs)){
  if(length(intersect(df.noSeliWInvQTNfinalGen$qtn_id, df.noSeliWMutsQTNs$mutID[i])) > 0){
    df.noSeliWMutsQTNs$inOut[i] <- "in"
  } else {
    df.noSeliWMutsQTNs$inOut[i] <- "out"
  }
}

for(i in 1:nrow(df.SeliWMutsQTNs)){
  if(length(intersect(df.SeliWInvQTNfinalGen$qtn_id, df.SeliWMutsQTNs$mutID[i])) > 0){
    df.SeliWMutsQTNs$inOut[i] <- "in"
  } else {
    df.SeliWMutsQTNs$inOut[i] <- "out"
  }
}

# Set FST column to numeric
df.noSelMutsNeut$FST <- as.numeric(as.character(df.noSelMutsNeut$FST))
df.SelMutsNeut$FST <- as.numeric(as.character(df.SelMutsNeut$FST))
df.noSelMutsQTNs$FST <- as.numeric(as.character(df.noSelMutsQTNs$FST))
df.SelMutsQTNs$FST <- as.numeric(as.character(df.SelMutsQTNs$FST))

df.noSeliWMutsNeut$FST <- as.numeric(as.character(df.noSeliWMutsNeut$FST))
df.SeliWMutsNeut$FST <- as.numeric(as.character(df.SeliWMutsNeut$FST))
df.noSeliWMutsQTNs$FST <- as.numeric(as.character(df.noSeliWMutsQTNs$FST))
df.SeliWMutsQTNs$FST <- as.numeric(as.character(df.SeliWMutsQTNs$FST))


# merge files
df.noSelMutsQTNsM <- rbind(df.noSelMutsQTNs, df.noSelMutsNeut)
df.SelMutsQTNsM <- rbind(df.SelMutsQTNs, df.SelMutsNeut)
df.SeliWQTNsM <- rbind(df.SeliWMutsQTNs, df.SeliWMutsNeut)
df.noSeliWQTNsM <- rbind(df.noSeliWMutsQTNs, df.noSeliWMutsNeut)

#noSelQTNs <- ggplot(df.noSelMutsQTNs, aes(x = FST, color = inOut)) +
                 # geom_histogram(alpha = 0.9, bins = 100) 
                  #geom_histogram(df.noSelMutsNeut, aes(x = FST))

#SelQTNs <- ggplot(df.SelMutsQTNs, aes(x = FST, color = inOut)) +
                  #geom_histogram(fill = "white", alpha = 0.5, bins = 100) 
                  #geom_histogram(df.SelMutsNeut, aes(x = FST))


noSelFSTpos <- ggplot(df.noSelMutsQTNsM, aes(x = position, y = FST, group = inOut)) + 
                 geom_point(aes(color = inOut)) + 
                 scale_color_manual(values=c("red", "goldenrod", "navy")) +
                 xlim(0, 2100000) +
                 ylim(0, 0.2) +
                 theme(legend.position = "none")

noSelFSTposNoylim <- ggplot(df.noSelMutsQTNsM, aes(x = position, y = FST, group = inOut)) + 
                        geom_point(aes(color = inOut)) + 
                        scale_color_manual(values=c("red", "goldenrod", "navy")) +
                        xlim(0, 2100000) 
              
                 
SeliWFSTpos <- ggplot(df.SeliWMutsQTNsM, aes(x = position, y = FST, group = inOut)) + 
                    geom_point(aes(color = inOut)) + 
                    scale_color_manual(values=c( "red", "goldenrod", "navy")) +
                    xlim(0, 2100000) + 
                    ylim(0, 0.2) +
                    theme(legend.position = "none")
                   # geom_vline()

noSeliWFSTpos <- ggplot(df.noSeliWMutsQTNsM, aes(x = position, y = FST, group = inOut)) + 
                  geom_point(aes(color = inOut)) + 
                  scale_color_manual(values=c("red", "goldenrod", "navy")) +
                  xlim(0, 2100000) +
                  ylim(0, 0.2) +
                  theme(legend.position = "none")

noSeliWFSTposNoylim <- ggplot(df.noSeliWMutsQTNsM, aes(x = position, y = FST, group = inOut)) + 
                        geom_point(aes(color = inOut)) + 
                        scale_color_manual(values=c("red", "goldenrod", "navy")) +
                        xlim(0, 2100000) 


SeliWFSTpos <- ggplot(df.SeliWMutsQTNsM, aes(x = position, y = FST, group = inOut)) + 
                      geom_point(aes(color = inOut)) + 
                      scale_color_manual(values=c( "red", "goldenrod", "navy")) +
                      xlim(0, 2100000) + 
                      ylim(0, 0.2) +
                      theme(legend.position = "none")
                      # geom_vline()
                    
              
ggarrange(SelFSTpos, noSelFSTpos, noSelFSTposNoylim, nrow = 1, ncol = 3)
ggarrange(SeliWFSTpos, noSeliWFSTpos, noSeliWFSTposNoylim, nrow = 1, ncol = 3)                  


highFST <- df.SelMutsQTNs[df.SelMutsQTNs$FST > 0.15, ]         
highFST <- highFST[!is.na(highFST$position),]  


# arrows(x0=df_sub$inv_pos[i], x1=df_sub$inv_end[i],
#        y0=-0.01, y1=-0.01, angle=90, length=0.05,
#        code=3, 
#        col=rgb(cr(df_sub$freq[i]), max=255))