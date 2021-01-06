### Code for analyzing inversion characteristics ###
options(scipen = 999)
## Download Data
folder <- "./results/Inversion/20201115_FullSet/"
df.invTime <- read.table(paste(folder, "FullSet_invTime.txt", sep = ""), header = TRUE)
df.invData <- read.table(paste(folder, "FullSet_invData.txt", sep = ""), header = TRUE)
df.params <- read.table(paste(folder, "invSimParams.txt", sep = ""), header = TRUE)
dim(df.params)
head(df.invTime)
head(df.invData)

#final.gen.invs <- subset(df.invData, subset = sim_gen == 50000)
final.gen.invsTime <- subset(df.invTime, subset = sim_gen == 50000)
df.allInv <- merge(df.invData, final.gen.invsTime, all.y = TRUE, by = c("seed", "inv_id"))
df.allInvData <- merge(df.allInvData, df.params, all.x = TRUE, by = "seed")
df.allInvData <- df.allInvData[, -which(names(df.allInvData) %in% "sim_gen.x")]

## Plotting
ggplot()

# merge inversion through time data with parameters
df.invTimeParam <- merge(df.invTime, df.params, all.x = TRUE, by = c("seed"))

# subset inversion qtn file for just the seeds of interest
df.invQTNTimeSub <- subset(df.invQTNTime, subset = seed == c(3383645, 3383657, 3383669, 3383681, 3383693))
colnames(df.invQTNData)[3] <- "qtn_id"

# merge inversion QTN data with parameters
movie.InvQTNSub.Data <- merge(df.invQTNTimeSub, df.params, all.x = TRUE, by = c("seed"))
head(movie.InvQTNSub.Data)

df.invQTNData[df.invQTNData$]
movie.InvQTN.Data <- merge(df.invQTNData, movie.InvQTNSub.Data, all.y = TRUE, by = c("seed", "qtn_id"))

