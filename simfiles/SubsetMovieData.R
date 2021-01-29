### Movie Data ###
#install.packages("dplyr")
library("dplyr")

all.dataLA <- read.table("results/Inversion/20201115_FullSet/FullSet_LAallData.txt", header = TRUE, stringsAsFactors = FALSE)
head(all.dataLA)
levels(as.factor(all.dataLA$mu_base))
head(subset(all.dataLA, subset = mig == 0.25 & sigmaK == 0.75 & mu_inv == 0.001 & alpha == 0.002 & enVar == 0 & mu_base == 0.0000001))

#### Inversiions through time data ####
movieDataInvTime <- read.table("results/Movie/3384725_outputInvTime.txt", header = TRUE, stringsAsFactors = FALSE)
movieDataInvData <- read.table("results/Movie/3384725_outputInvSumInfo.txt", header = TRUE, stringsAsFactors = FALSE)
movieDataInvTime2 <- read.table("results/Movie/Rerun/3384725_outputInvTime.txt", header = TRUE, stringsAsFactors = FALSE)
movieDataInvData2 <- read.table("results/Movie/Rerun/3384725_outputInvSumInfo.txt", header = TRUE, stringsAsFactors = FALSE)

head(movieDataInvData)
head(movieDataInvTime)

movie.all.data <- merge(movieDataInvTime, movieDataInvData, all.x = TRUE, by = "inv_id")
movie.all.data2 <- merge(movieDataInvTime2, movieDataInvData2, all.x = TRUE, by = "inv_id")
movie.all.data <- movie.all.data[,-12]
colnames(movie.all.data)[2] <- "sim_gen"
head(movie.all.data)
str(movieDataInvData)
str(movieDataInvTime)

write.table(movie.all.data, "results/Movie/movieData.txt", row.names = FALSE)

#### QTNs through time data ####
movieDataQTNTime <- read.table("results/Movie/3384725_outputInvQtn.txt", header = TRUE, stringsAsFactors = FALSE)
movieDataQTNData <- read.table("results/Movie/3384725_outputInvQtnSumInfo.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(movieDataQTNData)[3] <- "qtn_id"
movie.all.QTN.data <- merge(movieDataQTNTime, movieDataQTNData, all.x = TRUE, by = "qtn_id")
movie.all.QTN.data <- movie.all.QTN.data[,-c(8:9)]
colnames(movie.all.QTN.data)[2:3] <- c("sim_gen", "inv_id")


### QTNs inside vs. outside inversions ###

df.final.mutations <- read.table("results/Movie/3384725_outputMutations.txt", header = TRUE, stringsAsFactors = FALSE)
head(df.final.mutations)
head(movieDataQTNData)
head(movieDataInvData)
colnames(df.final.mutations)


df.QTNInvSummaryData <- merge(movieDataQTNData, movieDataInvData, by = "inv_id")
df.metaData <- merge(df.final.mutations, df.QTNInvSummaryData, by = "qtn_id")

# Inside Outside QTNs
for(i in 1:nrow(df.metaData)){
  if(is.na(inv_pos)){
    df.metaData$inOut <- "OUT"
  } else {
    df.metaData$inOut <- "IN"
  }
}

# Origin

df.noOutside <- drop_na(df.metaData)

for(i in 1:nrow(df.noOutside)){
  if()
}

