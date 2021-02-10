rerun.seeds <- read.table("./src/seedsToRerun.txt", header = TRUE)
df.params <- read.table("./src/InvSimParams.txt", header = TRUE)
colnames(rerun.seeds) <- "seed"
rerun.sims <- merge(rerun.seeds, df.params, by = "seed", all.x = TRUE)

write.table(rerun.sims, "20201214_RerunSims.txt")
