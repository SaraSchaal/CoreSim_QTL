ind.pheno <- read.table("results/Inversion/20200827_testForLA/1618803499285_outputIndPheno.txt", header = TRUE)

pops <- 2
n <- 30
inds.sub <- NULL
sampleLA <- NULL
for (i in 1:pops){
  subpop <- ind.pheno[ind.pheno$subpop==i,]
  inds <- sample(subpop$id, size = n,
                replace=FALSE, prob = subpop$fitness)
  sampleLA <- rbind(sampleLA, subpop[inds,])
}



for(i in 1:nrow(sampleLA)){
  if(sampleLA$subpop[i]==1){
    sampleLA$opPopFitness[i] <- dnorm(sampleLA$phenotype[i], -1, sd = 0.5)/dnorm(0,0,sd=0.5)
  } else {
    sampleLA$opPopFitness[i] <- dnorm(sampleLA$phenotype[i], 1, sd = 0.5)/dnorm(0,0,sd=0.5)
  }
}

mean(sampleLA$fitness)-mean(sampleLA$opPopFitness)

