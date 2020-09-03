ind.pheno <- read.table("results/PrelimSims/3383282_outputIndPheno.txt", header = TRUE)

pops <- 2
n <- 30
inds.sub <- NULL
for (i in 1:pops){
  inds <- sample(ind.pheno$id[ind.pheno$subpop==i], size = n,
                replace=FALSE, prob = ind.pheno$fitness[ind.pheno$subpop==i])
  inds.sub <- c(inds.sub, inds)
}


