## 20201104 - first parameter set run	

Ran 2160 simulations with the following parameter combinations:

	mu_base <- c(1e-07, 1e-06)
 	mu_inv <- c(0, 1e-03, 1e-06)
 	sigmaK <- c(0.75, 1.5, 3) 
 	alpha <- c(0.002, 0.2)
  	reps <- c(1:5)
  	N <- 1000 
  	envar <- c(0,0.1) 
  	mig1 <- c(0.001, 0.01, 0.1, 0.25, 0.4, 0.5)
  	mig2 <- c(0.001, 0.01, 0.1, 0.25, 0.4, 0.5)
  	chrom_num <- 21
  	recom <- 1e-6
  	burnin <- 10000
  	theta1 <- 1
  	theta2 <- -1
  	dom <- 2 	#2 is false so no dominance 

Each simulation runs between ~1.5-4 hrs 
Peak memory ~1.4 GB
