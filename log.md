## 20190213 KEL
- talked with Kevin about simparams.txt and which parameters will be input into Slim run
- also talked with Bodie about this file and how he will need to use mutation rate and recombination rate when he overlays the neutral loci

## 20190213 KSF
- developed first (working) version of callSlimWithParams.R, script to call SliM w/ different params using the command line. Pushed to the "cli_params" branch
- To do on Rscript: implement paralellization, add additional parameters (coordinate with Sara)
	- current parameters: mu, Ne, m, seed, r
        - to do: nqtl, alpha, envi_var, sel_strenght/sigma_k

## 20190213 SMS'
- finished adjusting simulation to calculate everything properly and included constants for all the parameters that we want to be able to call from the command line
- sent Kevin a tree seq output file to play around with for the machine learning
- simulation is running slow for 2000 QTLs have it running on my computer tonight to see how long it actually takes
- question: I noticed one of your comments for Bodie is that we need to consider the data either as a single population or as two populations with the same migration rate. Does this mean we can't have varying migration rates between the populations?

## 20190220 KSF

- merged cli_params branch with master on github
- added parallelization to Rscript callSlimWithParams.R
- Worked with Bodie to update SLiM on comp5 from v 2.6 -> 3.2.1
- to do: test script on small subset of simulations

## 20190228 KSF

- Tested script on cluster - took 8.5hrs to run a batch of sims with 2000 qtls in parallel
