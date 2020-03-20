 ## 20200302

 Ran test simulations to determine whether dynamics were working as expected. Ran two four migration rates 0.001, 0.01, 0.49, and 0.5 with and without inversions. 

 Finding: under moderate levels of migation inversions allow for local adaptation to occur. 

 Other parameters:
 
  sigmaK = 0.9  ## strong stabilizing selection 
  alpha = 0.002 ## low QTN effect size
  enVar = 0.01  ## low environmental variance added to phenotype
  mu = 1e-7     ## low qtn mutation rate
  muInv = 1e-3  ## high inversion mutation rate
