# Project Checklist

* Develop a list of simulation parameters - *Katie*
    * have 200K simulations that we want to run

* Develop a script to call a simulation from the command line with different parameters - *Kevin*

* Develop pseudoreplication script for cluster - *Kevin*

* SLiM simulation - *Sara*
    * need to coordinate with Kevin to take in parameters:
        * nqtl
        * alpha
        * Ne
        * mu
        * r
        * envi_var
        * sel.strength ($\sigma_K$ ?)
        * m
        * seed
    * need to coordinate with Bodie on output file names for tree processing
        * seed should be first 6 digits of the name of the `.trees` file
    * want to make sure we are outputting all relevant info that would allow us to process the files, like mean phenotype and mean fitness of each deme

* trees processing to vcf file - *Bodie*
    * we need to decide how we want to recapitate - assume 1 populations or two population with same migration? If latter, check with Aki about recapitation process for multiple populations
    * will use `grep` to get information on m and 

## Run on cluster
* decide on how results will be organized
* write a bash script to run simulation and processes trees output with pseudo-parallelization
* run pipeline on computer first
* run 1 replicate of pipeline on one processor of cluster second. Make sure all paths are set up correctly.
* run a subset of replicates with bash script on cluster to make sure psueodoparallelization is working properly
    * at this point we want to meet and make sure we have all the outputs we want and they are well organized
* let 'em go!

## For machine learning project
* Run test statistics (many pipelines set up from recombination paper)
* Run machine learning algorithm on test statistics



