## 20190213 KEL
- talked with Kevin about simparams.txt and which parameters will be input into Slim run
- also talked with Bodie about this file and how he will need to use mutation rate and recombination rate when he overlays the neutral loci

## 20190213 KSF
- developed first (working) version of callSlimWithParams.R, script to call SliM w/ different params using the command line. Pushed to the "cli_params" branch
- To do on Rscript: implement paralellization, add additional parameters (coordinate with Sara)
	- current parameters: mu, Ne, m, seed, r
        - to do: nqtl, alpha, envi_var, sel_strenght/sigma_k
