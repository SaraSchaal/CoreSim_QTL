library(parallel)
#########################################################################################
# File    : callSlimWithParams.R
# History : 2019-2-13  Created by Kevin Freeman (KF)
# --------------------------------------------------------------------------------------
# Usage: Rscript callSlimWithParams.R [output directory] [params file] [SliM Script]
# 
# Command line options must be in order. If called without inputs, it will default to 
# paths that work on my local environment. 
# ---------------------------------------------------------------------------------------
# Description:
# callSlimWithParams.R takes a file where each row is a set of parameters for 1 
# replication of the SliM simulation. It constructs calls to SliM for each row
# and runs the given script with the given params in parallel.
#
########################################################################################
burnin <- 4000   # use standard burnin

### Parse arguments or set to default
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  path       <- "/media/kevin/TOSHIBA_EXT/CoreSim_QTL/results"
  param_path <- "simparams.txt"
  script     <- "../simfiles/IndependentQTNsModel.slim"
} else if (length(args) == 3){
  path       <- args[1]
  path       <- normalizePath(path)     # strip trailing / and expand relative path
  param_path <- args[2]
  script     <- args[3]
} else {
  stop(paste("Wrong number of arguments. Script takes 3: (1) the desired output directory, (2) the params file, and \
             (3) the SliM script. You provided:", 
             length(args)), call.=TRUE)
}
usage <- "Rscript callSlimWithParams.R [output directory] [params file] [SliM Script]"

### some input checking
if (!dir.exists(path)){
  stop(paste("Specified output path is not a directory or does not exist. Make sure you provided the parameters in the right order\n" , usage), 
       call.=TRUE)
}
if (!file.exists(param_path)){
  stop(paste("Specified sim param file is not a valid file or does not exist. Make sure you provided the parameters in the right order\n" , usage), 
       call.=TRUE)
}
if (!file.exists(script)){
  stop(paste("Specified SLiM script is not a valid file or does not exist. Make sure you provided the parameters in the right order\n" , usage), 
       call.=TRUE)
}
## todo: implement try-catch
sim_params <- read.csv(param_path, header = TRUE, sep = " ")

#-------------------------------------------------------------------------
# callSlim(param_row, path, script)
# ------------------------------------------------------------------------
# Inputs  (3) : row of parameters, path for output files, SliM script
# Outputs (0) : none
#------------------------------------------------------------------------
# callSlim() does the bulk of the work in this script. It takes a row
# of parameters from the data frame of all parameters, a path for the output
# file, the name of the SliM script and constructs a slim command to run 
# the script using the given parameters
#-------------------------------------------------------------------------
callSlim <- function(param_row, path, script){
  # extract params from the given row
  seed    <- param_row$seed                # seed
  mu      <- param_row$mu                  # mutation rate
  Ne      <- param_row$Ne                  # sample size
  mig     <- param_row$m                   # migration rate
  alpha   <- param_row$alpha_levels        # effect size
  qtls    <- param_row$nqtl_levels         # number of qtl 
  envVar  <- param_row$envi_var            # environmental variance (0,1,2)
  sigma_k <- param_row$sel.strength        # strength of selection
  r       <- param_row$r                   # recombination rate
  
  ## construct the command and call it with system()
  cmd <- paste0("slim -l -d mu=", mu, " -d Ne=", Ne, " -d mig1=", mig," -d mig2=", mig,
                " -d sigma_K=", sigma_k, " -d alpha=", alpha, " -d burnin=", burnin,
                " -d my_seed=", seed, " -d \"path='", path, "'\" -d r=", r, " ", script)
  print(cmd)
  system(cmd)
}

callSlim(param_row = sim_params[1,], path = path, script = script)
