### Example command line script ###

nohup slim -d muBase=1e-05 -d N=1000 -d mig1=0.01 -d mig2=0.01 -d sigmaK=0.9 -d alpha=0.002 -d seed=3383282 -d envVar=0.01 -d path="/scratch/schaal.s/InversionSimulations/results" -d muInv=1e-3 -d dom=F -m -t simfiles/InversionModel.slim &> nohup3383282.out&