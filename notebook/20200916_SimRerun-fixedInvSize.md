## Re ran simulations this week and 
### Update 1 - fixed the inversion size issue  
Now we have inversions that evolve to be the mean of the size range that we are drawing from in the simulation.  
![Average Inversion Size Through Time](../figures/20200911/invLengthAverage.pdf)  
This makes sense but is this what we want? It is potentially just a figment of how we create inversions.  
  
### Update 2 - convertToSubstitution
Message from Katie:    
       
         I have a hypothesis about the weird inversion dynamics. Could you take the case where large inversions arise prior to selection, and for the QTNs set  convertToSubstitution=FALSE (for reasons we already talked about, we need to keep this set to TRUE, but just change it to test the hypothesis).
  
However, we have always had this as FALSE because we do NOT want SLiM to remove mutations when fixed because they contribute to the phenotypes. I reran all sims
with TRUE for m2 mutations, but this changed results and after I reread section in manual we should have this as FALSE. I should have reread before running with TRUE.  
![Local Adaptation convertToSubstitution Summary](../figures/20200911/summary_ConvToSub.pdf)

Going to confirm with Katie that this all checks out.
	
#### KEL notes

If the (m3) inversion dynamics depend on whether m2 mutations are kept in the sim when they are fixed, this is still a bug we need to run by Ben.

covertToSubstitution for m2 QTNs
* when it is FALSE this means we keep the QTNs in the genome
* when it is TRUE this means that when QTNs get fixed, they are "removed" from the sim and no longer contribute to the phenotype

We want to keep this as FALSE for m2 mutations, because they contribute to the absolute value of the mean phenotype when they are fixed.
