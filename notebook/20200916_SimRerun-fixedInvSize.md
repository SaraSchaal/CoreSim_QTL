## Re ran simulations this week and fixed the inversion size issue. 
Now we have inversions that evolve to be the mean of the size range that we are drawing from in the simulation.

This makes sense but is this what we want? It is potentially just a figment of how we create inversions.

Second update is convertToSubstitution. 
Message from Katie:
I have a hypothesis about the weird inversion dynamics. Could you take the case where large inversions arise prior to selection, and for the QTNs set convertToSubstitution=FALSE (for reasons we already talked about, we need to keep this set to TRUE, but just change it to test the hypothesis).

However, we have always had this as FALSE because we do NOT want SLiM to remove mutations when fixed because they contribute to the phenotypes. I reran all sims
with TRUE for m2 mutations, but this changed results and after I reread section in manual we should have this as FALSE. I should have reread before running with TRUE.
