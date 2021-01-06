### QTN summary data file missing data ###

For both inversion mutations and QTN mutations we output two different files. One with the frequencies through time and the other with the summary data about the individual mutations so it doesn't print out every generation that each specific mutation is present in the population. The second file with the summary data for the qtns doesn't include all the data.

*Solution*: found that the if statment was missing a set of parantheses that caused the if statement to be incorrect. 

	Was: if(sim.generation - InvQTNsGen[i] <= 200)

	Now: if((sim.generation - InvQTNsGen[i]) <= 200)