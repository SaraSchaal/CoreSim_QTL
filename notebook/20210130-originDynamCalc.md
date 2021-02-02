# Calculate origin dyanmics

To understand how inversions arose we need another calculation in the simulation. This is because we don't have any way currently of recording what the genetic background was in the individual at the time the inversion mutated.

1) It’d be great if a mutation could have two tag values because then we could store the number of QTNs in inversion at the time of mutation as the tag value.
This would be the easiest fix but I don’t think this is possible and we need the tag value for the length of the inversions.
2) Instead we can create a new output file that is a dataframe with all inversion mutations and each line represents a qtn in that inversion.
So the number of times an inversion appears in that data file will be tied to how many qtns it had. We will also output the selection coefficient for each QTN and 
the dominance value.
3) VCF files every generation would work too but that would be too memory intensive

We are going with the second option. 

Message to Katie about Inversion QTNs:
Hey Katie I have a quick question about the QTNs and the window description that we’ve been talking about. 
So the way I identify QTNs inside inversions is by subsetting for all genomes with the inversion in it and getting the IDs of the QTNs that are actually 
found in inversions. The frequencies, like we’ve talked about before, will include any individuals (with or without the inversion) that have the QTN in the 
genome because it takes the mutation ID of that QTN and calculates its frequency in the whole population no matter the genetic background. 
But what we aren’t getting with how its currently coded are QTNs that are within the “window” that aren’t in any inverted individuals. 
So this is your scenario 1 in your drawing from 20210127 that you did. None of those QTNs are considered inversion QTNs at the moment.
I am now not sure if that is what we want. I think it still is what we want because those QTNs although in the “window” never contribute to the dynamics of 
the inversion, but I just wanted to get your thoughts. I can 1) leave as is or 2) add any QTN thats in the window EVEN if its not in any inverted individuals.

