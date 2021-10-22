
Example of what we want slim to do when we have an inversion at [X,Y]. In this example, D is at Y + 1 and A is at X - 1 in the standard orientation and C is some locus in the middle of the inversion just as an example locus to follow:
	
	AX--C----YD standard

	AY----C--XD inverted


   Homozygote Inverted (* is proposed breakpoint):

	   *  
	AY----C--XD - chrom1

	..    .  ..
	AY----C--XD - chrom2


   After recombination:
```
      	      .  ..
	AY----C--XD

	..    
	AY----C--XD

```


Example of what SLiM ACTUALLY does with this proposed recombination break regardless of whether an individual is homozygous inverted or not beceause SLiM only recognizes the standard orientation in the genetic map:

   Homozygote whether inverted or not (* is proposed breakpoint):

	   *  
	AX--C----XD - chrom1

	..  .    ..
	AX--C----YD - chrom2


   After recombination:
   
       	    .    ..
	AX--C----YD

	..    
	AX--C----YD


So we manually have to make a recombination rate with the standard orientation that results in the same inheritance as what we showed in the "what we want slim to do" example.


Example of how to make the recombination for an individual homozygous for the inversion to match the first example by adding breakpoints at X and Y + 1 and again * represent proposed breakpoints. We want AY to be inherited together and X and D to be inherited together:

	 * *      *
	AX--C----YD - chrom1

	..  .    ..
	AX--C----YD - chrom2


	After recombination:

	 .        .
	AX--C----YD

        .   .    .
	AX--C----YD

																	   
Here by adding those new breakpoints at X and Y + 1, we successfully get AY and the dotted AY to be inherited together and the dotted XD and XD to be inherited together. What doesn't match is the haplotypes inside the inversion. You can see that C is now inherited with dotted XD and not XD. Is this a problem? It could change the number of segregating haplotypes. 


Now I checked the code that was proposed by Vince, Peter and Andrew and did simple checks in Eidos to see that it did this properly and I do not believe it does. In the pictured example, I have a breakpoints vector with one breakpoint that falls inside the inverted region and therefore sum(inInv) will be odd. There are no breakpoints at the ends of the inversion (that needs to be checked next). What happens in their currect code is that they check the breakpoint vector for any that are either at the start or the end + 1 (where we need to add breakpoints) and stores a logical vector that is the length of break points and stores them in the left or right variable. This will only have a T if one of those proposed breakpoints is either at the start (T at the position in the left variable) or end + 1 (T at the position in the right variable). In my example they are all F because there are no breakpoints at the ends and what SHOULD happen is that they are then added to the vector. BUT this doesn't happen because they add the start and end+1 locations in this way ```c(inv_start, inv_end+1)[c(sum(left)> 0, sum(right)>0)]``` and they will only be added then if there is a T in those logical vectors. I THINK it should be that only if they are all F they should be added.

<img src="../src/Inv_Issue/Current_Inv_Code.jpeg" width = "500">

In my proposed change that last segment of code changes to ```c(inv_start, inv_end+1)[c(sum(left)==0, sum(right)==0)]``` which adds the breakpoints at the start and end + 1 of the inversion if they aren't already present. 

<img src="../src/Inv_Issue/Proposed_Inv_Correction.jpeg" width = "500">
