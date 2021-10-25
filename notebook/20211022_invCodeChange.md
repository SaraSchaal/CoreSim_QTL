## Example of what we want slim to do when we have an inversion at [X,Y]   
In this example, A is at X - 1, B is at X + 1, C is at Y - 1 and D is at Y + 1 in the standard orientation and W is some locus in the middle of the inversion just as an example locus to follow: 
```	

# standard (uninverted):
   X  			   Y
---+-------+---------------+-----
  A B      W              C D

# inverted:

   X  			   Y
---+-------+---------------+-----
  A C      W              B D

```

Homozygote Inverted (* is proposed breakpoint):
```
	   *  
   X  			   		   Y
---+-------+---------------+-----
  A C      W              B D       - chrom1

   X  			   		   Y
---+-------+---------------+-----
  a c      w              b     	- chrom2
	
```

After recombination:
```

   X  			   		   Y
---+-------+---------------+-----
  A C      w              b d       

   X  			   		   Y
---+-------+---------------+-----
  a c      W              B D    	

```

## Example of what SLiM ACTUALLY does with this proposed recombination break regardless of whether an individual is homozygous inverted or not  
SLiM only recognizes the standard orientation in the genetic map. Below is an example of how SLiM handles recombination. 

Homozygote whether inverted or not (* is proposed breakpoint):
```
   		*
   X  			  		   Y
---+-------+---------------+-----
  A B      W              C D

   X  			  		   Y
---+-------+---------------+-----
  a b      w              c d

```
After recombination:
```   
    X  			  		   Y
---+-------+---------------+-----
  A B      w              c d

   X  			  		   Y
---+-------+---------------+-----
  a b      W              C D

```
Now A and B are inherited together and C and D are inherited together which is not what we want. So we manually have to make a recombination event with the standard orientation that results in the same inheritance as what we showed in the "what we want slim to do" example.

## Example of how to make the recombination for an individual homozygous for the inversion to match the first example  
We do this by adding breakpoints at X and Y + 1 and again * represent proposed breakpoints. We want AY to be inherited together and X and D to be inherited together:  
```
   *	*				    *
   X  			  		   Y
---+-------+---------------+-----
  A B      W              C D

   X  			  		   Y
---+-------+---------------+-----
  a b      w              c d

```

After recombination:
```
   		
   X  			  		   Y
---+-------+---------------+-----
  A b      W              C d

   X  			  		   Y
---+-------+---------------+-----
  a B      w              c D

```
																	   
Here by adding those new breakpoints at X and Y + 1, we successfully get A and C to be inherited together then B and D to be inherited together. What doesn't match is the haplotypes inside the inversion. You can see that C is now inherited with dotted XD and not XD. Is this a problem? It could change the number of segregating haplotypes. (UNRESOVELED QUSTION)

  
## Checking the new code to make sure this happens properly
Now I checked the code that was proposed by Vince Buffalo, Peter Ralph and Andrew Kern. I did simple checks in Eidos to see that it did this properly and I do not believe it does. In the pictured example, I have a breakpoints vector with one breakpoint that falls inside the inverted region and, therefore, sum(inInv) will be odd. There are no breakpoints at the ends of the inversion (that needs to be checked next). What happens in their currect code is that they check the breakpoint vector for any that are either at the start or the end + 1 (where we need to add breakpoints) and creates a logical vector that is the length of break points and stores them in the variables left or right. This will only have a T if one of those proposed breakpoints is either at the start (T at the position in the left variable) or end + 1 (T at the position in the right variable). In my example, they are all F because there are no breakpoints at the ends and what SHOULD happen is that they are then added to the vector. BUT this doesn't happen because they add the start and end+1 locations in this way ```c(inv_start, inv_end+1)[c(sum(left)> 0, sum(right)>0)]``` and they will only be added then if there is a T in those logical vectors. I THINK it should be that only if they are all F they should be added.

<img src="../src/Inv_Issue/Current_Inv_Code.png" width = "1000">

I've added a proposed change that I think works to first deal with whether there are no proposed breakpoints at the start or end + 1. This correctly adds those breakpoints in if they aren't in the proposed breakpoints vector. Then if there are is a breakpoint at either the start or end + 1, we need to remove it and add one at the opposite end. In my example, there is one present at the start so we remove it and add one to the end + 1 with this bit of code ```breakpoints = sort(c(breakpoints[!(left | right)],c(inv_start, inv_end + 1)[c(sum(right)>0,sum(left)>0)]))```that is similar to what they had except they had ```c(inv_start, inv_end + 1)[c(sum(left)>0, sum(right)>0)]``` which doesn't add the opposite end, it just re-ads the one that was present already. Finally, if there are proposed breakpoints both at the start and at end+1 then we remove them. Here are three Eidos examples that address all three examples:

No breakpoints at ends but one inside:
<img src="../src/Inv_Issue/Proposed_Inv_Correction_ExampleNoEnds.png" width = "1000">
  


One breakpoint at the start of the inversion:
<img src="../src/Inv_Issue/Proposed_Inv_Correction_ExampleOneEnd.png" width = "1000">

  

Breakpoints proposed at both ends:
<img src="../src/Inv_Issue/Proposed_Inv_Correction_ExampleBothEnds.png" width = "1000">
  
## Katie's Questions

#### some questions I have - how does SliM model breakpoints? Are they stored as a location, and if so, does the breakpoint occur before or after the mutation? 
They are stored as a position in the genome. I'm not sure I know what you mean by mutation. If you are asking, what happens when a proposed breakpoint is exactly on a mutation. That would mean it flips to the opposite strand at that mutation position. This is confirmed by Ben in the thred of the Github issue: "A breakpoint at n means that [0, n-1] and [n, L] are inherited separately". 
Example (proposed breakpoint at the position 4 which is C):

```
	   *	
	A--C--B-D

	a--c--b-d
```	
After recombination:

```
	A--c--b-d

	a--C--B-D
```

#### If we are going to ignore double crossovers, shouldn’t we ignore triple crossovers? 
We aren't ignoring all double crossovers. We are just ignoring all crossovers in individuals heterozygous for an inversion. All crossovers happen in homozygous individuals. If it is an odd number of crossovers (so a single or triple.. although I doubt we get triple crossovers in our inversions with our length and recombination rate), thats when we have to think through all of these scenarios and adjust for having the correct inheritence. If it is even (so double crossovers), recombination works properly and we don't have to worry about adjust the recombination events in inverted homozygous individuals.

#### How are the m2 and m3 locations that mark the inversion start and end converted into inv_start and inv_end, and are inv_start and inv_end modeled as inclusive of the location they are on, or exclusive? Where are m2 and m3 on your conceptual graphs?

We didn't do this in our simulation. This is the part that I modeled differently. Our m2 mutations are our QTNs and m3 are our inversions. The m3 position is the inversion starting position and we draw a value from the uniform distribution to get the length of the inversion and store that as the mutations tag value. I grab the tag value add it to the start to get the end point. I store each of those values as inv_start and inv_end. ```start = invS.position; end = start + invS.tag; ```. Here is an in depth explaination from Peter about exactly this:
https://github.com/MesserLab/SLiM/issues/203#issuecomment-884635953
but essentially the inversion is from [X,Y] inclusive. So we need to subset for breakpoints that are more than X or less than or equal to Y because of how the breakpoints work: breakpoint at n means [0, n-1] is inherited together and [n, L] are inherited together. 
 
####  I’m getting confused between the different problems of the single crossover within the inversion, a crossover at one inversion end,  a crossover at two inversion ends, a triple crossover inside the inversion… it would help me to see a conceptual example of each one.
  
Fair the end dynamics are still confusing me. I do an example of how I think this should work below. I'm going to use the same syntax as the other group to make things more streamlined with the second chromosome having lower case letters:
  
#### double cross over
  
In this example, inside the inversion is between positions [X,Y] and A is just to the left of X and B is just to the right of position X in the standard orientation.   

```

# standard (uninverted):
   X  			   Y
---+-------+---------------+-----
  A B      W              C D

# inverted:

   X  			   Y
---+-------+---------------+-----
  A C      W              B D


```
two break points inside the inversion **how SLiM does it**:
  
```

   X  	*	*	   Y
---+-------+---------------+-----
  A B      W              C D


   X  		 	   Y
---+-------+---------------+-----
  a b      w              c d

After recombination:

   X  			   Y
---+-------+---------------+-----
  A B      w              C D

   X  			   Y
---+-------+---------------+-----
  a b      W              c d
```
  
**how we want SLiM to do it**:
  
```

   X  	*	*	   Y
---+-------+---------------+-----
  A C      W              B D


   X  			   Y
---+-------+---------------+-----
  a c      w              b d

After recombination:

   X  			   Y
---+-------+---------------+-----
  A C      w              B D

   X  			   Y
---+-------+---------------+-----
  a c      W              b d
```
  
As you can see the end dynamics stay the same when there are two break points inside the inversion no matter what orientation the inversion is in. But again this doesn't solve the haplotypes inside the inversion. In this examples, I'm just focused on what they are talking about in the post which is the end dynamics, but you can imagine the orientation of the inversions means if there is variation inside the inversion (multiple haplotypes of the inversion segregating in the population) then those won't be inherited properly. 
   


#### breakpoint at exactly the start position of the inversion and one inside (i.e., odd number of breakpoints inside the inversion)
  
For this example, we have a breakpoint at the location where X falls in the standard orientation which is the starting location of the inversion and one in the inversion. X is a locus position that falls between A and C. Remembering that a breakpoint at X means [0, X-1] and [X-L] are inherited together.
  
**how SLiM does it**

```
   *              *
   X  			   Y
---+-------+---------------+-----
  A B      W              C D


   X  			   Y
---+-------+---------------+-----
  a b      w              c d

After recombination:

   X  			  Y
---+-------+---------------+-----
  A b      w              C D

   X  			   Y
---+-------+---------------+-----
  a B      W              c d
```
   
**how we want SLiM to do it**

```
   *              *
   X  			   Y
---+-------+---------------+-----
  A C      W              B D


   X  			   Y
---+-------+---------------+-----
  a c      w              b d

After recombination:

   X  			   Y
---+-------+---------------+-----
  A c      w              B D

   X  			   Y
---+-------+---------------+-----
  a C      W              b d
```
Here A is inherited with B and D, but not C which is what we want for an odd number of breakpoints inside the inversion and a breakpoint at exactly one end. Now lets force SLiM to do this from the standard orientation.
    
**flip breakpoint to the other side to fix it**

```
                  *         *
   X  			   Y
---+-------+---------------+-----
  A B      W              C D


   X  			   Y
---+-------+---------------+-----
  a b      w              c d

After recombination:

   X  			   Y
---+-------+---------------+-----
  A B      W              c D

   X  			   Y
---+-------+---------------+-----
  a b      w              C d


```
  
Now as you can see the A, B, and D are inherited together and C isn't which is what we want. Again this is another example where the W isn't inherited properly and I'm not sure how to handle that. 

#### breakpoints at exactly X and Y + 1 plus one in the middle (this would be extremely rare in our sims)

**how SLiM does it**

```
   *		*	    *
   X  			   Y
---+-------+---------------+-----
  A B      W              C D


   X  			   Y
---+-------+---------------+-----
  a b      w              c d

After recombination:

   X  			   Y
---+-------+---------------+-----
  A b      w              C d

   X  			   Y
---+-------+---------------+-----
  a B      W              c D
```

**how we want SLiM to do it**

```
   *		*	    *
   X  			   Y
---+-------+---------------+-----
  A C      W              B D


   X  			   Y
---+-------+---------------+-----
  a c      w              b d

After recombination:

   X  			   Y
---+-------+---------------+-----
  A c      w              B d

   X  			   Y
---+-------+---------------+-----
  a C      W              b D
```
Here A and B and C and D are inherited together. Now lets force SLiM to do this from the standard orientation.  

**remove breakpoints from the two ends** 

```
   		*		    
   X  			   Y
---+-------+---------------+-----
  A B      W              C D


   X  			   Y
---+-------+---------------+-----
  a b      w              c d

After recombination:

   X  			   Y
---+-------+---------------+-----
  A B      W              c d

   X  			   Y
---+-------+---------------+-----
  a b      w              C D
```

Now A and B then C and D are each inherited together which is how we want it to work.
