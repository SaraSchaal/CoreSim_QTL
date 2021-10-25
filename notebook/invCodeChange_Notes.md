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


## Example of what we want slim to do when we have an inversion at [X,Y]   
In this example, D is at Y + 1 and A is at X - 1 in the standard orientation and C is some locus in the middle of the inversion just as an example locus to follow: 
```	
	AX--C----YD standard

	AY----C--XD inverted
```

Homozygote Inverted (* is proposed breakpoint):
```
	   *  
	AY----C--XD - chrom1

	..    .  ..
	AY----C--XD - chrom2
```

After recombination:
```
      	      .  ..
	AY----C--XD

	..    
	AY----C--XD

```


## Example of what SLiM ACTUALLY does with this proposed recombination break regardless of whether an individual is homozygous inverted or not  
SLiM only recognizes the standard orientation in the genetic map. Below is an example of how SLiM handles recombination. 

Homozygote whether inverted or not (* is proposed breakpoint):
```
	   *  
	AX--C----XD - chrom1

	..  .    ..
	AX--C----YD - chrom2

```
After recombination:
```   
       	    .    ..
	AX--C----YD

	..    
	AX--C----YD

```
So we manually have to make a recombination event with the standard orientation that results in the same inheritance as what we showed in the "what we want slim to do" example.


## Example of how to make the recombination for an individual homozygous for the inversion to match the first example  
We do this by adding breakpoints at X and Y + 1 and again * represent proposed breakpoints. We want AY to be inherited together and X and D to be inherited together:  
```
	 * *      *
	AX--C----YD - chrom1

	..  .    ..
	AX--C----YD - chrom2
```

After recombination:
```
	 .        .
	AX--C----YD

        .   .    .
	AX--C----YD
```
																	   
Here by adding those new breakpoints at X and Y + 1, we successfully get AY and the dotted AY to be inherited together and the dotted XD and XD to be inherited together. What doesn't match is the haplotypes inside the inversion. You can see that C is now inherited with dotted XD and not XD. Is this a problem? It could change the number of segregating haplotypes. (UNRESOVELED QUSTION)
