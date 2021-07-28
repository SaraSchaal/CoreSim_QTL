Currently inversion characteristics are calculated by averaging the characteristic (i.e., age, size, or scaled num QTNs) across all inversions that are found 
in any replicate for a parameter combination:
```mean(c(inv1, inv2, ..., invN))```
so that is the second way you were describing yesterday. For total LA and percent VA I get a single value that I then average across the 5 reps, 
but for all the characteristics I use separate dataframe with a line per inversion and I aggregate over all replicates to get those points and 
standard deviations in the point plot. 

head of inversion characteristics file:
```
     seed  inv_id inv_age  inv_length num_qtns_Lscaled     adaptInv muBase muInv sigmaK alpha rep   n  enVar mig1 
1  3383674 1263852   38059      22607     0.0009289158     Adaptive  1e-07 0.001   0.75 0.002   3 1000   0.1  0.1
2  3383674 1991588   31194      46681     0.0010496776     Adaptive  1e-07 0.001   0.75 0.002   3 1000   0.1  0.1
3  3383674 3873076   13379      19996     0.0009501900     Adaptive  1e-07 0.001   0.75 0.002   3 1000   0.1  0.1
4  3383674 4050760   11698      44368     0.0009917057     Adaptive  1e-07 0.001   0.75 0.002   3 1000   0.1  0.1
5  3383674 4120319   11037      49761     0.0009646108     Adaptive  1e-07 0.001   0.75 0.002   3 1000   0.1  0.1
6  3383674 4876035    3883      33785     0.0010063638     Adaptive  1e-07 0.001   0.75 0.002   3 1000   0.1  0.1
7  3383674 5086344    1898      28474     0.0006672754     Adaptive  1e-07 0.001   0.75 0.002   3 1000   0.1  0.1
8  3383674 3835642   13735      34414     0.0007264485  Nonadaptive  1e-07 0.001   0.75 0.002   3 1000   0.1  0.1
9  3383674 4110800   11126      27987     0.0007503484  Nonadaptive  1e-07 0.001   0.75 0.002   3 1000   0.1  0.1
10 3383674 4418551    8218        311     0.0000000000  Nonadaptive  1e-07 0.001   0.75 0.002   3 1000   0.1  0.1
11 3383674 4654764    5981      26065     0.0009591406  Nonadaptive  1e-07 0.001   0.75 0.002   3 1000   0.1  0.1
12 3383674 4873573    3906      26111     0.0005744705  Nonadaptive  1e-07 0.001   0.75 0.002   3 1000   0.1  0.1

```

aggregate line:
```
df.invChar.muInv3.av <- aggregate(cbind(inv_age, inv_length, num_qtns_Lscaled)~adaptInv + muBase + muInv + sigmaK + alpha + enVar + mig1 + mig2, data = df.invChar.muInv3, FUN = mean)

```

The boxplots are then just the raw data across all replicates. So for each parameter combination thats 
the distribution of all the sizes, ages, and scaled number of QTNs across all 5 replicates for the 3 inversion categories (i.e., adaptive, nonadaptive, 
no-selection). After talking yesterday though, I thought I should be doing it at the level of the replicate.  Now looking at the number of adaptive 
inversions and how big those SD are I think the results have the potential to be biased int he simulations where we had a lot of inversions evolving. 
That bias will also vary depending on the parameters because of how different some of the SD are. 

