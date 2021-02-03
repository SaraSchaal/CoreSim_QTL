## Trying to find the best way to plot origin dynamics of inversions. 
For example, do inversions first arise neutrally then gain QTNs through later mutation or do inversions capture QTNs at the time that it mutates. In my example simulation, all inversions that were present in the final generation captured QTNs at the time of mutation but then gained more through time. 

In this plot I have the generation the inversion arose on the x axis, the sum of the effect sizes of the qtns found in the inversion when it mutates with the size being FST of the inversion by the end of the simulation and color being which population it has the highest frequency in.
  
![Inversion Origin with FST and QTN effects](../figures/OriginDynamics/20210203_origin_FSTandQTNeffect.pdf)

In this plot, I have the origin generation of the inversion on the x-axis with the number of QTNs captured when the inversion mutates on the y-axis in panel 1 and the number of QTNs in the final generation on the y-axis of panel 2
  
![Inversion Origin with number of QTNs](../figures/OriginDynamics/20210203_origin_numQTNs.pdf)

What we can conclude:
  
At least with the simulation we’ve been working with. All inversions that persist through the simulation capture some number of qtns at the time it first mutates then gains more throughout the simulation.
 
Interestingly some inversions that persist through the whole simulation arose during the burnin period when selection wasn't occuring. 
  
Follow up questions:
  
1. Do we want to just look at the inversions in the final generation or do we want to try to get a plot through time? 
2. Do we want to calculate the sum of the QTN effects on the phenotype through time instead of the mean and sd?
