Hi Vince, Andrew, & Peter,

I was trying to make sure your code was implemented in the way that I expect it to based on the logic you show here: https://github.com/MesserLab/SLiM/issues/203#issuecomment-884635953. That logic makes sense to me. I've written out all the different scenarios on paper and they all check out. When I tested the code, however, I was not getting breakpoint additions how I expected for an inversion homozygote. Keep in mind I may be overlooking something so please correct me if I am way off base! 

The way the logic works is that if there is a single breakpoint inside the inversion, we want to add breakpoints on both the inv_start and the inv_end + 1. To handle this scenario, the current code is written like this: ```breakpoints = sort(c(breakpoints[!(left | right)], c(inv_start, inv_end + 1)[c(sum(left) > 0, sum(right) > 0)]));```. 

Take the example here where we have a breakpoints vector that includes one breakpoint inside the inversion, but no breakpoints on either the start or end + 1:

```
breakpoints = c(6, 7, 12, 23, 49, 51);
inv_start = 10;
inv_end = 21;
left = (breakpoints == inv_start); # F F F F F F
right = (breakpoints == inv_end + 1); # F F F F F F
breakpoints = sort(c(breakpoints[!(left | right)], c(inv_start, inv_end + 1)[c(sum(left) > 0, sum(right) > 0)])); # 6 7 12 23 49 51
```

What I would expect here is that the code would add a breakpoint at start and at end + 1, but what it outputs is just the same breakpoints vector. I think this comes down to the fact that it only adds inv_start and inv_end + 1 conditioned on whether those logical vectors have trues in them. I think what we want here instead is that they do NOT have trues in them. It became a little more complicated when thinking through the end dynamics and what needs to be added to those. But here is my proposed change that 1) adds the breakpoints at the start and end + 1 if there is a breakpoint inside but not at either end, 2) flips the breakpoint to the opposite end if there is one inside the inversion and one at either the start or end + 1, but not both, and 3) removes both breakpoints at the start and end + 1 if they were both drawn as breakpoints and there is one inside the inversion. 

Scenario 1 - add breakpoints at the start and end + 1 if there is a breakpoint inside but not at either end 
```
breakpoints = c(6, 7, 12, 23, 49, 51);
inv_start = 10;
inv_end = 21;
left = (breakpoints == inv_start); # F F F F F F 
right = (breakpoints == inv_end + 1); # F F F F F F 

if(sum(c(left,right)) == 0){ 
	# SLiM does not call breakpoints at the inversion ends, we add them
	breakpoints = sort(c(breakpoints,c(inv_start, inv_end + 1)));
} else if(sum(left) > 0 & sum(right) > 0){
	# SLiM call breakpoints at both inversion ends, we remove them
	breakpoints = sort(c(breakpoints[!(left | right)]));
} else if (sum(left) > 0 | sum(right) > 0){
	# SLiM calls a breakpoint at one inversion end, we remove the one that SLiM calls and add the one that slim doesn't call
	breakpoints = sort(c(breakpoints[!(left | right)],c(inv_start, inv_end + 1)[c(sum(right)>0,sum(left)>0)]));
}
breakpoints # 6 7 10 12 22 23 49 51

```

Scenario 2a - flip the breakpoint to the opposite end if there is one inside the inversion and one at either the start or end + 1, but not both. Here is an example with a breakpoint at the start.
```
breakpoints = c(6, 10, 12, 23, 49, 51);
inv_start = 10;
inv_end = 21;
left = (breakpoints == inv_start); # F T F F F F 
right = (breakpoints == inv_end + 1); # F F F F F F 

if(sum(c(left,right)) == 0){
	# SLiM does not call breakpoints at the inversion ends, we add them
	breakpoints = sort(c(breakpoints,c(inv_start, inv_end + 1)));
} else if(sum(left) > 0 & sum(right) > 0){
	# SLiM call breakpoints at both inversion ends, we remove them
	breakpoints = sort(c(breakpoints[!(left | right)]));
} else if (sum(left) > 0 | sum(right) > 0){ 
	# SLiM calls a breakpoint at one inversion end, we remove the one that SLiM calls and add the one that slim doesn't call
	breakpoints = sort(c(breakpoints[!(left | right)],c(inv_start, inv_end + 1)[c(sum(right)>0,sum(left)>0)]));
}
breakpoints # 6 12 22 23 49 51

```

Scenario 2b - Here is an example with a breakpoint at the end + 1.
```
breakpoints = c(6, 7, 12, 22, 49, 51);
inv_start = 10;
inv_end = 21;
left = (breakpoints == inv_start); # F F F F F F 
right = (breakpoints == inv_end + 1); # F F F T F F 

if(sum(c(left,right)) == 0){
	# SLiM does not call breakpoints at the inversion ends, we add them
	breakpoints = sort(c(breakpoints,c(inv_start, inv_end + 1)));
} else if(sum(left) > 0 & sum(right) > 0){ 
	# SLiM call breakpoints at both inversion ends, we remove them
	breakpoints = sort(c(breakpoints[!(left | right)]));
} else if (sum(left) > 0 | sum(right) > 0){
	# SLiM calls a breakpoint at one inversion end, we remove the one that SLiM calls and add the one that slim doesn't call
	breakpoints = sort(c(breakpoints[!(left | right)],c(inv_start, inv_end + 1)[c(sum(right)>0,sum(left)>0)]));
}
breakpoints # 6 7 10 12 49 51

```


Scenario 3 -  remove both breakpoints at the start and end + 1 if they were both drawn as breakpoints and there is additional breakpoint inside the inversion
```
breakpoints = c(6, 10, 12, 22, 49, 51);
inv_start = 10;
inv_end = 21;
left = (breakpoints == inv_start); # F T F F F F 
right = (breakpoints == inv_end + 1); # F F F T F F 

if(sum(c(left,right)) == 0){
	# SLiM does not call breakpoints at the inversion ends, we add them
	breakpoints = sort(c(breakpoints,c(inv_start, inv_end + 1)));
} else if(sum(left) > 0 & sum(right) > 0){ 
	# SLiM call breakpoints at both inversion ends, we remove them
	breakpoints = sort(c(breakpoints[!(left | right)]));
} else if (sum(left) > 0 | sum(right) > 0){
	# SLiM calls a breakpoint at one inversion end, we remove the one that SLiM calls and add the one that slim doesn't call
	breakpoints = sort(c(breakpoints[!(left | right)],c(inv_start, inv_end + 1)[c(sum(right)>0,sum(left)>0)]));
}
breakpoints # 6 12 49 51

```

I may definitely be overlooking something that would make the original code work, so please correct me if I am wrong. But I would appreciate your thoughts on this! 
