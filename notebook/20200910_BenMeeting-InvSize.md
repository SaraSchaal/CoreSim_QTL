## 20200910 - NAs when I merge files together in R when mutation ID is > 1,000,000 ##
# Summary of issue #

After I merge the summary data file with the inversion dynamics through time file, there are many rows that have NA for the summary data about that inversion. It happens after the inversion mutation id is > 1,000,000 and then seems to randomly add in inversion data after that. After inspecting the metadata file, these inversions do not have data in this file. This means it is an issue during the simulation not in R afterwards. 

# Questions
1) What happens when it is greater than 1,000,000? My first thought was that it has to do with scientific notation, but I don't see that anywhere in my file.
*ANSWER* it was in the file just buried in the through time file not the metadata file
2) Are the mutation IDs also unique across different mutation types? So this would explain why we only see this under high QTN mutation scenarios because we are reaching that 1,000,000 id much sooner under the higher mutation rate. 
*ANSWER* yes they are all given a unique number across every type of mutation (e.g., m1, m2, m3, etc.)

# Explanation 
Okay - the issue IS scientific notation! when you use the c() function it will coerce all items within the c() function to a single type. In the outputInvSumInfo.txt file everything is converted to an integer which is no problem. However, the bigger file which is the inversions through time file,  outputInvTime.txt, they get converted to a float which converts the mutation IDs to scientific notation. They no longer have any match in the metadata file. So the metadata file wasn't the issue, the inversions through time file was. 

Response from Ben: 

	BH: Note that some of the mutation IDs are being output in scientific notation; that loses significant digits
	doing paste(c(…)) is not really a good way to construct output; much better is something like: "" + a + " " + b + " " + …
	That way each value a b etc. gets converted to a string on its own, in the best way available

	When you do c() you tell Eidos to throw all the values together into one vector, which has to be of a single type; so the integer mutation ID gets coerced upward to be a float, and thus gets output as a float

	SMS: Oh goodness. my first inkling of it being scientific notation was correct. I had checked the metadata file thinking it was that one which doesn’t have this issue, but not the InvTime output. Why would it do it in one but not the other?

	BH: It’s due to the way you use c().  In the second file’s output, one of the values you are concatenating together with c() must be a float, so everybody gets coerced up to float.  In the first file’s output, all the values must be integers.  Implicit type coercion is problematic…

	Yeah, the frequencies and the means are float, so everybody goes up to float because of c()

	No worries.  It’s a good bug.  I’m thinking about whether there’s a way to improve Eidos to avoid it.

	R does the same thing with c()

	Perhaps the root of the problem is that my paste() can’t take multiple arguments; you wouldn’t be using c() if it could.  I should upgrade the priority of that work.The other possible fix would be to disallow implicit type conversion with c() in the first place.  It really is kind of an evil language feature.  R does lots of hidden/magic/implicit things that can lead to hard-to-find bugs.  In many cases I chose not to follow that precedent in Eidos, but this is a case where I did, and perhaps it was a mistake.

	I could add a new c_coerce() that provides the old behavior, and make c() give an error if mixed types are passed to it.  I wonder how much end-user code I would break with that change?

	SMS: Yeah I think it would be a good thing to get rid of the implicit conversion. Especially given that I didn’t know this was even an issue in my code until we looked at weird dynamics in results plots and tried to dive deep into why these were happening. I’m trying to think of a scenario of when I would want it to convert everything…nothing off the top of my head, but your suggestion of having two different functions makes sense to me. Would it be less of an issue to end-user code if you left c() as is and instead of giving an error message give a warning that mixed types were passed to c() and then make a second function c_noCoerce() .. or something like that.

	BH: Yeah, perhaps.  I’ll ponder it.  :->  Anyway, thanks for the interesting bug, good luck!



