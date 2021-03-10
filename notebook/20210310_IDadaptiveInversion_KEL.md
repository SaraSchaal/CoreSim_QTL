# Categorization of inversions harboring genetic basis of adaptation:

In order to answer our research questions about the characteristics of adaptive inversions, we had to determine from our simulations which inversions were underlying the genetic basis of local adaptation. Because of the dynamic nature of our simulations and the inference of local adaptation as a metapopulation phenomenon, this was not immediately straightforward for a few reasons. First, an inversion could arise in frequency in one set of individuals at a genetic map location that partially overlapped with a common inversion in a different set of individuals on the genetic map (even though overlapping inversions within an individual was not allowed in the code) - see Supp Figure for example(?). To address this first point, we identified inversion windows on the genetic map as contiguous sets of inversion QTNs. A QTN was determined to be an inversion QTN if it was found in any individual that had an inversion (MAF > 0.01) encompassing it. 

Additional reasons needed to be carefully considered when identifying inversion windows as adaptive. Inversions could arise and contain QTNs, but not contribute to phenotypic divergence (e.g., inversion haplotypes not differentiated among populations). Conversely, inversions could arise and drift to different frequencies in each population (especially under the low gene flow scenarios), but the amount of differentiation in the inversions was not more than observed in control simulations where individuals had equal fitness (and the inversion mutation rate was the same). In certain areas of the parameter space, comparing inversion to no inversion simulations clearly led to local adaptation and adaptive inversions clearly had signals of selection (i.e., inversions increased the amount of local adaptation and harbored the genetic basis). In other areas of the parameter space, inversion and no-inversion simulations led to the same amount of local adaptation, but some inversions had signals of selection (i.e., inversions did not increase the amount of local adaptation, but still harbored the genetic basis of adaptation).

Therefore, we developed specific criteria to determine if an inversion was adaptive (i.e., harbored the genetic basis of trait differentiation): 

(i) the focal inversion window contained more FST outlier QTNs per unit map distance when compared to a null distribution generated from QTNs in inversion windows from a paired control simulation with the same inversion mutation rate but no selection (this identifies inversion outliers above that expected by genetic drift in a control scenario); 

(ii) the focal inversion window contained more FST outlier QTNs per unit map distance when compared to a null distribution generated from QTNs in non-inversion windows from the same simulation (this identifies outliers above that expected by genetic drift in the same genome, and is a typical approach used in population genomic studies). The issue with using this criteria without criteria (i) is that neutrally drifting inversions could be enriched with FST outliers due to low recombination; therefore it is of interest to determine how reliable this criteria is for identifying adaptive inversions. (The potential issue with relying entirely on criteria (i) is that selection in the highly polygenic simulations was affecting the entire genome, so it isn’t clear to me yet if drift operates the same way in the selection vs. control sims).

(iii) the QTNs in the focal inversion window collectively harbor a greater percentage of the additive genetic variation (VA) in the trait when compared to a null distribution generated from QTNs in non-inversion windows from the same genome.

(iv) the focal inversion window contained more SNPs per map distance when compared to a null distribution generated from non-inversion windows from the same genome. We used this criteria because (idk?) theory predicts that the establishment probability of like-effect QTNs is higher when they are in close proximity (see Yeaman papers) -  as adaptive inversions started to diverge, this process led to a concentration of SNPs within adaptive inversions.

(v) ? the focal inversion window contained more SNPs per map distance when compared to a null distribution generated from inversion windows from the paired control simulation.

(I”M NOT SURE WHICH OF THE LAST TWO MAKE MORE SENSE. (v) MIGHT BE TOUGH SINCE THERE AREN;T THAT MANY INVERSIONS).

To generate the null distribution, we cut the genome into 1000 chunks (500 bases each? need to check) and labeled them as inversion or non-inversion windows. Chunks that spanned window boundaries were removed. Depending on the criteria, a focal inversion window was labeled as meeting that criteria if it’s observed test statistic was outside the 99% quantile for the relevant null distribution.


* Cut genome into 1000 chunks, remove anything that overlaps with inversion, label as inversion and non-inversion windows
  * Va_out_scaled = sum(Va_perc[“OUT”]))/map_distance (chunk i)
  * SNPs_out_scaled = # snps / map distance
  * num outliers compared to control sims / map_distance
  * num outliers OutFLANK / map_distance

For inversion window i, 
  * get Va_perc_inv_i = sum(Va_perc[inv_i])) - all SNPS in that inversion window
  * Va_inv_i_scaled = Va_inv_i/(map distance inversion window)
  * Get # snps / map distance
  * num outliers compared to control sims / map_distance
  * num outliers OutFLANK / map_distance

compare obsreved to > null dsitribution
