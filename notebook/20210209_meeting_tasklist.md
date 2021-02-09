#KEL and SMS

- Sara outputting inversion capture stats

Action items: 

* inversions within inversions - sara will check how often this happens only in last generation
* output smaller VCF files - to 100 individuals per pop and filter to MAF 0.01
* muts dataframe: sara will add inversion mutations
* using MAF 0.01
* analyze control sims first
* ID inversion windows for analysis - FIX CURRENT CODE FOR THIS
* MERGE KATIE AND SARA'S R CODES
* write up the methods and results
* THEN, RERUN ALL THE SIMS

How are we going to label important inversions?
* FST outlier maybe - first we need to show control sims do not give outliers
  * Way 1: run outflank on all controls sims and show no outliers ---> then we're good
  * Way 2: compare inversion sims to control sims --> empirical outlier
* VA? We discussed this, but kept coming back to FST
  * We can get VA for each inversion window --> 100 SNPs in the window
  * We can get VA for the non-inverted genome scaled to inersion window size - randomly sample 100 SNPs from genome and VA 1000 times to create a null
  * Call outlier

https://htmlpreview.github.io/?https://github.com/SaraSchaal/CoreSim_QTL/blob/master/notebook/2021_02_02InversionMovie_vcf_KEL.html
