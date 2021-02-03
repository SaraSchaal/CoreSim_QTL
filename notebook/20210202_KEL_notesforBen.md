# Notes for meeting with Ben and Sara tomorrow

- SLiM output
  - the muts dataframe is missing m1 and m3 mutations
    -  The FSTs for the neutral loci aren't outputting with the muts file. So haven't looked at FSTs of neutral loci yet

- outputting VCF file - filter for muts to match Sara's output and filter for 100 individuals from each subpop
  - when I filter the VCF file for MAF 0.05 with VCFtools, there is not a perfect line-up of mutations with Sara's output in her muts file (rounding error)
  - the VCF file is huge - we only need 100 individuals from each subpop to show dynamics
  
- "overlapping" inversions
  - I don't think it's serious, but we haven't evaluated how often this occurs
  - There's also the potential for "inversion within an inversion", which is troublesome because it should have recombination. Again, don't know how often, but looked like it was happening in one video.
  
  ## Ask Ben about evolving seas RCN SLiM Workshop in exchange for stipend
  
  # TO DO
  
  -  FST distributions within vs. outside inversions vs. neutral loci
    - need to write more code
  
  - Proportion of FST outliers within inversions vs. outside
    - Add OutFLANK analysis to Katie's R script
  
  - How to identify 'important' inversion windows? Do we really need this?
    - mean FST within inversion window
    - proportion of Va within inversion window (to do this, need to add inversion ID to each mutation ID in Katie's R script - straightforward)
    - direction of effect on the genotype (+ or minus) * # copies * FST
