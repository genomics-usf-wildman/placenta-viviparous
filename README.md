Placenta Viviparous Paper
-------------------------

This repository contains all of the code necessary to generate the
analyses presented in this paper, including the manuscript itself.

Subdirectories
--------------

* data -- basic analyses and combining of aligned/called sequences
* fastq -- fastq files and makefiles to assemble, align, call, and
  annotate all of the sequences in all of the species examined
* manuscript -- R and LaTeX code to produce the published manuscript
* figures -- various complicated figures which are produced separately
  from the manuscript
* alignments -- alignments of genomic regions of interest in species
  studied with known genomes
* promoters -- promoter peak data for placenta samples in humans


Precomputed Data (Missing Supplemental Material)
------------------------------------------------

* [all_species_mean_fpkm.txt.xz](https://www.donarmstrong.com/ld/pv2016/all_species_mean_fpkm.txt.xz) contains
  all of the mean FPKM for each species. This was supposed to be
  included in the supplemental material of the paper, but was missing
  for some reason.
* [all_species_per_sample_fpkm.txt.xz](https://www.donarmstrong.com/ld/pv2016/all_species_per_sample_fpkm.txt.xz) contains
  the FPKM per sample across all samples in each species. This file
  was also missing from supplemental material for some reason.

* [placenta_core_transcriptome.txt.xz](https://www.donarmstrong.com/ld/pv2016/placenta_core_transcriptome.txt.xz)
  contains the FPKM of all genes in the core placenta transcriptome in
  long format. [One row per gene in each species, including the "all
  but one", "all", "all but metatherian" minimum FPKM in the set of
  species as depicted in the graphs included in the paper.]
* [placenta_classification_shape_p_fdr.txt.xz](https://www.donarmstrong.com/ld/pv2016/placenta_classification_shape_p_fdr.txt.xz)
  contains the p, fdr, and overall fdr values of the glm of the
  association of gene expression with placenta shape, interdigitation,
  and intimacy which are shown in supplemental figures 57,58, and 60.
  [I have excluded the ordered analysis for clarity.]
  
