Data Directory
==============

This directory contains the rules to make various data files

+ `gene_trees` -- contains the ensembl gene trees parsed from
   Compara.80.protein.nhx.emf.gz. Run `make gene_trees` to build it.
    -  trees -- ensembl gene trees
    -  prot.to.tree -- mapping from ensembl protein id to gene trees
    -  genes.to.tree -- mapping from ensembl gene id to gene trees
+ `combined_fpkm` -- contains the combined fpkm counts from the fastq
   files which are aligned in `../fastq`
+ `homology_table` -- contains the fpkm counts of all human 1:1
   orthologs in all species

Examine [Makefile](Makefile) for details about them.

### diamond_against_human subdirectory

This directory contains rules to generate a blast of all species
against human, primarily useful for annotating un-annotated proteins
to get an idea of what they might actually do.

The important file is `diamond_against_human/gene_to_human_gene.txt`
which is combined into combined_fpkm.

