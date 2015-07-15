#!/usr/bin/make -f

SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data


all: hormone_gene_expression.pdf

hormone_gene_expression.pdf: hormone_genes.R \
	../../data/combined_fpkm ../../data/go_gene_association \
	../species_ordering.R
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

