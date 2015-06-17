#!/usr/bin/make -f

SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data


all: top_gene_expression.pdf

top_gene_expression.pdf: top_gene_expression.R ../../data/combined_fpkm ../../data/go_gene_association
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

