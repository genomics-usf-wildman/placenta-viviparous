#!/usr/bin/make -f

SSCONVERT ?=ssconvert
SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data
PERL?=perl

all: go_tissue_analysis

# include species.mk

# species.mk: species_list.R combined_fpkm
# 	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@


go_tissue_analysis: go_tissue_analysis.R combined_fpkm housekeeping_genes_superset
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

