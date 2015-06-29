#!/usr/bin/make -f

SSCONVERT ?=ssconvert
SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data


all: go_functional_grouping heatmap_gene_by_functional_grouping.pdf


go_functional_grouping: get_functional_grouping.R ../../data/combined_fpkm
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

heatmap_gene_by_functional_grouping.pdf: heatmap_gene_by_functional_grouping_plot.R \
	../../data/combined_fpkm \
	go_functional_grouping
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@
