#!/usr/bin/make -f

SSCONVERT ?=ssconvert
SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data
PERL?=perl

all: $(GO_ANALYSES)

include species.mk

species.mk: species_list.R combined_fpkm
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

GO_ANALYSES=$(patsubst %,go_%_mf,$(SPECIES)) $(patsubst %,go_%_bp,$(SPECIES))

go_tissue_analysis: go_tissue_analysis.R combined_fpkm housekeeping_genes_superset
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

