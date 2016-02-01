#!/usr/bin/make -f

SSCONVERT ?=ssconvert
SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data
PERL?=perl


all: combined_ancestral_states combined_ancestral_plots.pdf

combined_ancestral_states: ancestral_state.R combined_fpkm_per_sample gene_trees gene_list.txt
	$(R) $(ROPTS) -f $< --args $(wordlist 2,3,$^) $(shell cat gene_list.txt) $@


combined_ancestral_plots.pdf: plot_combined_states.R combined_ancestral_states
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

clean:
	rm -f combined_ancestral_states combined_ancestral_plots.pdf

.PHONY: clean
