#!/usr/bin/make -f

SSCONVERT ?=ssconvert
SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data
PERL?=perl

all: all_species_mean_fpkm.txt.xz all_species_per_sample_fpkm.txt.xz

%.txt.xz: %.txt
	xz $<

all_species_mean_fpkm.txt: species_fpkm.R ../combined_fpkm
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $* $@

all_species_per_sample_fpkm.txt: species_per_sample_fpkm.R ../combined_fpkm_per_sample
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $* $@

