#!/usr/bin/make -f

SSCONVERT ?=ssconvert
SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data
PERL?=perl

SPECIES=ateles_fusciceps bos_taurus canis_familiaris				\
dasypus_novemcinctus equus_caballus homo_sapiens loxodonta_africana	\
monodelphis_domestica mus_musculus spalax_galili ovis_aries	\
pan_paniscus spalax_carmeli sus_scrofa

SPECIES_FPKM=$(foreach specie,$(SPECIES),$(specie)_mean_fpkm.txt)

all: $(SPECIES_FPKM)


$(SPECIES_FPKM): %_mean_fpkm.txt: species_fpkm.R ../combined_fpkm
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $* $@
