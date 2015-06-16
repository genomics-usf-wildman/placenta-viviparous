#!/usr/bin/make -f

R=R
ROPTS=-q --no-save --no-restore-data


GENE_FPKM_FILES=$(wildcard *_genes.fpkm_tracking)
ISOFORM_FPKM_FILES=$(wildcard *_isoforms.fpkm_tracking)


all: combined_fpkms

combined_fpkms: combine_fpkms.R $(GENE_FPKM_FILES) $(ISOFORM_FPKM_FILES)
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@
