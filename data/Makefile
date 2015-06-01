#!/usr/bin/make -f

SSCONVERT ?=ssconvert
SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data

FPKM_FILES=Supplemental_Table_2_Homo_NEW_UPDATED.txt $(patsubst %.xlsx,%.txt,$(filter-out $(wildcard *Homo*), $(wildcard *_fpkm.xlsx)))


all: $(FPKM_FILES)


HOM_MouseHumanSequence.rpt: 
	wget -O $@ "ftp://ftp.informatics.jax.org/pub/reports/HOM_MouseHumanSequence.rpt"

HOM_AllOrganism.rpt:
	wget -O $@ "ftp://ftp.informatics.jax.org/pub/reports/HOM_AllOrganism.rpt"

combined_fpkm: combine_fpkm.R HOM_AllOrganism.rpt $(FPKM_FILES)
	$(R) $(ROPTS) -f $< --args $@ $(wordlist 2,$(words $^),$^)

%.txt: %.xlsx
	$(SSCONVERT) --export-type Gnumeric_stf:stf_assistant -O 'separator="	" quote="" eol="unix"' $< $@
