#!/usr/bin/make -f

SSCONVERT ?=ssconvert
SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data

FPKM_FILES=$(patsubst %.xlsx,%.txt,$(wildcard *_fpkm.xlsx))


all: $(FPKM_FILES)


HOM_MouseHumanSequence.rpt: 
	wget -O $@ "ftp://ftp.informatics.jax.org/pub/reports/HOM_MouseHumanSequence.rpt"

combined_fpkm: combine_fpkm.R HOM_MouseHumanSequence.rpt $(FPKM_FILES)
	$(R) $(ROPTS) -f $< --args $@ HOM_MouseHumanSequence.rpt $(FPKM_FILES)

%.xls: %.gnumeric
	$(SSCONVERT) --export-type Gnumeric_Excel:excel_dsf $< $@

%.txt: %.gnumeric
	$(SSCONVERT) --export-type Gnumeric_stf:stf_assistant -O 'separator="	" quote="" eol="unix"' $< $@

%.txt: %.xls
	$(SSCONVERT) --export-type Gnumeric_stf:stf_assistant -O 'separator="	" quote="" eol="unix"' $< $@

%.txt: %.xlsx
	$(SSCONVERT) --export-type Gnumeric_stf:stf_assistant -O 'separator="	" quote="" eol="unix"' $< $@

%.xls: %.txt
	txt2xls $< > $@

