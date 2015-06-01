#!/usr/bin/make -f

SSCONVERT ?=ssconvert
SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data
PERL?=perl


FPKM_FILES=Supplemental_Table_2_Homo_NEW_UPDATED.txt $(patsubst %.xlsx,%.txt,$(filter-out $(wildcard *Homo*), $(wildcard *_fpkm.xlsx)))


all: $(FPKM_FILES)


go_gene_association.txt: go_gene_association.pl go-basic.obo gene_association.goa_human
	$(PERL) $^ |grep -v '^!' > $@

go_gene_association: go_gene_association.R go_gene_association.txt
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

go-basic.obo:
	wget -O $@ "http://purl.obolibrary.org/obo/go/go-basic.obo"

gene_association.goa_human: gene_association.goa_human.gz
	gzip -dc $^ > $@

gene_association.goa_human.gz: 
	wget -O $@ "http://geneontology.org/gene-associations/gene_association.goa_human.gz"

HOM_MouseHumanSequence.rpt: 
	wget -O $@ "ftp://ftp.informatics.jax.org/pub/reports/HOM_MouseHumanSequence.rpt"

HOM_AllOrganism.rpt:
	wget -O $@ "ftp://ftp.informatics.jax.org/pub/reports/HOM_AllOrganism.rpt"

combined_fpkm: combine_fpkm.R HOM_AllOrganism.rpt $(FPKM_FILES)
	$(R) $(ROPTS) -f $< --args $@ $(wordlist 2,$(words $^),$^)

%.txt: %.xlsx
	$(SSCONVERT) --export-type Gnumeric_stf:stf_assistant -O 'separator="	" quote="" eol="unix"' $< $@
