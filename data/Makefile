#!/usr/bin/make -f

SSCONVERT ?=ssconvert
SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data
PERL?=perl


FPKM_FILES=Supplemental_Table_2_Homo_NEW_UPDATED.txt $(patsubst %.xlsx,%.txt,$(filter-out $(wildcard *Homo*), $(wildcard *_fpkm.xlsx)))


all: combined_fpkm


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

oma-ensembl.txt.gz:
	wget -O $@ "http://omabrowser.org/All/oma-ensembl.txt.gz"

oma-groups.txt.gz:
	wget -O $@ "http://omabrowser.org/All/oma-groups.txt.gz"

oma-groups.txt: oma-groups.txt.gz
	gzip -dc $^ > $@

oma-ensembl.txt: oma-ensembl.txt.gz
	gzip -dc $^ > $@

oma-groups_long.txt: oma-group-to-long.pl oma-groups.txt oma-ensembl.txt
	$(PERL) $^ > $@

combined_fpkm: combine_fpkm.R HOM_AllOrganism.rpt oma-groups_long.txt fastq_fpkm/combined_fpkms ensembl_orthologs
	$(R) $(ROPTS) -f $< --args $@ $(wordlist 2,$(words $^),$^)

homology_table: homology_table.R combined_fpkm
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

ensembl_orthologs: ensembl_orthologs.R
	$(R) $(ROPTS) -f $< --args $@

ensembl_orthologs_long: ensembl_orthologs_long.R
	$(R) $(ROPTS) -f $< --args $@

Compara.80.protein.nhx.emf.gz:
	rsync -avP "rsync://ftp.ensembl.org/ensembl/pub/current_emf/ensembl-compara/homologies/$@" $@

gene_trees: parse_tree.R Compara.80.protein.nhx.emf.gz
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

%.txt: %.xlsx
	$(SSCONVERT) --export-type Gnumeric_stf:stf_assistant -O 'separator="	" quote="" eol="unix"' $< $@


### placenta specific transcripts and go analysis of them
placenta_core_transcriptome: placenta_core_transcriptome.R combined_fpkm housekeeping_genes_superset
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

placenta_go_results_mf: go_analysis.R placenta_core_transcriptome
	$(R) $(ROPTS) -f $< --args $(wordlist 2,2,$^) placenta.transcriptome.list MF $@

placenta_go_results_bp: go_analysis.R placenta_core_transcriptome
	$(R) $(ROPTS) -f $< --args $(wordlist 2,2,$^) placenta.transcriptome.list BP $@
