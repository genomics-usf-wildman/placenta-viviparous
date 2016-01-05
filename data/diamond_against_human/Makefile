#!/usr/bin/make -f

SHELL=/bin/bash

### module is how the biocluster loads specific versions; if we're not
### running there, we'll assume the correct version is installed and
### just echo what we're loading
ifdef MODULEPATH
MODULE=module
else
MODULE=echo
endif

ifdef PBS_NUM_PPN
CORES=$(PBS_NUM_PPN)
else
CORES=8
endif


SPECIES_PEP_FASTA=Dasypus_novemcinctus.Dasnov3.0.pep.all.fa \
	Pan_troglodytes.CHIMP2.1.4.pep.all.fa \
	Bos_taurus.UMD3.1.pep.all.fa \
	Canis_familiaris.CanFam3.1.pep.all.fa \
	Loxodonta_africana.loxAfr3.pep.all.fa \
	Equus_caballus.EquCab2.pep.all.fa \
	Mus_musculus.GRCm38.pep.all.fa \
	Monodelphis_domestica.BROADO5.pep.all.fa \
	Sus_scrofa.Sscrofa10.2.pep.all.fa \
	Ovis_aries.Oar_v3.1.pep.all.fa

ALIGNMENT_TARGET=Homo_sapiens.GRCh38.pep.all.fa

SPECIES_ALIGNMENT_FILES=$(patsubst %.pep.all.fa,%_diamond.txt,$(SPECIES_PEP_FASTA))


all: gene_to_human_gene.txt $(SPECIES_ALIGNMENT_FILES)

$(SPECIES_PEP_FASTA) $(ALIGNMENT_TARGET): %.pep.all.fa: %.pep.all.fa.gz
	gzip -dc $< |perl -pe 's/^(>\S+)(.+gene:)(\S+)/$${1}_$${3} $${2}$${3}/' > $@

%.pep.all.fa.gz:
	rsync -avP "rsync://ftp.ensembl.org/ensembl/pub/release-80/fasta/$(shell echo $@|sed 's/\([^.]*\).*/\L\1/')/pep/$@" $@

homo_sapiens.dmnd: $(ALIGNMENT_TARGET)
	$(MODULE) load diamond/0.7.9 ; \
	diamond makedb --threads $(CORES) \
		--in $@ \
		--db $(patsubst %.dmnd,%,$<)


%.daa: %.pep.all.fa homo_sapiens.dmnd
	$(MODULE) load diamond/0.7.9 ; \
	diamond blastp --threads $(CORES) \
		--query $< \
		--daa $@ \
		--db $(patsubst %.dmnd,%,$(wordlist 2,2,$^))

%_diamond.txt: %.daa
	$(MODULE) load diamond/0.7.9 ; \
	diamond view --threads $(CORES) -o $@ -a $<

gene_to_human_gene.txt: gene_to_human_gene.pl protein_id_to_gene_id.txt $(SPECIES_ALIGNMENT_FILES)
	perl $^ > $@
