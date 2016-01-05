#!/usr/bin/make

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

GTF_FILES=Dasypus_novemcinctus.Dasnov3.0.80.gtf \
	Callithrix_jacchus.C_jacchus3.2.1.80.gtf \
	Pan_troglodytes.CHIMP2.1.4.80.gtf \
	Bos_taurus.UMD3.1.80.gtf \
	Canis_familiaris.CanFam3.1.80.gtf \
	Loxodonta_africana.loxAfr3.80.gtf \
	Equus_caballus.EquCab2.80.gtf \
	Homo_sapiens.GRCh38.80.gtf \
	Mus_musculus.GRCm38.80.gtf \
	Monodelphis_domestica.BROADO5.80.gtf \
	Sus_scrofa.Sscrofa10.2.80.gtf \
	Ovis_aries.Oar_v3.1.80.gtf

all: $(GTF_FILES)

%.gtf: %.gtf.gz
	gzip -dc $< > $@

%.gtf.gz:
	rsync -avP "rsync://ftp.ensembl.org/ensembl/pub/release-80/gtf/$(shell echo $@|sed 's/\([^.]*\).*/\L\1/')/$@" $@

gene_id_to_type.txt: gene_id_to_type.pl $(GTF_FILES)
	perl $^ > $@
