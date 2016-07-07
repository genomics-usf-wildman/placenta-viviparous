#!/usr/bin/make -f

SPECIES=nannospalax_galili
ALIGNMENT_SPECIES=nannospalax_galili
GTF=nannospalax_galili.gtf
FASTA=nannospalax_galili.toplevel.fa

COMMON_MAKEFILE=../common_makefile
# this is not an ensembl assembly
NOTENSEMBL=1

NREADS=2

# we need to use trinity to do the alignment
USE_TRINITY=1
PROTEIN_ALIGNMENT_SPECIES=mus_musculus
PROTEIN_ALIGNMENT_FASTA=Mus_musculus.GRCm38.pep.all.fa

BOWTIE_INDEX_DIR=./

include ../srr_makefile


# get the fasta files; the perl regex is to just keep the ref name
# which is used in the gtf file
1026970_ref_S.galili_v1.0_chrUn.fa.gz: 
	rsync -avP \
	"rsync://ftp.ncbi.nih.gov/genomes/Nannospalax_galili/CHR_Un/$@" $@
1026970_ref_S.galili_v1.0_chrMT.fa.gz:
	rsync -avP \
	"rsync://ftp.ncbi.nih.gov/genomes/Nannospalax_galili/CHR_MT/$@" $@ 

nannospalax_galili.toplevel.fa: 1026970_ref_S.galili_v1.0_chrUn.fa.gz \
	1026970_ref_S.galili_v1.0_chrUn.fa.gz
	gzip -dc $^ | \
		perl -pe 's/^>gi\|\d+\|ref\|([^\|]+)\|.+/>$$1/; $$_;' > \
		$@
ref_S.galili_v1.0_top_level.gff3.gz:
	rsync -avP "rsync://ftp.ncbi.nih.gov/genomes/Nannospalax_galili/GFF/$@" $@

ref_S.galili_v1.0_top_level.gff3: ref_S.galili_v1.0_top_level.gff3.gz
	gzip -dc $< > $@

# build the gtf file from the gff3 file
nannospalax_galili.gtf: nannospalax_galili.toplevel.fa \
	ref_S.galili_v1.0_top_level.gff3
	$(MODULE) load cufflinks/2.2.1; \
	gffread -o $@ -T -g $< $(wordlist 2,2,$^)
